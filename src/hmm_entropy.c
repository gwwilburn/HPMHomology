/* hmm_entropy.c */
#include "hmm_entropy.h"
#include "p7_config.h"
#include "hmmer.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"
#include "potts.h"


float
hmm_entropy_Calculate(ESL_SQ *sq, P7_PROFILE *gm, float *ret_H)
{

	return eslOK;
}


#ifdef HPM_3MER
#include "p7_config.h"
#include "esl_alphabet.h"
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static ESL_OPTIONS options[] = {
	/* name           type      default  env  range  toggles reqs incomp  help                                    docgroup*/
	{ "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",         0 },
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of calculating posterior path entropy for a given hmm and sequence.";

int main(int argc, char *argv[])
{
	ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
	char          *hmmfile = esl_opt_GetArg(go, 1);
	char          *seqfile = esl_opt_GetArg(go, 2);
	ESL_ALPHABET  *abc     = NULL;
	P7_HMMFILE    *hfp     = NULL;
	P7_HMM        *hmm     = NULL;
	ESL_SQ       **sq      = NULL;                      /* array of sequences */
	ESL_SQFILE    *sqfp    = NULL;
	P7_BG         *bg      = NULL;
	P7_PROFILE    *gm      = NULL;
	P7_REFMX      *fwd     = p7_refmx_Create(100, 100);
	int            format  = eslSQFILE_UNKNOWN;
	int            totseq  = 0;
	int            nseq;
	int            status;
	int            i;
	float          fsc;
	float          nullsc;
	float          H       = 0;

	/* Read in one HMM */
	if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
	p7_hmmfile_Close(hfp);


	/* open sequence file */
	status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
	if      (status == eslENOTFOUND) p7_Fail("No such file.");
	else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
	else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
	else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

	/* read sequences into array */
	ESL_REALLOC(sq, sizeof(ESL_SQ *) * (totseq + 1));
	sq[totseq] = esl_sq_CreateDigital(abc);
	nseq = 0;
	while ((status = esl_sqio_Read(sqfp, sq[totseq+nseq])) == eslOK) {
		nseq++;
		ESL_REALLOC(sq, sizeof(ESL_SQ *) * (totseq+nseq+1));
		sq[totseq+nseq] = esl_sq_CreateDigital(abc);
	}
	/* error handling and cleanup */
	if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
			                                   sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
	else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);
	totseq = nseq;


	/* configure background model */
	bg = p7_bg_Create(abc);

	/* Configure a profile from the HMM in uniglocal mode*/
	gm = p7_profile_Create(hmm->M, abc);
	if (p7_profile_ConfigUniglocal(gm, hmm, bg, 400) != eslOK) esl_fatal("failed to configure profile");

	/* calculate posterior entropy for each sequence */
	for (i=0; i < totseq; i++) {
		if (i % 1000 == 0) fprintf(stdout, "%d\n", i);

		/* Set the profile and null model's target length models */
		p7_bg_SetLength     (bg, sq[i]->n);
		p7_profile_SetLength(gm, sq[i]->n);

		/* get fwd matrix */
		p7_ReferenceForward (sq[i]->dsq, sq[i]->n, gm, fwd, &fsc);

		p7_bg_NullOne(bg, sq[i]->dsq, sq[i]->n, &nullsc);

		/* calculate H(pi | x) */
		hmm_entropy_Calculate(sq[i], gm, &H);

		/* reuse fwd matrix for next iteration */
		p7_refmx_Reuse(fwd);
	}

	fprintf(stdout, "hello world, I am entropy!\n");

	/* clean up and return */
	esl_getopts_Destroy(go);
	esl_alphabet_Destroy(abc);
	esl_sqfile_Close(sqfp);
	p7_hmm_Destroy(hmm);
	free(sq);
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);

	return 0;

	ERROR:
		return status;

}

#endif
