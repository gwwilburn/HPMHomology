/* hmmalign_uniglocal.c */
/* based on example section of h4's dp_reference/reference_viterbi.c */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_msa.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "simple program to produce a uniglocal alignment to an hmm";

int
main(int argc, char **argv)
{
	ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
	char           *hmmfile = esl_opt_GetArg(go, 1);
 	char           *seqfile = esl_opt_GetArg(go, 2);
 	ESL_ALPHABET   *abc     = NULL;
 	P7_HMMFILE     *hfp     = NULL;
	P7_HMM         *hmm     = NULL;
	ESL_SQ        **sq      = NULL;                      /* array of sequences                    */
	ESL_SQFILE     *sqfp    = NULL;
	P7_BG          *bg      = NULL;                      /* background model                      */
	P7_PROFILE     *gm      = NULL;
	P7_REFMX       *vit     = p7_refmx_Create(200, 400); /* will grow as needed                   */
	P7_TRACE      **tr      = NULL;
	ESL_MSA        *msa     = NULL;	                    /* resulting multiple alignment          */
	int             outfmt  = eslMSAFILE_STOCKHOLM;
	float           vsc;                                 /* viterbi score                         */
	int             totseq  = 0;
	int             msaopts = 0;                         /* option flags for p7_tracealign_Seqs() */
	int             nseq;
	int             i;
	int             format  = eslSQFILE_UNKNOWN;
	int				 status;

   /* Read in one HMM */
   if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
   if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
   p7_hmmfile_Close(hfp);

	/* Open sequence file */
 	status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
	if      (status == eslENOTFOUND) p7_Fail("No such file.");
	else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
	else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
	else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

	/* read in sequences into array */
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
	esl_sqfile_Close(sqfp);
	totseq = nseq;


	/* allocate memory for traces */
	ESL_REALLOC(tr, sizeof(P7_TRACE *) * totseq);

	for (i = 0; i < totseq; i++) {
		tr[i] = p7_trace_Create();
	}

	/* Configure a profile from the HMM */
	bg = p7_bg_Create(abc);
	gm = p7_profile_Create(hmm->M, abc);
	/* Configure in uniglocal mode */
	p7_profile_ConfigUniglocal(gm, hmm, bg, 400);

	/* loop through sequence array and align seqs */
	for (i = 0; i < totseq; i++) {

		/* Set the profile and null model's target length models */
		p7_bg_SetLength     (bg, sq[i]->n);
		p7_profile_SetLength(gm, sq[i]->n);

		/* Run Viterbi - get raw score and optimal trace */
		p7_ReferenceViterbi(sq[i]->dsq, sq[i]->n, gm, vit, tr[i], &vsc);

		/* prepare dp matrix for next iteration */
		p7_refmx_Reuse(vit);
	}

	/* create MSA from traces and seqs */
	msaopts |= p7_ALL_CONSENSUS_COLS; /* include all consensus columns in alignment */
	p7_tracealign_Seqs(sq, tr, totseq, hmm->M, msaopts, hmm, &msa);

	/* write MSA to stdout */
	esl_msafile_Write(stdout, msa, outfmt);

	/* clean up and return */
 	p7_hmm_Destroy(hmm);
 	esl_alphabet_Destroy(abc);
	esl_msa_Destroy(msa);
 	esl_getopts_Destroy(go);
	free(sq);
	free(tr);
	//esl_sq_Destroy(sq);
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	return 0;

	ERROR:
		return status;

}
