/* hmmscore_uniglocal.c */
/* inspired by example in HMMER's generic_fwdback.c */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "hmm_scoreset.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile> <scorefile>";
static char banner[] = "calculate forward and backward scores with a glocal alignment method";


int main(int argc, char **argv) {
	ESL_GETOPTS    *go        = p7_CreateDefaultApp(options, 3, argc, argv, banner, usage);
   char           *hmmfile   = esl_opt_GetArg(go, 1);
   char           *seqfile   = esl_opt_GetArg(go, 2);
	char           *scorefile = esl_opt_GetArg(go, 3);
	ESL_ALPHABET   *abc       = NULL;
	P7_HMMFILE     *hfp       = NULL;
	P7_HMM         *hmm       = NULL;
	ESL_SQ        **sq        = NULL;                      /* array of sequences */
	ESL_SQFILE     *sqfp      = NULL;
	P7_BG          *bg        = NULL;
	P7_PROFILE     *gm        = NULL;
	P7_REFMX       *fwd       = p7_refmx_Create(100, 100);
	P7_REFMX       *bck       = p7_refmx_Create(100, 100);
	P7_REFMX       *vit       = p7_refmx_Create(200, 400); /* will grow as needed                   */
	P7_TRACE       *tr        = p7_trace_Create();
	HMM_SCORESET   *hmm_ss    = NULL;
	FILE				*hmm_ss_fp = NULL;
	float           fsc;
	float           bsc;
	float           vsc;
	float           nullsc;
	int             totseq = 0;
	int             nseq;
	int             i;
	int             format  = eslSQFILE_UNKNOWN;
	int             status;



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
	else if (status != eslEOF) esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);
	totseq = nseq;

	/* Configure a profile from the HMM */
	bg = p7_bg_Create(abc);
	gm = p7_profile_Create(hmm->M, abc);



	/* Configure model in unihit glocal mode */
	//p7_ProfileConfig(hmm, bg, gm, sq->n, p7_UNIGLOCAL);  \\ h3 function
	if (p7_profile_ConfigUniglocal(gm, hmm, bg, 400) != eslOK) esl_fatal("failed to configure profile");

	/* set up hmm scoreset object */
	hmm_ss = hmm_scoreset_Create(nseq);

	for (i = 0; i < totseq; i++) {

		if (i % 1000 == 0) fprintf(stdout, "%d\n", i);

		/* Set the profile and null model's target length models */
		p7_bg_SetLength     (bg, sq[i]->n);
		p7_profile_SetLength(gm, sq[i]->n);


		/* get fwd and bck scores */
		p7_ReferenceForward (sq[i]->dsq, sq[i]->n, gm, fwd, &fsc);
		p7_ReferenceBackward(sq[i]->dsq, sq[i]->n, gm, bck, &bsc);


		/* get viterbi score */
		p7_ReferenceViterbi(sq[i]->dsq, sq[i]->n, gm, vit, tr, &vsc);

		/* get null score */
		p7_bg_NullOne(bg, sq[i]->dsq, sq[i]->n, &nullsc);

		/* add scores to score file */
		hmm_ss->sqname[i] = sq[i]->name;
		hmm_ss->fsc[i]    = fsc;
		hmm_ss->bsc[i]    = bsc;
		hmm_ss->vsc[i]    = vsc;
		hmm_ss->nullsc[i] = nullsc;

		/* prepare objects for next iteration */
		p7_refmx_Reuse(fwd);
		p7_refmx_Reuse(bck);
		p7_refmx_Reuse(vit);
		p7_trace_Reuse(tr);
	}

	/* write scorefile */
	if ((hmm_ss_fp = fopen(scorefile, "w")) == NULL) esl_fatal("Failed to open output hmm score setfile %s for writing", scorefile);
	hmm_scoreset_Write(hmm_ss_fp, hmm_ss);
	fclose(hmm_ss_fp);

	/* Clean up and return */
	esl_sqfile_Close(sqfp);
	free(sq);
	p7_hmm_Destroy(hmm);
	esl_alphabet_Destroy(abc);
   p7_refmx_Destroy(fwd);
   p7_refmx_Destroy(bck);
   p7_refmx_Destroy(vit);
   p7_profile_Destroy(gm);
	free(tr);
	free(hmm_ss);
   p7_bg_Destroy(bg);
	esl_getopts_Destroy(go);

	ERROR:
		return status;



	return 0;
}

