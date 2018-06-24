/* hmmscore_uniglocal.c */
/* inspired by example in HMMER's generic_fwdback.c */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "calculate forward and backward scores with a glocal alignment method";


int main(int argc, char **argv) {
	ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
   char           *hmmfile = esl_opt_GetArg(go, 1);
   char           *seqfile = esl_opt_GetArg(go, 2);
	ESL_ALPHABET   *abc     = NULL;
	P7_HMMFILE     *hfp     = NULL;
	P7_HMM         *hmm     = NULL;
	ESL_SQ         *sq      = NULL;
	ESL_SQFILE     *sqfp    = NULL;
	P7_BG          *bg      = NULL;
	P7_PROFILE     *gm      = NULL;
	P7_REFMX       *fwd     = p7_refmx_Create(100, 100);
	P7_REFMX       *bck     = p7_refmx_Create(100, 100);
	P7_REFMX       *pp      = p7_refmx_Create(100, 100);
	float           fsc;
	float           bsc;
	float           nullsc;
	int             format  = eslSQFILE_UNKNOWN;
	int             status;



   /* Read in one HMM */
   if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
   if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
   p7_hmmfile_Close(hfp);

	/* Read in one sequence */
   sq     = esl_sq_CreateDigital(abc);
	status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
	/* handle possible errors w/ seq file reading */
	if      (status == eslENOTFOUND) p7_Fail("No such file.");
	else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
	else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
	else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

	/* Configure a profile from the HMM */
	bg = p7_bg_Create(abc);
	gm = p7_profile_Create(hmm->M, abc);

	/* Configure model in unihit glocal mode */
	//p7_ProfileConfig(hmm, bg, gm, sq->n, p7_UNIGLOCAL);  \\ h3 function
	if (p7_profile_ConfigUniglocal(gm, hmm, bg, sq->n) != eslOK) esl_fatal("failed to configure profile");

	/* print header to .csv output */
	fprintf(stdout, "id, fwd (raw), bck (raw), fwd (bits), bck (bits)\n");

	/* loop through sequences */

	while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF) {
		if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
		else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);



		 /* Set the profile and null model's target length models */
		 p7_bg_SetLength     (bg, sq->n);
		 p7_profile_SetLength(gm, sq->n);

		  /* Run Forward, Backward, Decoding;
       * after decoding, <bck> becomes the pp matrix;
       * after alignment, <fwd> becomes the OA alignment DP matrix
       */
		 p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, &fsc);
		 p7_ReferenceBackward(sq->dsq, sq->n, gm, bck, &bsc);
		 p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);

		 /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

		/* print out scores */
		fprintf(stdout, "%s,%10.4f,%10.4f,%10.4f,%10.4f\n",
			  	  sq->name,
			  	  fsc,
				  bsc,
			     (fsc - nullsc) / eslCONST_LOG2,
			     (bsc - nullsc) / eslCONST_LOG2);

		/* prepare objects for next iteration */
		p7_refmx_Reuse(fwd);
		p7_refmx_Reuse(bck);
		p7_refmx_Reuse(pp);
		esl_sq_Reuse(sq);
	}


	/* Clean up and return */
	esl_sqfile_Close(sqfp);
	esl_sq_Destroy(sq);
	p7_hmm_Destroy(hmm);
	esl_alphabet_Destroy(abc);
   p7_refmx_Destroy(pp);
   p7_refmx_Destroy(fwd);
   p7_refmx_Destroy(bck);
   p7_profile_Destroy(gm);
   p7_bg_Destroy(bg);
	esl_getopts_Destroy(go);



	return 0;
}

