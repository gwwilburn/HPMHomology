/* hpmscoreIS.c */
/* infant stages of scoring sequences with Potts models using importance samplign and HMM's */

#include "p7_config.h"
#include "easel.h"
#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"
#include "hpm.h"
#include "hpmfile.h"
#include "hmm_entropy.h"
#include "hpm_scoreset.h"
#include "hpm_scoreops.h"

#include "hmmer.h"

/* declaration of internal functions */
int Calculate_IS_scores(HPM *hpm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int totseq, HPM_IS_SCORESET *hpm_is_ss, int verbose);


static ESL_OPTIONS options[] = {
	/* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
	{ "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL,           "help; show brief info on version and usage",                 0 },
	{ "-s",         eslARG_INT,     "0", NULL, NULL, NULL, NULL, NULL,           "set random number seed to <n>",                              0 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<seqfile> contains protein sequences",                       1 },
   { "--rna",      eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--amino",  "<seqfile> contains RNA sequences",                           1 },
   { "--dna",      eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--rna,--amino",  "<seqfile> contains DNA sequences",                           1 },

	/* debugging tools */
   { "--v",        eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,NULL,              "Verbose mode: print info on intermediate scoring steps",    2 },
{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hpmfile> <hmmfile> <seqfile> <scoreoutfile>";
static char banner[] = "Score unaligned sequences with an HPM using importance sampling";

int main(int argc, char **argv){

	ESL_GETOPTS      *go            = p7_CreateDefaultApp(options, 4, argc, argv, banner, usage);
	ESL_RANDOMNESS   *rng           = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET     *abc;
	char             *hpmfile       = esl_opt_GetArg(go, 1);
	char             *hmmfile       = esl_opt_GetArg(go, 2);
	char             *seqfile       = esl_opt_GetArg(go, 3);
	char             *scorefile     = esl_opt_GetArg(go, 4);
	HPM              *hpm;
	P7_HMMFILE       *hmmfp         = NULL;
	P7_HMM           *hmm           = NULL;
	ESL_SQ          **sq            = NULL;                      /* array of sequences */
	ESL_SQFILE       *sqfp          = NULL;
	HPM_IS_SCORESET  *hpm_is_ss     = NULL;
	FILE				  *hpm_is_ss_fp  = NULL;
	int               format        = eslSQFILE_UNKNOWN;
	int               totseq        = 0;
	int               nseq          = 0;
	int               v             = 0;                         /* Boolean for verbose mode */
	int               status;
	char              errbuf[eslERRBUFSIZE];

	/* if user has defined an alphabet we define it here */
	abc = NULL;
	if        (esl_opt_GetBoolean(go, "--amino"))		abc = esl_alphabet_Create(eslAMINO);
	else if   (esl_opt_GetBoolean(go, "--dna"))        abc = esl_alphabet_Create(eslDNA);
	else if   (esl_opt_GetBoolean(go, "--rna"))        abc = esl_alphabet_Create(eslRNA);

	/* check for verbose mode */
	if (esl_opt_GetBoolean(go, "--v")) v=1;

	/* read in hmm, autodetect alphabet if not manually set by user  */
	if (p7_hmmfile_OpenE(hmmfile, NULL, &hmmfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hmmfp, &abc, &hmm)        != eslOK) p7_Fail("Failed to read HMM");
	p7_hmmfile_Close(hmmfp);

	/* open hpm */
	hpm = hpmfile_Read(hpmfile, abc, errbuf);

	/* open sequence file */
	status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
	if      (status == eslENOTFOUND) p7_Fail("No such file.");
	else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
	else if (status ==eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
	else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

	/* read sequences into array */
	ESL_REALLOC(sq, sizeof(ESL_SQ *) * (totseq + 1));
	sq[totseq] = esl_sq_CreateDigital(abc);
	while ((status = esl_sqio_Read(sqfp, sq[totseq+nseq])) == eslOK) {
		nseq++;
		ESL_REALLOC(sq, sizeof(ESL_SQ *) * (totseq+nseq+1));
		sq[totseq+nseq] = esl_sq_CreateDigital(abc);

	}

	/* error handling and cleanup */
	if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
			                                   sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
	else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
			                                    status, sqfp->filename);

	totseq = nseq;

	/* set up score set object */
	hpm_is_ss = hpm_is_scoreset_Create(nseq);


	/* calculate importance sampling score for all seqs */
	Calculate_IS_scores(hpm, hmm, sq, rng, totseq, hpm_is_ss, v);

	/* write hpm is scores to output csv */
	if ((hpm_is_ss_fp = fopen(scorefile, "w")) == NULL) esl_fatal("Failed to open output hpm IS scoreset file %s for writing", scorefile);
	hpm_is_scoreset_Write(hpm_is_ss_fp, hpm_is_ss);
	fclose(hpm_is_ss_fp);

	/* clean up and return */
	esl_sqfile_Close(sqfp);
	free(sq);
	p7_hmm_Destroy(hmm);
	free(hpm);
	esl_alphabet_Destroy(abc);
	esl_getopts_Destroy(go);


	ERROR:
		return status;

	return 0;
}


/* outer loop over all sequences in seqfile */
int Calculate_IS_scores(HPM *hpm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int totseq, HPM_IS_SCORESET *hpm_is_ss, int verbose) {


	P7_BG          *bg        = NULL;
	P7_PROFILE     *gm        = NULL;

	P7_REFMX       *fwd       = p7_refmx_Create(100, 100);
	P7_TRACE       *tr        = p7_trace_Create();
	float           fsc;                                     /* partial forward log odds score; output *
																				* calculated in p7_ReferenceForward()    */
	float          *wrk       = NULL;

	int             i;
	int             r;                                      /* importance sample index                */
	int             R         = 1;                          /* total number of samples per sequence   */
	float           sc_ld;   				                    /* ln( Q( x, \pi) ) under an hpm          */
	float           hsc;
	float           esc;
	float           nesc;
	float           nmesc;
	float           ntsc;
	float           tsc;
	double          pr[R];
	float           H;
	int             status;


	/* create a profile from the HMM */
	gm = p7_profile_Create(hmm->M, hmm->abc);
	/* create null model*/
	bg = p7_bg_Create(hmm->abc);
	/* configure model in uniglocal mode */
	if (p7_profile_ConfigUniglocal(gm, hmm, bg, 400) != eslOK) esl_fatal("failed to configure profile");

	/* outer loop over all sequences in seqfile */
	for (i = 0; i < totseq; i++) {
	//for (i = 0; i < 1; i++) {
		if (i % 100 == 0) fprintf(stdout, "%d\n", i);

		/* Set the profile and null model's target length models */
		p7_profile_SetLength(gm, sq[i]->n);

		/* calculate forward dp matrix, forward score */
		p7_ReferenceForward(sq[i]->dsq, sq[i]->n, gm, fwd, &fsc);

		/* calculate null emission scores */
		hpm_scoreops_ScoreNullEmissions(hpm, sq[i], &nesc);

		/* calculate null transition scores */
		hpm_scoreops_ScoreNullTransitions(bg, sq[i], &ntsc);

		/* calculate posterior entropy H(pi | x) */
		hmm_entropy_Calculate(gm, fwd, &H, 0);

		/* inner loop over sampled paths */
		for (r = 0; r < R; r++) {
			if (verbose) fprintf(stdout, "\nsequence %s, stochastic trace %d\n", sq[i]->name, r);

			/* perform stochastic traceback */
			p7_reference_trace_Stochastic(rng, &wrk, gm, fwd, tr);
			/* Calculate Q(x,\pi) under hmm */
			p7_trace_Score(tr, sq[i]->dsq, gm, &sc_ld);
			/* calculate Potts Hamiltonian */
			hpm_scoreops_CalculateHamiltonian(hpm, tr, sq[i]->dsq, &hsc, &esc);
			/* score match state insert emissions */
			hpm_scoreops_ScoreNullMatchEmissions(hpm, tr, sq[i]->dsq, &nmesc);
			/* score transitions */
			hpm_scoreops_ScoreTransitions(hpm, tr, sq[i]->dsq, sq[i]->L, &tsc);

			/* calculate log of this samples contribution to importance sampling sum */
			pr[r] = hsc + esc + tsc - sc_ld - nmesc;

			if (verbose) {
				fprintf(stdout, "# Q(x,pi) details\n");
				p7_trace_DumpAnnotated(stdout, tr, gm, sq[i]->dsq);
				fprintf(stdout, "\n");
				fprintf(stdout, "# length of traceback: %d\n", tr->N);
				fprintf(stdout, "# partial log odds score of x, pi: %.4f\n", sc_ld);
				fprintf(stdout, "# local field hamiltonian sum: %.4f\n", hsc);
				fprintf(stdout, "# coupling hamiltonian sum: %.4f\n", esc);
				fprintf(stdout, "# HPM log odds transitions: %.4f\n", tsc);
				fprintf(stdout, "# Null match state emission log prob: %.4f\n", nmesc);
				fprintf(stdout, "# Contribution to importance sampling sum: %.4f\n", pr[r]);
				fprintf(stdout, "\n");
				//fprintf(stdout, "\t tr->N hpm->ld, hsc, esc, nmesc, tsc\n");
				//fprintf(stdout, "\t%d, %.2f, %.2f, %.2f, %2f, %.2f %.2f,\n", tr->N, sc_ld, hsc, esc, nmesc, tsc);
			}

			p7_trace_Reuse(tr);

		}

		float ls = esl_vec_DLogSum(pr, R);
		p7_refmx_Reuse(fwd);
		float ld = (fsc - logf(R) - ntsc + ls) / eslCONST_LOG2;
		if (verbose) {
			fprintf(stdout, "\n\n## final details for sequence %s\n", sq[i]->name);
			fprintf(stdout, "## id, fsc, ntsc, fwd_ld, H(pi | x), hpm_ld\n");
			fprintf(stdout, "## %s, %.2f, %.2f, %.2f,  %.2f, %.2f\n", sq[i]->name, fsc, ntsc, (fsc-ntsc) / eslCONST_LOG2, H, ld);
		}
		/* add scoring info to scoreset object */
		hpm_is_ss->sqname[i]  = sq[i]->name;
		hpm_is_ss->R[i]       = R;
		hpm_is_ss->H[i]       = H;
		hpm_is_ss->fwd[i]     = fsc-ntsc;
		hpm_is_ss->is_ld[i]   = ld;
	}


	/* clean up */
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	p7_refmx_Destroy(fwd);
	free(tr);
	free(wrk);

	return eslOK;

	ERROR:
		return status;

}

