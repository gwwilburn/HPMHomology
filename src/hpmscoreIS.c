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
#include "hpm_trace.h"

#include "hmmer.h"

#include <omp.h>

/* declaration of internal functions */
int Calculate_IS_scores(HPM *hpm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int totseq, HPM_IS_SCORESET *hpm_is_ss, int start, int end, ESL_MSA **msa, int A, int verbose);


static ESL_OPTIONS options[] = {
	/* name         type         default  env range  togs   reqs  incomp         help                                                   docgroup */
	{ "-h",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL,            "help; show brief info on version and usage",              0 },
	{ "-s",         eslARG_INT,      "0", NULL, NULL, NULL, NULL, NULL,            "set random number seed to <n>",                           0 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,"--dna,--rna",    "<seqfile> contains protein sequences",                    1 },
   { "--rna",      eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,"--dna,--amino",  "<seqfile> contains RNA sequences",                        1 },
   { "--dna",      eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL,"--rna,--amino",  "<seqfile> contains DNA sequences",                        1 },

	/* options for bounding sequences we score */
   { "--msastart", eslARG_INT,     "-1", NULL, NULL, NULL, NULL, NULL,            "Start sequence",                                          2 },
   { "--msaend",   eslARG_INT,     "-2", NULL, NULL, NULL, NULL, NULL,            "End sequence",                                            2 },


	/* control of output */
	{ "-A",         eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, NULL,            "save multiple alignment of all hits to file <s>",         3 },

	/* debugging tools */
   { "--v",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL,            "Verbose mode: print info on intermediate scoring steps",  4 },
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
	FILE             *afp           = NULL;                      /* alignment output file (-A)   */
	ESL_MSA          *msa           = NULL;                      /* alignment output object (-A) */
	int               outfmt        = eslMSAFILE_STOCKHOLM;      /* alignment output format (-A) */
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
	int               A             = 0;                         /* Boolean for output alignment */
	int               start;
	int               end;
	int               status;
	char              errbuf[eslERRBUFSIZE];

	/* if user has defined an alphabet we define it here */
	abc = NULL;
	if        (esl_opt_GetBoolean(go, "--amino"))		abc = esl_alphabet_Create(eslAMINO);
	else if   (esl_opt_GetBoolean(go, "--dna"))        abc = esl_alphabet_Create(eslDNA);
	else if   (esl_opt_GetBoolean(go, "--rna"))        abc = esl_alphabet_Create(eslRNA);

	/* check for verbose mode */
	if (esl_opt_GetBoolean(go, "--v")) v=1;

	/* if output msa requested by user, try to open it */
	if (esl_opt_IsOn(go, "-A")) {
		if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A"));
		A = 1;
  	}

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

	/* if user specified bounds, set them here */
	start = 0;
	end   = nseq;

	if (esl_opt_GetInteger(go, "--msastart") > -1) start = esl_opt_GetInteger(go, "--msastart");
	if (esl_opt_GetInteger(go, "--msaend") > -1)     end = esl_opt_GetInteger(go, "--msaend");
	if (end < start)  esl_fatal("msa start > msa end: %d, %d\n", start, end);

	/* total number of sequences we will be scoring */
	totseq = end-start;

	/* set up score set object */
	hpm_is_ss = hpm_is_scoreset_Create(totseq);


	/* calculate importance sampling score for all seqs */
	Calculate_IS_scores(hpm, hmm, sq, rng, totseq, hpm_is_ss, start, end, &msa, A, v);

	/* write hpm is scores to output csv */
	if ((hpm_is_ss_fp = fopen(scorefile, "w")) == NULL) esl_fatal("Failed to open output hpm IS scoreset file %s for writing", scorefile);
	hpm_is_scoreset_Write(hpm_is_ss_fp, hpm_is_ss);
	fclose(hpm_is_ss_fp);

	/* if output msa requested, write it to otput file */
	if (A) {
		esl_msafile_Write(afp, msa, outfmt);
		fclose(afp);
		esl_msa_Destroy(msa);
	}

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
int Calculate_IS_scores(HPM *hpm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int totseq, HPM_IS_SCORESET *hpm_is_ss, int start, int end, ESL_MSA **msa, int A, int verbose) {


	P7_BG          *bg        = NULL;
	P7_PROFILE     *gm        = NULL;

	P7_REFMX       *fwd       = p7_refmx_Create(100, 100);
	P7_TRACE       *tr        = p7_trace_Create();
	float           fsc;                                     /* partial forward log odds score; output *
																				* calculated in p7_ReferenceForward()    */
	float          *wrk       = NULL;

	P7_TRACE      **out_tr    = NULL;                       /* traces used to construct output msa  (-A) */
	ESL_SQ        **out_sq    = NULL;                       /* traces used to construct output msa  (-A) */
	int             msaopts   = 0;
	int             i,j;
	int             r;                                      /* importance sample index                */
	int             R         = 1;                          /* total number of samples per sequence   */
	float           sc_ld;   				                    /* ln( Q( x, \pi) ) under an hpm          */
	float           hsc;
	float           esc;
	float           nesc;
	float           nmesc;
	float           ntsc;
	float           tsc;
	float           isc;
	float           hpmsc;
	float           hpmsc_max;
	float           H;
	int             status;


	/* create a profile from the HMM */
	gm = p7_profile_Create(hmm->M, hmm->abc);
	/* create null model*/
	bg = p7_bg_Create(hmm->abc);
	/* configure model in uniglocal mode */
	if (p7_profile_ConfigUniglocal(gm, hmm, bg, 400) != eslOK) esl_fatal("failed to configure profile");

	if (A) {
		/* allocate memory for traces  and initialize */
		ESL_ALLOC(out_tr, sizeof(P7_TRACE *) * totseq);
		ESL_ALLOC(out_sq, sizeof(ESL_SQ *) *   (totseq+1));
	}

	/* outer loop over all sequences in seqfile */
	for (i = start; i < end; i++) {

	   /* index for score set: must start at 0 */
		j=i-start;

		/* for getting optimal path under hpm */
		if (A) {
			hpmsc_max = -eslINFINITY;
			out_sq[j] = esl_sq_CreateDigital(hmm->abc);
			esl_sq_Copy(sq[i], out_sq[j]);
		}


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
		R = pow(2,ceil(H));
		if (H > 17.0)
		{
			R = 130000;  /* approx 2^17 */
		}
		else if (R < 1000) {
			R = 1000;
		}
		if (verbose) fprintf(stdout, "seq: %s, H: %.4f, R: %d\n", sq[i]->name, H, R);

		double pr[R];


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

			/* for output alignment */
			if (A) {

				/* score insert state  emissions */
				hpm_scoreops_ScoreInsertEmissions(hpm, tr, sq[i]->dsq, &isc);

				/* calculate numerator of log-odds score */

				hpmsc = hsc + esc + isc + tsc;

				/* see if we've got a better path */
				if (hpmsc > hpmsc_max) {
					hpmsc_max = hpmsc;
					out_tr[j] = hpm_trace_Clone(tr);
				}
			}


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
		hpm_is_ss->sqname[j]  = sq[i]->name;
		hpm_is_ss->R[j]       = R;
		hpm_is_ss->H[j]       = H;
		hpm_is_ss->fwd[j]     = fsc-ntsc;
		hpm_is_ss->is_ld[j]   = ld;
	}

	/* if output msa requesed, create MSA from trace */
	if (A) {
		msaopts |= p7_ALL_CONSENSUS_COLS; /* include all consensus columns in alignment */
		p7_tracealign_Seqs(out_sq, out_tr, totseq, hmm->M, msaopts, hmm, msa);
	}

	/* clean up */
	if (A) {
		for (j=0; j<totseq; j++) {
			p7_trace_Destroy(out_tr[j]);
			esl_sq_Destroy(out_sq[j]);
		}
		free(out_tr);
		free(out_sq);
	}
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	p7_refmx_Destroy(fwd);
	p7_trace_Destroy(tr);
	free(wrk);

	return eslOK;

	ERROR:
		return status;

}

