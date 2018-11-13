/* hpmscoreIS.c */
/* infant stages of scoring sequences with Potts models using importance samplign and HMM's */

#include "p7_config.h"
#include "easel.h"
#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "hpm.h"
#include "hpmfile.h"
#include "hpmscoreIS.h"

#include "hmmer.h"


static ESL_OPTIONS options[] = {
	/* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
	{ "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL,           "help; show brief info on version and usage",           0 },
	{ "-s",         eslARG_INT,     "0", NULL, NULL, NULL, NULL, NULL,           "set random number seed to <n>",                        0 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<seqfile> contains protein sequences",                 1 },
   { "--rna",      eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--amino",  "<seqfile> contains RNA sequences",                     1 },
   { "--dna",      eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--rna,--amino",  "<seqfile> contains DNA sequences",                     1 },

{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hpmfile> <hmmfile> <seqfile> <scoreoutfile>";
static char banner[] = "Score unaligned sequences with an HPM using importance sampling";

int main(int argc, char **argv){

	ESL_GETOPTS    *go        = p7_CreateDefaultApp(options, 4, argc, argv, banner, usage);
	ESL_RANDOMNESS *rng       = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET   *abc;

	char           *hpmfile   = esl_opt_GetArg(go, 1);
	char           *hmmfile   = esl_opt_GetArg(go, 2);
	char           *seqfile   = esl_opt_GetArg(go, 3);
	char           *scorefile = esl_opt_GetArg(go, 3);

	HPM            *hpm;

	P7_HMMFILE     *hmmfp     = NULL;
	P7_HMM         *hmm       = NULL;

	ESL_SQ        **sq        = NULL;                      /* array of sequences */
	ESL_SQFILE     *sqfp      = NULL;

	int             format    = eslSQFILE_UNKNOWN;
	int             totseq    = 0;
	int             nseq      = 0;
	int             i;

	int             status;
	char            errbuf[eslERRBUFSIZE];

	/* if user has defined an alphabet we define it here */
	abc = NULL;
	if        (esl_opt_GetBoolean(go, "--amino"))		abc = esl_alphabet_Create(eslAMINO);
	else if   (esl_opt_GetBoolean(go, "--dna"))        abc = esl_alphabet_Create(eslDNA);
	else if   (esl_opt_GetBoolean(go, "--rna"))        abc = esl_alphabet_Create(eslRNA);

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
	else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
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

	/* calculate importance sampling score for all seqs */
	Calculate_IS_scores(hpm, hmm, sq, rng, totseq);

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
int Calculate_IS_scores(HPM *hpm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int totseq) {


	P7_BG          *bg        = NULL;
	P7_PROFILE     *gm        = NULL;

	P7_REFMX       *fwd       = p7_refmx_Create(100, 100);
	P7_TRACE       *tr        = p7_trace_Create();
	float           fsc;
	float          *wrk       = NULL;

	int             i;
	int             r;                                      /* importance sample index                */
	int             R         = 1;                          /* total number of samples per sequence   */
	float           sc;   				                       /* ln( Q( x, \pi) ) under an hpm          */
	float           hsc;
	float           esc;
	float           isc;
	float           tsc;

	/* create a profile from the HMM */
	gm = p7_profile_Create(hmm->M, hmm->abc);
	/* create null model*/
	bg = p7_bg_Create(hmm->abc);
	/* configure model in uniglocal mode */
	if (p7_profile_ConfigUniglocal(gm, hmm, bg, 400) != eslOK) esl_fatal("failed to configure profile");


	fprintf(stdout, "tot seq: %d\n", totseq);
	/* outer loop over all sequences in seqfile */
	//for (i = 0; i < totseq; i++) {
	for (i = 0; i < 1; i++) {
		if (i % 100 == 0) fprintf(stdout, "%d\n", i);


		/* Set the profile and null model's target length models */
		p7_bg_SetLength     (bg, sq[i]->n);
		p7_profile_SetLength(gm, sq[i]->n);

		/* calculate forward dp matrix, forward score */
		p7_ReferenceForward(sq[i]->dsq, sq[i]->n, gm, fwd, &fsc);
 		fprintf(stdout, "%s: %.2f\n", sq[i]->name, fsc);


		/* inner loop over sampled paths */
		for (r = 0; r < R; r++) {
			fprintf(stdout, "%d %d\n", i, r);
			/* perform stochastic traceback */
			p7_reference_trace_Stochastic(rng, &wrk, gm, fwd, tr);
			/* Calculate Q(x,\pi) under hmm */
			p7_trace_Score(tr, sq[i]->dsq, gm, &sc);
			/* calculate Potts Hamiltonian */
			hpmscore_CalculateHamiltonian(hpm, tr, sq[i]->dsq, &hsc, &esc);
			/* score insert emissions */
			hpmscore_ScoreInsertions(hpm, tr, sq[i]->dsq, &isc);
			/* score transitions */
			hpmscore_ScoreTransitions(hpm, tr, sq[i]->dsq, sq[i]->L, &tsc);

			fprintf(stdout, "%d, %f\n", tr->N, sc);
			p7_trace_Reuse(tr);

		}

		p7_refmx_Reuse(fwd);
	}


	/* clean up */
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	p7_refmx_Destroy(fwd);
	free(tr);
	free(wrk);

	return eslOK;
}


int hpmscore_CalculateHamiltonian(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_hsc, float *ret_esc) {
	int   z;         /* index for trace elements 	               */
	int   y;         /* index for trace elements 	               */
	int   i;         /* match state index                          */
	int   j;         /* match state index                          */
	int   a;         /* residue index            	               */
	int   b;         /* residue index                              */
	int   idx;       /* index for potts parameters                 */
	float hsc = 0;   /* Hamiltonian contribution from h_i's        */
	float esc = 0;   /* Hamiltonian contribution from e_ij's       */


	int K = hpm->abc->K;

	for (z = 0; z < tr->N; z++) {

		fprintf(stdout, "%d\n", z);
		/* check if we are in a match or delete state */
		if (tr->st[z] == 2 || tr->st[z] == 6) {
			i = tr->k[z];

			/* we have a match position */
			if (tr->st[z] == 2) {
				a = dsq[tr->i[z]];
				if (a > K) a = K;

				//fprintf(stdout, "%d\n", a);
			}

			/* we have a delete position */
			else if (tr->st[z] == 6) {
				a = K;
			}
			hsc += hpm->h[i][a];

			/* now add e_ij terms to pseudo-energy */
			for (y = z+1; y < tr->N; y++) {

				/* check if we are in a match position */
				if (tr->st[y] == 2 || tr->st[y] == 6) {
					j = tr->k[y];

					/* we have a match state */
					if (tr->st[y] == 2) {
						b = dsq[tr->i[y]];
						if (b > K) b = K;
					}

					/* we have a delete state */
					else if (tr->st[y] == 6) {
						b = K;
					}

					idx = IDX(a,b,K+1);
					//fprintf(stdout, "i=%d, j=%d, a=%d, b=%d, idx=%d\n", i, j, a, b, idx);
					esc += hpm->e[i][j][idx];
				}
			}
		}
	}

	*ret_hsc = hsc;
	*ret_esc = esc;
	return eslOK;
}

int hpmscore_ScoreInsertions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_isc) {
	int     z;          /* index for trace elements 	     */
	int     i;          /* match state index			        */
	int     a;          /* residue index            	     */
	float   isc = 0.0;  /* log prob of insert emissions     */

	int K = hpm->abc->K;

	for (z = 0; z < tr->N; z++) {
		/* we have a match state */
		if (tr->st[z] == 4 || (tr->st[z] == 8 && tr->i[z] > 0 ) || ((tr->st[z] == 13 && tr->i[z] > 0 ))) {
			/* previous match state */
			i = tr->k[z];
			a = dsq[tr->i[z]];
			/* treat all degenerate residues as first letter in abc  */
			/* insert emissions cancel in log odds score anyways     */
			/* this is possibly a half baked idea, 11/13/2018        */
			if (a > K) a = 0;
		  	isc += log(hpm->ins[i][a]);
		}
	}

	*ret_isc = isc;

	return eslOK;
}

int hpmscore_ScoreTransitions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, int L, float *ret_tsc) {
	int   z;              /* index for trace elements 	             */
	int   st;
	int   stprev  = -1;   /* state id index           	             */
	int   i;              /* match state index			             */
	int   iprev;          /* match state index			             */
	float tsc     = 0.0;  /* log of transition probs                */
	float tsc_NN  = 0.0;  /* log of N->N and C->C transition prob   */
	float tsc_NB  = 0.0;  /* log of N->B and C->T transition prob   */

	/* calculate N->N and C->C transition probability, in log space */
	tsc_NN = log(L) - log(L+2);
	/* calculate N->B and C->T transition probability, in log space */
	tsc_NB = log(2) - log(L+2);

	/* loop over trace positions for this seq */
	for (z = 0; z < tr->N; z++) {
		/* state type */
		st = tr->st[z];
		/* node in model */
		i = tr->k[z];

		/* handle transitions into N state */
		if (st == 8) {
			/*ignore S->N transitions, they have prob 1 */

			/* handle N->N transitions */
			if (stprev == 8)                     tsc = tsc + tsc_NN;
		}

		/* handle transitions into B state */
		else if (st == 9) {
			/* handle N->B transitons */
			if (stprev == 8)                     tsc += tsc_NB;
			/* no other transitions into B state possible */
		}

		/* handle transitions into hybrid match-del states */
		else if (st == 2 || st == 6) {
			/* handle M->M transitions */
			if      (stprev == 2 || stprev == 6) tsc += log(hpm->t[iprev][HPM_MM]);
			/* handle I->M transitions */
			else if (stprev == 4)                tsc += log(hpm->t[iprev][HPM_IM]);
			/* ignore B->M transitions, they have prob 1 */
		}

		/* handle transitions into I states */
		else if (st == 4) {
			/* handle M->I transitions */
			if      (stprev == 2 || stprev == 6) tsc += log(hpm->t[iprev][HPM_MI]);
			/* handle I->I transitions */
			else if (stprev == 4)                tsc += log(hpm->t[iprev][HPM_II]);
		}

		/* handle transitions into C state */
		else if (st == 13) {
			if (stprev == 13)                    tsc += tsc_NN;
			/* ignore E->C transitions, they have prob 1 */
		}

		/*handle transitions into T state */
		else if (st == 15) {
			if (stprev == 13)                    tsc += tsc_NB;
			/* no other possible transitions to T state */
		}

		stprev = st;
		iprev = i;
	}

	*ret_tsc = tsc;

	return eslOK;
}
