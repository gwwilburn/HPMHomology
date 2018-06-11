#include "hpm.h"

#include "p7_config.h"
#include "hmmer.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"
#include "potts.h"

HPM *
hpm_Create(int M, ESL_ALPHABET *abc)
{
	int    Kg  			 = abc->K+1;
	int 	 Kg2 			 = Kg*Kg;
	int 	 nTransition = 4;

	HPM 	*hpm = NULL;
	int 	 status;
	int 	 i;
	int 	 j;
	int 	 k;
	int 	 l;

	ESL_ALLOC(hpm, sizeof(HPM));

	hpm->M   = M;
	hpm->abc = abc;
	hpm->nTransition = nTransition;

	/* allocate memory for transition and insert emission params */

	/* level 1 */
	ESL_ALLOC(hpm->t,    (M+1) * sizeof(double *));
	ESL_ALLOC(hpm->ins,  (M+1) * sizeof(double *));
	hpm->t[0]   = NULL;
	hpm->ins[0] = NULL;

	/* level 2 */
	ESL_ALLOC(hpm->t[0],   (nTransition*(M+1)) * sizeof(double));
	ESL_ALLOC(hpm->ins[0], ((abc->K)*(M+1))      * sizeof(double));
	for (k = 0; k < M+1; k++) {
		hpm->t[k]   = hpm->t[0] 	+ (k * nTransition);
		hpm->ins[k] = hpm->ins[0]  + (k * hpm->abc->K);
		if (abc->type == eslAMINO) {
			esl_composition_SW50(hpm->ins[k]);
		}
		else {
			for (l=0; l < abc->K; l++){
				hpm->ins[k][l] = 0.25;
			}
		}

	}


	/* allocate memory for potts params */
	ESL_ALLOC(hpm->h,  sizeof(double  *) * (M+1));
	ESL_ALLOC(hpm->e,  sizeof(double **) * (M+1));

	for (i=0; i<M+1; i++){

		ESL_ALLOC(hpm->h[i], sizeof(double ) * Kg);
		ESL_ALLOC(hpm->e[i], sizeof(double  *) * (M+1));
		for (j=0; j<M+1; j++) {
			ESL_ALLOC(hpm->e[i][j], sizeof(double  ) * Kg2);


		}
	}

	return hpm;

  ERROR:
	return NULL;
}

HPM *
hpm_Create_hmm_potts(P7_HMM *hmm, POTTS *potts, ESL_ALPHABET *abc) {

	HPM   *hpm      = NULL;   /* hpm object to return */
	int 	 status;
	int    i;                 /* position index 		    */
	int    j;				     /* position index         */
	int    a;				     /* character index        */
	int    b;                 /* character index        */
	int    idx;               /* character combo-index  */
	int    idx_rev;           /* character combo-index  */
	float *mocc    = NULL;    /* hmm match state probs  */

	/* check to make sure number of match states are equal */
	if      (potts->L != hmm->M)     p7_Fail("Length of potts model is %d. HMM has %d match states. These must match!\n", potts->L, hmm->M);
	/* check to make sure alphabets are equal */
	if      (potts->abc != hmm->abc) p7_Fail("Potts alphabet and hmm alphabet do not match!\n", potts->L, hmm->M);

	ESL_ALLOC(hpm, sizeof(HPM));
	ESL_ALLOC(mocc, sizeof(float)*(hmm->M+1));

	/* create hpm with appropriate  number of match states */
	hpm = hpm_Create(hmm->M, potts->abc);

	/* copy insert emission from hmm as-is */
	/* this throws an error for now */
	/* but the swissprot freqs are assigned in hpm_create, will fix later */
	//hpm->ins = hmm->ins;


	/* copy over potts parameters */
	for (i=0; i<hpm->M; i++) {
		for (a = 0; a < (hpm->abc->K+1); a++){
			/* assign h_i's */
			hpm->h[i+1][a] = potts->h[i][a];

			for (j = i+1; j<hpm->M; j++) {
				for (b = 0; b < (hpm->abc->K+1); b++) {

					/* assign e_ij's */
					idx = IDX(a,b,hpm->abc->K+1);
					idx_rev = IDX(b,a,hpm->abc->K+1);
					hpm->e[i+1][j+1][idx]     = potts->e[i][j][idx];
					hpm->e[j+1][i+1][idx_rev] = potts->e[j][i][idx_rev];
				}
			}
		}
	}

	/* calculate hpm transition probahilities */

	/* indices for hmm transitions */
	int hmm_MM = 0;
	int hmm_MI = 1;
	int hmm_MD = 2;
	int hmm_IM = 3;
	int hmm_II = 4;

	/* obtain match state probabilities from hmm */
	p7_hmm_CalculateOccupancy(hmm, mocc, NULL);

	/* transitions from begin state */
	hpm->t[0][HPM_MM] = hmm->t[0][hmm_MM] + hmm->t[0][hmm_MD];
	hpm->t[0][HPM_MI] = hmm->t[0][hmm_MI];
	hpm->t[0][HPM_IM] = hmm->t[0][hmm_IM];
	hpm->t[0][HPM_II] = hmm->t[0][hmm_II];

	/* intermediate transitions */
	for (i = 1; i < hpm->M; i++) {
		hpm->t[i][HPM_MM] = ((hmm->t[i][hmm_MM] + hmm->t[i][hmm_MD]) * mocc[i]) + (1.0 - mocc[i]);
		hpm->t[i][HPM_MI] = hmm->t[i][hmm_MI] * mocc[i];
		hpm->t[i][HPM_IM] = hmm->t[i][hmm_IM];
		hpm->t[i][HPM_II] = hmm->t[i][hmm_II];
	}

	/* transitions into end state */
	hpm->t[hpm->M][HPM_MM] = (hmm->t[hpm->M][hmm_MM] * mocc[i]) + (1.0 - mocc[i]);
	hpm->t[hpm->M][HPM_MI] = mocc[i]*hmm->t[0][hmm_MI];
	hpm->t[hpm->M][HPM_IM] = hmm->t[0][hmm_IM];
	hpm->t[hpm->M][HPM_II] = hmm->t[0][hmm_II];



	return hpm;

	ERROR:
		return NULL;


}

HPM *
hpm_Create_3mer(ESL_ALPHABET *abc) {
	HPM    *hpm   = NULL;       /* hpm object to return      */
	int     M     = 3;          /* number of nodes           */
	int     i;                  /* site index                */
	int 	  j;                  /* site index                */
	int     a;                  /* alphabet index            */
	int     b;                  /* alphabet index            */
	int     idx;                /* alphabet index for eij's  */
	int     idx2;               /* alphabet index for eij's  */
	int     Kg    = abc->K+1;   /* alphabet size w/ gap      */
	int     status;

	/* allocate memory for hpm */
   ESL_ALLOC(hpm, sizeof(HPM));

   /* create hpm with 3 match states */
   hpm = hpm_Create(M, abc);

	/* set potts parameters */
	for (i = 1; i < (M+1); i++) {
		for(a = 0; a < Kg; a++) {
			/* set all h_i's to 0 */
			hpm->h[i][a] = 0;

			/* set specific set of e_13's to be 1, else 0 */
			for (j=i+1; j < (M+1); j++) {
				for (b=0; b < Kg; b++) {
					idx  = IDX(a,b,Kg);
					idx2 = IDX(b,a,Kg);
					if ( j == (i+1) && (a+b) == abc->K) {
						hpm->e[i][j][idx]  = 1.0;
						hpm->e[j][i][idx2] = 1.0;
					}
					else {
						hpm->e[i][j][idx]  = 0.0;
						hpm->e[j][i][idx2] = 0.0;
					}
				}
			}
		}
	}

	/* set transition parameters */
	for (i = 0; i < M+1; i++) {
		hpm->t[i][HPM_MM] = 1.0 - (0.1)*(i+1);
		hpm->t[i][HPM_MI] = (0.1)*(i+1);
		hpm->t[i][HPM_IM] = 1.0 - (0.15)*(i+1);
		hpm->t[i][HPM_II] = (0.15)*(i+1);
	}

	return hpm;

	ERROR:
		return NULL;
}

HPM_SCORESET *
hpm_scoreset_Create(int nseq)
{
	HPM_SCORESET *hpm_ss = NULL;
	int           n;
	int           status;

	ESL_ALLOC(hpm_ss, sizeof(HPM_SCORESET));

	/* set number of sequences */
	hpm_ss->nseq = nseq;

	/* allocate memory for sequence name array */
	ESL_ALLOC(hpm_ss->sqname, sizeof(char *) * nseq);

	/* allocate memory for probabilities/scores */
	ESL_ALLOC(hpm_ss->E_potts,      sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lp_ins,       sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lp_trans,     sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lpnull_match, sizeof(float) * nseq);

	for (n = 0; n < nseq; n++) {
		hpm_ss->E_potts[n]      = 0.0;
		hpm_ss->lp_ins[n]       = 0.0;
		hpm_ss->lp_trans[n]     = 0.0;
		hpm_ss->lpnull_match[n] = 0.0;
	}

	return hpm_ss;

	ERROR:
		return NULL;
}

/* for now, write hpm scores to a csv */
int
hpm_scoreset_Write(FILE *fp, HPM_SCORESET *hpm_ss){
	int n;

	/* write csv header line */

	fprintf(fp,"id,hamiltonian,lp_ins,lp_trans,lpnull_match\n");

	/* loop over sequences, print id and scores */
	for (n = 0; n < hpm_ss->nseq; n++) {
		fprintf(fp, "%s,%.4f,%.4f,%.4f,%.4f\n",
				  hpm_ss->sqname[n],
				  hpm_ss->E_potts[n],
				  hpm_ss->lp_ins[n],
				  hpm_ss->lp_trans[n],
				  hpm_ss->lpnull_match[n]);

	}

	return eslOK;
}

int IDX(int i, int j, int K) {

	return (i*K) + j;

}


