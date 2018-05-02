#include "hpm.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"


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

int IDX(int i, int j, int K) {

	return (i*K) + j;

}


