#include "hpm.h"

#include "easel.h"
#include "esl_alphabet.h"


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

	ESL_ALLOC(hpm, sizeof(HPM));

	hpm->M   = M;
	hpm->abc = abc;

	/* allocate memory for transition and insert emission params */

	/* level 1 */
	ESL_ALLOC(hpm->t,    (M+1) * sizeof(float *));
	ESL_ALLOC(hpm->ins,  (M+1) * sizeof(float *));
	hpm->t[0]   = NULL;
	hpm->ins[0] = NULL;

	/* level 2 */
	ESL_ALLOC(hpm->t[0],   (nTransition*(M+1)) * sizeof(float));
	ESL_ALLOC(hpm->ins[0], (abc->K*(M+1))      * sizeof(float));
	for (k = 0; k <= M; k++) {
		hpm->t[k]   = hpm->t[0] 	+ (k * nTransition);
		hpm->ins[k] = hpm->ins[0]  + (k * hpm->abc->K);
	}

	/* allocate memory for potts params */
	ESL_ALLOC(hpm->h,  sizeof(double  *) * M);
	ESL_ALLOC(hpm->e,  sizeof(double **) * M);

	for (i=0; i<M; i++){

		ESL_ALLOC(hpm->h[i], sizeof(double ) * Kg);
		ESL_ALLOC(hpm->e[i], sizeof(double  *) * M);
		for (j=0; j<M; j++) {
			ESL_ALLOC(hpm->e[i][j], sizeof(double  ) * Kg2);
	/* set up insert state emission probs to be swissprot freqs for proteins */
	/* right now this is the only way I'm going to deal with insert emissions */
	/* what to do for RNA??? */
	double *bg[] = NULL;


		}
	}

	return hpm;

  ERROR:
	return NULL;
}

int IDX(int i, int j, int K) {

	return (i*K) + j;

}


