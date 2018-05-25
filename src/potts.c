/* potts.c */
#include <potts.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "hpm.h"

POTTS *
potts_Create(int L, ESL_ALPHABET *abc)

{
	POTTS 	*potts     = NULL;      /* potts model object to return */
	int 		 status;

	int       Kg        = abc->K+1;  /* alphabet size w/ gap char     */
	int       Kg2       = Kg*Kg;

	int       i;
	int       j;


	ESL_ALLOC(potts, sizeof(POTTS));

	potts->L   = L;
	potts->abc = abc;

	/* allocate memory for potts params */
	ESL_ALLOC(potts->h,  sizeof(double  *) * (L));
	ESL_ALLOC(potts->e,  sizeof(double **) * (L));

	for (i=0; i<L; i++){
		ESL_ALLOC(potts->h[i], sizeof(double ) * Kg);
		ESL_ALLOC(potts->e[i], sizeof(double  *) * (L));

		for (j=0; j<L; j++) {
			ESL_ALLOC(potts->e[i][j], sizeof(double  ) * Kg2);
		}
	}

	return potts;

	ERROR:
		return NULL;

}
