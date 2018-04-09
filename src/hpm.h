/* hpm.h */

#ifndef HPM_INCLUDED
#define HPM_INCLUDED

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"

enum hpm_transitions {
	HPM_MM = 0,
	HPM_MI = 1,
	HPM_IM = 2,
	HPM_II = 3,
};

#define HPM_NTRANSITIONS 4



typedef struct hpm_s{

	int		  		 M;				/* number of nodes in model 	*/
	float 		  **t;				/* transition probs 				*/
	float	 		**ins;				/* insert emission probs		*/
	float    	  **h;				/* site-specific h_i params	*/
	float   		 ***e;				/* coupling e_ij params			*/

	ESL_ALPHABET *abc;            /* alphabet							*/

} HPM;

/* hpm.c	*/
extern HPM	*hpm_Create(int M, ESL_ALPHABET *abc);
extern int   IDX(int i, int j, int K);


#endif
