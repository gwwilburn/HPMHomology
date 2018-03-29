/* hpm.h */

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"

#ifndef HPM_INCLUDED
#define HPM_INCLUDED

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










#endif
