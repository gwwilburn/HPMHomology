/* hpm.h */

#ifndef HPM_INCLUDED
#define HPM_INCLUDED

#include <stdio.h>

#include "hmmer.h"
#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "potts.h"

enum hpm_transitions {
	HPM_MM = 0,
	HPM_MI = 1,
	HPM_IM = 2,
	HPM_II = 3,
};

#define HPM_NTRANSITIONS 4



typedef struct hpm_s{

	int		  		 M;				/* number of nodes in model 					*/
	double 		  **t;				/* transition probs 								*/
	double	 		**ins;			/* insert emission probs					*/
	float    	  **h;				/* site-specific h_i params					*/
	float   		 ***e;				/* coupling e_ij params							*/
	int				 nTransition;  /* number of transitions, should be 4 		*/

	ESL_ALPHABET *abc;            /* alphabet											*/

} HPM;

typedef struct hpmscoreset_s{
	char   **sqname;   /* sequence names [0..nseq-1][], \0-terminated      */
	int 		nseq;     /* number of sequences                              */
	float   *E_potts;  /* Potts pseudo-energies [0...nseq-1]               */
	float   *p_ins;    /* product of insertion probabilities [0...nseq-1]  */
	float   *p_trans;  /* product of transition probabilities [0...nseq-1] */

} HPM_SCORESET;

/* hpm.c	*/
extern HPM	         *hpm_Create(int M, ESL_ALPHABET *abc);
extern HPM 				*hpm_Create_hmm_potts(P7_HMM *hmm, POTTS *potts, ESL_ALPHABET *abc);
extern HPM_SCORESET	*hpm_scoreset_Create(int nseq);
extern int				 hpm_scoreset_Write(FILE *fp, HPM_SCORESET *hpm_SS);
extern int            IDX(int i, int j, int K);


#endif
