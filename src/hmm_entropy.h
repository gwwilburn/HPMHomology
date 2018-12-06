/* hmm_entropy.h */

#ifndef HMM_ENTROPY_INCLUDED
#define HMM_ENTROPY_INCLUDED

#include <stdio.h>

#include "hmmer.h"
#include "p7_config.h"
#include "easel.h"


/* hmm_entropy.c */
extern float     hmm_entropy_Calculate(ESL_SQ *sq, P7_PROFILE *gm, float *ret_H);

#endif
