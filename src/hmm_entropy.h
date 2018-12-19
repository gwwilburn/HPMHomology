/* hmm_entropy.h */

#ifndef HMM_ENTROPY_INCLUDED
#define HMM_ENTROPY_INCLUDED

#include <stdio.h>


#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_config.h"
#include "dp_reference/p7_refmx.h"
#include "base/p7_profile.h"

/* hmm_entropy.c */
extern float     hmm_entropy_Calculate(P7_PROFILE *gm, P7_REFMX *fwd, float *ret_H, int verbose);
#endif
