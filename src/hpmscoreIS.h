/* hpmscoreIS.h */

#ifndef HPMSCOREIS_INCLUDED
#define HPMSCOREIS_INCLUDED

#include "p7_config.h"
#include "easel.h"
#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "hpm.h"
#include "hpmfile.h"

#include "hmmer.h"

/* hpmscoreIS.c	*/
extern int hpmscore_CalculateHamiltonian(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_hsc, float *ret_esc);
extern int hpmscore_ScoreInsertions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_isc);
extern int hpmscore_ScoreTransitions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, int L, float *ret_tsc);

#endif
