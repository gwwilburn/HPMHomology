/* hpm_scoreops.h */

#ifndef HPM_SCOREOPS_INCLUDED
#define HPM_SCOREOPS_INCLUDED

#include "easel.h"
#include "hmmer.h"
#include "p7_config.h"
#include "esl_sq.h"
#include "hpm.h"

/* hpm_scoreops.c */
extern int hpm_scoreops_CalculateHamiltonian(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_hsc, float *ret_esc);
extern int hpm_scoreops_ScoreInsertEmissions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_isc);
extern int hpm_scoreops_ScoreNullEmissions(HPM *hpm, ESL_SQ *sq, float *ret_nesc);
extern int hpm_scoreops_ScoreNullMatchEmissions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, float *ret_nmesc);
extern int hpm_scoreops_ScoreTransitions(HPM *hpm, P7_TRACE *tr, ESL_DSQ *dsq, int L, float *ret_tsc);
extern int hpm_scoreops_ScoreNullTransitions(P7_BG  *bg, ESL_SQ *sq, float *ret_ntsc);
#endif
