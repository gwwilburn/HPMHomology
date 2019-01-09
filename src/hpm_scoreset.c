#include "hpm_scoreset.h"

#include "p7_config.h"
#include "hmmer.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"


/* Function hpm_scoreset_Create()
 *
 * Synopsis: Allocate score set object for original (viterbi path) hpm scoring method
 *
 * Args:     nseq: number of sequences to be scored
 *
 * Returns:  pointer to the new hpm scoreset object
 *
 * Throws:   <NULL> if any allocation or initialization fails
 */

HPM_SCORESET *
hpm_scoreset_Create(int nseq)
{
	HPM_SCORESET *hpm_ss = NULL;
	int           n;
	int           status;

	ESL_ALLOC(hpm_ss, sizeof(HPM_SCORESET));

	/* set number of sequences */
	hpm_ss->nseq = nseq;

	/* allocate memory for sequence name array */
	ESL_ALLOC(hpm_ss->sqname, sizeof(char *) * nseq);

	/* allocate memory for probabilities/scores */
	ESL_ALLOC(hpm_ss->E_hi,         sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->E_eij,        sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lp_ins,       sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lp_trans,     sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lpnull_match, sizeof(float) * nseq);
	ESL_ALLOC(hpm_ss->lpnull_trans, sizeof(float) * nseq);

	for (n = 0; n < nseq; n++) {
		hpm_ss->E_hi[n]         = 0.0;
		hpm_ss->E_eij[n]        = 0.0;
		hpm_ss->lp_ins[n]       = 0.0;
		hpm_ss->lp_trans[n]     = 0.0;
		hpm_ss->lpnull_match[n] = 0.0;
		hpm_ss->lpnull_trans[n] = 0.0;
	}

	return hpm_ss;

	ERROR:
		return NULL;
}



/* Function hpm_scoreset_Create()
 *
 * Synopsis: Allocate score set object for v2 (importance sampling, or "is")
 *           hpm scoring method
 *
 * Args:     nseq: number of sequences to be scored
 *
 * Returns:  pointer to the new hpm_is scoreset object
 *
 * Throws:   <NULL> if any allocation or initialization fails
 */

HPM_IS_SCORESET *
hpm_is_scoreset_Create(int nseq)
{
	HPM_IS_SCORESET *hpm_is_ss = NULL;
	int              n;
	int              status;

	ESL_ALLOC(hpm_is_ss, sizeof(HPM_IS_SCORESET));

	/* set number of sequences */
	hpm_is_ss->nseq = nseq;

	/* allocate memory for sequence name array */
	ESL_ALLOC(hpm_is_ss->sqname, sizeof(char *) * nseq);

	/* allocate memory for probabilities/scores */
	ESL_ALLOC(hpm_is_ss->R,           sizeof(int)   * nseq);
	ESL_ALLOC(hpm_is_ss->fwd,         sizeof(float) * nseq);
	ESL_ALLOC(hpm_is_ss->is_ld,       sizeof(float) * nseq);


	for (n = 0; n < nseq; n++) {
		hpm_is_ss->R[n]         = 0;
		hpm_is_ss->fwd[n]       = 0.0;
		hpm_is_ss->is_ld[n]     = 0.0;
	}

	return hpm_is_ss;

	ERROR:
		return NULL;
}



/* Function hpm_scoreset_Write()
 *
 * Synopsis: Write original (Viterbi-path) hpm scores to an output .csv file
 *
 * Args:     fp:     output file
 *           hpm_ss: original hpm scoreset object with scores
 *
 * Returns:  <eslOK> on success
 */

int
hpm_scoreset_Write(FILE *fp, HPM_SCORESET *hpm_ss){
	int n;

	/* write csv header line */

	fprintf(fp,"id,E_hi,E_eij,hamiltonian,lp_ins,lp_trans,lpnull_match,lpnull_trans,logodds_unnormalized\n");

	/* loop over sequences, print id and scores */
	for (n = 0; n < hpm_ss->nseq; n++) {
		fprintf(fp, "%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
				  hpm_ss->sqname[n],
				  hpm_ss->E_hi[n],
				  hpm_ss->E_eij[n],
				  hpm_ss->E_hi[n] + hpm_ss->E_eij[n],
				  hpm_ss->lp_ins[n],
				  hpm_ss->lp_trans[n],
				  hpm_ss->lpnull_match[n],
				  hpm_ss->lpnull_trans[n],
				  hpm_ss->E_hi[n] + hpm_ss->E_eij[n] + hpm_ss->lp_trans[n] - hpm_ss->lpnull_match[n] - hpm_ss->lpnull_trans[n]);

	}

	return eslOK;
}


/* Function hpm_is_scoreset_Write()
 *
 * Synopsis: Write v2 (importance sampling) hpm scores to an output .csv file
 *
 * Args:     fp:        output file
 *           hpm_is_ss: v2 hpm_is scoreset object with scores
 *
 * Returns:  <eslOK> on success
 */

int
hpm_is_scoreset_Write(FILE *fp, HPM_IS_SCORESET *hpm_is_ss){
	int n;

	/* write csv header line */

	fprintf(fp,"id,N_sample,fwd_logodds,HPM_IS_logodds\n");

	/* loop over sequences, print id and scores */
	for (n = 0; n < hpm_is_ss->nseq; n++) {
		fprintf(fp, "%s,%d,%.4f,%.4f\n",
				  hpm_is_ss->sqname[n],
				  hpm_is_ss->R[n],
				  hpm_is_ss->fwd[n],
				  hpm_is_ss->is_ld[n]);

	}

	return eslOK;
}
