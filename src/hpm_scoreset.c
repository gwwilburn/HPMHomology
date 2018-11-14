#include "hpm_scoreset.h"

#include "p7_config.h"
#include "hmmer.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"

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
	ESL_ALLOC(hpm_ss->E_eij,       sizeof(float) * nseq);
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

/* for now, write hpm scores to a csv */
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
