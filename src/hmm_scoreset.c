#include "hmm_scoreset.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"

HMM_SCORESET *
hmm_scoreset_Create(int nseq)
{
	HMM_SCORESET *hmm_ss = NULL;
	int           n;
	int           status;

	ESL_ALLOC(hmm_ss, sizeof(HMM_SCORESET));

	/* set number of sequences */
	hmm_ss->nseq = nseq;

	/* allocate memory for sequence name array */
	ESL_ALLOC(hmm_ss->sqname, sizeof(char *) * nseq);

	/* allocate memory for probabilities/scores */
	ESL_ALLOC(hmm_ss->vsc,       sizeof(float) * nseq);
	ESL_ALLOC(hmm_ss->fsc,       sizeof(float) * nseq);
	ESL_ALLOC(hmm_ss->bsc,       sizeof(float) * nseq);
	ESL_ALLOC(hmm_ss->nullsc,    sizeof(float) * nseq);

	for (n = 0; n < nseq; n++) {
		hmm_ss->vsc[n]     = 0.0;
		hmm_ss->fsc[n]     = 0.0;
		hmm_ss->bsc[n]     = 0.0;
		hmm_ss->nullsc[n]  = 0.0;
	}

	return hmm_ss;

	ERROR:
		return NULL;
}

/* for now, write hmm scores to a csv */
int
hmm_scoreset_Write(FILE *fp, HMM_SCORESET *hmm_ss){
	int n;

	/* write csv header line */

	fprintf(fp,"id,raw_viterbi,raw_forward,raw_backward,logodds_viterbi,logodds_forward,logodds_backward\n");

	/* loop over sequences, print id and scores */
	for (n = 0; n < hmm_ss->nseq; n++) {
		fprintf(fp, "%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
				  hmm_ss->sqname[n],
				  hmm_ss->vsc[n],
				  hmm_ss->fsc[n],
				  hmm_ss->bsc[n],
				  (hmm_ss->vsc[n] - hmm_ss->nullsc[n])/ eslCONST_LOG2,
				  (hmm_ss->fsc[n] - hmm_ss->nullsc[n])/ eslCONST_LOG2,
				  (hmm_ss->bsc[n] - hmm_ss->nullsc[n])/ eslCONST_LOG2);

	}

	return eslOK;
}

