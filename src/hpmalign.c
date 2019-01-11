/* hpmalign.c */
/* align sequences to an hpm via importance sampling w/ an hmm */

#include "hmmer.h"
#include "easel.h"

#include "p7_config.h"

#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hpm.h"
#include "hpmfile.h"
#include "hmm_entropy.h"
#include "hpm_scoreops.h"

/* declaration of internal functions */



static ESL_OPTIONS options[] = {
	/* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
	{ "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL,           "help; show brief info on version and usage",                 0 },
	{ "-s",         eslARG_INT,     "0", NULL, NULL, NULL, NULL, NULL,           "set random number seed to <n>",                              0 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
	{ "--amino",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<seqfile> contains protein sequences",                       1 },
	{ "--rna",      eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--amino",  "<seqfile> contains RNA sequences",                           1 },
	{ "--dna",      eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--rna,--amino",  "<seqfile> contains DNA sequences",                           1 },

	/* debugging tools */
	{ "--v",        eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,NULL,              "Verbose mode: print info on intermediate alignment steps",   2 },
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hpmfile> <hmmfile> <seqfile> <msa_outfile>";
static char banner[] = "align sequences to a profile HPM";

int main(int argc, char **argv){

	ESL_GETOPTS      *go            =  p7_CreateDefaultApp(options, 4, argc, argv, banner, usage);
	ESL_RANDOMNESS   *rng           = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	int               v             = 0;
	char             *hpmfile       = esl_opt_GetArg(go, 1);
	char             *hmmfile       = esl_opt_GetArg(go, 2);
	char             *seqfile       = esl_opt_GetArg(go, 3);
	char             *msafile       = esl_opt_GetArg(go, 4);
	ESL_ALPHABET     *abc           = NULL;
	HPM              *hpm           = NULL;
	P7_HMMFILE       *hmmfp         = NULL;
	P7_HMM           *hmm           = NULL;
	ESL_SQ          **sq            = NULL;                      /* array of sequences */
	ESL_SQFILE       *sqfp          = NULL;
	FILE             *afp           = NULL;                      /* output alignment file */
	ESL_MSA          *msa           = NULL;                      /* resulting msa */
	P7_TRACE        **tr            = NULL;
	int               outfmt        = eslMSAFILE_STOCKHOLM;
	int               format        = eslSQFILE_UNKNOWN;
	int               totseq        = 0;
	int               nseq          = 0;
	int status;
	char              errbuf[eslERRBUFSIZE];

	/* check for verbose mode */
	if (esl_opt_GetBoolean(go, "--v")) v=1;

	/* if user has defined an alphabet we define it here */
	if        (esl_opt_GetBoolean(go, "--amino"))       abc = esl_alphabet_Create(eslAMINO);
	else if   (esl_opt_GetBoolean(go, "--rna"))         abc = esl_alphabet_Create(eslRNA);
	else if   (esl_opt_GetBoolean(go, "--dna"))         abc = esl_alphabet_Create(eslDNA);

	/* read in hmm, autodetect alphabet if not manually set by user  */
	if (p7_hmmfile_OpenE(hmmfile, NULL, &hmmfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hmmfp, &abc, &hmm)        != eslOK) p7_Fail("Failed to read HMM");
	p7_hmmfile_Close(hmmfp);

	/* open hpm */
	hpm = hpmfile_Read(hpmfile, abc, errbuf);

	/* open sequence file */
	status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
	if      (status == eslENOTFOUND) p7_Fail("No such file.");
	else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
	else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
	else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

	ESL_REALLOC(sq, sizeof(ESL_SQ *) * (totseq + 1));
	sq[totseq] = esl_sq_CreateDigital(abc);
	while ((status = esl_sqio_Read(sqfp, sq[totseq+nseq])) == eslOK) {
		nseq++;
		ESL_REALLOC(sq, sizeof(ESL_SQ *) * (totseq+nseq+1));
		sq[totseq+nseq] = esl_sq_CreateDigital(abc);
	}

	/* clean up and return */
	esl_sqfile_Close(sqfp);
	free(sq);
	p7_hmm_Destroy(hmm);
	free(hpm);
	esl_alphabet_Destroy(abc);
	esl_getopts_Destroy(go);
	esl_msa_Destroy(msa);
	free(tr);

	return 0;

	ERROR:
		return status;
}
