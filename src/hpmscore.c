#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_getopts.h"

#include "p7_config.h"
#include "hmmer.h"

#include "hpm.h"
#include "hpmfile.h"

struct cfg_s { /* shared configuration in master and workers */
	ESL_ALPHABET	 *abc;
	ESL_MSAFILE   	 *afp;     	/* input ali file						*/
	ESL_MSA			 *msa;     	/* input msa							*/
	P7_HMMFILE		 *hfp;		/* input hmm file						*/
	P7_HMM			 *hmm;		/* input hmm							*/
	HPM				 *hpm;		/* hidden potts model				*/

};


static ESL_OPTIONS options[] = {
	/* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
	{ "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
	{ "--amino",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<msafile> contains protein alignments",                   2 },
	{ "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<msafile> contains DNA alignments",                       2 },
	{ "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<msafile> contains RNA alignments",                       2 },
	{ 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <msafile>";

static void
cmdline_failure(char *argv0, char *format, ...)
{
	va_list argp;
	va_start(argp, format);
	vfprintf(stderr, format, argp);
	va_end(argp);
	esl_usage(stdout, argv0, usage);
	printf("\nTo see more help on available options, do %s -h\n\n", argv0);
	exit(1);
}


int main(int argc, char *argv[])
{
	ESL_GETOPTS     *go;											/* cmd line configuration 			*/
	struct cfg_s    cfg;  										/* application configuration 		*/
	char            *alifile 					= NULL;		/* alignment file path 				*/
	char				 *hmmfile 					= NULL;		/* hmm file path						*/
	int				 status;										/* easel return code 				*/
	char            errbuf[eslERRBUFSIZE];

	/* parse command line */
	go = esl_getopts_Create(options);
	if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
	if (esl_opt_VerifyConfig(go) 					 != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
	if (esl_opt_ArgNumber(go)						 != 2) 	  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

	alifile  = esl_opt_GetArg(go, 1);
	hmmfile  = esl_opt_GetArg(go, 2);

	/* Set up the configuration structure shared amongst functions here */
	cfg.abc = NULL; /* until we read the msa file below */

	/* Open the MSA file, digital mode; determine alphabet */
	if 		(esl_opt_GetBoolean(go, "--amino"))		cfg.abc = esl_alphabet_Create(eslAMINO);
	else if  (esl_opt_GetBoolean(go, "--dna"))		cfg.abc = esl_alphabet_Create(eslDNA);
	else if  (esl_opt_GetBoolean(go, "--rna"))		cfg.abc = esl_alphabet_Create(eslRNA);

	status = esl_msafile_Open(&(cfg.abc), alifile, NULL, eslMSAFILE_UNKNOWN, NULL, &cfg.afp);
	if (status != eslOK) esl_msafile_OpenFailure(cfg.afp, status);

	esl_alphabet_SetEquiv(cfg.abc, '.', '-'); /* allow '.' as gap char too */

	/* Read in the MSA */
	status = esl_msafile_Read(cfg.afp, &cfg.msa);
	if (status != eslOK) esl_msafile_ReadFailure(cfg.afp, status);

	/* Open the .hmm file */
	status = p7_hmmfile_OpenE(hmmfile, NULL, &cfg.hfp, errbuf);
	/* error handling below copied from hmmsearch */
	if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
	else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
	else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);


	/* read the .hmm file */
	status = p7_hmmfile_Read(cfg.hfp, &cfg.abc, &cfg.hmm);

	/* temporary 3-mer hpm for testing, to be deleted */
	int M = 3;
	cfg.hpm = hpm_Create(M, cfg.abc);

	int a;
	int b;
	int i;
	int j;


	/* fill in e_ij's */
	for (a = 0; a < cfg.abc->K; a++) {
		for (i = 0; i < M; i++) {
			cfg.hpm->h[i][a] = (0.1 * i)  + (0.0001 * a);
			for (b = 0; b < cfg.abc->K; b++) {
				for (j = 0; j < i; j++) {
					cfg.hpm->e[i][j][IDX(a,b,cfg.hpm->abc->K+1)] = (0.1 * i) + (0.01*j) + (0.0001 * a) + (0.000001 * b);
					cfg.hpm->e[j][i][IDX(b,a,cfg.hpm->abc->K+1)] = (0.1 * i) + (0.01*j) + (0.0001 * a) + (0.000001 * b);
				}
			}

		}

	}



	int k;
	for (k = 0; k < M+1; k++) {
		cfg.hpm->t[k][HPM_MM] = 1.0 - (0.1*(k+1));
		cfg.hpm->t[k][HPM_MI] = (0.1*(k+1));
		if (k < M+1) {
			cfg.hpm->t[k][HPM_IM] = 1.0 - (0.1*(k+1));
			cfg.hpm->t[k][HPM_II] = (0.1*(k+1));
		}

	}


	char *hpm_outfile = "test.hpm";
	FILE *hpm_outfp 	= 	NULL;
	if ((hpm_outfp = fopen(hpm_outfile, "w")) == NULL) esl_fatal("Failed to open output hpm file %s for writing", hpm_outfile);
	hpmfile_Write(hpm_outfp, cfg.hpm);
	fprintf(stdout, "hello world!\n");

	/* clean up */
	esl_alphabet_Destroy(cfg.abc);
	esl_msa_Destroy(cfg.msa);
	esl_msafile_Close(cfg.afp);
	p7_hmm_Destroy(cfg.hmm);
	p7_hmmfile_Close(cfg.hfp);
	esl_getopts_Destroy(go);
	fclose(hpm_outfp);

	return 0;
}

