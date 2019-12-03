#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"

#include "p7_config.h"
#include "hmmer.h"

#include "hpm.h"
#include "hpmfile.h"
#include "potts.h"

static int trim_msa_rf(ESL_MSA *msa, ESL_ALPHABET *abc, char *errbuf);

static ESL_OPTIONS options[] = {
   /* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
   { "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },

   /* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",  eslARG_NONE,  FALSE,  NULL, NULL, NULL,NULL,"--dna,--rna",   "<pottsfile> and <hmmfile> are for protein sequences",                    2 },
   { "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<pottsfile> and <hmmfile> are for dna sequences",                        2 },
   { "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<pottsfile> and <hmmfile> are for rna sequences",                        2 },
   { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <hpmfile> <hmmfile> <msafile>";

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

int main(int argc, char *argv[]) {
	ESL_GETOPTS  *go          = NULL;                /* application configuration       */
	ESL_ALPHABET *abc         = NULL;                /* biological alphabet             */
	char         *alifile     = NULL;	             /* alignment file name             */
	ESL_MSAFILE  *afp         = NULL;	             /* open msa file                   */
	ESL_MSA      *msa         = NULL;	             /* one multiple sequence alignment */
	int           fmt         = eslMSAFILE_UNKNOWN;  /* format code for alifile         */
	int64_t       alen;		                         /* alignment length                */
	int           nseq;                              /* number of sequences in the msa  */
	char         *hmmfile     = NULL;                /* hmm input file path             */
	P7_HMMFILE   *hfp;                               /* input hmm file                  */
   P7_HMM       *hmm;                               /* input hmm                       */
	char         *hpm_outfile;                       /* hpm output file path            */
   HPM          *hpm;                               /* hidden potts model              */
	POTTS        *potts;                             /* potts model                     */
	int           status;                            /* easel return code               */
	char          errbuf[eslERRBUFSIZE];

	/* parse command line */
	go = esl_getopts_Create(options);
	if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
	if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);

	if (esl_opt_ArgNumber(go)                  != 3)     cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);


	/* read arguments */
	hpm_outfile    = esl_opt_GetArg(go, 1);
	hmmfile        = esl_opt_GetArg(go, 2);
	alifile        = esl_opt_GetArg(go, 3);


	/* if user has defined an alphabet we define it here */
	if        (esl_opt_GetBoolean(go, "--amino"))       abc = esl_alphabet_Create(eslAMINO);
	else if   (esl_opt_GetBoolean(go, "--rna"))         abc = esl_alphabet_Create(eslRNA);
	else if   (esl_opt_GetBoolean(go, "--dna"))         abc = esl_alphabet_Create(eslDNA);

	/* Open the .hmm file */
	status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
	if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
	else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
	else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);

	/* read first hmm  from hmm file */
	if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

	/* open the alignment file */
	if ( (status = esl_msafile_Open(&abc, alifile, NULL, fmt, NULL, &afp)) != eslOK) esl_msafile_OpenFailure(afp, status);

	/* read first msa from alignment file */
	if ( (status = esl_msafile_Read(afp, &msa)) != eslOK) esl_msafile_ReadFailure(afp, status);

	/* calculate Henikoff PB weights from MSA */
	status = esl_msaweight_PB(msa);

	/* trim msa to consensus */
	trim_msa_rf(msa, abc, errbuf);

	fprintf(stdout, "hello world!\n");

	/* clean up and return */
	esl_msafile_Close(afp);
	esl_msa_Destroy(msa);
	p7_hmmfile_Close(hfp);
	p7_hmm_Destroy(hmm);
	esl_alphabet_Destroy(abc);
	esl_getopts_Destroy(go);

	return 0;
}

static int trim_msa_rf(ESL_MSA *msa, ESL_ALPHABET *abc, char *errbuf) {
	int  status;
	int  rflen = 0;
	int  apos  = 0;
	int *useme = NULL;

	ESL_ALLOC(useme, sizeof(int) * msa->alen);

	/* contract check */
	if(msa->rf == NULL) ESL_FAIL(eslEINVAL, errbuf, "Error, trying to map RF positions to alignment positions, but msa->rf is NULL.");

	/* loop over alignment columns, see which are annotated match columns */
	for(apos = 0; apos < msa->alen; apos++) {
		if( (! esl_abc_CIsGap(abc, msa->rf[apos]))         &&
			 (! esl_abc_CIsMissing(abc, msa->rf[apos]))     &&
			 (! esl_abc_CIsNonresidue(abc, msa->rf[apos]))) {

				rflen++;
				useme[apos] = TRUE;
		}
		else{
			useme[apos] = FALSE;
		}
	}


	/* apply rf mask to MSA */
	if ((status = esl_msa_ColumnSubset         (msa, errbuf, useme)) != eslOK) esl_fatal(errbuf);

	/* clean up and return */
	free(useme);
	return eslOK;

	ERROR:
		return status;
}



