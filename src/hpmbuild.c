#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "p7_config.h"
#include "hmmer.h"

#include "hpm.h"
#include "hpmfile.h"
#include "potts.h"
#include "pottsfile.h"

struct cfg_s { /* shared configuration in master and workers */
   ESL_ALPHABET    *abc;
	P7_HMMFILE      *hfp;      /* input hmm file                */
   P7_HMM          *hmm;      /* input hmm                     */
   HPM             *hpm;      /* hidden potts model            */
	POTTS           *potts;    /* potts model                   */

};

static ESL_OPTIONS options[] = {
   /* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
   { "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },

   /* Options forcing which alphabet we're working in (normally autodetected) */
   { "--amino",  eslARG_NONE,  FALSE,  NULL, NULL, NULL,NULL,"--dna,--rna",   "<pottsfile> and <hmmfile> are for protein sequences",                    2 },
   { "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<pottsfile> and <hmmfile> are for dna sequences",                        2 },
   { "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<pottsfile> and <hmmfile> are for rna sequences",                        2 },
	/* Option to create 3-mer dummy protein hpm */
	{ "--3mer",  eslARG_NONE,  FALSE,  NULL, NULL, NULL,NULL,NULL,             "create dummy 3-mer hpm, usage = [-options] <hpmfile>",                   3 },
   { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <hpmfile> <hmmfile> <pottsfile>";

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

	ESL_GETOPTS    *go;                                            /* cmd line configuration        */
   struct cfg_s    cfg;                                           /* application configuration     */
   char           *hpm_outfile            = NULL;                 /* hpm output file path          */
   char           *hmmfile                = NULL;                 /* hmm file path                 */
   char           *pottsfile              = NULL;                 /* hmm file path                 */
 	FILE           *hpm_outfp              = NULL;						/* hpm output file object        */
 	int             status;                                        /* easel return code             */
   char            errbuf[eslERRBUFSIZE];

	/* parse command line */
	go = esl_getopts_Create(options);
   if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
   if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);

	/* check and see if we are making a dummy 3-mer model */
	if       (esl_opt_GetBoolean(go, "--3mer")) {
		/* we only want 1 argument */
		if (esl_opt_ArgNumber(go)   != 1)     cmdline_failure(argv[0], "Incorrect number of command line arguments w/ dummy 3-mer hpm.\n", go->errbuf);
		hpm_outfile    = esl_opt_GetArg(go, 1);

		/* set alphabet to amino, forget what user says */
		cfg.abc = esl_alphabet_Create(eslAMINO);

		/* make 3-mer hpm*/
		cfg.hpm = hpm_Create_3mer(cfg.abc);

		/* write to output file */
		/* open hpm outfile for writing, write hpm to it */
   	if ((hpm_outfp = fopen(hpm_outfile, "w")) == NULL) esl_fatal("Failed to open output hpm file %s for writing", hpm_outfile);
   	hpmfile_Write(hpm_outfp, cfg.hpm);

		return 0;
	}


	if (esl_opt_ArgNumber(go)                  != 3)     cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

	hpm_outfile    = esl_opt_GetArg(go, 1);
	hmmfile        = esl_opt_GetArg(go, 2);
	pottsfile      = esl_opt_GetArg(go, 3);

   /* Open the .hmm file */
   status = p7_hmmfile_OpenE(hmmfile, NULL, &cfg.hfp, errbuf);

	/* Set up the configuration structure shared amongst functions here */
	/* Determine alphabet: do what user says, else default to amino */
   if       (esl_opt_GetBoolean(go, "--amino"))    cfg.abc = esl_alphabet_Create(eslAMINO);
   else if  (esl_opt_GetBoolean(go, "--dna"))      cfg.abc = esl_alphabet_Create(eslDNA);
   else if  (esl_opt_GetBoolean(go, "--rna"))      cfg.abc = esl_alphabet_Create(eslRNA);
	else                                            cfg.abc = esl_alphabet_Create(eslAMINO);

   /* Open the .hmm file */
   status = p7_hmmfile_OpenE(hmmfile, NULL, &cfg.hfp, errbuf);
   /* error handling below copied from hmmsearch */
   if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
   else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
   else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);
	else
	/* read the .hmm file */
   status = p7_hmmfile_Read(cfg.hfp, &cfg.abc, &cfg.hmm);

	/* read the potts file */
	cfg.potts = pottsfile_Read(pottsfile, cfg.abc, errbuf);

	/* combine hmm and potts object to create hpm object */
	cfg.hpm = hpm_Create_hmm_potts(cfg.hmm, cfg.potts, cfg.abc);

	/* open hpm outfile for writing, write hpm to it */
	if ((hpm_outfp = fopen(hpm_outfile, "w")) == NULL) esl_fatal("Failed to open output hpm file %s for writing", hpm_outfile);
	hpmfile_Write(hpm_outfp, cfg.hpm);

	/* close hpm outfile */
	fclose(hpm_outfp);
	return 0;

	//ERROR:
   //   return status;

}
