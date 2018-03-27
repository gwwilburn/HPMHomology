#include <stdio.h>

#include "esl_msa.h"
#include "esl_getopts.h"

struct cfg_s { /* shared configuration in master and workers */
	int		argc;
	char		**argv;
};


static ESL_OPTIONS options[] = {
	/* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
	{ "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
	{ "--amino",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<msafile> contains protein alignments",                   3 },
	{ "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<msafile> contains DNA alignments",                       3 },
	{ "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<msafile> contains RNA alignments",                       3 },
};


int main(int argc, char *argv[])
{
	ESL_GETOPTS     *go;

	/* parse command line */
	go = esl_getopts_Create(options);

	fprintf(stdout, "hello world!\n");
}

