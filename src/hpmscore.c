#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_msafile.h"
#include "esl_getopts.h"
#include "esl_composition.h"

#include "p7_config.h"
#include "hmmer.h"

#include "hpm.h"
#include "hpmfile.h"

struct cfg_s { /* shared configuration in master and workers */
	ESL_ALPHABET	 *abc;
	ESL_MSAFILE   	 *afp;     	/* input ali file						*/
	ESL_MSA			 *msa;     	/* input msa							*/
	ESL_SQFILE		 *sqfp;     /* input sequence file           */
	ESL_SQ			 *sq;       /* sequence object					*/
	P7_HMMFILE		 *hfp;		/* input hmm file						*/
	P7_HMM			 *hmm;		/* input hmm							*/
	HPM				 *hpm;		/* hidden potts model				*/
	HPM_SCORESET    *hpm_ss;   /* for storing hpm sores         */

};


static ESL_OPTIONS options[] = {
	/* name         type        default env   range togs  reqs  incomp      help                                                   docgroup */
	{ "-h",         eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },

	/* Options forcing which alphabet we're working in (normally autodetected) */
	{ "--amino",  eslARG_NONE,  FALSE,  NULL, NULL, NULL,NULL,"--dna,--rna",    "<seqfile> contains protein sequences",                    2 },
	{ "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<seqfile> contains DNA sequences",                        2 },
	{ "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<seqfile> contains RNA sequences",                        2 },
	{ 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <seqfile>";

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

int CalculateHamiltonian(HPM *hpm, P7_TRACE **tr, ESL_MSA *msa, HPM_SCORESET *hpm_ss) {
	int   z;    /* index for trace elements 	     */
	int   y;    /* index for trace elements 	     */
	int   n;    /* sequence index          	     */
	int   i;    /* match state index			        */
	int   j;    /* match state index	    		     */
	int   a;    /* residue index            	     */
	int   b;    /* residue index           	     */
	int   idx;  /* index for potts parameters      */
	float E;    /* pseudo-energy, aka Hamiltonian  */

	/* loop over sequences */
	for (n=0; n < msa->nseq; n++) {

		E = 0.0;

		/* print out sequence info */
		fprintf(stdout, "%d: %d %d %d %s\n", z, tr[n]->L, tr[n]->M, tr[n]->N, msa->sqname[n]);

		/* loop over trace positions for this seq */
		for (z = 0; z < tr[n]->N; z++) {

			/* check if we are in a match or delete state */
			if (tr[n]->st[z] == 2 || tr[n]->st[z] == 6) {

				i = tr[n]->k[z];

				//fprintf(stdout, "\t %d\n", i);

				/* we have a match position */
				//if (tr[n]->st[z] == P7T_M ) {
				if (tr[n]->st[z] == 2) {
					a = msa->ax[n][tr[n]->i[z]];
				}

				/* we have a delete position */
				//else if ( tr[n]->st[z] == p7T_D ) {
				else if (tr[n]->st[z] == 6) {
					a = 20;
				}

				E = E + hpm->h[i][a];

				/* now add e_ij terms to pseudo-energy */
				for (y = z+1; y < tr[n]->N; y++) {

					/* check if we are in a match position */
					 if (tr[n]->st[y] == 2 || tr[n]->st[y] == 6) {

						 j = tr[n]->k[y];

						/* we have a match state */
						if (tr[n]->st[z] == 2) {
							b =  msa->ax[n][tr[n]->i[y]];
						}

						/* we have a match position */
						else if (tr[n]->st[z] == 6) {
							b = 20;
						}

						idx = IDX(a,b,msa->abc->K+1);
						//fprintf(stdout, "\t\t %d %d, %d, %d, %d \n", i, j, idx, a, b);
						E += hpm->e[i][j][idx];

					 }
				}
			}

		}
		fprintf(stdout, "E: %.6f\n", E);
	}


	return eslOK;
}

int main(int argc, char *argv[])
{
	ESL_GETOPTS    *go;															/* cmd line configuration 			*/
	struct cfg_s    cfg;  														/* application configuration 		*/
	char           *seqfile 					= NULL;						/* sequence file path		 		*/
	char           *alifile 					= NULL;						/* alignment file path		 		*/
	int             infmt               	= eslSQFILE_UNKNOWN;
	int				 alphatype    				= eslUNKNOWN;         	/*	for guessing alphabet			*/
	char				*hmmfile 					= NULL;						/* hmm file path						*/
	int				 status;														/* easel return code 				*/
	char            errbuf[eslERRBUFSIZE];

	/* parse command line */
	go = esl_getopts_Create(options);
	if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
	if (esl_opt_VerifyConfig(go) 					 != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
	//if (esl_opt_ArgNumber(go)						 != 2) 	  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);
	if (esl_opt_ArgNumber(go)						 != 3) 	  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

	hmmfile  = esl_opt_GetArg(go, 1);
	seqfile  = esl_opt_GetArg(go, 2);
	alifile  = esl_opt_GetArg(go, 3);

	/* Set up the configuration structure shared amongst functions here */
	cfg.abc = NULL; /* until we read the seq file below */

	/* open the sequence file */
	status = esl_sqfile_Open(seqfile, infmt, NULL, &cfg.sqfp);

	/* Determine alphabet */
	if 		(esl_opt_GetBoolean(go, "--amino"))		cfg.abc = esl_alphabet_Create(eslAMINO);
	else if  (esl_opt_GetBoolean(go, "--dna"))		cfg.abc = esl_alphabet_Create(eslDNA);
	else if  (esl_opt_GetBoolean(go, "--rna"))		cfg.abc = esl_alphabet_Create(eslRNA);
	else {
		status = esl_sqfile_GuessAlphabet(cfg.sqfp, &alphatype);
		cfg.abc = esl_alphabet_Create(alphatype);
	}

	esl_alphabet_SetEquiv(cfg.abc, '.', '-'); /* allow '.' as gap char too */

	/* prepare to read in sequences */
	cfg.sq  = esl_sq_CreateDigital(cfg.abc);
	esl_sqfile_SetDigital(cfg.sqfp, cfg.abc);

	/* loop through seqs in file */

	while ((status = esl_sqio_Read(cfg.sqfp, cfg.sq)) == eslOK) {
		//esl_sqio_Write(stdout, cfg.sq, eslSQFILE_FASTA, /*update=*/FALSE);
		esl_sq_Reuse(cfg.sq);

	}

	/* common error handling */
	if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(cfg.sqfp));
	else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);


	/* Open the .hmm file */
	status = p7_hmmfile_OpenE(hmmfile, NULL, &cfg.hfp, errbuf);
	/* error handling below copied from hmmsearch */
	if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
	else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
	else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);


	/* read the .hmm file */
	status = p7_hmmfile_Read(cfg.hfp, &cfg.abc, &cfg.hmm);

	/* open the msa file */
	status = esl_msafile_Open(&(cfg.abc), alifile, NULL, eslMSAFILE_UNKNOWN, NULL, &cfg.afp);
	if (status != eslOK) esl_msafile_OpenFailure(cfg.afp, status);

	/* Read in the MSA */
	status = esl_msafile_Read(cfg.afp, &cfg.msa);
	if (status != eslOK) esl_msafile_ReadFailure(cfg.afp, status);


	/* Set up the score set object */
	cfg.hpm_ss = hpm_scoreset_Create(cfg.msa->nseq);

	/* temporary N-mer hpm for testing, to be deleted */
	int M = 99; // for PF00042 seed
	cfg.hpm = hpm_Create(M, cfg.abc);

	int a;
	int b;
	int i;
	int j;


	/* fill in e_ij's */
	for (a = 0; a < cfg.abc->K+1; a++) {
		for (i = 1; i < M+1; i++) {
			cfg.hpm->h[i][a] = (0.1 * (i))  + (0.0001 * a);
			for (b = 0; b < cfg.abc->K+1; b++) {
				for (j = 1; j < i; j++) {
					cfg.hpm->e[i][j][IDX(a,b,cfg.hpm->abc->K+1)] = (1.0 * (i)) + (0.01*(j)) + (0.0001 * a) + (0.000001 * b);
					cfg.hpm->e[j][i][IDX(b,a,cfg.hpm->abc->K+1)] = (1.0 * (i)) + (0.01*(j)) + (0.0001 * a) + (0.000001 * b);
				}
			}

		}

	}



	int k;
	for (k = 0; k < M+1; k++) {
		cfg.hpm->t[k][HPM_MM] = 1.0 - (0.001*(k+1));
		cfg.hpm->t[k][HPM_MI] = (0.001*(k+1));
		if (k < M+1) {
			cfg.hpm->t[k][HPM_IM] = 1.0 - (0.001*(k+1));
			cfg.hpm->t[k][HPM_II] = (0.001*(k+1));
		}

	}



	char *hpm_outfile = "test.hpm";
	FILE *hpm_outfp 	= 	NULL;
	if ((hpm_outfp = fopen(hpm_outfile, "w")) == NULL) esl_fatal("Failed to open output hpm file %s for writing", hpm_outfile);
	hpmfile_Write(hpm_outfp, cfg.hpm);
	fclose(hpm_outfp);
	HPM *hpmr = NULL;
	hpmr = hpmfile_Read("test.hpm", cfg.abc, errbuf);

	char *hpm_outfile2 = "test2.hpm";
	FILE *hpm_outfp2   = NULL;
	if ((hpm_outfp2 = fopen(hpm_outfile2, "w")) == NULL) esl_fatal("Failed to open output hpm file %s for writing", hpm_outfile2);
	hpmfile_Write(hpm_outfp2, hpmr);
	fclose(hpm_outfp2);

	/* check if file is being read, print line numbers */

	/*
	HPM *hpmDummy = NULL;
	hpmDummy = hpmfile_ReadDummy("test.hpm", cfg.abc, errbuf);
	*/


	/* get traces for seqs in msa */

	/* initiate empty trace array and match column array*/
	P7_TRACE **tr = NULL;
	int *matassign = NULL;

	/* allocate memory for trace and match arrays */
	ESL_ALLOC(tr, sizeof(P7_TRACE *) * cfg.msa->nseq);
	ESL_ALLOC(matassign, sizeof(int) * (cfg.msa->alen + 1));

	for (i=0; i < cfg.msa->alen; i++){
		if (cfg.msa->rf[i] != '.') {
			matassign[i+1] = 1;
		}
		else {
			matassign[i+1] = 0;
		}
		//fprintf(stdout, "%d: %c, %d\n", i, cfg.msa->rf[i], matassign[i+1]);
	}

	/* assign traces to msa seqs */
	p7_trace_FauxFromMSA(cfg.msa, matassign, p7_MSA_COORDS, tr);

	/* now explore the newly created trace object */
	/*
	int z;
	for (z=0; z < cfg.msa->nseq; z++) {
		fprintf(stdout, "%d: %d %d %d %s\n", i, tr[z]->L, tr[z]->M, tr[z]->N, cfg.msa->sqname[z]);
		fprintf(stdout, "\tj,i,k\n");
		for (j = 0; j < tr[z]->N; j++) {
			fprintf(stdout, "\t%d,%d,%d,%d\n", j, tr[z]->i[j],tr[z]->k[j],tr[z]->st[j]);
			if (tr[z]->st[j] == 2) fprintf(stdout, "Match state!\n");
		}
	}
	*/

	/* Calculate hamiltonian for all sequences */
	CalculateHamiltonian(cfg.hpm, tr, cfg.msa, cfg.hpm_ss);

	fprintf(stdout, "hello world!\n");

	/* clean up */
	esl_alphabet_Destroy(cfg.abc);
	esl_sqfile_Close(cfg.sqfp);
	esl_sq_Destroy(cfg.sq);
	esl_msa_Destroy(cfg.msa);
	esl_msafile_Close(cfg.afp);
	p7_hmm_Destroy(cfg.hmm);

	free(matassign);

	p7_hmmfile_Close(cfg.hfp);
	esl_getopts_Destroy(go);

	return 0;

	ERROR:
		return status;

}
