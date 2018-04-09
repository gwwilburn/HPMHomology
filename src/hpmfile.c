#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"

#include "hpm.h"
#include "hpmfile.h"


/* Functions for writing .hpm files */
static int printprob(FILE *fp, int fieldwidth, float p);


int hpmfile_Write(FILE *fp, HPM *hpm) {
	int x;
	int status;

	fprintf(fp, ".hpm file\n");
	fprintf(fp, "LENG 	%d\n", hpm->M);

	/* model parameters */
	fprintf(fp, "HPM     ");
	for (x = 0; x < hpm->abc->K; x++) {
		fprintf(fp, "     %c   ", hpm->abc->sym[x]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "        %8s %8s %8s %8s\n",
				"m->m", "m->i", "i->m", "i->i");

	/* need to add avg h_i values in compo line? */
	fprintf(fp, "  COMPO ");
	/* temporary fix to make this look like a .hmm file */
	float compo = 0.5;
	for (x = 0; x < hpm->abc->K; x++) {
		if ( (status = printprob(fp, 8, compo)) != eslOK) return status;
	}





	return eslOK;
}


/* function for writing probabilities in log2 */
static int
printprob(FILE *fp, int fieldwidth, float p) {
	if      (p == 0.0) { if (fprintf(fp, " %*s",   fieldwidth, "*")      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
	else if (p == 1.0) { if (fprintf(fp, " %*.5f", fieldwidth, 0.0)      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
	else               { if (fprintf(fp, " %*.5f", fieldwidth, -logf(p)) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
	return eslOK;
}



