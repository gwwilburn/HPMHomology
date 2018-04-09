/* hpmfile.h */

#ifndef HPMFILE_INCLUDED
#define HPMFILE_INCLUDED

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"

#include "hpm.h"




typedef struct hmmfile {
	FILE		*f; 		/* pointer to stream for reading models 		*/
	char		*fname;	/* name of the hpm file								*/

}	HPMFILE;

/* hpmfile.c */
extern int hpmfile_Write(FILE *f, HPM *hpm);
#endif
