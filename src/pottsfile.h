/* pottsfile.h */

#ifndef POTTSFILE_INCLUDED
#define POTTSFILE_INCLUDED

#include <stdio.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"

#include "potts.h"

typedef struct pottsfile{
	FILE     *f;      /* pointer to stream for reading models      */
	char     *fname;  /* name of the hpm file                      */
} POTTSFILE;

/* pottsfile.c */
POTTS *pottsfile_Read(char *f, ESL_ALPHABET *abc, char *errbuff);


#endif
