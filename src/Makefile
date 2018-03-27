#usefull directories
ESLDIR = ../lib/hmmer/easel
HMMERDIR = ../lib/hmmer/src

#variables that I will run
INCLUDES = -I${HMMERDIR} -I${ESLDIR}
IFLAGS   = -I${ESLDIR} -I${HMMERDIR}
LFLAGS   = -L${ESLDIR} -L${HMMERDIR}
LIBS     = -leasel     -lm -lhmmer

#compiler information
CC = gcc
CFLAGS = -Wall -g

# targest for make

all: ScoreSequenceHPM

ScoreSequenceHPM: ScoreSequenceHPM.c
	${CC} ${CFLAGS} -o ScoreSequenceHPM ${IFLAGS} ${LFLAGS} ${LIBS} ScoreSequenceHPM.c

clean:
	rm -f ScoreSequenceHPM

.PHONY: all clean
