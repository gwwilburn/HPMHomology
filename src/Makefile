#usefull directories
ESLDIR = ../lib/hmmer/easel
HMMERDIR = ../lib/hmmer/src

#variables that I will run
INCLUDES = -I${HMMERDIR} -I${ESLDIR}
IFLAGS   = -I${ESLDIR} -I${HMMERDIR} -I../src
LFLAGS   = -L${ESLDIR} -L${HMMERDIR} -L../src
LIBS     = -leasel     -lm -lhmmer

#compiler information
CC = gcc
CFLAGS = -Wall -g

# targest for make

all: hpm hpmfile potts pottsfile hpmscore

hpmscore: hpmscore.c
	${CC} ${CFLAGS} -o hpmscore ${IFLAGS} ${LFLAGS} ${LIBS} hpmscore.c hpm.c hpmfile.c potts.c pottsfile.c

hpm: hpm.c
	${CC} ${CFLAGS} -c -o hpm.o ${IFLAGS} hpm.c

hpmfile: hpmfile.c
	${CC} ${CFLAGS} -c -o hpmfile.o ${IFLAGS} hpmfile.c

potts: potts.c
	${CC} ${CFLAGS} -c -o potts.o ${IFLAGS} potts.c

pottsfile: pottsfile.c
	${CC} ${CFLAGS} -c -o pottsfile.o ${IFLAGS} pottsfile.c



clean:
	rm -f hpmscore hpm.o hpmfile.o potts.o pottsfile.o

.PHONY: all clean
