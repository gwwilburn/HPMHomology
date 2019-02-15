#usefull directories
ESLDIR = ../lib/hmmer/easel
HMMERDIR = ../lib/hmmer/src

#variables that I will run
INCLUDES = -I${HMMERDIR} -I${ESLDIR}
IFLAGS   = -I${ESLDIR} -I${HMMERDIR} -I../src
LFLAGS   = -L${ESLDIR} -L${HMMERDIR} -L../src
DFLAGS   = -DHPM_3MER
LIBS     = -lm
SLIBS    = ${HMMERDIR}/libhmmer.a ${ESLDIR}/libeasel.a

#compiler information
CC = gcc
CFLAGS = -Wall -g

# targest for make

all: hpm hpmfile hpm_scoreset hpm_trace potts pottsfile hpmscore hpmbuild hmm_scoreset hpm_scoreops hmm_entropy hmmscore_uniglocal hmmalign_uniglocal hpmscoreIS hpmalign

hpmscore: hpmscore.c
	${CC} ${CFLAGS} -o hpmscore ${IFLAGS} ${LFLAGS} ${LIBS} hpmscore.c hpm.o hpmfile.o hpm_scoreset.o potts.o pottsfile.o ${SLIBS}

hpmbuild: hpmbuild.c
	${CC} ${CFLAGS} -o hpmbuild ${IFLAGS} ${LFLAGS} ${LIBS} hpmbuild.c hpm.o hpmfile.o potts.o pottsfile.o ${SLIBS}

hmmscore_uniglocal: hmmscore_uniglocal.c
	${CC} ${CFLAGS} -o hmmscore_uniglocal ${IFLAGS} ${LFLAGS} ${LIBS} hmmscore_uniglocal.c hmm_scoreset.o ${SLIBS}

hmmalign_uniglocal: hmmalign_uniglocal.c
	${CC} ${CFLAGS} -o hmmalign_uniglocal ${IFLAGS} ${LFLAGS} ${LIBS} hmmalign_uniglocal.c ${SLIBS}

hpmscoreIS: hpmscoreIS.c
	${CC} ${CFLAGS} -o hpmscoreIS -fopenmp ${IFLAGS} ${LFLAGS} ${LIBS} hpmscoreIS.c hpm.o hpmfile.o potts.o pottsfile.o hmm_entropy.o hpm_scoreops.o hpm_scoreset.o ${SLIBS}

hpmalign: hpmalign.c
	${CC} ${CFLAGS} -o hpmalign ${IFLAGS} ${LFLAGS} ${LIBS} hpmalign.c hpm.o hpmfile.o potts.o pottsfile.o hmm_entropy.o hpm_scoreops.o hpm_trace.o ${SLIBS}

hpm: hpm.c
	${CC} ${CFLAGS} -c -o hpm.o ${IFLAGS} hpm.c

hpmfile: hpmfile.c
	${CC} ${CFLAGS} -c -o hpmfile.o ${IFLAGS} hpmfile.c

hpm_scoreset: hpm_scoreset.c
	${CC} ${CFLAGS} -c -o hpm_scoreset.o ${IFLAGS} hpm_scoreset.c

hpm_scoreops: hpm_scoreops.c
	${CC} ${CFLAGS} -c -o hpm_scoreops.o ${IFLAGS} hpm_scoreops.c

hpm_trace: hpm_trace.c
	${CC} ${CFLAGS} -c -o hpm_trace.o ${IFLAGS} hpm_trace.c

potts: potts.c
	${CC} ${CFLAGS} -c -o potts.o ${IFLAGS} potts.c

pottsfile: pottsfile.c
	${CC} ${CFLAGS} -c -o pottsfile.o ${IFLAGS} pottsfile.c

hmm_scoreset: hmm_scoreset.c
	${CC} ${CFLAGS} -c -o hmm_scoreset.o ${IFLAGS} hmm_scoreset.c

hmm_entropy: hmm_entropy.c
	${CC} ${CFLAGS} -c -o hmm_entropy.o ${IFLAGS} hmm_entropy.c


3mer:
	 ${CC} ${CFLAGS} -o hpm_3mer ${IFLAGS} ${LFLAGS} ${LIBS} hpm.c ${DFLAGS} hpmfile.o hpm_scoreset.o potts.o pottsfile.o ${SLIBS}
	 ${CC} ${CFLAGS} -o hmm_entropy_3mer ${IFLAGS} ${LFLAGS} ${LIBS} hmm_entropy.c ${DFLAGS} ${SLIBS}


clean:
	rm -f hpmscore hpmbuild hmmscore_uniglocal hmmalign_uniglocal hpm.o hpmfile.o hpm_scoreset.o potts.o pottsfile.o hmm_scoreset.o hpm_scoreops.o hpm_trace.o hmm_entropy.o hpm_3mer hpmalign

.PHONY: all clean
