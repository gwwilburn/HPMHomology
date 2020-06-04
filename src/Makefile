#useful directories
ESLDIR = ../../hmmer4/easel
HMMERDIR = ../../hmmer4/src
H4DIR = ../../hmmer4/nwo

#variables that I will run
INCLUDES = -I${HMMERDIR} -I${H4DIR} -I${ESLDIR}
IFLAGS   = -I${ESLDIR} -I${HMMERDIR} -I${H4DIR} -I../src
LFLAGS   = -L${ESLDIR} -L${HMMERDIR} -L${H4DIR} -L../src
DFLAGS   = -DHPM_3MER
LIBS     = -lm
SLIBS    = ${HMMERDIR}/libhmmer.a ${H4DIR}/libhmmer.a ${ESLDIR}/libeasel.a

#compiler information
CC = gcc
CFLAGS = -Wall -g

# targets for make

all: hpm hpmfile hpm_scoreset hpm_trace potts pottsfile hpmscore hpmbuild hmm_scoreset hpm_scoreops hmm_entropy h4_path_hpm h4_pathalign hmmscore_uniglocal hmmalign_uniglocal hpmscoreIS hpmscore_path hpmalign hpmemit hpmbuild_masked

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

h4_path_hpm: h4_path_hpm.c
	${CC} ${CFLAGS} -c -o h4_path_hpm.o ${IFLAGS} h4_path_hpm.c

h4_pathalign: h4_pathalign.c
	${CC} ${CFLAGS} -c -o h4_pathalign.o ${IFLAGS} h4_pathalign.c

hpmscore: hpmscore.c
	${CC} ${CFLAGS} -o hpmscore ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmscore.c hpm.o hpmfile.o hpm_scoreset.o potts.o pottsfile.o

hpmbuild: hpmbuild.c
	${CC} ${CFLAGS} -o hpmbuild ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmbuild.c hpm.o hpmfile.o potts.o pottsfile.o

hmmscore_uniglocal: hmmscore_uniglocal.c
	${CC} ${CFLAGS} -o hmmscore_uniglocal ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hmmscore_uniglocal.c hmm_scoreset.o

hmmalign_uniglocal: hmmalign_uniglocal.c
	${CC} ${CFLAGS} -o hmmalign_uniglocal ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hmmalign_uniglocal.c

hpmscoreIS: hpmscoreIS.c
	${CC} ${CFLAGS} -o hpmscoreIS ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmscoreIS.c hpm.o hpmfile.o potts.o pottsfile.o hpm_scoreops.o hpm_scoreset.o h4_path_hpm.o h4_pathalign.o

hpmscore_path: hpmscore_path.c
	${CC} ${CFLAGS} -o hpmscore_path ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmscore_path.c hpm.o hpmfile.o hpm_scoreops.o h4_path_hpm.o

hpmalign: hpmalign.c
	${CC} ${CFLAGS} -o hpmalign ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmalign.c hpm.o hpmfile.o potts.o pottsfile.o hmm_entropy.o hpm_scoreops.o hpm_trace.o

hpmemit: hpmemit.c
	${CC} ${CFLAGS} -o hpmemit ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmemit.c hpm.o hpmfile.o hpm_scoreops.o

hpmbuild_masked: hpmbuild_masked.c
	${CC} ${CFLAGS} -o hpmbuild_masked ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpmbuild_masked.c hpm.o hpmfile.o potts.o

3mer:
	 ${CC} ${CFLAGS} -o hpm_3mer ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hpm.c ${DFLAGS} hpmfile.o hpm_scoreset.o potts.o pottsfile.o
	 ${CC} ${CFLAGS} -o hmm_entropy_3mer ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hmm_entropy.c ${DFLAGS} ${SLIBS}

entropy_msa:
	${CC} ${CFLAGS} -o hmm_entropy_msa ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hmm_entropy.c  -DENTROPY_MSA

utest:
	 ${CC} ${CFLAGS} -o hmm_entropy_utest ${IFLAGS} ${LFLAGS} ${SLIBS} ${LIBS} hmm_entropy.c -DHMMENTROPY_TESTDRIVE

clean:
	rm -f hpm.o hpmfile.o hpm_scoreset.o hpm_scoreops.o  hpm_trace.o potts.o pottsfile.o hmm_scoreset.o hmm_entropy.o h4_path_hpm.o h4_pathalign.o hpmscore hpmbuild  hmmscore_uniglocal  hmmalign_uniglocal hpmscoreIS hpmscore_path hpmalign  hpmemit hpmbuild_masked hpm_3mer hmm_entropy_3mer hmm_entropy_msa hmm_entropy_utest




.PHONY: all clean
