BASEDIR = /home/jitu/btp/iitm/preprocessing
SRCDIR = ../src
DSPLIB = ../src
OBJDIR = ../obj
BINDIR = ../bin

all:diarization

diarization: preprocessing.o VQ_Modified.o GMM.o viterbi_realign.o
	gcc preprocessing.o VQ_Modified.o GMM.o viterbi_realign.o -o diarization -lm

viterbi_realign.o: viterbi_realign.c
	gcc -c viterbi_realign.c

GMM.o: GMM.c
	gcc -c GMM.c

VQ_Modified.o: VQ_Modified.c
	gcc -c VQ_Modified.c

preprocessing.o: preprocessing.c
	gcc -c preprocessing.c
clean: rm *.o
	rm diarization
