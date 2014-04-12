BASEDIR = /home/jitu/btp/c_code
SRCDIR = src
DSPLIB = src
OBJDIR = obj
BINDIR = bin

all:bin/diar

bin/diar: $(OBJDIR)/preprocessing.o $(OBJDIR)/VQ_Modified.o $(OBJDIR)/GMM.o 
	gcc -g $(OBJDIR)/preprocessing.o $(OBJDIR)/VQ_Modified.o $(OBJDIR)/GMM.o -o $(BINDIR)/diar -lm

$(OBJDIR)/GMM.o: $(SRCDIR)/GMM.c
	gcc -c $(SRCDIR)/GMM.c -o $(OBJDIR)/GMM.o

$(OBJDIR)/VQ_Modified.o: $(SRCDIR)/VQ_Modified.c
	gcc -c $(SRCDIR)/VQ_Modified.c -o $(OBJDIR)/VQ_Modified.o

$(OBJDIR)/preprocessing.o: $(SRCDIR)/preprocessing.c
	gcc -c $(SRCDIR)/preprocessing.c -o $(OBJDIR)/preprocessing.o

clean:  $(OBJDIR)/*.o $(BINDIR)/diar
	rm -rf $(BINDIR)/diar $(OBJDIR)/*.o
