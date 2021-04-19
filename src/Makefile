#COPT = -O3 -Wall -ansi -pedantic -I$(HOME)/include
COPT = -g -Wall -ansi -pedantic -I$(HOME)/include
LOPT = -L$(HOME)/lib
LIBS = -lbiop -lgen -lm -lxml2
INCFILES = properties.h

all : matchpatch matchpatchsurface

matchpatch.o : matchpatch.c $(INCFILES)
	$(CC) $(COPT) -c -o $@ $<

matchpatchsurface.o : matchpatchsurface.c $(INCFILES)
	$(CC) $(COPT) -c -o $@ $<

matchpatch : matchpatch.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

matchpatchsurface : matchpatchsurface.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

clean : 
	\rm *.o

install :
	mkdir -p $(HOME)/bin
	cp matchpatch matchpatchsurface $(HOME)/bin
