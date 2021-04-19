#COPT = -O3 -Wall -ansi -pedantic -I$(HOME)/include
COPT = -g -Wall -ansi -pedantic -I$(HOME)/include
LOPT = -L$(HOME)/lib
LIBS = -lbiop -lgen -lm -lxml2
INCFILES = properties.h

all : match surface

match.o : match.c $(INCFILES)
	$(CC) $(COPT) -c -o $@ $<

surface.o : surface.c $(INCFILES)
	$(CC) $(COPT) -c -o $@ $<

match : match.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

surface : surface.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

clean : 
	\rm *.o

install :
	mkdir -p $(HOME)/bin
	cp match surface $(HOME)/bin
