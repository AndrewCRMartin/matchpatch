#COPT = -O3 -Wall -ansi -pedantic -I$(HOME)/include
COPT = -g -Wall -ansi -pedantic -I$(HOME)/include
LOPT = -L$(HOME)/lib
LIBS = -lbiop -lgen -lm -lxml2

all : match surface


.c.o :
	$(CC) $(COPT) -c -o $@ $<

match : match.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

surface : surface.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

clean : 
	\rm *.o

