CC = cc
CFLAGS = -O2 
#CFLAGS = -g
LFLAGS = -L$(HOME)/lib
IFLAGS = -I$(HOME)/include
DBMLIB = -lgdbm

all : makecadb searchcadb


makecadb : makecadb.c
	$(CC) $(IFLAGS) $(LFLAGS) $(CFLAGS) -o makecadb makecadb.c -lbiop -lgen -lm

searchcadb : searchcadb.c
	$(CC) $(IFLAGS) $(LFLAGS) $(CFLAGS) -o searchcadb searchcadb.c -lgen $(DBMLIB)
