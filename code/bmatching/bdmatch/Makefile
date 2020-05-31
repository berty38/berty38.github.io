CC=gcc
CCFLAGS= -Wall -g -O3
EXECUTABLE_LIST=bdmatch

## Targets ##

all:
	make bdmatch

bdmatch: bdmatch.c bp.o quickselect.o
	$(CC) $(CCFLAGS) -o bdmatch bdmatch.c bp.o quickselect.o -lm -lpthread

bp.o: bp.c
	$(CC) $(CCFLAGS) -c bp.c

quickselect.o: quickselect.c bdmatch.h
	$(CC) $(CCFLAGS) -c quickselect.c

clean:
	@/bin/rm -vrf $(EXECUTABLE_LIST) *~ *.o *.s *.dSYM a.out
