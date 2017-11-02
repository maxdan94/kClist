CC=gcc
CFLAGS=-O9

all: Arboricity kClist kClistCore kClistDens kClistEdgeParallel kClistNodeParallel

Arboricity : Arboricity.c
	$(CC) $(CFLAGS) Arboricity.c -o Arboricity

kClist : kClist.c
	$(CC) $(CFLAGS) kClist.c -o kClist

kClistCore : kClistCore.c
	$(CC) $(CFLAGS) kClistCore.c -o kClistCore

kClistDens : kClistDens.c
	$(CC) $(CFLAGS) kClistDens.c -o kClistDens -fopenmp

kClistEdgeParallel : kClistEdgeParallel.c
	$(CC) $(CFLAGS) kClistEdgeParallel.c -o kClistEdgeParallel -fopenmp

kClistNodeParallel : kClistNodeParallel.c
	$(CC) $(CFLAGS) kClistNodeParallel.c -o kClistNodeParallel -fopenmp


clean:
	rm Arboricity kClist kClistCore kClistDens kClistEdgeParallel kClistNodeParallel
