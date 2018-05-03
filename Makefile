all: clcg4.h clcg4.c assignment4-5.c
	gcc -I. -O3 -c clcg4.c -o clcg4.o
	mpicc -I. -O3 assignment4-5.c clcg4.o -o assignment4-5 
	
bgqall:clcg4.h clcg4.c assignment4-5.c	
	gcc -I. -O3 -c clcg4.c -o clcg4.o
	mpicc -I. -O3 assignment4-5.c clcg4.o -o assignment4-5 -lpthread

dbg: clcg4.h clcg4.c assignment4-5.c
	gcc -I. -Wall -g -c clcg4.c -o clcg4.o
	mpicc -I. -Wall -g assignment4-5.c clcg4.o -o assignment4-5 -lpthread

clean:
	rm assignment4-5 *.dat out*

