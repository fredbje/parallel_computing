CFLAGS+= -std=c99 -g -O1 -Wall

# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
CFLAGS+= -DDEBUG
util_objects = CA.o bitmap.o
parallel_objects = RPS_MPI.o

SHELL := /bin/bash

all: $(util_objects) $(parallel_objects)
	mpicc $(util_objects) $(parallel_objects) -o parallelRPS

RPS_MPI.o: RPS_MPI.c RPS_MPI.h
	mpicc -c RPS_MPI.c RPS_MPI.h

remake : clean all

.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4 && rm -f paralellRPS

