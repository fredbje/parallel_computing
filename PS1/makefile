CC = mpicc
CFLAGS+= -Wall -std=c99 -O3 -I.
# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
# CFLAGS+= -DDEBUG
OBJECTS = julia_mpi.o bitmap.o
SHELL := /bin/bash
TARGET = julia_mpi

# When running "make" in your terminal, the first target will be chosen.
# The syntax of make is basically:
# name : [stuff I need to make this target]

# In this case the dependencies are easy to figure out, so I will not elaborate further
all : $(OBJECTS)
	$(CC) $(CLFLAGS) $(OBJECTS) -o $(TARGET)

# In this target [stuff I need to make this target] is two other targets, namely clean and all
# This command simply runs the clean target, then the all target, thus recompiling everything.
remake : clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f $(TARGET)

# Finally, the test target. Builds the 'all' target, then runs the test script on the output
#.PHONY : test
#test : all
#	./test.sh myProgram

# julia_mpi: julia_mpi.c bitmap.c
#	gcc -std=c99 -o julia_mpi julia_mpi.c bitmap.c -I.
