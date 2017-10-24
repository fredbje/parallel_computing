SHELL := /bin/bash


omp :
	gcc -std=c99 -g -O3 -fopenmp CA.c bitmap.c RPS.c -o RPS_omp -lm
pth :
	gcc -std=c99 -g -O3 CA.c bitmap.c p_RPS.c -o RPS_pthread -lm -pthread

remake : clean all

.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4
