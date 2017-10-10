SHELL := /bin/bash

pth :
	gcc -std=c99 -g -O1 -Wall CA.c bitmap.c p_RPS.c -o pthreadRPS -lpthread

remake : clean all

.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4
