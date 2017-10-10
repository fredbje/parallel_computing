CFLAGS+= -std=c99 -g -O3 -Wall

# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
CFLAGS+= -DDEBUG
util_objects = CA.o bitmap.o
serial_objects = RPS.o

SHELL := /bin/bash

all : $(util_objects) $(serial_objects)
	gcc $(util_objects) $(serial_objects) -o serialRPS

remake : clean all

.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4 && rm -f serialRPS

.PHONY : video
video : clean all
	./serialRPS && ffmpeg -framerate 60 -i data/CA-%000d.bmp -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
