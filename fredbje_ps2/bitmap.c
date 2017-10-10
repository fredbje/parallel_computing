#include "bitmap.h"

/* save 24-bits bmp file, buffer must be in bmp format: upside down */
void savebmp(char *name, uchar *buffer, int x, int y){
	FILE *f = fopen(name, "wb");
	if(!f){
		printf("Error writing image to disk.\n");
		return;
	}
	unsigned int size = x*y*3 + 54;
	uchar header[54]={'B','M',size&255,(size>>8)&255,(size>>16)&255,size>>24,0,
                    0,0,0,54,0,0,0,40,0,0,0,x&255,x>>8,0,0,y&255,y>>8,0,0,1,0,24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	fwrite(header, 1, 54, f);
	fwrite(buffer, 1, IMG_X*IMG_Y*3, f);
	fclose(f);
}

void colorize(uchar* p, cell my_cell) {
	if(my_cell.color == WHITE){
		p[0] = 255;
		p[1] = 255;
		p[2] = 255;
	 }
		else if(my_cell.color == ROCK) {
	    p[0] = 0;
	    p[1] = 0;
	    p[2] = 255;
	 }
		else if(my_cell.color == SCISSOR) {
	    p[0] = 0;
	    p[1] = 255;
	    p[2] = 0;
	 }
		else if(my_cell.color == PAPER) {
	    p[0] = 255;
	    p[1] = 0;
	    p[2] = 0;
	}
}

void make_bmp(cell** image){
	unsigned char* buffer = malloc(IMG_X*IMG_Y*3);
	for(int r = 0; r < IMG_Y; r++){
		for(int c = 0; c < IMG_X; c++){
			// Make image upside down because bmp
			int p = ((IMG_Y-r-1) * IMG_X + c) * 3;
		  	colorize(buffer+p, image[r][c]);
		}
	}

	char filename[50];
	sprintf(filename, "MPI_petri.bmp");
	/* write image to disk */
	printf("Saving BMP\n");
	savebmp(filename,buffer, IMG_X, IMG_Y);
	free(buffer);
}


