#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "julia_mpi.h"
#include <string.h>
#include "bitmap.h"

complex_t square_complex(complex_t a){
    complex_t b;
    b.re = a.re*a.re - a.im*a.im;
    b.im = 2*a.re*a.im;
    return b;
}

complex_t add_complex(complex_t a, complex_t b){
    complex_t c;
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    return c;
}

complex_t add_real(complex_t a, int b){
    complex_t c;
    c.re = a.re + b;
    c.im = a.im;
    return c;
}

int main(int argc,char **argv) {
	if(argc==1) {
		puts("Usage: JULIA\n");
		puts("Input real and imaginary part. ex: ./julia 0.0 -0.8");
		return 0;
	}

	int rank, numprocs;
	MPI_Init(NULL,NULL);

 	double startT, endT;
	startT=MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	/* Calculate the range in the y-axis such that we preserve the
	   aspect ratio */
	double x_start=-2.01;
   	double x_end=1;
	double step=(x_end-x_start)/XSIZE;
	double ycenter = step*YSIZE*rank/(numprocs) - step*YSIZE*(numprocs-1)/(2*numprocs);
	double ylower = ycenter - (step*YSIZE)/(2*numprocs);

    complex_t julia_C;
    julia_C.re = strtod(argv[1], NULL);
    julia_C.im = strtod(argv[2], NULL);

	complex_t z;
	int iter;
	int* local_pixel = (int*)malloc(sizeof(int)*YSIZE*XSIZE/numprocs);
	for(int y = 0; y < YSIZE/numprocs; y++){
		for(int x = 0; x < XSIZE; x++){
			z.re = (x_start + step*x);
			z.im = (ylower + step*y);
			iter=0;
	        while(z.re*z.re + z.im*z.im < 4) {
	            if(++iter==MAXITER) break;
	            else z=add_complex(square_complex(z),julia_C);
	        }
			local_pixel[PIXEL(x,y)]=iter;
		}
	}

	int* pixel = NULL;
	if(rank==0){
		pixel = (int*)malloc(XSIZE*YSIZE*sizeof(int));
	}

	MPI_Gather(local_pixel, XSIZE*YSIZE/numprocs, MPI_INT, pixel, XSIZE*YSIZE/numprocs, MPI_INT, 0, MPI_COMM_WORLD);

  	/* create nice image from iteration counts. take care to create it upside
     down (bmp format) */
	if(rank==0){
	    unsigned char *buffer=calloc(XSIZE*YSIZE*3,1);
	    for(int i=0;i<XSIZE;i++) {
	        for(int j=0;j<YSIZE;j++) {
	            int p=((YSIZE-j-1)*XSIZE+i)*3;
	            fancycolour(buffer+p,pixel[PIXEL(i,j)]);
	        }
	    }
	    savebmp("julia.bmp",buffer,XSIZE,YSIZE);
		free(pixel);
	}
	free(local_pixel);
	endT=MPI_Wtime();
	printf("Time elapsed for process %d: %f [s]\n", rank, endT-startT);
	MPI_Finalize();
    return 0;
}
