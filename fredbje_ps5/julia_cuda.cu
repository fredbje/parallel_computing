#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define XSIZE 2560
#define YSIZE 2048
#define MAXITER 255
#define PIXEL(i,j) ((i)+(j)*XSIZE)

typedef struct {
    double re;
    double im;
} complex_t;

typedef unsigned char uchar;

double walltime(){
	static struct timeval t;
	gettimeofday(&t, NULL);
	return (t.tv_sec + 1e-6 * t.tv_usec);
}

/* save 24-bits bmp file, buffer must be in bmp format: upside-down */
void savebmp(const char *name, uchar *buffer, int x, int y) {
	FILE *f = fopen(name, "w");
	if(!f) {
		printf("Error writing image to disk.\n");
		return;
	}
	unsigned int size = x*y*3 + 54;
	uchar header[54]={'B','M',size&255,(size>>8)&255,(size>>16)&255,size>>24,0,
                    0,0,0,54,0,0,0,40,0,0,0,x&255,x>>8,0,0,y&255,y>>8,0,0,1,0,24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	fwrite(header, 1, 54, f);
	fwrite(buffer, 1, XSIZE*YSIZE*3, f);
	fclose(f);
}

/* given iteration number, set a colour */
void fancycolour(uchar *p, int iter) {
	if(iter == MAXITER);
	else if(iter < 8) { p[0] = 128 + iter*16; p[1] = p[2] = 0; }
	else if(iter < 24) { p[0] = 255; p[1] = p[2] = (iter - 8)*16; }
	else if(iter < 160) { p[0] = p[1] = 255 - (iter - 24)*2; p[2] = 255; }
	else { p[0] = p[1] = (iter - 160)*2; p[2] = 255 - (iter - 160)*2; }
}

__device__ void square_complex(complex_t a, complex_t *res){
    (*res).re = a.re*a.re - a.im*a.im;
    (*res).im = 2*a.re*a.im;
}

__device__ void add_complex(complex_t a, complex_t b, complex_t *res){
    (*res).re = a.re + b.re;
    (*res).im = a.im + b.im;
}

__device__ void add_real(complex_t a, int b, complex_t *res){
    (*res).re = a.re + b;
	(*res).im = a.im;
}

__global__ void iterate_pixel(int *pixel, complex_t julia_C){

	// Calculate the range in the y-axis such that we preserve the aspect ratio
	double x_start = -2.01;
	double x_end = 1.0;
	double y_center = 1e-6;
	double step = (x_end - x_start)/XSIZE;
	double y_start = y_center - (step*YSIZE)/2;

	int x =  blockIdx.x * blockDim.x + threadIdx.x;
	int y =  blockIdx.y * blockDim.y + threadIdx.y;

	complex_t z;
	z.re = (x_start + step*x);
	z.im = (y_start + step*y);

	int iter = 0;
	complex_t z_squared;
	while(z.re*z.re + z.im*z.im < 4) {
		if(++iter==MAXITER){
			break;
		} else {
			square_complex(z, &z_squared);
			add_complex(z_squared, julia_C, &z);
		}
	}
	pixel[PIXEL(x,y)]=iter;
}

int main(int argc,char **argv) {
	if(argc==1) {
		puts("Usage: JULIA\n");
		puts("Input real and imaginary part. ex: ./julia 0.0 -0.8");
		return 0;
	}

	complex_t julia_C;
    julia_C.re = strtod(argv[1], NULL);
    julia_C.im = strtod(argv[2], NULL);


/*  // Runs much slower 513ms as opposed to 0.088ms without mallocManaged
	double start_gpu = walltime();

	int *pixeldouble start_gpu = walltime();;
	cudaMallocManaged(&pixel, XSIZE*YSIZE*sizeof(int));
	dim3 grid(XSIZE, YSIZE);
	iterate_pixel<<<grid, 1>>>(pixel, julia_C);
	cudaDeviceSynchronize();

	double end_gpu = walltime();
	printf("Computation complete. It took %7.3f ms\n", end_gpu - start_gpu);

  	//create nice image from iteration counts. take care to create it upside down (bmp format)
 	uchar *buffer = (uchar*)calloc(XSIZE*YSIZE*3, sizeof(uchar));
    for(int i = 0; i < XSIZE; i++) {
        for(int j = 0; j < YSIZE; j++) {
            int p = ((YSIZE - j - 1)*XSIZE + i)*3;
			fancycolour(buffer + p, pixel[PIXEL(i, j)]);
        }
    }
    savebmp("julia_cuda.bmp", buffer, XSIZE, YSIZE);

	cudaFree(pixel);
	free(buffer);

*/



	double start_gpu = walltime();

	int *device_pixel;
	int *host_pixel;
	host_pixel = (int*)malloc(XSIZE*YSIZE*sizeof(int));
	cudaMalloc((void**)&device_pixel, XSIZE*YSIZE*sizeof(int));
	int num_threads = 32;
	dim3 grid_dim(XSIZE/num_threads, YSIZE/num_threads);
	dim3 block_dim(num_threads, num_threads);
	iterate_pixel<<<grid_dim, block_dim>>>(host_pixel, julia_C);
	cudaMemcpy(host_pixel, device_pixel, XSIZE*YSIZE*sizeof(int), cudaMemcpyDeviceToHost);

	double end_gpu = walltime();
	printf("Computation complete. It took %7.3f ms\n", end_gpu - start_gpu);

  	// create nice image from iteration counts. take care to create it upside down (bmp format)
 	uchar *buffer = (uchar*)calloc(XSIZE*YSIZE*3, sizeof(uchar));
    for(int i = 0; i < XSIZE; i++) {
        for(int j = 0; j < YSIZE; j++) {
            int p = ((YSIZE - j - 1)*XSIZE + i)*3;
			fancycolour(buffer + p, host_pixel[PIXEL(i, j)]);
        }
    }
    savebmp("julia_cuda.bmp", buffer, XSIZE, YSIZE);

	cudaFree(device_pixel);
	free(host_pixel);
	free(buffer);

/*
	Number of SMs per card; CL_DEVICE_MAX_COMPUTE_UNITS: 30.
	Each SM can run 8 blocs simultaneously.
	Thread blocks have maximum dimensions CL_DEVICE_MAX_WORK_ITEM_SIZES: 512 / 512 / 64
	Maximum number of threads, 512, that can be grouped together in a block. CL_DEVICE_MAX_WORK_GROUP_SIZE: 512. F.ex. 512 x 1 x 1
	specifies the warp size: CL_NV_DEVICE_WARP_SIZE: 32
*/
    return 0;
}
