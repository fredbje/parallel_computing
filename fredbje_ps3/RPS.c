#include "RPS.h"
#include <time.h>
#include <omp.h>
#include <math.h>

void swap_petris();
void init_petri();

cell* petri_A;
cell* petri_B;

int num_threads, thread_id;

int main(int argc, char** argv){

    if(argc==1) {
        puts("Input number of threads. ex: ./julia 4");
        return 0;
    }
    num_threads = strtod(argv[1], NULL);
    int temp = sqrt(num_threads);
    if (num_threads != temp*temp){
        puts("Number of threads must be a square integer");
        return 0;
    }

    printf("running %d iterations..\n",ITERATIONS);

    srand(time(NULL));
    petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
    petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

    init_petri(petri_A);

    #pragma omp parallel private(num_threads, thread_id){

        int bound = IMG_Y / num_threads;
        int start = thread_id * bound;
        int finish = start + bound;
       

       for(int i = 0; i < ITERATIONS; i++){
            if(i % 100 == 0 && *(thread_id)==0){printf("Progress %d %%\n", i*100 / ITERATIONS);}
            if(i % 2 == 0){
                iterate_image(petri_A, petri_B, start, finish);
            }else{
                iterate_image(petri_B, petri_A, start, finish);
            #pragma omp barrier
        }
        
    }

    if(ITERATIONS % 2 == 0){
        make_bmp(petri_A, 0);
    } else {
        make_bmp(petri_B, 0);
    }

    free(petri_A);
    free(petri_B);
    return 0;
}

void init_petri(cell* petri){
    for(int ii = 0; ii < 100; ii++){
        int rx = rand() % (IMG_X - 1);
        int ry = rand() % (IMG_Y - 1);
        int rt = rand() % 4;

        petri[TRANS(rx,ry)].color = rt;
        petri[TRANS(rx,ry)].strength = 1;
    }
}






