#include "RPS.h"
#include <time.h>
#include <omp.h>
#include <math.h>

void init_petri();

cell* petri_A;
cell* petri_B;

int main(int argc, char** argv){

    if(argc==1) {
        puts("Input number of threads. ex: ./RPS_omp 4");
        return 0;
    }
    int n = (int)strtod(argv[1], NULL);
    int temp = sqrt(n);
    if (n != temp*temp){
        puts("Number of threads must be a square integer");
        return 0;
    }

    printf("running %d iterations..\n",ITERATIONS);

    srand(time(NULL));

    petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
    petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

    init_petri(petri_A);

    clock_t start_clock = clock();

    #pragma omp parallel num_threads(n)
    {
        int thread_id = omp_get_thread_num();
        int bound = IMG_X / n;
        int start = thread_id * bound;
        if(thread_id == 0) {start++;}
        int stop = start + bound;
        if(thread_id == n-1) {stop--;}
       

       for(int i = 0; i < ITERATIONS; i++){
            if(i % 2 == 0){
                iterate_image(petri_A, petri_B, start, stop);
            }else{
                iterate_image(petri_B, petri_A, start, stop);
            }
            #pragma omp barrier
        }
    }

    clock_t end_clock = clock();
    double time_spent = (double)(end_clock-start_clock) / CLOCKS_PER_SEC;
    printf("Time spent on main section: %f\n", time_spent);

    if(ITERATIONS % 2 == 0){
        make_bmp(petri_A, "RPS_omp");
    } else {
        make_bmp(petri_B, "RPS_omp");
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


