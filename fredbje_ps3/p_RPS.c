#include "RPS.h"
#include <time.h>
#include <pthread.h>
#include <math.h>

void* entry_function(void *threadid);
void init_petri();

cell* petri_A;
cell* petri_B;

int num_threads;
int* thread_ids;
pthread_t* threads;
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int count = 0;

int main(int argc, char** argv){

    if(argc==1) {
        puts("Input number of threads. ex: ./RPS_pthread 4");
        return 0;
    }
    num_threads = (int)strtod(argv[1], NULL);
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

    clock_t start_clock = clock();

    thread_ids = (int*)malloc(num_threads*sizeof(int));
    threads = (pthread_t*)malloc(num_threads*sizeof(pthread_t));
    int err;
    for(int i = 0; i < num_threads; i++){
        thread_ids[i] = i;
        err = pthread_create(&threads[i], NULL, entry_function, (void*)&thread_ids[i]);
        if (err){
            printf("ERROR; return code from pthread_create() is %d\n", err);
            exit(-1);
        }
    }

    for (int i=0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_t end_clock = clock();
    double time_spent = (double)(end_clock-start_clock) / CLOCKS_PER_SEC;
    printf("Time spent on main section: %f\n", time_spent);

    if(ITERATIONS % 2 == 0){
        make_bmp(petri_A, "RPS_pthread");
    } else {
        make_bmp(petri_B, "RPS_pthread");
    }

    free(threads);
    free(thread_ids);
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

void* entry_function(void *thread_id){

    int bound = IMG_Y / num_threads;
    int start = *((int*)thread_id) * bound;
    int finish = start + bound;
   
   for(int i = 0; i < ITERATIONS; i++){
        if(i % 2 == 0){
            iterate_image(petri_A, petri_B, start, finish);
            pthread_mutex_lock(&mut);
            count++;
            if(count == num_threads){
                count = 0;
                for(int i = 0; i < num_threads; i++){pthread_cond_signal(&cond);}
            }else{
                pthread_cond_wait(&cond, &mut);
            }
        }else{
            iterate_image(petri_B, petri_A, start, finish);
            pthread_mutex_lock(&mut);
            count++;
            if(count == num_threads){
                count = 0;
                for(int i = 0; i < num_threads; i++){pthread_cond_signal(&cond);}
            }else{
                pthread_cond_wait(&cond, &mut);
            }
        }
        pthread_mutex_unlock(&mut);
    }
    
   return NULL;
}






