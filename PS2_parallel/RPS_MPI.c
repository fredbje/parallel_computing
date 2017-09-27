#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stddef.h>
#include "RPS_MPI.h"

void initialize();
void iterate_CA();
void gather_petri(cell** local_petri);
void create_types();
void initialize_petri();
void exchange_borders(cell** local_petri);

int rank;
int size;

// The dimensions of the processor grid. Same for every process
int p_row_dims;
int p_col_dims;

// The location of a process in the process grid. Unique for every process
int p_my_row_coord;
int p_my_col_coord;

int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int local_petri_row_dim;
int local_petri_col_dim;
int local_petri_size;
//#define LOCAL_PIXEL(i,j) ((i)+(j)*local_petri_row_dim)
int local_petri_exp_row_dim;
int local_petri_exp_col_dim;
int local_petri_exp_size;
//#define LOCAL_EXP_PIXEL(i,j) ((i)+(j)*local_petri_exp_row_dim)

// The global petri
cell** petri;
int petri_row_dim = IMG_Y;
int petri_col_dim = IMG_X;
int petri_size = IMG_X*IMG_Y;

MPI_Comm cart_comm;

MPI_Datatype border_row_t;  
MPI_Datatype border_col_t;  
MPI_Datatype local_petri_t; 
MPI_Datatype mpi_cell_t;    
MPI_Datatype receive_t;

// Each process has two local petris to iterate the CA in a lockstep fashion
cell** local_petri_A;
cell** local_petri_B;

int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	srand(rank);

	int dims[2];
	dims[0] = p_row_dims;
	dims[1] = p_col_dims;

	int periods[2]; // 0 for no wrap-around
	periods[0] = 0;
	periods[1] = 0;

	int coords[2];
	coords[0] = p_my_row_coord;
	coords[1] = p_my_col_coord;

	MPI_Dims_create(size, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, 2, coords);

	MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
	MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

	p_row_dims = dims[0];
	p_col_dims = dims[1];

	p_my_row_coord = coords[0];
	p_my_col_coord = coords[1];

// printf("rank: %d p_my_row_coord: %d p_my_col_coord: %d p_row_dims: %d p_col_dims: %d\n", rank, p_my_row_coord, p_my_col_coord, p_row_dims, p_col_dims);
	// printf("rank: %d p_north: %d p_south: %d p_west: %d p_east: %d \n", rank, p_north, p_south, p_west, p_east);
	////////////////////////////////
	////////////////////////////////

	initialize();
	create_types();

	for(int i = 0; i < N_ITERATIONS; i++){
					//old			new
		iterate_CA(local_petri_A, local_petri_B, local_petri_row_dim, local_petri_col_dim);
		MPI_Barrier(MPI_COMM_WORLD);
		exchange_borders(local_petri_B);
		MPI_Barrier(MPI_COMM_WORLD);
		iterate_CA(local_petri_B, local_petri_A, local_petri_row_dim, local_petri_col_dim);
		MPI_Barrier(MPI_COMM_WORLD);
		exchange_borders(local_petri_A);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(rank==0){
		petri = allocate_petri(petri_row_dim, petri_col_dim);
	}

	gather_petri(local_petri_A);
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){
		printf("Writing to file\n");
		make_bmp(petri);
	}
	printf("Rank %d about to free stuff\n", rank);

	// free_stuff()
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		free_petri(petri);
	}
	free_petri(local_petri_A);
	free_petri(local_petri_B);

	MPI_Finalize();
  	return 0;
}


void create_types(){
	const int    nitems=2;
	int          blocklengths[2] = {1,1};
	MPI_Aint     offsets[2];
	MPI_Datatype types[2] = {MPI_INT, MPI_INT};

	offsets[0] = offsetof(cell, color);
	offsets[1] = offsetof(cell, strength);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
	MPI_Type_commit(&mpi_cell_t);

	MPI_Type_vector(local_petri_col_dim, 1, 1, mpi_cell_t, &border_row_t); 
	MPI_Type_commit(&border_row_t);

	MPI_Type_vector(local_petri_row_dim, 1, local_petri_col_dim, mpi_cell_t, &border_col_t); 
	MPI_Type_commit(&border_col_t);

	MPI_Type_vector(local_petri_row_dim-2, local_petri_col_dim-2, local_petri_col_dim, mpi_cell_t, &local_petri_t);
	MPI_Type_commit(&local_petri_t);

	MPI_Type_vector(local_petri_row_dim-2, local_petri_col_dim-2, petri_col_dim, mpi_cell_t, &receive_t);
	MPI_Type_commit(&receive_t);
}


void initialize(){
	local_petri_row_dim = petri_row_dim/p_row_dims+2;
	local_petri_col_dim = petri_col_dim/p_col_dims+2;
	local_petri_size = local_petri_row_dim*local_petri_col_dim;

	// Allocating extra space for borders
	local_petri_A = allocate_petri(local_petri_row_dim, local_petri_col_dim);
	local_petri_B = allocate_petri(local_petri_row_dim, local_petri_col_dim);

	 for(int ii = 0; ii < N_START_CELLS / size; ii++){
    	int r_row = rand() % (local_petri_row_dim - 1);
    	int r_col = rand() % (local_petri_col_dim - 1);
    	int r_color = rand() % 4;

    	local_petri_A[r_row][r_col].color = r_color;
    	local_petri_A[r_row][r_col].strength = 1;
  	}

	// for(int r = 0; r < local_petri_row_dim; r++){
	// 	for(int c = 0; c < local_petri_col_dim; c++){
	// 		if(rank == 0){
	// 			local_petri_A[r][c].color = ROCK; // Red
	// 		}else if(rank == 1){
	// 			local_petri_A[r][c].color = PAPER; // Blue
	// 		}else if(rank == 2){
	// 			local_petri_A[r][c].color = SCISSOR; // Green
	// 		}
	// 	}
	// }	
}


void exchange_borders(cell** local_petri){
				//sendbuf														dest		TAG	receivebuf														source		TAG	
	MPI_Sendrecv(&local_petri[1][0],						1, border_row_t, 	p_north, 	0, 	&local_petri[local_petri_row_dim-1][0], 	1, border_row_t,	p_south, 	0, cart_comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&local_petri[local_petri_row_dim-2][0],	1, border_row_t, 	p_south, 	0, 	&local_petri[0][0], 						1, border_row_t, 	p_north, 	0, cart_comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&local_petri[0][local_petri_col_dim-2], 	1, border_col_t, 	p_east, 	0, 	&local_petri[0][0], 						1, border_col_t, 	p_west, 	0, cart_comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&local_petri[0][1], 						1, border_col_t, 	p_west, 	0, 	&local_petri[0][local_petri_col_dim-1], 	1, border_col_t, 	p_east, 	0, cart_comm, MPI_STATUS_IGNORE);
}

void gather_petri(cell** local_petri){
	MPI_Request request;
	MPI_Status status;
	MPI_Isend(&local_petri[1][1], 1, local_petri_t, 0, 0, MPI_COMM_WORLD, &request);
  	if(rank == 0){   		
  		for(int i = 0; i < size; i++){
  			MPI_Irecv(&petri[0][(i/p_col_dims)*petri_col_dim*(local_petri_row_dim-2) + (i%p_col_dims)*(local_petri_col_dim-2)], 
  					1, receive_t, i, 0, MPI_COMM_WORLD, &request);
  		}
  	}
  	MPI_Wait(&request, &status);
	// TODO find out what test and wait does. Free memory.
}








