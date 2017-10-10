#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stddef.h>
#include "RPS_MPI.h"

void create_cart_comm();
void create_mpi_types();
void init_local_petri();
void exchange_borders(cell** local_petri);
void gather_petri(cell** local_petri);
void free_mpi_types();

int rank;
int world_size;

// The dimensions of the processor grid. Same for every process
int p_row_dims;
int p_col_dims;

// The location of a process in the process grid. Unique for every process
int p_my_row_coord;
int p_my_col_coord;

// Rank of adjacent process in cardinal direction
int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int local_petri_row_dim;
int local_petri_col_dim;
int local_petri_size;

// The dimensions for the global petri
int petri_row_dim = IMG_Y;
int petri_col_dim = IMG_X;
int petri_size = IMG_X*IMG_Y;

// Cartesian communicator
MPI_Comm cart_comm;

// Each process has two local petris to iterate the CA in a lockstep fashion
cell** local_petri_A;
cell** local_petri_B;

// The root process gathers the local petris in petri and writes it to file
cell** petri;

// Petris are arrays of cell_t
MPI_Datatype cell_t; 

// Border exchange types
MPI_Datatype border_row_t;  
MPI_Datatype border_col_t;  

// Used in gather_petri() to gather local petris in global petri in correct order 
MPI_Datatype local_petri_send_t;    
MPI_Datatype local_petri_receive_t;

/////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double start_time = MPI_Wtime();
	
	srand(rank);

	create_cart_comm();
	init_local_petri();
	create_mpi_types();

	for(int i = 0; i < N_ITERATIONS; i++){
		if(i % 2 == 0){
			iterate_CA(local_petri_A, local_petri_B, local_petri_row_dim, local_petri_col_dim);
			exchange_borders(local_petri_B);
		}
		else{
			iterate_CA(local_petri_B, local_petri_A, local_petri_row_dim, local_petri_col_dim);
			exchange_borders(local_petri_A);
		}
	}

	if(rank == 0){
		petri = allocate_petri(petri_row_dim, petri_col_dim);
	}

	if(N_ITERATIONS % 2 == 0){
		gather_petri(local_petri_A);
	}
	else{
		gather_petri(local_petri_B);
	}	

	if(rank == 0){
		printf("Writing petri to file\n");
		make_bmp(petri);
	}

	if(rank == 0){
		free_petri(petri);
	}
	free_petri(local_petri_A);
	free_petri(local_petri_B);
	free_mpi_types();
	
	double end_time = MPI_Wtime();
	printf("Time elapsed for process %d: %f [s]\n", rank, end_time-start_time);

	MPI_Finalize();
  	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void create_cart_comm(){
	int dims[2] = {0, 0}; // Init to 0 so MPI_Dims_create will change it
	int periods[2] = {0, 0}; // Init to 0 for no wrap-around
	int coords[2];

	MPI_Dims_create(world_size, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, 2, coords);

	MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
	MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

	p_row_dims = dims[0];
	p_col_dims = dims[1];

	p_my_row_coord = coords[0];
	p_my_col_coord = coords[1];
}

void create_mpi_types(){
	const int    nitems=2;
	int          blocklengths[2] = {1,1};
	MPI_Aint     offsets[2];
	MPI_Datatype types[2] = {MPI_INT, MPI_INT};

	offsets[0] = offsetof(cell, color);
	offsets[1] = offsetof(cell, strength);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &cell_t);
	MPI_Type_commit(&cell_t);

	MPI_Type_vector(local_petri_col_dim, 1, 1, cell_t, &border_row_t); 
	MPI_Type_commit(&border_row_t);

	MPI_Type_vector(local_petri_row_dim, 1, local_petri_col_dim, cell_t, &border_col_t); 
	MPI_Type_commit(&border_col_t);

	MPI_Type_vector(local_petri_row_dim-2, local_petri_col_dim-2, local_petri_col_dim, cell_t, &local_petri_send_t);
	MPI_Type_commit(&local_petri_send_t);

	MPI_Type_vector(local_petri_row_dim-2, local_petri_col_dim-2, petri_col_dim, cell_t, &local_petri_receive_t);
	MPI_Type_commit(&local_petri_receive_t);
}

void init_local_petri(){
	// Local petris has extra space for borders
	local_petri_row_dim = petri_row_dim/p_row_dims + 2;
	local_petri_col_dim = petri_col_dim/p_col_dims + 2;
	local_petri_size = local_petri_row_dim*local_petri_col_dim;

	local_petri_A = allocate_petri(local_petri_row_dim, local_petri_col_dim);
	local_petri_B = allocate_petri(local_petri_row_dim, local_petri_col_dim);

	// Some cells may be chosen twice - not important
	for(int i = 0; i < N_START_CELLS / world_size; i++){
    	int r_row = rand() % (local_petri_row_dim - 1);
    	int r_col = rand() % (local_petri_col_dim - 1);
    	int r_color = rand() % 4;

    	local_petri_A[r_row][r_col].color = r_color;
    	local_petri_A[r_row][r_col].strength = 1;
  	}
}

void exchange_borders(cell** local_petri){
				//sendbuf														dest			receivebuf														source			
	MPI_Sendrecv(&local_petri[1][0],						1, border_row_t, 	p_north, 	0, 	&local_petri[local_petri_row_dim-1][0], 	1, border_row_t,	p_south, 	0, cart_comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&local_petri[local_petri_row_dim-2][0],	1, border_row_t, 	p_south, 	0, 	&local_petri[0][0], 						1, border_row_t, 	p_north, 	0, cart_comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&local_petri[0][local_petri_col_dim-2], 	1, border_col_t, 	p_east, 	0, 	&local_petri[0][0], 						1, border_col_t, 	p_west, 	0, cart_comm, MPI_STATUS_IGNORE);
	MPI_Sendrecv(&local_petri[0][1], 						1, border_col_t, 	p_west, 	0, 	&local_petri[0][local_petri_col_dim-1], 	1, border_col_t, 	p_east, 	0, cart_comm, MPI_STATUS_IGNORE);
}

void gather_petri(cell** local_petri){
	MPI_Request send_request, recv_request;
	MPI_Status status;
	// Don't include borders when sending local petri to global petri
	MPI_Isend(&local_petri[1][1], 1, local_petri_send_t, 0, 0, MPI_COMM_WORLD, &send_request);
	int start_row, start_col;
  	if(rank == 0){   		
  		for(int i = 0; i < world_size; i++){
  			// Global petri coordinate of first element in local petri 
  			start_row = (i/p_col_dims)*(local_petri_row_dim-2);
  			start_col = (i%p_col_dims)*(local_petri_col_dim-2);

  			// Utilize MPI types to arrange local petri correctly in global petri
  			MPI_Irecv(&petri[start_row][start_col], 1, local_petri_receive_t, i, 0, MPI_COMM_WORLD, &recv_request);
  			MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
  		}
  	}
  	MPI_Wait(&send_request, MPI_STATUS_IGNORE);
}

void free_mpi_types(){
	MPI_Type_free(&border_row_t);
	MPI_Type_free(&border_col_t);
	MPI_Type_free(&local_petri_receive_t);
	MPI_Type_free(&local_petri_send_t);
	MPI_Type_free(&cell_t);
	MPI_Comm_free(&cart_comm);
}

