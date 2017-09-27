#include "CA.h"

cell pick_neighbor(int row, int col, cell** petri);
bool neighborhood_contains(int row, int col, cell** petri, int color);
int next_color(int row, int col, cell** petri);

cell** allocate_petri(int row_dim, int col_dim){
    cell** new_petri = (cell**)calloc(row_dim, sizeof(cell**));
    new_petri[0] = (cell*)calloc(row_dim*col_dim, sizeof(cell*));
    for(int r = 0; r < row_dim; r++){
        new_petri[r] = &new_petri[0][r * col_dim];
    }
    return new_petri;
}
  
void free_petri(cell** petri){
    free(petri[0]);
    free(petri);
}

// void swap_petris(cell*** petri_A, cell*** petri_B){
//     cell** temp = *petri_A;
//     *petri_A = *petri_B;
//     *petri_B = temp;
// }

bool beats(cell me, cell other){
  return
    (((me.color == SCISSOR) && (other.color == PAPER)) ||
     ((me.color == PAPER) && (other.color == ROCK))    ||
     ((me.color == ROCK) && (other.color == SCISSOR))  ||
     (me.color == other.color));
}

cell next_cell(int row, int col, cell** petri){

  	cell neighbor_cell = pick_neighbor(row, col, petri);
  	cell my_cell = petri[row][col];
  	if(neighbor_cell.color == WHITE){
    	return my_cell;
  	}
  	cell next_cell = my_cell;

  	if(my_cell.color == WHITE){
    	next_cell.strength = 1;
    	next_cell.color = neighbor_cell.color;
    	return next_cell;
  	}
  	else {
    	if(beats(my_cell, neighbor_cell)){
      		next_cell.strength++;
    	}
    	else{
      	next_cell.strength--;
    	}
  	}

	if(next_cell.strength == 0){
        next_cell.color = neighbor_cell.color;
   	    next_cell.strength = 1;
    }

  	if(next_cell.strength > 4){
        next_cell.strength = 4;
  	}

  	return next_cell;
}

cell pick_neighbor(int row, int col, cell** petri){
  	int chosen = rand() % 8;

  	if(chosen == 4){ chosen++; } // a cell cant choose itself
  	int delta_row = chosen / 3;
  	int delta_col = chosen % 3;

  	return petri[row + delta_row - 1][col + delta_col - 1];
}

void iterate_CA(cell** old_petri, cell** next_petri, int row_dim, int col_dim){
  for(int r = 1; r < row_dim - 1; r++){
    for(int c = 1; c < col_dim - 1; c++){
      // printf("%d %d\n", xx, yy);
      next_petri[r][c] = next_cell(r, c, old_petri);
    }
  }
}
