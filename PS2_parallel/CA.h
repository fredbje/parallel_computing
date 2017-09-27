#ifndef CA_H
#define CA_H

#include "RPS_MPI.h"

cell** allocate_petri(int row_dim, int col_dim);
void free_petri(cell** petri);
void iterate_CA(cell** old_petri, cell** next_petri, int row_dim, int col_dim);

#endif
