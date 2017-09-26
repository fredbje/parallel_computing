#ifndef CA_H
#define CA_H

#include "RPS_MPI.h"


void free_img(cell** image);
bool neighborhood_contains(int x, int y, cell** image, int color);
int next_color(int x, int y, cell** image);
void iterate_image(cell** old_image, cell** next_image);
cell** alloc_img(cell* buffer, int index);

#endif
