#ifndef BITMAP_H
#define BITMAP_H

#include "CA.h"
#include "RPS_MPI.h"

// note that we are extra careful with preprocessor macros. Adding parenthesises is never the
// wrong choice.
#define PIXEL(i,j) ((i)+(j)*IMG_X)
#define LOCAL_PIXEL(i,j) ((i)+(j)*p_local_petri_x_dim)

typedef unsigned char uchar;
void make_bmp(cell* image, int index);

#endif
