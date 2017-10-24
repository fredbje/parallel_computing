#include <stdio.h>
#include <stdlib.h>

int main(void){
	int *w = (int*)malloc(2*sizeof(int));
	w[0] = 300; 
	w[1] = 500;
	printf("Value of *w++: %d\n", *w ++);
	printf("Value of *w: %d\n", *w);
	return 0;
}