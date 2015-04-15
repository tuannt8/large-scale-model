#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>
#include <time.h>
#include <stdlib.h>


void init_state(int *state, int row, int col){
	// Dead boundary
	// And random inside
	time_t t;
	srand(time(&t));
	for(int i = 0; i < row; i++){
		for(int j = 0; j < col; j++){
			if(i == 0 || i == row - 1
			   || j == 0 || j == col - 1){
				state[i*col + j] = 0;
			}
			else{
				
				state[i*col + j] = (rand() > RAND_MAX/2.0)? 0:1;
			}
		}
	}
};

void write_to_file(int *mat, int nrow, int ncol, char *file_name){
	FILE *f = fopen(file_name, "w");
	if(!f)
		printf("File %s can not be open\n", file_name);
		
	// Write to file	
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			fprintf(f, "%d ", mat[i*ncol + j]);
		}
		fprintf(f, "\n");
	}
		
	fclose(f);
}


#endif
