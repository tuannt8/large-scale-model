#pragma once

#define VAR(a) #a
#define PRINT_MAT(a,b,c) print_mat(a, b, c, #a)

void print_mat(double *mat, int row, int col, char *name){
	printf("%s = \n", name);
	for(int i = 0; i < row; i++){
		for(int j = 0; j < col; j++){
			printf("%f\t", mat[i*col + j]);
		}
		printf("\n");
	}
}
