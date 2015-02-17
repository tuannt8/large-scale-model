
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"

#define LOG
#define PROFILE
#define VALIDATE

#define PRINT(a, ...) 	\
	if(rank == root){ 	\
		printf(a, ##__VA_ARGS__);	\
	} 					\

int main(int argc, char *argv[]){
	int N = atoi(argv[1]);

	// MPI rank and size
	int root = 0;
	int rank,size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	
	// 1. Initialize
	int block = ceil((double)N/size);
	PRINT("N = %d; block = %d\n", N, block);
	double *A_sub = (double*)malloc(block*N*sizeof(double));
	double *x_sub = (double*)malloc(N*sizeof(double));
	double *y_sub = (double*)malloc(block*sizeof(double));
	for(int i = 0; i < block; i++){
		y_sub[i] = 1.0;
	}
	for(int j = 0; j < N; j++){
		for(int i = 0; i < block; i++){
			A_sub[j*block + i] = 1.0;
		}
	}
	
	// 2. Multiply
#ifdef PROFILE
	double t;
	if(rank == root){
		t = MPI_Wtime();
	}
#endif

	
	if(rank == root){
		double **buf = (double**)malloc(size*sizeof(double*));
		MPI_Request request[size-1];
		for(int i = 1; i < size; i++){	
			buf[i-1] = (double*)malloc(N*sizeof(double));
			MPI_Irecv(buf[i-1], N, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &request[i-1]);
		}		
		
		for(int j = 0; j < N; j++){
			double sum = 0.0;
			for(int i = 0; i < block; i++){
				if(rank*block + i >= N)
					break;
				
				sum += A_sub[j*block + i] * y_sub[i];
			}
			x_sub[j] = sum;
		}
		
		MPI_Status stat[size - 1];
		MPI_Waitall(size-1, request, stat);
		
		for(int i = 0; i < N; i++){
			for(int j = 0; j < size - 1; j++){
				x_sub[i] += buf[j][i];
			}
		}		
	}
	else{
		for(int j = 0; j < N; j++){
			double sum = 0.0;
			for(int i = 0; i < block; i++){
				if(rank*block + i >= N)
					break;
				
				sum += A_sub[j*block + i] * y_sub[i];
			}
			x_sub[j] = sum;
		}
		MPI_Request request;
		MPI_Isend(x_sub, N, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &request);
	}

#ifdef PROFILE
	if(rank == root){
		t = MPI_Wtime() - t;
		printf("Walltime: %f \n", t);
	}
#endif	

#ifdef VALIDATE
	if(rank == root){
		double err = 0.0;
		for(int i = 0; i < N; i++){
			err += N - x_sub[i];
		}
		printf("Error: %f\n", err);
	}
#endif

	MPI_Finalize();	
		
	free(A_sub);
	free(x_sub);
	free(y_sub);
	
	return 0;
}
