#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"

// #define LOG 
#define PROFILE

int main(int argc, char *argv[]){

	// MPI rank and size
	int rank,size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	
	int root = 0;
	
	// 1. Read input size
	int N;
	double alpha;
	if(rank == root){//Read size from command line
		N = atoi(argv[1]);
		alpha = atof(argv[2]);
	}
	MPI_Bcast(&N, 1, MPI_INT, root, MPI_COMM_WORLD);	
	MPI_Bcast(&alpha, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	// 2. Initialize array
	double *r, *v, *w;
	int vec_size = size * N;
	r = (double*)malloc(N*sizeof(double));
	v = (double*)malloc(N*sizeof(double));
	w = (double*)malloc(N*sizeof(double));
	for(int i = 0; i < N; i++){
		int index = rank*N + i;
		v[i] = index;
		w[i] = index + 1;
	}

#ifdef PROFILE
	double t;
	if(rank == root)
		t = MPI_Wtime();
#endif
	
	// 3. Sum r = v + alpha.w
	for(int i = 0; i < N; i++){
		r[i] = v[i] + alpha * w[i];
	}
	// Send to process root
	double *result;
	if(rank == root){
		result = (double*)malloc(N*size*sizeof(double));
	}
	MPI_Gather(r, N, MPI_DOUBLE, result, N, MPI_DOUBLE, root, MPI_COMM_WORLD);

#ifdef PROFILE
	if(rank == root){
		t = MPI_Wtime() - t;
		printf("%f ", t);
	}
#endif

#ifdef LOG
	if(rank == root){
		PRINT_MAT(result, 1, size*N);
	}
#endif
	
	// 4. product s = v(t)*w
#ifdef PROFILE
	if(rank == root)
		t = MPI_Wtime();
#endif

	double sum = 0.0;
	for(int i = 0; i < N; i++){
		sum += v[i] * w[i];
	}
	// Reduction
	double s = 0.0;
	MPI_Reduce(&sum, &s, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	
#ifdef PROFILE
	if(rank == root){
		t = MPI_Wtime() - t;
		printf("%f ", t);
	}
#endif	
	
#ifdef PROFILE
	if(rank == root){
		printf("# sum; product \n");
	}
#endif

	if(rank == root){
		printf("s = %f\n", s);
	}
	
	
	MPI_Finalize();		
	
	free(r);
	free(v);
	free(w);
	if(rank == root)
		free(result);

	return 0;
}
