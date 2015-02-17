#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int main(int argc, char* argv[])
{
	int N = 100000000;
	
	if(argc == 2)
		N = atoi(argv[1]);
	
	//double pi = pi_mpi(N);
	int rank,size;
	double pi_;
	
	double time = MPI_Wtime();
	double h = 1.0/(double)N;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	
	// There are [size] block, each block sum [block] elements
	double sum = 0.0;
	int block = ceil((double)N/size);
	
	// 1. Compute the block

	
	// 2. Send and receive the result
	if(rank == 0){
		MPI_Request request[size - 1];
		double pp[size-1];
		for(int i = 1; i < size; i++){	
			MPI_Irecv(&pp[i-1], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &request[i-1]);
		}
		for(int i = 0; i < block; i++){
			double x = (rank*block + i - 0.5)*h;
			sum += 4.0/(1 + x*x);
		}
	
		MPI_Status stat[size - 1];
		MPI_Waitall(size-1, request, stat);
		
		for(int i = 0; i < size - 1; i ++){
			sum += pp[i];
		}
		
		pi_ = sum*h;
		printf("pi = %f\n", pi_);		
	}
	else{
		for(int i = 0; i < block; i++){
			double x = (rank*block + i - 0.5)*h;
			sum += 4.0/(1 + x*x);
		}
		// All send to proc 0
		MPI_Request request;
		MPI_Isend(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
	}
	
	MPI_Finalize();

	return 0;
}
