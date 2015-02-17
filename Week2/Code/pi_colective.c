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
	
	double sum = 0.0;
	int block = ceil((double)N/size);
		
	for(int i = 0; i < block; i++){
		double x = (rank*block + i - 0.5)*h;
		sum += 4.0/(1 + x*x);
	}
	
	double *buf;
	if(rank == 0){
		buf = (double*)malloc(size*sizeof(double));
		for(int i = 0; i < size; i++){
			pi_ += buf[i];
			printf("Process %d: %f\n", i, buf[i]);			
		}
	}
	MPI_Gather(&sum, 1, MPI_DOUBLE, buf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	
	if(rank == 0){
		pi_ = 0.0;
		for(int i = 0; i < size; i++){
			pi_ += buf[i];
			printf("Process %d: %f\n", i, buf[i]);			
		}
	
		pi_ = pi_*h;
		printf("pi = %f\n", pi_);	
	}
	
	MPI_Finalize();

	return 0;
}
