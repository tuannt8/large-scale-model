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
	MPI_Status stat;
	double pi_;
	
	double time = MPI_Wtime();
	
	double h = 1.0/(double)N;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	
	double sum = 0.0;
	int block = ceil((double)N/size);
	
	int i;
	for(i = 0; i < block; i++){
		double x = (rank*block + i - 0.5)*h;
		sum += 4.0/(1 + x*x);
	}
	
	double pp;
	if(rank == 0){
		for(i = 1; i < size; i++){
			MPI_Recv(&pp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stat); 
			sum += pp;
		}
		pi_ = sum*h;
		printf("pi = %f\n", pi_);		
	}
	else{
		MPI_Send(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();

	return 0;


}
