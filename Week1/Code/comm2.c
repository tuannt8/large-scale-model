#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	int rank,size;
	MPI_Status stat;
	double time = MPI_Wtime();


	MPI_Init(&argc,&argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	
	int pre = -1, after = -1;
	if(rank > 0){ // Send and receive from it neighbor
		MPI_Recv(&pre, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &stat); 	
		MPI_Send(&rank, 1, MPI_INT, rank -1, 1, MPI_COMM_WORLD);
	}

	if(rank < size - 1){ // Send and receive from next neighbor
		MPI_Send(&rank, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
		MPI_Recv(&after, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &stat);
	}
	
	int total = pre + after;
	
	if(rank % 10){
		time = MPI_Wtime() - time;
		printf("Proc %d: %d + %d = %d in %f s\n", rank, pre, after, total, time);
	}
	MPI_Finalize();



	return 0;
}
