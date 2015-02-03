#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	int rank,size;



	MPI_Init(&argc,&argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	

	if(rank == 0){ // Send to 1
		int mess = 10;
		MPI_Status status;
		MPI_Send(&mess, 1, MPI_INT, 1, 1, MPI_COMM_WORLD); // Send to 1
		printf("Process 0 sent the message [%d]\n", mess);

		MPI_Recv(&mess, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &status);
		printf("Process 0 received final message: %d\n", mess);
	}

	int pre = 0, after = 0;
	if(rank > 0){ // Send and receive from it neighbor
		MPI_Send(&rank, 1, MPI_INT, rank -1, 1, MPI_COMM_WORLD);
		MPI_Recv(&pre, 1, MPI_INT, rank - 1, 1, MPI_COMM_WOLRLD, MPI_Status);
	}

	if(rank < size - 1){ // Send and receive from next neighbor
		MPI_Send(&rank, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
		MPI_Recv(&pre, 1, MPI_INT, rank + 1, 1, MPI_COMM_WOLRLD, MPI_Status);
	}
	
	int total = pre + after;
	printf("Proc %d: %d + %d = %d", rank, pre, after, total);
	
	MPI_Finalize();



	return 0;
}
