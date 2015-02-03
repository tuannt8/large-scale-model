#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	int rank,size;



	MPI_Init(&argc,&argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	

	double time = MPI_Wtime();

	if(rank == 0){ // Send to 1
		int mess = 10;
		MPI_Status status;
		MPI_Send(&mess, 1, MPI_INT, 1, 1, MPI_COMM_WORLD); // Send to 1
		printf("Process 0 sent the message [%d]\n", mess);

		MPI_Recv(&mess, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &status);
		printf("Process 0 received final message: %d\n", mess);
	}else if (rank < size - 1){ // Send to next process
		int mess;
		MPI_Status status;
		MPI_Recv(&mess, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
		mess ++;
		printf("Process %d increased the message to %d\n", rank, mess);
		MPI_Send(&mess, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD); // Send to rank + 1
	}else if(rank == size - 1){ // Last process, send to 0
		int mess;
		MPI_Status status;
		MPI_Recv(&mess, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
		mess++;
		printf("Process %d increased the message to %d\n", rank, mess);

		printf("Process %d sent message to process %d\n", rank, 0);
		MPI_Send(&mess, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); // Send to 0
	}
	time = MPI_Wtime() - time;
	printf("Process %d Total time: %f s\n", rank, time);


	MPI_Finalize();



	return 0;
}
