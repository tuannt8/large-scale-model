#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	int rank;
	int size;
	short tag=0;
	int msg_buf;
	MPI_Status status;

	MPI_Init(&argc,&argv);

        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size (MPI_COMM_WORLD, &size);

	if(rank==0){
		//send my rank to all processes with rank > 0
		printf( "Proc %d sending messages. \n",rank);
		for(int i=1;i<size;i++){
			MPI_Send(&rank,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		}
	}
	else{
		//receive the message and print it
		MPI_Recv(&msg_buf, 1, MPI_INT, 0, tag, MPI_COMM_WORLD,&status);
		printf( "Proc %d receved message from %d \n",rank, msg_buf);
	}

	MPI_Finalize();
        return 0;
}
