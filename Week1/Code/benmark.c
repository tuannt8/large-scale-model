#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EXIT(msg) {printf("%s", msg); \
					return 1;}

int main(int argc, char* argv[])
{
	long long N = 1000;
	if(argc == 2)
		N = atoi(argv[1]);
		
	if(N < 0)
		EXIT("Input failed\n");
	 
	int rank,size;
	MPI_Status stat;
	
	int loop = 100;
	
	MPI_Init(NULL, NULL);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	


	if(rank == 0){
		for(int i = 0; i < loop; i++){
			double *mess = (double*)calloc(N, sizeof(double));
			double time = MPI_Wtime();
			mess[0] = time;
		
			MPI_Send(mess, N, MPI_DOUBLE, 1, i, MPI_COMM_WORLD);
			free(mess);
		}
	}else if(rank == 1){
		double speed = 0.0;
		for(int i = 0; i < loop; i++){
			double *mess = (double*)calloc(N, sizeof(double));
			MPI_Recv(mess, N, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &stat);
			double time = MPI_Wtime() - mess[0];
			speed += (double)N*sizeof(double)/time/1000.0;
			free(mess);
		}
		
		printf("%d %f\n", rank, speed/loop);
	}

	
	MPI_Finalize();

	return 0;


}
