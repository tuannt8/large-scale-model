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
	int nbMess = 100;
	
	MPI_Init(NULL, NULL);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);		

	if(rank == 0){
		double *mess = calloc(N, sizeof(double));
	
		double total_time;
		for(int i = 1; i < (size-1)*nbMess; i++){
			MPI_Recv(mess, N, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
			double time = MPI_Wtime() - mess[0];
			printf("I got mess from %d - %f\n", stat.MPI_SOURCE, time);
		
			total_time += time;
		}
		free(mess);
		
		printf("%f %f# Speed time\n", N*8.0/1000.0/(total_time/(N-1)), total_time/(N-1));
	}else{
		double *mess = calloc(N, sizeof(double));
		for(int i = 0; i < nbMess; i++){
			mess[i] = MPI_Wtime();
			MPI_Send(mess, N, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
		
		}
		free(mess);			

	}
	
	MPI_Finalize();

	return 0;


}
