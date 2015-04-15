
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include <stdbool.h>

#define MAIN_PROC 0
#define LOG(a, ...) \
	if(rank == MAIN_PROC){ \
		printf(a, ##__VA_ARGS__); \
	}
#define LOG_1(a, ...) \
		printf("[%d / %d: ]", rank, size); \
		printf(a, ##__VA_ARGS__); \
// Global array


// For each process
int rank,size;
int *g_state;
int *g_temp;
int n_row = 20;
int n_col = 20;


int sub_mat_size = 5;// Number of rows in each sub matrix

/************************************/
int *row(int r){
	return g_state + r*n_col;
}

void set_state(int r, int c, int v){
	r++;
	g_state[r*n_col + c] = v;
}

void set_temp(int r, int c, int v){
	r++;
	g_temp[r*n_col + c] = v;
}

int state_at(int r, int c){
	r ++;
	return g_state[r*n_col + c];
}

int temp_at(int r, int c){
	r ++;
	return g_temp[r*n_col + c];
}

void log_sub_mat(){
	printf("I am proc %d / %d\n", rank, size);
	for(int i = 0; i < sub_mat_size + 2; i++){
		for(int j = 0; j < n_col; j++){
			printf("%d ", g_state[i*n_col+j]);
		}
		printf("\n");
	}
}

void init_sub_mat(){
	int _size = (sub_mat_size + 2)*n_col;
	g_state = (int *)malloc(_size * sizeof(int));
	g_temp = (int *)malloc(_size * sizeof(int));
	
	time_t t;
	srand(time(&t));
	for(int i = 0; i < sub_mat_size; i++){
		for(int j = 0; j < n_col; j++){
			if(j == 0 || j == n_col - 1){
				set_state(i,j,0);
			}
			else{
				int v = (rand() > RAND_MAX/2.0)? 0:1;
				set_state(i,j,v);
			}
		}
	}
	
	// Special for first and last
	if(rank == 0){
		for(int i = 0; i < n_col; i++){
			set_state(0, i, 0);
		}
	}
	if(rank == size-1){
		for(int i = 0; i < n_col; i++){
			set_state(sub_mat_size - 1, i, 0);
		}
	}
	
//	log_sub_mat();
}

void syn_share_data_blocking(){
	MPI_Status status;
	if(rank != 0){ // get data from rank-1
		MPI_Send(row(1), n_col, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
		MPI_Recv(row(sub_mat_size), n_col, MPI_INT, rank-1, 0, 
					MPI_COMM_WORLD, &status);
	}
	
	if(rank != size - 1){ // Send - receive lower row
		MPI_Recv(row(0), n_col, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
		MPI_Send(row(sub_mat_size+1), n_col, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	}
}

void gather_to_main(char *name){
	int * big_mat;
	if(rank == MAIN_PROC){
		big_mat = (int*) malloc(n_row*n_col*sizeof(int));
	}
	MPI_Gather(row(1), sub_mat_size*n_col, MPI_INT, big_mat, 
		sub_mat_size*n_col, MPI_INT, MAIN_PROC, MPI_COMM_WORLD);
	if(rank == MAIN_PROC){
		write_to_file(big_mat, n_row, n_col, name);
		free(big_mat);
	}	
}

int get_sum(int r, int c){
	int sum = 0;
	for(int i = r-1; i <= r+1; i++)
		for(int j = c-1; j <= c+1; j++){
			
			if(i >= 0 && i < sub_mat_size
				&& j >= 0 && j < n_col){
					sum += state_at(i,j);
				}
		}
	return sum - state_at(r,c);
	
}

void update_value(){
	for(int i = 0; i < sub_mat_size; i++){
		for(int j = 0; j < n_col; j++){
			int sum_n = get_sum(i,j);
			if(sum_n == 3 || (state_at(i,j) == 1 && sum_n == 2)){
				set_temp(i,j,1);
			}else{
				set_temp(i,j,0);
			}
		}
	}
	
	int *d = g_temp;
	g_temp = g_state;
	g_state = d;
};

int is_unchange_(){
	int unchange = 1;
	for(int i = 0; i < sub_mat_size; i++){
		for(int j = 0; j < n_col; j++){
			if(state_at(i,j) != temp_at(i,j)){
				unchange = 0;
				break;
			}
		}
		if(unchange == 0)
			break;
	}
	
	return unchange;
}

void syn_share_data_non_blocking(){
	//TODO
}

int main(int argc, char* argv[]){
		
	if(argc > 1){
		n_row = atoi(argv[1]);
		n_col = n_row;
	}

	// MPI rank and size	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);	
	
	LOG("Started, nb threads = %d\n", size);
	
	sub_mat_size = n_row / size + 1;
	n_row = size*sub_mat_size;
	init_sub_mat();
	syn_share_data_non_blocking();
	gather_to_main("origin.txt");
	// test
	
	int count = 0;
	while(count++ < 40){
		// In non blocking,
		// 1. Send data
		// 2. Update known elements
		// 3. Recieve data
		// 4. Update the 2 rows
		
		syn_share_data_non_blocking();
		update_value();
		
		// Check if should return here
		int unchange = is_unchange_();
		int total;
		
		
		MPI_Allreduce(&unchange, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if(total == size)
			break;
	}
	
	LOG("Finish after %d iterations\n", count);
	gather_to_main("final.txt");
	MPI_Finalize();	
	
	free(g_state);
	free(g_temp);
	return 0;
}
