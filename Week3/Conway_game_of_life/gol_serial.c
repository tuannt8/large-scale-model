#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int row, col;
int *state = NULL;
int *temp = NULL;
int *sum;

#define dis_type 0

void init_state(){
	// Dead boundary
	// And random inside
	time_t t;
	srand(time(&t));
	for(int i = 0; i < row; i++){
		for(int j = 0; j < col; j++){
			if(i == 0 || i == row - 1
			   || j == 0 || j == col - 1){
				state[i*col + j] = 0;
			}
			else{
				
				state[i*col + j] = (rand() > RAND_MAX/2.0)? 0:1;
			}
		}
	}
};

void display(int *a){

	printf("------------------------------\n");
	for(int i = 0; i < row; i++){
	printf("|");
		for(int j = 0; j < col; j++){
			if(dis_type == 0){
			if(state[i*col + j] == 0){
				printf(" ");
			}else{
				printf("*");
			}
			}else if(dis_type == 1)
				printf("%d", a[i*col + j]);
		}
		printf("|\n");
	}
	printf("------------------------------\n");
};

int neighbor_sum(int r, int c){
	int sum = 0;
	
	for(int i = r-1; i <= r+1; i++)
		for(int j = c-1; j <= c+1; j++){
			
			if(i >= 0 && i < row
				&& j >= 0 && j < col){
					sum += state[i*col + j];
				}
		}
	return sum - state[r*col + c];
}

void update_state(){
	for(int i = 0; i < row; i++){
		for(int j = 0; j < col; j++){
			int sum_n = neighbor_sum(i,j);
			sum[i*col + j] = sum_n;
			if(sum_n == 3 || (state[i*col + j] == 1 && sum_n == 2)){
				temp[i*col + j] = 1;
			}else{
				temp[i*col + j] = 0;
			}
		}
	}
	
	int *d = temp;
	temp = state;
	state = d;	
}

void erase_old(){
	for(int i = 0; i < 5; i++){
		printf("\n");
	}
}

void test(){
	row = 20;
	col = 20;
	state = (int*)malloc(row*col*sizeof(int));
	temp = (int*)malloc(row*col*sizeof(int));
	sum = (int*)malloc(row*col*sizeof(int));
	
	init_state();
	
	display(state);
	
	char c;
	while(1){
		c = getchar();
		if(c=='\n'){
			update_state();

			erase_old();
			display(state);
		}
		else if(c == 'q'){
			break;
		}
	}	
	
	free(state);
	free(temp);	
}

void test_final(char * file_name){

}


int main(){
	test();
	return 0;
}
