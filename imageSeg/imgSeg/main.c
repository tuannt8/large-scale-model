//
//  main.c
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 3/5/15.
//
//

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include <stdbool.h>

// Process information
const int main_proc = 0;
int rank, size;

int main(int argc, char* argv[]){
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    
    print_info(); // Some information for debugging
    
    MPI_Finalize();
    return 0;
}
