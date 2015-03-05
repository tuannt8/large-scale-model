
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include <stdbool.h>

#define MAIN_PROC 0

// Image information
int n_row, n_col;
unsigned int * image = NULL; // 0-255

// Level set function
// Consider level set size same as image size
double *Phi0, *Phi;


int main(int argc, char* argv[]){
    
    // 1. Read input image
    
    image = read_my_image("../Data/square.my" , &n_row, &n_col);
    if (!image) {
        printf("Read image failed\n");
        return 1;
    }
    
    // 2. Initialize level set functions
    Phi0 = malloc(n_row * n_col * sizeof(double));
    Phi = malloc(n_row * n_col * sizeof(double));
    
    
    
    // Free memory
    free(image);
    free(Phi);
    free(Phi0);
    
    return 0;
}
