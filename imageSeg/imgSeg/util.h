//
//  util.h
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 3/5/15.
//
//

#ifndef imageSeg_util_h
#define imageSeg_util_h

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "define.h"



void write_to_file(int *mat, int nrow, int ncol, char *file_name){
    FILE *f = fopen(file_name, "w");
    if(!f)
        printf("File %s can not be open\n", file_name);
    
    // Write to file
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            fprintf(f, "%d ", mat[i*ncol + j]);
        }
        fprintf(f, "\n");
    }
    
    fclose(f);
}

unsigned int * read_my_image(const char * file_name, int * n_row, int * n_col){
    FILE *f = fopen(file_name, "r");
    if (f) {
        
        fscanf(f, "%d %d\n", n_row, n_col);
        assert(*n_row > 0 && *n_col > 0, "image file error");
        unsigned int* img = malloc( (*n_row) * (*n_col) * sizeof(unsigned int));
        
        for (int i = 0; i < *n_row; i++) {
            for (int j = 0; j < *n_col; j++) {
                fscanf(f, "%d", img + i* (*n_col) + j);
            }
            fscanf(f, "\n");
        }
        
        fclose(f);
        return img;
    }
    
    return NULL;
}


#endif
