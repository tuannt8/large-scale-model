//
//  util.h
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 3/5/15.
//  function in main.c
//
//

#ifndef imageSeg_util_h
#define imageSeg_util_h

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "define.h"
#include "cliio.h"
#include <mpi.h>

#define MAX_LEN_S_T 100
#define MAX_DIM_PROCS_ARRAY 100

#define x_dim 0
#define y_dim 1

#define BOUNDARY_TRANSFER 100
#define GATHER_PHI 101

#define GAP_SIZE 20 // pixel

typedef struct{
    int x;
    int y;
}vec2;


void print_help();

/*
 Local index: (-1, -1) -> (sub_size, sub_size)
 Access the main local domain: (0,0) -> ( sub_size - 1, sub_size -1 )
 */

typedef struct {
    num lamda1;
    num lamda2;
    num mu;
    num nu;
    num dt;
    int max_iter;
    num phi_tol;
    num phi_diff_norm;
}chan_vese_opt;

// Global variables
typedef struct {
    // Process information
    int main_proc;
    int rank, size;
    
    // Data information
    // Image information
    int image_width, image_height;
    
    // Block information
    num* sub_image;
    
    int bl_idx_x, bl_idx_y;
    int bl_dim_x, bl_dim_y;
    int block_size; // block_size x block_size is number of pixels in each process
    int sub_size; // = block_size + 2
    int active_size_x;
    int active_size_y;
    
    // Level set function
    num* phi;
    chan_vese_opt opt;
    
	int img_size;
    
    // Input informtaion
    // Image path
    char file_path[MAX_LEN_S_T];
    
}global_info;

extern global_info g;

/////////////////////////////////////////////////
// Parser argument
int parse_arguments(int argc, char* argv[]);

/////////////////////////////////////////////////
// Read image data
void generate_image();
int read_data();

int log_local_phi();
int log_local_image();

/////////////////////////////////////////////////
// Init level set function
int init_phi();
void gather_phi();
void gather_phi_p(int iter);

/////////////////////////////////////////////////
void chan_vese_loop();
void region_average(num *c1, num *c2);
void update_boundary();
void exchange_boundary();

/////////////////////////////////////////////////
inline int proc_index(int x, int y);

/////////////////////////////////////////////////
// Information of program settings
void print_info();

void print_g();
int max_(int a, int b);
void print_mat(num* data, int width, int height);
int get_global_pixel_index(vec2 const local, vec2 *global);
extern inline num get_sub_image_data(int x, int y);
extern inline void set_sub_image_data(int x, int y, int inten);
extern inline num get_phi_data(int x, int y);
extern inline void set_phi_data(int x, int y, int inten);

// local (-1 -> block) -> index
extern inline int local_array_idx(int x, int y); 
// vec2 global_idx_convert(vec2 local);

int block_idx(int x, int y);

/////////////////////////////////////////////////
// Memory allocation and release
void init_data();
void release_memory();

// Print image for debuging
void debug_print();

#endif
