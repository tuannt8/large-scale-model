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

global_info g;

void init_data(){
    g.main_proc = 0;
    g.sub_image = NULL;
    g.phi = NULL;
    
    // Chan-Vese option
    g.opt.dt = 0.5;
    g.opt.lamda1 = 1;
    g.opt.lamda2 = 1;
    g.opt.max_iter = 15;
    g.opt.mu = 0.5;
    g.opt.nu = 0.0;
    g.opt.phi_tol = 1.0e-5;
    g.opt.phi_diff_norm = 1.;
}

int main(int argc, char* argv[]){

    /*
        Memory allocation and variable initializtion
     */
    init_data();
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank (MPI_COMM_WORLD, &g.rank);
    MPI_Comm_size (MPI_COMM_WORLD, &g.size);

    /* 
        1. Parse arguments
            And print debug information
    */

    parse_arguments(argc, argv);
    print_info();

    /* 
        2. Read image
            Compute optimal block size and number of procs
            Load partial image to each proc
     */ 
    if(!read_data())
        goto catch;
    
    /*
        3. Init Phi
     */
    if(!init_phi())
        goto catch;
    
    /*
        4. Main loop
     */
    
    chan_vese_loop();
    
    
    /*
        Gather level set function and log
     */
    gather_phi();
    LOG_LINE;
    LOG("Done!!!!\n")
    
catch:
    MPI_Finalize();
    release_memory();
    
    return 0;
}
