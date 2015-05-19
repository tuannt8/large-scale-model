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
MPI_Comm grid_comm;

void init_data(){
    g.main_proc = 0;
    g.sub_image = NULL;
    g.phi = NULL;
    
    g.img_size = 0;
    
    // Chan-Vese option
    g.opt.dt = 0.5;
    g.opt.lamda1 = 1;
    g.opt.lamda2 = 1;
    g.opt.max_iter = 10;
    g.opt.mu = 0.5;
    g.opt.nu = 0.0;
    g.opt.phi_tol = 1.0e-5;
    g.opt.phi_diff_norm = 1.;
}

/*
image_seg
	file_path 
*/

// #define TIME_LOG

int main(int argc, char* argv[]){

    /*
        Memory allocation and variable initializtion
     */
    init_data(); // In main.c
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank (MPI_COMM_WORLD, &g.rank);
    MPI_Comm_size (MPI_COMM_WORLD, &g.size);

    // 1. Parse arguments
    //    And print debug information
    if(parse_arguments(argc, argv) != 0){
    	print_help();
    	goto catch;
    }
    print_info();
    
    if(g.img_size != 0){
    	strcpy(g.file_path, "LOG/dummy.bmp");
//    	generate_image();
    }
    
    num other_t = 0.0, loop_t = 0.0;
    

    // 2. Read image
    //    Compute optimal block size and number of procs
    //    Load partial image to each proc
//    MPI_Barrier(MPI_COMM_WORLD);
    double t = MPI_Wtime();
    
    if(read_data() != 0)
        goto catch;
    
    MPI_Barrier(MPI_COMM_WORLD);    
   	t = MPI_Wtime() - t;
   	
   	LOG("Read image time: %f sec \n", t);
   	other_t += t;
   	
   	// Log local image
   	// log_local_image();
 
    // 3. Init Phi
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    
    if(init_phi() != 0)
        goto catch;
        
  //  log_local_phi();
    
    MPI_Barrier(MPI_COMM_WORLD);    
    t = MPI_Wtime() - t;
   	LOG("Init level set time: %f sec \n", t);
   	other_t += t;
    
    // 4. Main loop
 
    LOG("Start Chan-Vese\n");
 
   
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    int count;
    if(grid_comm != MPI_COMM_NULL)
    {
		count = chan_vese_loop();
		
		MPI_Barrier(MPI_COMM_WORLD);
		t = MPI_Wtime() - t;
	   	LOG("Main loop: %f sec \n", t);
	   	loop_t += t;
    }
    // Gather level set function and log  

//    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    LOG("Start gathering phi\n");  
        
    gather_phi();
    LOG_LINE;
    
//    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
   	LOG("gathering level set time: %f sec \n", t);
   	other_t += t;
   	
   	// Output
   	if(g.rank == g.main_proc){
   		num mem = 2 
   					* g.bl_dim_x*g.bl_dim_y 
   					* g.sub_size * g.sub_size
   					* sizeof(num);
   		printf("%f %f %f %d %d # Mem other_time segment_time num_thread iters\n", 
   				mem, other_t, loop_t, g.bl_dim_x*g.bl_dim_y , count-1);
   	}
   	
  
    MPI_Finalize();
    release_memory();
    
    LOG("Done!!!!\n")
    
    return 0;
    
catch:
	LOG("Some errors occur \n");
    MPI_Finalize();
    release_memory();
    
    return 1;
}
