//
//  util.c
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 3/24/15.
//
//

#include "util.h"



#ifdef USE_WINDOWS
num * buffer; 
num * receive_buffer;
MPI_Win win;
// Buffer for other process to access
// Left - down - right - top
#define LEFT 0*g.block_size
#define BOTTOM g.block_size
#define RIGHT 2*g.block_size
#define TOP 3*g.block_size

#else
num *left, *right, *top, *bottom;
num *right_recv, *left_recv;
#endif

/////////////////////////////////////////////////
// Read image data
int read_data(){
    /*
    //
     1. Find optimal block size
     */
//    LOG("Start reading image data \n");
    
    if(g.rank == g.main_proc){
 
		int width, height;
		if(!ReadImageSize_BMP(&width, &height, g.file_path)){
		    LOG("==== Failed to read header BMP file");
		    goto Catch;
		} 
		
		// Find optimal block size
		int nb_proc_aray_x[MAX_DIM_PROCS_ARRAY] = {1};
		int nb_proc_aray_y[MAX_DIM_PROCS_ARRAY] = {1};
		int block_size[MAX_DIM_PROCS_ARRAY] = {max_(width, height)};
		
		int length = 1;
		while (length < MAX_DIM_PROCS_ARRAY) {
		    int x,y, s;
		    if (width > height) {
		        y = length + 1;
		        x = ceil(y * width/(num)height);
		        s = ceil((num)height / (num)y);
		    }else{
		        x = length + 1;
		        y = ceil(x * (num)height/width);
		        s = ceil((num)width / (num)x);
		    }
		    nb_proc_aray_x[length] = x;
		    nb_proc_aray_y[length] = y;
		    block_size[length] = s;
		    
		    length++;
		    if (x*y > g.size) {
		        break;
		    }
		}
		
		g.image_width = width;
		g.image_height = height;
		g.bl_dim_x = nb_proc_aray_x[length - 2];
		g.bl_dim_y = nb_proc_aray_y[length - 2];
		g.block_size = block_size[length - 2];

		
		    int num_proc = g.bl_dim_x * g.bl_dim_y;
    
		LOG("Image: [%d x %d]\n", g.image_width, g.image_height);
		LOG("(%d)Number of procs. (with optimal) %d x %d= %d ~ %d \n", g.rank,
		           g.bl_dim_x, g.bl_dim_y, num_proc, g.size);
		LOG("Next should be: %d = [%d x %d]\n", nb_proc_aray_x[length-1]*nb_proc_aray_y[length-1],
		    nb_proc_aray_x[length-1], nb_proc_aray_y[length-1]);
		LOG("Proc size: %d x %d\n", g.block_size, g.block_size);
		LOG_LINE
    }
    
    MPI_Bcast(&g.image_width, 1, MPI_INT, g.main_proc, MPI_COMM_WORLD);
    MPI_Bcast(&g.image_height, 1, MPI_INT, g.main_proc, MPI_COMM_WORLD);
    MPI_Bcast(&g.bl_dim_x, 1, MPI_INT, g.main_proc, MPI_COMM_WORLD);
    MPI_Bcast(&g.bl_dim_y, 1, MPI_INT, g.main_proc, MPI_COMM_WORLD);
    MPI_Bcast(&g.block_size, 1, MPI_INT, g.main_proc, MPI_COMM_WORLD);

//    g.bl_idx_x = g.rank % g.bl_dim_x;
//	g.bl_idx_y = g.rank / g.bl_dim_x;
	g.sub_size = g.block_size + 2;
	
	//////////
	// Create Cartersian communication
	//////////
	create_cartersian_comm();
    
    int length_y, length_x;
    if (g.bl_idx_x < g.bl_dim_x - 1) {
        length_x = g.block_size;
    }else
        length_x = g.image_width - (g.bl_dim_x-1)*g.block_size;
    
    if(g.bl_idx_y < g.bl_dim_y - 1)
        length_y = g.block_size;
    else
        length_y = g.image_height - (g.bl_dim_y-1)*g.block_size;
    
    g.active_size_x = length_x;
    g.active_size_y = length_y;
    
    // 2. Load partial image
    // TODO: Optimize later. Now load all image
    image orignal;
    if(ReadImageObjGrayscale(&orignal, g.file_path) == 0){
        LOG("===Failed to read image data");
        goto Catch;
    }
        
    num *partial_data = malloc(g.sub_size * g.sub_size * sizeof(num));
    for (int y = -1; y < g.sub_size -1; y++) {
        for(int x = -1; x < g.sub_size -1; x++){
            vec2 local = {x, y}, global ={0,0};
            
            if (get_global_pixel_index(local, &global))
            {
                partial_data[local_array_idx(x, y)] =
                    orignal.Data[global.y * g.image_width + global.x];
            }
            else 
            {
                partial_data[local_array_idx(x, y)] = 0.;
     		}
        }
    }
  
    
    g.sub_image = partial_data;
    
    FreeImageObj(orignal);

    LOG("Read data successfully \n");
    return 0;
    

    
Catch:
	LOG("Failed in reading image \n");
    return 1;
};

void create_cartersian_comm(){
    int dims [ ] = { g.bl_dim_x , g.bl_dim_y } , 
    				periodic [ ] = { 0 , 0 } , 
    				reorder = 1 ;
    int error ;
    // Create Cart communication
    error = MPI_Cart_create ( MPI_COMM_WORLD
                             , 2  // Cartersian grid dimension
                             , dims // Number of process in each dimension
                             , periodic // true if first item connect to last item
                             , reorder // Ranking may be reorder, give the system a chance to optimize
                             , &grid_comm ) ;
                             
   // Set block index
   if(grid_comm != MPI_COMM_NULL){
   		int coordinates[2]; int my_grid_rank;
        int err = MPI_Comm_rank(grid_comm, &my_grid_rank);
        MPI_Cart_coords( grid_comm, my_grid_rank, 2, coordinates);
        
        g.bl_idx_x = coordinates[0];
        g.bl_idx_y = coordinates[1];
        
        g.rank = my_grid_rank;
   }
}

int log_local_phi(){
	if(grid_comm != MPI_COMM_NULL){
		// Write
		image img_;
		img_.Data = g.phi;
		img_.Height = g.sub_size;
		img_.Width = g.sub_size;
		img_.NumChannels = 1;
		char name[MAX_LEN_S_T];
		sprintf(name, "LOG/phi_%d.bmp", g.rank);
		
		printf("Write to %s \n", name);
		

		int s = WriteImage(img_.Data, img_.Width, img_.Height, name,
		                     IMAGEIO_NUM | IMAGEIO_GRAYSCALE | IMAGEIO_PLANAR, 1);
	   	if(s == 0 )
	   		printf("Error writing %s\n", name);
		return s;
	}
	return 1;
}

int log_local_image(){
	if(grid_comm != MPI_COMM_NULL){
		char name[MAX_LEN_S_T];
		sprintf(name, "LOG/img_%d_%d_%d.bmp", g.rank, g.bl_idx_x, g.bl_idx_y);
		
		printf("Write to %s \n", name);
		

		int s = WriteImage(g.sub_image, g.sub_size, g.sub_size, name,
		                     IMAGEIO_NUM | IMAGEIO_GRAYSCALE | IMAGEIO_PLANAR, 1);
	   	if(s == 0 ) // error
	   		printf("Error writing %s\n", name);
		return s;
	}
	return 1;
}

void debug_print(){

    char name[MAX_LEN_S_T];
    
    sprintf(name, "LOG/sub_%d.bmp", g.rank);
    printf("Write image: %s\n", name);
    image a;
    a.Data = g.sub_image;
    a.Width = g.sub_size;
    a.Height = g.sub_size;
    a.NumChannels = 1;
    WriteImageObj(a, name, 1);
}

int init_phi(){
    num * phi_data = malloc(g.sub_size * g.sub_size * sizeof(num));
    
    for (int x = -1; x < g.sub_size-1; x++) {
        for (int y = -1; y< g.sub_size-1; y++) {
            vec2 local = {x, y}, global ={0,0};
            if (get_global_pixel_index(local, &global))
            {
                phi_data[local_array_idx(x, y)] =
                    (num)(sin(global.x*M_PI/5.0)*sin(global.y*M_PI/5.0));
            }
            else
                phi_data[local_array_idx(x, y)] = 0.;
        }
    }
    
    g.phi = phi_data;
    
    LOG("Init phi success \n");
    return 0; // sucess
    
catch:
    return 1;
}

void region_average(num *c1, num *c2){
    num sum1 = 0.0, sum2 = 0.0;
    long count1 = 0, count2 = 0;
    
    for (int y = 0; y < g.active_size_y; y++) {
    	for (int x = 0; x < g.active_size_x; x++) {
            if (get_phi_data(x, y) >= 0) {
                count1++;
                sum1 += get_sub_image_data(x, y);
            }else{
                count2++;
                sum2 += get_sub_image_data(x, y);
            }
        }
    }
    
    num total_sum1 = 0.0, total_sum2 = 0.0;
    long total_count1 = 0, total_count2 = 0;
    
    
    num total[4];
    num temp[4] = {sum1, sum2, count1, count2};
    
    MPI_Allreduce(temp, total, 4, MPI_NUM, MPI_SUM, grid_comm);
    *c1 = total[0] / total[2];
    *c2 = total[1] / total[3];
    
//    MPI_Allreduce(&sum1, &total_sum1, 1, MPI_NUM, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&sum2, &total_sum2, 1, MPI_NUM, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&count1, &total_count1, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&count2, &total_count2, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); 
//    *c1 = total_sum1 / (num)total_count1;
//    *c2 = total_sum2 / (num)total_count2;
}

void update_boundary(){
	// For the edge process
	// We copy data from neighbor pixel
	// For derivation
    if (g.bl_idx_x == 0) {
        for (int y = 0; y < g.active_size_y; y++) {
            set_phi_data(-1, y, get_phi_data(0, y));
        }
    }
    if (g.bl_idx_x == (g.bl_dim_x - 1)) {
        for (int y = 0; y < g.active_size_y; y++) {
            set_phi_data(g.active_size_x, y, get_phi_data(g.active_size_x-1, y));
        }
    }
    if (g.bl_idx_y == 0) {
        for (int x = 0; x < g.active_size_x; x++) {
            set_phi_data(x, -1, get_phi_data(x, 0));
        }
    }
    if (g.bl_idx_y == (g.bl_dim_y-1)) {
        for (int x = 0; x < g.active_size_x; x++) {
            set_phi_data(x, g.active_size_y, get_phi_data(x, g.active_size_y - 1));
        }
    }
}
#ifndef USE_WINDOWS

#define EXCHANGE_BOUND 10001
void exchange_boundary(){

#ifdef DEBUG
	//printf("%d - Exchange bound\n", g.rank);
#endif

	if(g.rank >= g.bl_dim_x*g.bl_dim_y){
		return;
	}
	
    if (g.bl_idx_x > 0) { // LEFT
        // exchange left boundary with g.bl_idx_x - 1
        // static num * left = malloc(g.block_size*sizeof(num));
        for (int j = 0; j < g.block_size; j++) {
            left[j] = get_phi_data(-1, j);
        }
        
 //       printf("Left: rank %d to %d \n", g.rank, block_idx(g.bl_idx_x-1, g.bl_idx_y));
        
        // Send
//        MPI_Send(left,
//                 g.block_size,
//                 MPI_NUM,
//                 block_idx(g.bl_idx_x-1, g.bl_idx_y),
//                 EXCHANGE_BOUND,
//                 MPI_COMM_WORLD);
       
       int source, dest;
       MPI_Cart_shift(grid_comm, 0, -1, &source, &dest);
       	MPI_Request req;
        MPI_Isend(left,
                 g.block_size,
                 MPI_NUM,
                 dest,
                 EXCHANGE_BOUND,
                 grid_comm,
                 &req);
        

 //       printf("Left: rank %d to %d ; start receive\n", g.rank, block_idx(g.bl_idx_x-1, g.bl_idx_y));
        
        // Receive
        MPI_Status stat;
        num * right_recv = malloc(g.block_size*sizeof(num));
        MPI_Recv(right_recv, g.block_size, MPI_NUM, 
        		dest,
                 EXCHANGE_BOUND, grid_comm, &stat);
        for (int j = 0; j < g.block_size; j++) {
            set_phi_data(-1, j, right_recv[j]);
        }
        
    }
    if (g.bl_idx_x < (g.bl_dim_x-1)) { // RIGHT
        // Exchange right boundary with g.bl_idx_x  + 1
       //  static num * right = malloc(g.block_size*sizeof(num));
        for (int j = 0; j < g.block_size; j++) {
            right[j] = get_phi_data(g.block_size, j);
        }
        // Send
 //       printf("Right: rank %d to %d \n", g.rank, block_idx(g.bl_idx_x+1, g.bl_idx_y));
//        MPI_Send(right,
//                 g.block_size,
//                 MPI_NUM,
//                 block_idx(g.bl_idx_x+1, g.bl_idx_y),
//                 EXCHANGE_BOUND,
//                 MPI_COMM_WORLD);

       int source, dest;
       int left = MPI_Cart_shift(grid_comm, 0, +1, &source, &dest);
       
        MPI_Request req;         
        MPI_Isend(right,
                 g.block_size,
                 MPI_NUM,
                 dest,
                 EXCHANGE_BOUND,
                 grid_comm,
                 &req);
                                  
  //      printf("Right: rank %d to %d ; start revceive\n", g.rank, block_idx(g.bl_idx_x+1, g.bl_idx_y));
 
        // Receive
        MPI_Status stat;
        num *left_recv = malloc(g.block_size*sizeof(num));
        MPI_Recv(left_recv,
                 g.block_size,
                 MPI_NUM,
                 dest,
                 EXCHANGE_BOUND,
                 grid_comm,
                 &stat);
        for (int j = 0; j < g.block_size; j++) {
            set_phi_data(g.block_size, j, left_recv[j]);
        }
    }
    if (g.bl_idx_y < (g.bl_dim_y-1)) { // BOTTOM
        // Exchange bottom boundary with g.bl_idx_y+1
        
//        MPI_Send(g.phi+local_array_idx(g.bl_idx_x, g.bl_idx_y+1), g.block_size, MPI_NUM,
//                 block_idx(g.bl_idx_x, g.bl_idx_y+1), EXCHANGE_BOUND, MPI_COMM_WORLD);
        
       int source, dest;
       int left = MPI_Cart_shift(grid_comm, 1, +1, &source, &dest);
        MPI_Request req; 
        MPI_Isend(g.phi+local_array_idx(g.bl_idx_x, g.bl_idx_y+1), 
        			g.block_size, MPI_NUM,
                 dest, 
                 EXCHANGE_BOUND, grid_comm,
                 &req);
               
       // static num* top_recv = malloc(g.block_size*sizeof(num));
        MPI_Status stat;
        MPI_Recv(top, g.block_size, MPI_NUM,
                 dest, 
                 EXCHANGE_BOUND, grid_comm, &stat);
        for (int i = 0; i < g.block_size; i++) {
            set_phi_data(i, g.block_size, top[i]);
        }
    } 
    if (g.bl_idx_y > 0) { // TOP
        // Exchagne top boundary with g.bl_idx_y-1
//        MPI_Send(g.phi+1, g.block_size, MPI_NUM,
//                 block_idx(g.bl_idx_x, g.bl_idx_y-1), EXCHANGE_BOUND, MPI_COMM_WORLD);
        
               int source, dest;
       int left = MPI_Cart_shift(grid_comm, 1, -1, &source, &dest);
       
        MPI_Request req; 
        MPI_Isend(g.phi+1, g.block_size, MPI_NUM,
                 dest, EXCHANGE_BOUND, 
                 grid_comm, &req);
               
        
      //  static num* bottom_recv = malloc(g.block_size*sizeof(num));
        MPI_Status stat;
        MPI_Recv(bottom, g.block_size, MPI_NUM,
                 dest, 
                 EXCHANGE_BOUND, grid_comm, &stat);
        for (int i = 0; i < g.block_size; i++) {
            set_phi_data(i, -1, bottom[i]);
        }
    }

}
#endif

int block_idx(int x, int y){
    return y*g.bl_dim_x + x;
}

int chan_vese_loop(double *time){
    
    // Temp var
    const long NumPixels = ((long)g.image_width) * ((long)g.image_height);
    const long NumEl = NumPixels;
    const num *fPtr;
    num PhiDiffNorm, PhiDiff;
    num *PhiPtr;
    num c1Scalar=0., c2Scalar=0., Mu, Nu, Lambda1, Lambda2, dt;
    num *c1 = &c1Scalar, *c2 = &c2Scalar;
    num PhiLast, Delta, PhiX, PhiY, IDivU, IDivD, IDivL, IDivR;
    num Dist1, Dist2, PhiTol;
    int MaxIter;
    int iu, id, il, ir;
    
    Mu = g.opt.mu;
    Nu = g.opt.nu;
    Lambda1 = g.opt.lamda1;
    Lambda2 = g.opt.lamda2;
    dt = g.opt.dt;
    MaxIter = g.opt.max_iter;
    PhiTol = g.opt.phi_tol;
    PhiDiffNorm = (PhiTol > 0) ? PhiTol*1000 : 1000;
    
    int count = 0;
    num compute_time = 0.0;
    num commute_time = 0.0;
    

    init_window_buffer();

    double t;
    double total_t = 0.0;
    while (1) {
        // Update level set function

		t = MPI_Wtime();
        
        region_average(c1, c2);
#ifndef USE_WINDOWS        
        update_boundary();
#endif
        
        PhiPtr = g.phi + local_array_idx(0, 0);
        fPtr = g.sub_image + local_array_idx(0, 0);
        PhiDiffNorm = 0;
 
        
        for (int j = 0; j < g.active_size_y ; j++)
        {
            
            for (int i = 0; i < g.active_size_x; i++)
            {
                PhiPtr = g.phi + local_array_idx(i, j);
                fPtr = g.sub_image + local_array_idx(i, j);
                
                iu = -g.sub_size;
                id = g.sub_size;

                
                il = -1;
                ir = 1;
                
                Delta = dt/(M_PI*(1 + PhiPtr[0]*PhiPtr[0])); // Delta
                PhiX = PhiPtr[ir] - PhiPtr[0]; // grad plus x
                PhiY = (PhiPtr[id] - PhiPtr[iu])/2; // Grad mid y
                IDivR = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY)); // A(i,j)
                PhiX = PhiPtr[0] - PhiPtr[il];
                IDivL = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                PhiX = (PhiPtr[ir] - PhiPtr[il])/2;
                PhiY =  PhiPtr[id] - PhiPtr[0];
                IDivD = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                PhiY = PhiPtr[0] - PhiPtr[iu];
                IDivU = (num)(1/sqrt(DIVIDE_EPS + PhiX*PhiX + PhiY*PhiY));
                
                Dist1 = fPtr[0] - c1Scalar;
                Dist2 = fPtr[0] - c2Scalar;
                Dist1 *= Dist1;
                Dist2 *= Dist2;
                
                /* Semi-implicit update of phi at the current point */
                PhiLast = PhiPtr[0];
                PhiPtr[0] = (PhiPtr[0] + Delta*(
                                                Mu*(PhiPtr[ir]*IDivR + PhiPtr[il]*IDivL
                                                    + PhiPtr[id]*IDivD + PhiPtr[iu]*IDivU)
                                                - Nu - Lambda1*Dist1 + Lambda2*Dist2) ) /
                (1 + Delta*Mu*(IDivR + IDivL + IDivD + IDivU));
                PhiDiff = (PhiPtr[0] - PhiLast);
                PhiDiffNorm += PhiDiff * PhiDiff;
                
                if (PhiDiffNorm == NAN) {
                    
                }
            }
        }
        
         LOG("Loop 1\n");
        
#ifdef USE_WINDOWS    
 LOG("Loop 2\n");    
        update_windows_buffer();
#endif
        
        MPI_Barrier(grid_comm);
        
      //  printf("(%d) finish main loop %d\n", g.rank, count);
        
        // Error
       
        num total_e2 = 0.0;
        MPI_Allreduce(&PhiDiffNorm, &total_e2, 1, 
        				MPI_NUM, MPI_SUM, grid_comm);
      //  total_e2 = sqrt(total_e2/NumEl);

	//	printf("(%d) finish Allreduce loop %d \n", g.rank, count);

#ifdef DEBUG        
        LOG("Iter: %d - error: %f \n", count, total_e2); 
#endif
        
        // Exchange boundary
 //       
#ifdef USE_WINDOWS
 		exchange_boundary_windows();
 		copy_buffer_to_boundary();
#else
        LOG("Loop 1\n");
		exchange_boundary();
		// MPI_Barrier(grid_comm);
#endif

		count++;

		MPI_Barrier(grid_comm);
		t = MPI_Wtime() - t;
		total_t += t;
		if(total_t > 300){
			break;
		}
        
        
        
        if (count >= g.opt.max_iter
      //      || total_e2 <= PhiTol
            ) {
            break;
        }
    }
  
    delete_window_buffer();
    
    if(g.rank == 0){
 //   	printf("- %f - ", total_t);
    }
   
    *time = total_t;
    return count;
}

#ifdef USE_WINDOWS
void copy_buffer_to_boundary(){
	
	// Left
	if(g.bl_idx_x > 0)
	{
		num * lptr = &receive_buffer[LEFT];
	    for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(-1, i, *lptr++);
	    }
	}
	else
	{
		for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(-1, i, get_phi_data(0, i));
	    }
	}
	
	// Right
	if(g.bl_idx_x < (g.bl_dim_x-1))
	{
     	num * rptr = &receive_buffer[RIGHT];
     	for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(g.block_size, i, *rptr++);
	    }
	}
	else
	{
		for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(g.block_size, i, get_phi_data(g.block_size-1, i));
	    }
	}
	
	// Bottom
	if (g.bl_idx_y < (g.bl_dim_y-1))
	{
      	num * bptr = &receive_buffer[BOTTOM];
     	for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(i, g.block_size, *bptr++);
	    }
	}
	else
	{
		for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(i, g.block_size, get_phi_data(i, g.block_size-1));
	    }
	}
	
	// top
	if (g.bl_idx_y > 0)
	{
 		num * tptr = &receive_buffer[TOP];
     	for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(i, -1, *tptr++);
	    }
	}
	else
	{
		for(int i = 0; i < g.block_size; i++)
	    {
	    	set_phi_data(i, -1, get_phi_data(i, 0));
	    }
	}
}
#endif

void init_window_buffer(){
#ifdef USE_WINDOWS
	buffer = malloc(4 * g.block_size * sizeof(num));
	receive_buffer = malloc(4 * g.block_size * sizeof(num));
	// Create windows
	MPI_Win_create(buffer
					, 4 * g.block_size // Window size
					, sizeof(num) // element size
					, MPI_INFO_NULL
					, grid_comm
					, &win);
#else
	left = malloc(g.block_size * sizeof(num));
	right = malloc(g.block_size * sizeof(num));
	top = malloc(g.block_size * sizeof(num));
	bottom = malloc(g.block_size * sizeof(num));
	right_recv = malloc(g.block_size * sizeof(num));
	left_recv = malloc(g.block_size * sizeof(num));
#endif
}

void delete_window_buffer(){
#ifdef USE_WINDOWS
	free(buffer);
	free(receive_buffer);
	MPI_Win_free(&win);
#else
	free(left);
	free(right);
	free(bottom);
	free(top);
	free(right_recv);
	free(left_recv);
#endif
}

#ifdef USE_WINDOWS
void update_windows_buffer(){
	// Left
	num * lptr = &buffer[LEFT];
	for(int i = 0; i < g.block_size; i++){
		*(lptr++) = get_phi_data(0, i);
	}
	
	// Right
	num * rptr = &buffer[RIGHT];
	for(int i = 0; i < g.block_size; i++){
		*(rptr++) = get_phi_data(g.block_size-1, i);
	}
	
	// Top
	num * tptr = &buffer[TOP];
	for(int i = 0; i < g.block_size; i++){
		*(tptr++) = get_phi_data(i, 0);
	}
	
	// Bottom
	num * bptr = &buffer[BOTTOM];
	for(int i = 0; i < g.block_size; i++){
		*(bptr++) = get_phi_data(i, g.block_size-1);
	}
}

void exchange_boundary_windows(){
	MPI_Win_fence(0, win);
	
	// Left
	// exchange left boundary with right of g.bl_idx_x - 1
	if(g.bl_idx_x > 0)
	{
	    int source, dest;
        MPI_Cart_shift(grid_comm, 0, -1, &source, &dest);
		MPI_Get(&receive_buffer[LEFT]	// copy to left
				, g.block_size			// size
				, MPI_NUM				// type
				, dest					// rank of target
				, RIGHT					// displace to right
				, g.block_size			// size
				, MPI_NUM				// type
				, win
				);
	}
	// Right
	// Copy left of x+1 to right
	if(g.bl_idx_x < (g.bl_dim_x-1))
	{
	   	int source, dest;
       	int left = MPI_Cart_shift(grid_comm, 0, +1, &source, &dest);
 		MPI_Get(&receive_buffer[RIGHT]	// copy to left
				, g.block_size			// size
				, MPI_NUM				// type
				, dest					// rank of target
				, LEFT					// displace to right
				, g.block_size			// size
				, MPI_NUM				// type
				, win
				);      
	}
	
	// Bottom
	if (g.bl_idx_y < (g.bl_dim_y-1))
	{
		int source, dest;
      	int left = MPI_Cart_shift(grid_comm, 1, +1, &source, &dest);
  		MPI_Get(&receive_buffer[BOTTOM]	// copy to bottom
				, g.block_size			// size
				, MPI_NUM				// type
				, dest					// rank of target
				, TOP					// displace to top
				, g.block_size			// size
				, MPI_NUM				// type
				, win
				);      	
	}
	
	// top
	if (g.bl_idx_y > 0)
	{
		int source, dest;
       	int left = MPI_Cart_shift(grid_comm, 1, -1, &source, &dest);
  		MPI_Get(&receive_buffer[TOP]	// copy to bottom
				, g.block_size			// size
				, MPI_NUM				// type
				, dest					// rank of target
				, BOTTOM					// displace to top
				, g.block_size			// size
				, MPI_NUM				// type
				, win
				);  
	}
	
	MPI_Win_fence(0, win);
}
#endif
int WriteBinary(image Phi, const char *File)
{
    unsigned char *Temp = NULL;
    const int NumPixels = Phi.Width*Phi.Height;
    int i, Success;
    
    if(!(Temp = (unsigned char *)malloc(Phi.Width*Phi.Height)))
        return 0;
    
    for(i = 0; i < NumPixels; i++)
        Temp[i] = (Phi.Data[i] >= 0) ? 255 : 0;
    
    Success = WriteImage(Temp, Phi.Width, Phi.Height, File,
                         IMAGEIO_U8 | IMAGEIO_GRAYSCALE, 0);
    
    free(Temp);
    return Success;
}

void gather_phi_p(int iter){
    num * main_domain = malloc(g.block_size * g.block_size * sizeof(num));
    
    // Copy main main part to main_domain
    for (int i = 0; i < g.active_size_y; i++) {
        memcpy(main_domain + i*g.block_size,		// Dest
               g.phi  + local_array_idx(0, i) ,		// Source
               g.active_size_x*sizeof(num));		// Size
    }
    
    if (g.size == 1) // Only one proc
    {
    	// Write
        image phi_;
        phi_.Data = main_domain;
        phi_.Height = g.image_height;
        phi_.Width = g.image_width;
        phi_.NumChannels = 1;
        char name[MAX_LEN_S_T];
        sprintf(name, "LOG/phi_%d.bmp", iter);
        
        WriteBinary(phi_, name);
        LOG("Level set function written to LOG/phi.bmp\n");
        
        free(main_domain);
    	return;
    }
    
    if (g.rank != g.main_proc) { 
    // Send domain part from all to main proc
        MPI_Send(main_domain,
                 g.block_size * g.block_size,
                 MPI_NUM,
                 g.main_proc,
                 GATHER_PHI,
                 grid_comm);
    }
    else { // receive
        num * phi_total = malloc(g.image_width*g.image_height*sizeof(num));
        for (int blx = 0; blx < g.bl_dim_x; blx ++) {
            for (int bly = 0; bly < g.bl_dim_y; bly++) {
                int idx = bly*g.bl_dim_x + blx;
                num * receive_domain = malloc(g.block_size * g.block_size * sizeof(num));
                
                
                
                if (idx == g.main_proc) {
                    memcpy(receive_domain, main_domain, g.block_size * g.block_size * sizeof(num));
                }
                else{
                    MPI_Status stat;
                    MPI_Recv(receive_domain,
                             g.block_size * g.block_size,
                             MPI_NUM,
                             idx,
                             GATHER_PHI,
                             grid_comm,
                             &stat);
                }
                
                // Copy
                int length_y, length_x;
                if (blx < g.bl_dim_x - 1) {
                    length_x = g.block_size;
                }else
                    length_x = g.image_width - (g.bl_dim_x-1)*g.block_size;
                
                if(bly < g.bl_dim_y - 1)
                    length_y = g.block_size;
                else
                    length_y = g.image_height - (g.bl_dim_y-1)*g.block_size;
                
                for (int y = 0; y < length_y; y++) {
                    int gx = blx*g.block_size;
                    int gy = bly*g.block_size + y;
                    memcpy(phi_total + gy*g.image_width+gx,
                           receive_domain + y*g.block_size,
                           length_x*sizeof(num));
                }
                
                free(receive_domain);
                
            }
            
            
        }
        
        // Write
        image phi_;
        phi_.Data = phi_total;
        phi_.Height = g.image_height;
        phi_.Width = g.image_width;
        phi_.NumChannels = 1;
        char name[MAX_LEN_S_T];
        sprintf(name, "LOG/phi_%d.bmp", iter);
        
        WriteBinary(phi_, name);
        LOG("Level set function written to LOG/phi.bmp\n");

        free(phi_total);
    } /* if (g.rank == g.main_proc) */
    
    free(main_domain);
}

void gather_phi(){
	

    num * main_domain = malloc(g.block_size * g.block_size * sizeof(num));
    
    // Copy main main part to main_domain
    for (int i = 0; i < g.active_size_y; i++) {
        memcpy(main_domain + i*g.block_size,		// Dest
               g.phi  + local_array_idx(0, i) ,		// Source
               g.active_size_x*sizeof(num));		// Size
    }
    
    if(g.rank == g.main_proc) 
    { // receive
    
        num * phi_total = malloc(g.image_width*g.image_height*sizeof(num));
        for (int blx = 0; blx < g.bl_dim_x; blx ++) {
            for (int bly = 0; bly < g.bl_dim_y; bly++) {
               int idx;
               int coord[2] = {blx, bly};
               MPI_Cart_rank(grid_comm, coord, &idx);
                
                num * receive_domain = malloc(g.block_size * g.block_size 
                							* sizeof(num));
                
                
                
                if (idx == g.main_proc) {
                    memcpy(receive_domain, main_domain, 
                    	g.block_size * g.block_size * sizeof(num));
                }
                else
                {
                    MPI_Status stat;
                    MPI_Recv(receive_domain,
                             g.block_size * g.block_size,
                             MPI_NUM,
                             idx,
                             GATHER_PHI,
                             MPI_COMM_WORLD,
                             &stat);
                    
                }
                
                // Copy
                int length_y, length_x;
                if (blx < g.bl_dim_x - 1) {
                    length_x = g.block_size;
                }else
                    length_x = g.image_width - (g.bl_dim_x-1)*g.block_size;
                
                if(bly < g.bl_dim_y - 1)
                    length_y = g.block_size;
                else
                    length_y = g.image_height - (g.bl_dim_y-1)*g.block_size;
                
                for (int y = 0; y < length_y; y++) {
                    int gx = blx*g.block_size;
                    int gy = bly*g.block_size + y;
                    memcpy(phi_total + gy*g.image_width+gx,
                           receive_domain + y*g.block_size,
                           length_x*sizeof(num));
                }
                
                free(receive_domain);
                
            }
            
            
        }
        
        // Write
        image phi_;
        phi_.Data = phi_total;
        phi_.Height = g.image_height;
        phi_.Width = g.image_width;
        phi_.NumChannels = 1;
        char name[MAX_LEN_S_T];
        sprintf(name, "LOG/phi.bmp");
        
        WriteBinary(phi_, name);
        LOG("Level set function written to LOG/phi.bmp\n");

        free(phi_total);
    } /* if (g.rank == g.main_proc) */
    
    if (g.rank != g.main_proc 
    	&& g.rank < g.bl_dim_x*g.bl_dim_y) 
    { 
    // Send domain part from all to main proc
    //	printf("Send phi from %d \n", g.rank);
    	
        MPI_Send(main_domain,
                 g.block_size * g.block_size,
                 MPI_NUM,
                 g.main_proc,
                 GATHER_PHI,
                 MPI_COMM_WORLD);
        
    }
    
    free(main_domain);
}

int get_global_pixel_index(vec2 const local, vec2 *global){
	
    global->x = g.bl_idx_x * g.block_size + local.x;
    global->y = g.bl_idx_y * g.block_size + local.y;
    
    if (global->x >= g.image_width || global->x < 0
        || global->y >= g.image_height || global->y < 0) {
        return 0; // Fail
    }
    
    return 1; // success
}
void print_g(){
    printf("Proc: %d / %d \n", g.rank, g.size);

}

void print_mat(num* data, int width, int height){
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width ; x++) {
            printf("%f ", data[y*width + x]);
        }
        printf("\n");
    }
}

int max_(int a, int b){return a>b? a:b;};

/////////////////////////////////////////////////
// Parser argument
#define FILE_OPTION "-file"
#define MAX_ITER 	"-iters"
#define IMG_SIZE	"-size"

void print_help(){
	if(g.rank == g.main_proc){
		printf("ARGUMENTS: \n");
		printf(" %s []: image file \n", FILE_OPTION);
		printf(" %s	[]: Max iteration \n", MAX_ITER);
		printf(" %s	[]: Image size \n", IMG_SIZE);
	}
}

int parse_arguments(int argc, char* argv[]){
	if(argc < 2)
		return 1;
		
	int k = 1;
	while( k < argc){
		char* option = argv[k];
		if(k + 1 > argc) return 1;
		char* value = argv[k+1]; // Should checlk also
		
		
		
		if(strcmp(option, FILE_OPTION) == 0){
			strcpy(g.file_path, value);
		}
		else if(strcmp(option, MAX_ITER) == 0){
			g.opt.max_iter = atoi(value);
		}
		else if (strcmp(option, IMG_SIZE) == 0){
			g.img_size = atoi(value);
		}
		else{
			LOG("Wrong input: %s \n", option);
			return 1;
		}
		
		k+=2;
	}
   	return 0;
}

void generate_image(){
	
	if(g.rank == g.main_proc){
		num * data = malloc(g.img_size * g.img_size * sizeof(num));
		num center = g.img_size/2;
		num radius2 = g.img_size / 4; radius2 = radius2 * radius2;
		for(int i = 0; i < g.img_size; i++){
			for(int j = 0; j < g.img_size; j++){
				if((i-center)*(i-center) + (j-center)*(j-center) < radius2){
					data[i * g.img_size + j] = 1.0;
				}else{
					data[i * g.img_size + j] = 0.0;
				}
			}
		}
	
		WriteImage(data, g.img_size, g.img_size, g.file_path,
		                     IMAGEIO_NUM | IMAGEIO_GRAYSCALE | IMAGEIO_PLANAR, 1);
		LOG("Dummy image generated: %s\n", g.file_path);
    
    }	
}

/////////////////////////////////////////////////
// Information of program settings
void print_info(){
    
    LOG_LINE
    LOG("%d proc\n", g.size);
    
    // 1. Precision
#ifdef NUM_SINGLE
    LOG("Single precision\n");
#else
    LOG("Double precision\n");
#endif /* SINGLE_PRECISION */
    
    // 2. File path
    LOG("Image file: %s \n", g.file_path);
    
    LOG_LINE
}


void release_memory(){
    if (g.sub_image) {
        free(g.sub_image);
    }
    
    if (g.phi) {
        free(g.phi);
    }
}
