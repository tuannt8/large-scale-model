//
//  define.h
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 3/5/15.
//
//

#ifndef imageSeg_define_h
#define imageSeg_define_h
///////////////////////////////////////////////////////////////
// PROGRAM SETTING

 #define NUM_SINGLE

// Comment for release mode
// #define DEBUG

///////////////////////////////////////////////////////////////
// GLOBAL VARIABLES


///////////////////////////////////////////////////////////////
#ifdef NUM_SINGLE
typedef float num;
#define MPI_NUM MPI_FLOAT
#else
typedef double num;
#define MPI_NUM MPI_DOUBLE
#endif /* SINGLE_PRECISION */

#define DIVIDE_EPS       ((num)1e-16)
#define BOUND_O_LAP         2

///////////////////////////////////////////////////////////////
// MARCO

// LOG(char*)
// Log by main proc
#ifdef DEBUG
#define LOG(a, ...)     if(g.rank == g.main_proc){ \
printf(a, ##__VA_ARGS__); \
}

#define LOG_LINE LOG("====================================\n");

#else
#define LOG(a, ...)
#define LOG_LINE
#endif



// LOG_ALL(char*)               |
//-------------------------------
// Log by all procs
#define LOG_ALL(a, ...)   printf("[%d / %d: ]", g.rank, g.size); \
printf(a, ##__VA_ARGS__); \

// assert(bool)
#ifdef DEBUG
    #define assert(a, b) if(!(a)) { printf("%s", b); }
#else
    #define assert(a, b)
#endif /* DEBUG */

#endif /* imageSeg_define_h */


///////////////////////////////////////////////////////////////
// Other
//#define MAX_STRING_T 500;
