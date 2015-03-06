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

// Comment for double precision
#define SINGLE_PRECISION

// Comment for release mode
#define DEBUG

///////////////////////////////////////////////////////////////
// GLOBAL VARIABLES

// Process parameters
extern int size, rank;
extern const int main_proc;


///////////////////////////////////////////////////////////////
#ifdef SINGLE_PRECISION
typedef float num;
#else
typedef double num;
#endif /* SINGLE_PRECISION */

///////////////////////////////////////////////////////////////
// MARCO

// LOG(char*)
// Log by main proc
#define LOG(a, ...)     if(rank == main_proc){ \
printf(a, ##__VA_ARGS__); \
}

// LOG_ALL(char*)
// Log by all procs
#define LOG_ALL(a, ...)   printf("[%d / %d: ]", rank, size); \
printf(a, ##__VA_ARGS__); \

// assert(bool)
#ifdef DEBUG
    #define assert(a, b) if(!(a)) { printf("%s", b); }
#else
    #define assert(a, b)
#endif /* DEBUG */

#endif /* imageSeg_define_h */
