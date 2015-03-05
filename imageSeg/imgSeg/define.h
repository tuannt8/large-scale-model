//
//  define.h
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 3/5/15.
//
//

#ifndef imageSeg_define_h
#define imageSeg_define_h

// Log by main thread
#define LOG(a, ...)     if(rank == MAIN_PROC){ \
printf(a, ##__VA_ARGS__); \
}

// Log performed by all threads
#define LOG_ALL(a, ...)   printf("[%d / %d: ]", rank, size); \
printf(a, ##__VA_ARGS__); \

// Assertion
#ifdef DEBUG
    #define assert(a, b) if(!(a)) { printf("%s", b); }
#else
    #define assert(a, b) 
#endif

#endif
