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
#include <time.h>
#include <stdlib.h>
#include "define.h"

// Information of program setting
void print_info(){
#ifdef DEBUG
#ifdef SINGLE_PRECISION
    LOG("Single precision\n");
#else
    LOG("Double precision\n");
#endif /* SINGLE_PRECISION */
#endif /* DEBUG */
}


#endif
