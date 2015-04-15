//
//  sub_image.c
//  imageSeg
//
//  Created by Tuan Nguyen Trung on 4/10/15.
//
//

#include "util.h"

inline num get_sub_image_data(int x, int y){
    return g.sub_image[(y+1)*g.sub_size + x + 1];
}

inline void set_sub_image_data(int x, int y, int inten){
    g.sub_image[(y+1)*g.sub_size + x + 1] = inten;
}

inline num get_phi_data(int x, int y){
    return g.phi[(y+1)*g.sub_size + x + 1];
}

inline void set_phi_data(int x, int y, int inten){
    g.phi[(y+1)*g.sub_size + x + 1] = inten;
}

inline int proc_index(int x, int y){
    return y*g.bl_dim_x + x;
}

inline int local_array_idx(int x, int y)
{
    return (y+1)*g.sub_size + x + 1;
}