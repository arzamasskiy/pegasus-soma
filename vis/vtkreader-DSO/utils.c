#include <stdlib.h>
#include "assert.h"

float **new_contiguous_2dArray(int nx, int ny){
    float *newarr;
    float **ret;
    int i;

    newarr = malloc(nx*ny*sizeof(float));
    ret= malloc(nx*sizeof(float));
    assert(newarr != NULL && ret != NULL);

    for(i = 0; i < nx; i++){
        ret[i] = &(newarr[i*ny]);
    }
    return ret;
}

void free_contiguous_2dArray(float **array){
    if(array == NULL) return;
    free(array[0]);
    free(array);
}


float ***new_contiguous_3dArray(int nx, int ny, int nz){
    float ***ret;
    float **nextarr;
    float *newarr;
    int i,j;

    newarr = (float*) malloc(nx*ny*nz*sizeof(float));
    nextarr = (float**) malloc(nx*ny*sizeof(float*));
    ret = (float***) malloc(nx*sizeof(float**));
    assert(newarr != NULL && nextarr != NULL && ret != NULL);

    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
        nextarr[ny*i + j] = &newarr[i*(ny*nz)+j*nz];
      }
    }
    for(i = 0; i < nx; i++){
    ret[i] = &nextarr[i*ny];
    }
    return ret;
}

void free_contiguous_3dArray(float ***array){
    if(array == NULL) return;
    free(array[0][0]);
    free(array[0]);
    free(array);
}
