/*
 * Utils : Helpful utilities.
 *
 *      Specifically: constructor of arrays that's guaranteed to be contiguous
 *                    in memory.
 */
#ifndef UTILS_H
#define UTILS_H

//Creates an allocated, contiguous-in-memory 2d array of size nx*nz
float **new_contiguous_2dArray(int nx, int ny);

//Frees 2d array (all frees the pointer passed to it!)
void free_contiguous_2dArray(float **array);


// See above, just in 3d
float ***new_contiguous_3dArray(int nx, int ny, int nz);
void free_contiguous_3dArray(float ***array);

#endif

