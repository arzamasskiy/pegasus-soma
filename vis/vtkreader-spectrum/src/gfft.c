/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "common.h"
#include "debug.h"

#ifdef MPI_SUPPORT
// That's a long MPI if def...

#include "transpose.h"

const int n_size2D[2] = {NX, NZ};
const int n_size1D[1] = {NY_COMPLEX};

fftw_plan	r2c_2d, c2r_2d, r2c_1d, c2r_1d;

double complex *wi1;
double *wir1;

/* GFFT (Like Geo's FFT as you might have guessed...) is an FFT wrapper for FFTW>=3.2
It takes care of the MPI part of the FFT while FFTW deals with the FFT themselves
GFFT can handle threaded FFTS (just ask...).
One concludes It's FAAAARRRR better than FFTW 2.1.5
*/

	
// This is an inplace real 2 complex transform
// Assumes wrin has the logical dimensions [NY/PROC, NX, NZ] of real positions
// physical dimensions [NY/NPROC, NX, NZ+2];
void gfft_r2c_t(double *wrin) {
	int i;
	double complex *win = (double complex *) wrin;
	
	//start transforming in 2D wrin
	fftw_execute_dft_r2c(r2c_2d, wrin, win);
	
	// The logical dimensions of win are [NY_COMPLEX/NPROC, NX_COMPLEX, NZ_COMPLEX]
	// transpose it
	transpose_complex_YX(win, win);
	
	// We now have an array with logical dimensions[NX_COMPLEX/NPROC, NY_COMPLEX, NZ_COMPLEX]
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(r2c_1d, &win[i*NY_COMPLEX*NZ_COMPLEX],&win[i*NY_COMPLEX*NZ_COMPLEX]);

	// done...
	return;
}


void gfft_c2r_t(double complex *win) {
	int i;
	double *wrin = (double *) win;
	// We now have an array with logical dimensions[NX_COMPLEX/NPROC, NY_COMPLEX, NZ_COMPLEX]
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif	
	for(i=0 ; i < NX_COMPLEX/NPROC ; i++) 
		fftw_execute_dft(c2r_1d, &win[i*NY_COMPLEX*NZ_COMPLEX],&win[i*NY_COMPLEX*NZ_COMPLEX]);
		
	// The logical dimensions of win are [NX_COMPLEX/NPROC, NY_COMPLEX, NZ_COMPLEX]
	// transpose it
	transpose_complex_XY(win, win);
	
	// The final 2D transform
	fftw_execute_dft_c2r(c2r_2d, win, wrin);
	 // and we're done !
	 return;
}
														  
	
// In place double transpose transforms
// Not Fast, but convenient...!

void gfft_r2c(double *wrin) {
	transpose_real(NX, NY, NZ+2, NPROC, wrin, wrin);
	gfft_r2c_t(wrin);
	return;
}

void gfft_c2r(double complex *win) {
	double *wrin = (double *) win;
	gfft_c2r_t(win);
	transpose_real(NY,NX,NZ+2,NPROC,wrin,wrin);
	return;
}	

void init_gfft() {
	// This will init the plans needed by gfft
	// Transform of NY/NPROC arrays of (logical) size [NX, NZ]
	// The physical size is [NX, NZ+2]
	// We use in-place transforms
	int i;
	
	DEBUG_START_FUNC;
	
	wi1 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (wi1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wi1 allocation");

	wir1 = (double *) wi1;
	
	for(i = 0 ; i < NTOTAL_COMPLEX; i++) {
		wi1[i]=1.0;
	}
	
#ifdef _OPENMP
	fftw_plan_with_nthreads( nthreads );
#endif
	r2c_2d = fftw_plan_many_dft_r2c(2, n_size2D, NY / NPROC, wir1, NULL, 1, (NZ+2)*NX,
															 wi1,  NULL, 1, (NZ+2)*NX/2, FFT_PLANNING);
	if (r2c_2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW R2C_2D plan creation failed");
														   
	c2r_2d = fftw_plan_many_dft_c2r(2, n_size2D, NY / NPROC, wi1,  NULL, 1, (NZ+2)*NX/2,
														    wir1, NULL, 1, (NZ+2)*NX  , FFT_PLANNING);
	if (c2r_2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW C2R_2D plan creation failed");
	
															 
	// 1D transforms: This are actually c2c transforms, but are used for global 3D transforms.
	// We will transform forward and backward an array of logical size [NX/NPROC, NY, (NZ+2)/2] along the 2nd dimension
	// We will do NZ_COMPLEX transforms along Y. Will need a loop on NX/NPROC
	// We use &w1[NZ_COMPLEX] so that alignement check is done properly (see SIMD in fftw Documentation)
	
#ifdef _OPENMP	
	fftw_plan_with_nthreads( 1 );
#endif	
	r2c_1d = fftw_plan_many_dft(1, n_size1D, NZ_COMPLEX, &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1,
														 &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1, FFTW_FORWARD, FFT_PLANNING);
	if (r2c_1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW R2C_1D plan creation failed");
																			  
	c2r_1d = fftw_plan_many_dft(1, n_size1D, NZ_COMPLEX, &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1,
														 &wi1[NZ_COMPLEX], NULL, NZ_COMPLEX, 1, FFTW_BACKWARD, FFT_PLANNING);
	if (c2r_1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW C2R_1D plan creation failed");

	// init transpose routines
	init_transpose();
	// Let's see which method is faster (with our without threads)
		
	fftw_free(wi1);

	DEBUG_END_FUNC;
	
	return;
}

void finish_gfft() {
	DEBUG_START_FUNC;
	
	fftw_destroy_plan(r2c_2d);
	fftw_destroy_plan(c2r_2d);
	fftw_destroy_plan(c2r_1d);
	fftw_destroy_plan(r2c_1d);
	
	finish_transpose();
	
	DEBUG_END_FUNC;
	
	return;
}

/**********************************************************************
***********************************************************************
*********     N O    M P I    R O U T I N E S        ******************
***********************************************************************
***********************************************************************/
// These routines are essentially wrappers for fftw3 routines

#else

// Here, we assume we don't have MPI

fftw_plan	r2cfft, c2rfft;

double complex *wi1;
double *wir1;

void gfft_r2c_t(double *wrin) {
	double complex *win = (double complex *) wrin;
	fftw_execute_dft_r2c(r2cfft, wrin, win);
	return;
}

void gfft_c2r_t(double complex *win){
	double *wrin = (double *) win;
	fftw_execute_dft_c2r(c2rfft, win, wrin);
	return;
}

void gfft_r2c(double *wrin) {
	double complex *win = (double complex *) wrin;
	fftw_execute_dft_r2c(r2cfft, wrin, win);
	return;
}

void gfft_c2r(double complex *win){
	double *wrin = (double *) win;
	fftw_execute_dft_c2r(c2rfft, win, wrin);
	return;
}

void init_gfft() {
	
	DEBUG_START_FUNC;
	
	wi1 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (wi1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for wi1 allocation");
	
	wir1 = (double *) wi1;
	
#ifdef _OPENMP
	fftw_plan_with_nthreads( nthreads );
#endif
	
	r2cfft = fftw_plan_dft_r2c_3d( NX, NY, NZ, wr1, w1,  FFT_PLANNING);
	if (r2cfft == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW R2C plan creation failed");

	c2rfft = fftw_plan_dft_c2r_3d( NX, NY, NZ, w1,  wr1, FFT_PLANNING);
	if (c2rfft == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW C2R plan creation failed");

	
	fftw_free(wi1);
	
	DEBUG_END_FUNC;
	
	return;
}

void finish_gfft() {
	DEBUG_START_FUNC;
	
	fftw_destroy_plan(r2cfft);
	fftw_destroy_plan(c2rfft);
	
	DEBUG_END_FUNC;
	
	return;
}

#endif

														   
	