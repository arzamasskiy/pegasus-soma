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



#include "snoopy.h"

#include <stdlib.h>

// Timing support
#ifndef MPI_SUPPORT
#ifndef _OPENMP
#include <stdio.h>
#include <time.h>
#endif
#endif

#include "error.h"
#include "gfft.h"
#include "debug.h"

// This are global variables used throughout the code
// Wave number pointers
double	*kx;	/**< x Wavevector */
double	*ky;	/**< y Wavevector */
double	*kz;	/**< z Wavevector */
double	*kxt;	/**< Time dependant x Wavevector. Different from kx only when SHEAR is present.*/
double	*k2t;	/**<  k squared Wavevector, function of time when SHEAR is present.*/
double	*ik2t;  /**< inverse of k2t Wavevector, function of time when SHEAR is present. set to 0 wheh k2t=0 to avoid singularity*/
double  *x;
double  *y;
double  *z;

double	kxmax,	kymax,  kzmax,	kmax;	/**< Maximum wavevectors */


// Mask for dealiasing
double   *mask;	/**< Deasliasing Mask*/

double	*wr1,	*wr2,	*wr3;		/** Temporary real array (alias of complex w**) */
double	*wr4,	*wr5,	*wr6;		/** Temporary real array (alias of complex w**) */
double	*wr7,	*wr8,	*wr9;		/** Temporary real array (alias of complex w**) */
double	*wr10,	*wr11,	*wr12;		/** Temporary real array (alias of complex w**) */
double	*wr13,	*wr14,	*wr15;		/** Temporary real array (alias of complex w**) */

struct Field			fld;

double complex		*w1,	*w2,	*w3;	/**< Temporary complex array (alias of real wr**) */
double complex		*w4,	*w5,	*w6;	/**< Temporary complex array (alias of real wr**) */
double complex		*w7,	*w8,	*w9;	/**< Temporary complex array (alias of real wr**) */
double complex		*w10,	*w11,	*w12;	/**< Temporary complex array (alias of real wr**) */
double complex		*w13,	*w14,	*w15;	/**< Temporary complex array (alias of real wr**) */
												

// Parameters
struct Parameters			param;

int		nthreads;								/**< Number of OpenMP threads */

/** Init all global variables, aligning them in memory */
void init_common(void) {
	/* This routine will initialize everything */
	int i,j,k;
	
	DEBUG_START_FUNC;
	
#ifdef _OPENMP
	if( !(fftw_init_threads()) ) ERROR_HANDLER( ERROR_CRITICAL, "Threads initialisation failed");
#endif
	
	/* We start with the coordinate system */
	kx = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kx allocation");
	
	ky = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (ky == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ky allocation");
	
	kz = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kz allocation");
	
	kxt = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kxt == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kxt allocation");
	
	k2t = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (k2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for k2t allocation");
	
	ik2t = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (ik2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ik2t allocation");
	
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX );
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX );
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX );
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");
	
	for( i = 0; i < NX_COMPLEX ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kx[ IDX3D ] = (2.0 * M_PI) / param.lx *
						(fmod( NX_COMPLEX + i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 );
					 
				ky[ IDX3D ] = (2.0 * M_PI) / param.ly *
						(fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 );
					 
				kz[ IDX3D ] = (2.0 * M_PI) / param.lz * k;

				kxt[ IDX3D ]= kx[IDX3D];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
								ky[IDX3D] * ky[IDX3D] +
								kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}
	
	kxmax = 2.0 * M_PI/ param.lx * ( (NX / 2) - 1);
	kymax = 2.0 * M_PI/ param.ly * ( (NY / 2) - 1);
	kzmax = 2.0 * M_PI/ param.lz * ( (NZ / 2) - 1);
	
	kmax=pow(kxmax*kxmax+kymax*kymax+kzmax*kzmax,0.5);
	
	/* Initialize the dealiazing mask Or the nyquist frequency mask (in case dealiasing is not required) */
	
	mask = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (mask == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for mask allocation");
	
	for( i = 0; i < NX_COMPLEX ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {

				mask[ IDX3D ] = 1.0;
				if(param.antialiasing) {
					if( fabs( kx[ IDX3D ] ) > 2.0/3.0 * kxmax)
						mask[ IDX3D ] = 0.0;
				
					if( fabs( ky[ IDX3D ] ) > 2.0/3.0 * kymax)
						mask[ IDX3D ] = 0.0;
					
					if( fabs( kz[ IDX3D ] ) > 2.0/3.0 * kzmax)
						mask[ IDX3D ] = 0.0;
				}
				else {			
					if ( i == NX_COMPLEX / 2 ) 
						mask[ IDX3D ] = 0.0;
					if ( j == NY_COMPLEX / 2 )  
						mask[ IDX3D ] = 0.0;
					if ( k == NZ_COMPLEX ) 
						mask[ IDX3D ] = 0.0;
				}
			}
		}
	}

	if(param.antialiasing) {
		kxmax = kxmax * 2.0 / 3.0;
		kymax = kymax * 2.0 / 3.0;
		kzmax = kzmax * 2.0 / 3.0;
		kmax = kmax * 2.0 / 3.0;
	}

// Allocate fields
// Complex fields
	
	fld.vx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vx allocation");
	
	fld.vy = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vy allocation");
	
	fld.vz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vz allocation");
	
#ifdef BOUSSINESQ
	fld.th = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.th allocation");
#endif
#ifdef MHD
	fld.bx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vx allocation");
	
	fld.by = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vy allocation");
	
	fld.bz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fld.bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fld.vz allocation");
#endif

	w1 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1 allocation");
	
	w2 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w2 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2 allocation");
	
	w3 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w3 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w3 allocation");
	
	w4 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w4 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w4 allocation");
	
	w5 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w5 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w5 allocation");
	
	w6 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w6 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w6 allocation");
	
	w7 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w7 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w7 allocation");
	
	w8 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w8 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w8 allocation");
	
	w9 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w9 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w9 allocation");
	
	w10 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w10 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w10 allocation");
	
	w11 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w11 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w11 allocation");
	
	w12 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w12 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w12 allocation");
	
	w13 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w13 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w13 allocation");
	
	w14 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w14 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w14 allocation");
	
	w15 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w15 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w15 allocation");
	
	/* Will use the same memory space for real and complex fields */
	
	wr1 = (double *) w1;
	wr2 = (double *) w2;
	wr3 = (double *) w3;
	wr4 = (double *) w4;
	wr5 = (double *) w5;
	wr6 = (double *) w6;
	wr7 = (double *) w7;
	wr8 = (double *) w8;
	wr9 = (double *) w9;
	wr10 = (double *) w10;
	wr11 = (double *) w11;
	wr12 = (double *) w12;
	wr13 = (double *) w13;
	wr14 = (double *) w14;
	wr15 = (double *) w15;

	for(i = 0 ; i < NX ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lx / 2 + (param.lx * i ) / NX;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.ly / 2 + (param.ly * j ) / NY;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lz / 2 + (param.lz * k ) / NZ;
			}
		}
	}

	DEBUG_END_FUNC;
	return;
}

void finish_common(void) {
	free(kx);
	free(ky);
	free(kz);
	free(kxt);
	free(k2t);
	free(ik2t);
	free(mask);
	free(fld.vx);
	free(fld.vy);
	free(fld.vz);
#ifdef BOUSSINESQ
	free(fld.th);
#endif
#ifdef MHD
	free(fld.bx);
	free(fld.by);
	free(fld.bz);
#endif
	free(w1);
	free(w2);
	free(w3);
	free(w4);
	free(w5);
	free(w6);
	free(w7);
	free(w8);
	free(w9);
	free(w10);
	free(w11);
	free(w12);
	free(w13);
	free(w14);
	free(w15);
	free(x);
	free(y);
	free(z);

	return;
}


/****************************************************/
/**
	Remove the divergence from a 3D field using
	the projector operator:
	q=q-k.q/k^2
	
	@param qx: x component of the field
	@param qy: y component of the field
	@param qz: z component of the field
*/
/****************************************************/
	
void projector( double complex qx[],
			    double complex qy[],
			    double complex qz[]) {
				
	int i;
	double complex q0;
	
	DEBUG_START_FUNC;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		q0 = kxt[i] * qx[i] + ky[i] * qy[i] + kz[i] * qz[i];
		qx[i] = qx[i] - kxt[i] * q0 * ik2t[i];
		qy[i] = qy[i] - ky[i] * q0 * ik2t[i];
		qz[i] = qz[i] - kz[i] * q0 * ik2t[i];
	}
	
	DEBUG_END_FUNC;
	
	return;
}

/*********************************************/
/** Compute the energy of a given field.
	@param q complex array containing the field for which
				  we want the total energy
*/
/*********************************************/

double energy(const double complex q[]) {
	
	int i,j,k;
	double energ_tot;
	
	energ_tot=0.0;
	
	for( i = 0; i < NX_COMPLEX; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k=0; k < NZ_COMPLEX; k++) {
				if( k == 0) 
					// k=0, we have all the modes.
					energ_tot = energ_tot + creal( 0.5 * q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				else
					// k>0, only half of the complex plane is represented.
					energ_tot = energ_tot + creal( q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
			}
		}
	}
//	energ_tot = 0;
	return(energ_tot);
}

/********************************************/
/**
Return the localtime in seconds. Use different
implementation depending on the avaiable
libraries  **/
/********************************************/
double get_c_time(void) {
#ifdef _OPENMP
	// We don't have MPI, but we have OpenMP
	return(omp_get_wtime());
#else
	// We really have nothing...
	clock_t now;
	now = clock();
	
	return( (double) now / ( (double) CLOCKS_PER_SEC));

#endif
}

void chknan_r(double *wri) {
	int i;
	DEBUG_START_FUNC;
	for(i=0; i < 2*NTOTAL_COMPLEX; i++) {
		if(wri[i] != wri[i]) ERROR_HANDLER(ERROR_CRITICAL,"NAN in real array");
	}
	DEBUG_END_FUNC;
	return;
}

void chknan_c(double complex *wri) {
	int i;
	DEBUG_START_FUNC;
	for(i=0; i < NTOTAL_COMPLEX; i++) {
		if(wri[i] != wri[i]) ERROR_HANDLER(ERROR_CRITICAL,"NAN in complex array");
	}
	DEBUG_END_FUNC;
	return;
}
