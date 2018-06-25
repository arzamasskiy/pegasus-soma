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



#include <string.h>

#include "common.h"
#include <fftw3.h>

// Global ifdef. No transpose if no MPI.

#ifdef MPI_SUPPORT

#ifdef FFTW3_MPI_SUPPORT
#include <fftw3-mpi.h>
#endif

double complex * temp1;
double complex * temp2;
#ifdef FFTW3_MPI_SUPPORT
fftw_plan	plan_t_XY, plan_t_YX;
#endif

double	transpose_timer;

void init_transpose() {
	temp1 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (temp1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for temp1 allocation");
	temp2 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (temp2 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for temp2 allocation");

#ifdef FFTW3_MPI_SUPPORT
#ifdef _OPENMP
	fftw_plan_with_nthreads( nthreads );
#endif
	plan_t_XY = fftw_mpi_plan_many_transpose(NX, NY, (NZ+2), NX/NPROC, NY/NPROC, wr1, wr1, MPI_COMM_WORLD, FFT_PLANNING);
	if (plan_t_XY == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW plan_t_XY plan creation failed");
	
	plan_t_YX = fftw_mpi_plan_many_transpose(NY, NX, (NZ+2), NY/NPROC, NX/NPROC, wr1, wr1, MPI_COMM_WORLD, FFT_PLANNING);
	if (plan_t_YX == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW plan_t_YX plan creation failed");
#endif	

	transpose_timer = 0.0;
	return;
}

void finish_transpose() {
	fftw_free(temp1);
	fftw_free(temp2);
	
#ifdef FFTW3_MPI_SUPPORT
	fftw_destroy_plan(plan_t_XY);
	fftw_destroy_plan(plan_t_YX);
#endif

	return;
}

double read_transpose_timer() {
	return(transpose_timer);
}

#ifdef FFTW3_MPI_SUPPORT
void transpose_complex_XY(double complex *qin, double complex *qout) {
	transpose_timer = transpose_timer - get_c_time();
	fftw_execute_r2r( plan_t_XY, (double *) qin, (double *) qout);
	transpose_timer = transpose_timer + get_c_time();
	return;
}

void transpose_complex_YX(double complex *qin, double complex *qout) {
	transpose_timer = transpose_timer - get_c_time();
	fftw_execute_r2r( plan_t_YX, (double *) qin, (double *) qout);
	transpose_timer = transpose_timer + get_c_time();
	return;
}

#else	
// transpose complex routines are optimized since they are going to be called by ffts routine
// the real transpose is not, since it's more a "convenient" routine. It can however be optimized easely...
void transpose_complex_XY(double complex *qin, double complex *qout) {

// Will transpose qin into qout, 
// qin have dimensions nxin/nproc, nyin, nzin
// Total qin array (all processors) nxin, nyin, nzin
// qout have dimensions nyin/nproc, nxin, nzin
// Total qout array (all processors) nyin, nxin, nzin

	int i,j,k,n;
	const int nxin = NX_COMPLEX;
	const int nyin = NY_COMPLEX;
	const int nzin = NZ_COMPLEX;
	int nproc = NPROC;
	
	int local_nxin = nxin / nproc;
	int local_nyin = nyin / nproc;
	
// First, transpose locally the array in qout (will be erased anyway...)

	transpose_timer = transpose_timer - get_c_time();
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif	
	for(i=0 ; i < local_nxin ; i++) {
		for(j=0 ; j < nyin ; j++) {
			for(k=0 ; k < nzin ; k++) {
				temp1[j*local_nxin*nzin + i*nzin + k] = qin[i*nyin*nzin + j*nzin + k];
			}
		}
	}

			
// Next, MPI the whole thing... Have to be out of place
// Here we could use qin as destination, if qin could be destroyed (might be an interesting optimisation...)
// This step corresponds to an exchange of chuncks of size (local_nyin,local_nxin,nzin)
	
	MPI_Alltoall(temp1, local_nxin*local_nyin*nzin*sizeof(double complex), MPI_BYTE,
				 temp2,   local_nxin*local_nyin*nzin*sizeof(double complex), MPI_BYTE, MPI_COMM_WORLD);
				 
// From here, temp is made of a contiguous array of chunks of size (local_nyin,local_nxin,nzin)
// Which can be seen as a 4D Array of size (nproc,local_nyin,local_nxin,nzin);
// One have to reorder the chunks to get the array right


#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k,n) schedule(static)	
#endif	
	for(i=0 ; i < local_nyin ; i++) {
		for(n=0 ; n < nproc ; n++) {
			for(j=0 ; j < local_nxin ; j++) {
				for(k=0 ; k < nzin ; k++) {
					qout[i*nxin*nzin + (j+n*local_nxin)*nzin + k] = temp2[n*local_nyin*local_nxin*nzin + i*local_nxin*nzin + j*nzin + k];
				}
			}
		}
	}
	transpose_timer = transpose_timer + get_c_time();
    return;
}

void transpose_complex_YX(double complex *qin, double complex *qout) {

// Will transpose qin into qout, 
// qin have dimensions nxin/nproc, nyin, nzin
// Total qin array (all processors) nxin, nyin, nzin
// qout have dimensions nyin/nproc, nxin, nzin
// Total qout array (all processors) nyin, nxin, nzin

	int i,j,k,n;
	const int nxin = NY_COMPLEX;
	const int nyin = NX_COMPLEX;
	const int nzin = NZ_COMPLEX;
	int nproc = NPROC;
	
	int local_nxin = nxin / nproc;
	int local_nyin = nyin / nproc;
	
// First, transpose locally the array in qout (will be erased anyway...)
	
	transpose_timer = transpose_timer - get_c_time();
	
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif	
	for(i=0 ; i < local_nxin ; i++) {
		for(j=0 ; j < nyin ; j++) {
			for(k=0 ; k < nzin ; k++) {
				temp1[j*local_nxin*nzin + i*nzin + k] = qin[i*nyin*nzin + j*nzin + k];
			}
		}
	}

				
// Next, MPI the whole thing... Have to be out of place
// Here we could use qin as destination, if qin could be destroyed (might be an interesting optimisation...)
// This step corresponds to an exchange of chuncks of size (local_nyin,local_nxin,nzin)
	
	MPI_Alltoall(temp1, local_nxin*local_nyin*nzin*sizeof(double complex), MPI_BYTE,
				 temp2,   local_nxin*local_nyin*nzin*sizeof(double complex), MPI_BYTE, MPI_COMM_WORLD);
	
// From here, temp is made of a contiguous array of chunks of size (local_nyin,local_nxin,nzin)
// Which can be seen as a 4D Array of size (nproc,local_nyin,local_nxin,nzin);
// One have to reorder the chunks to get the array right


#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k,n) schedule(static)	
#endif	
	for(i=0 ; i < local_nyin ; i++) {
		for(n=0 ; n < nproc ; n++) {
			for(j=0 ; j < local_nxin ; j++) {
				for(k=0 ; k < nzin ; k++) {
					qout[i*nxin*nzin + (j+n*local_nxin)*nzin + k] = temp2[n*local_nyin*local_nxin*nzin + i*local_nxin*nzin + j*nzin + k];
				}
			}
		}
	}

	transpose_timer = transpose_timer + get_c_time();
	
    return;
}

#endif
void transpose_complex_YZ(double complex *qin, double complex *qout) {
	// this transposition is **out of place**
	int i,j,k;
	const int nxin = NX_COMPLEX;
	const int nyin = NY_COMPLEX;
	const int nzin = NZ_COMPLEX;
	int nproc = NPROC;
	
	int local_nxin = nxin / nproc;

	for(i=0 ; i < local_nxin ; i++) {
		for(j=0 ; j < nyin ; j++) {
			for(k=0 ; k < nzin ; k++) {
				qout[i*nzin*nyin + k*nyin + j] = qin[i*nyin*nzin + j*nzin + k];
			}
		}
	}
	return;
}

void transpose_complex_ZY(double complex *qin, double complex *qout) {
	// this transposition is **out of place**
	int i,j,k;
	const int nxin = NX_COMPLEX;
	const int nyin = NZ_COMPLEX;
	const int nzin = NY_COMPLEX;
	int nproc = NPROC;
	
	int local_nxin = nxin / nproc;
	
// First, transpose locally the array in qout (will be erased anyway...)

	for(i=0 ; i < local_nxin ; i++) {
		for(j=0 ; j < nyin ; j++) {
			for(k=0 ; k < nzin ; k++) {
				qout[i*nzin*nyin + k*nyin + j] = qin[i*nyin*nzin + j*nzin + k];
			}
		}
	}
	return;
}	

void transpose_real(const int nxin, const int nyin, const int nzin, const int nproc, double *qin, double *qout) {

// Will transpose qin into qout, 
// qin have dimensions nxin/nproc, nyin, nzin
// Total qin array (all processors) nxin, nyin, nzin
// qout have dimensions nyin/nproc, nxin, nzin
// Total qout array (all processors) nyin, nxin, nzin

	int i,j,k,n;
	const int local_nxin = nxin / nproc;
	const int local_nyin = nyin / nproc;
	
// Typecast for compatibility
	double *tempc1 = (double *) temp1;
	double *tempc2 = (double *) temp2;
	
// First, transpose locally the array in qout (will be erased anyway...)

	for(i=0 ; i < local_nxin ; i++) {
		for(j=0 ; j < nyin ; j++) {
			for(k=0 ; k < nzin ; k++) {
				// tempc is seen as a (nyin, local_nxin, nzin) array here
				tempc1[j*local_nxin*nzin + i*nzin + k] = qin[i*nyin*nzin + j*nzin + k];
			}
		}
	}
		
// Next, MPI the whole thing... Have to be out of place
// Here we could use qin as destination, if qin could be destroyed (might be an interesting optimisation...)
// This step corresponds to an exchange of chuncks of size (local_nyin,local_nxin,nzin)
	MPI_Alltoall(tempc1, local_nxin*local_nyin*nzin*sizeof(double), MPI_BYTE,
				 tempc2,local_nxin*local_nyin*nzin*sizeof(double), MPI_BYTE, MPI_COMM_WORLD);
				 
// From here, temp is made of a contiguous array of chunks of size (local_nyin,local_nxin,nzin)
// Which can be seen as a 4D Array of size (nproc,local_nyin,local_nxin,nzin);
// One have to reorder the chunks to get the array right
//	for(i=0 ; i < NTOTAL_COMPLEX ; i++) {
//		qout[i]=w1[i];
//	}

	for(n=0 ; n < nproc ; n++) {
		for(i=0 ; i < local_nyin ; i++) {
			for(j=0 ; j < local_nxin ; j++) {
				for(k=0 ; k < nzin ; k++) {
					qout[i*nxin*nzin + (j+n*local_nxin)*nzin + k] = tempc2[n*local_nyin*local_nxin*nzin + i*local_nxin*nzin + j*nzin + k];
				}
			}
		}
	}

    return;
}

#endif
