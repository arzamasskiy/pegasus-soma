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



#ifndef __SNOOPY_H__
#define __SNOOPY_H__

#include <math.h>
#include <complex.h>
#include <fftw3.h>


#ifdef _OPENMP
#include <omp.h>
#endif

#include "gvars.h"

#define		NTOTAL			NX * NY * NZ	/**< Total number of grid points over all the processors */

#define		NX_COMPLEX		NX				/**< Number of complex point in X over all the processes (Fourier space) */
#define		NY_COMPLEX		NY				/**< Number of complex point in Y over all the processes (Fourier space) */
#define		NZ_COMPLEX		(NZ / 2 + 1)	/**< Number of complex point in Z over all the processes (Fourier space) */

#define		NTOTAL_COMPLEX	(NX_COMPLEX * NY_COMPLEX * NZ_COMPLEX)	/**< Number of complex points in one MPI process */

#define		IDX3D			(k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i)  /**< General wrapper for 3D Arrays */

#define		MPI_Printf		printf

// Structures

struct Field {
	double complex *vx;
	double complex *vy;
	double complex *vz;
#ifdef BOUSSINESQ
	double complex *th;
#endif
#ifdef MHD
	double complex *bx;
	double complex *by;
	double complex *bz;
#endif
};

struct Parameters {
	// Physics Parameters
	double lx;				/**< Box length in X*/
	double ly;				/**< Box length in Y*/
	double lz;				/**< Box length in Z*/

	double shear;			/**< Shear rate (only when WITH_SHEAR is on) */
	
	double omega_shear;		/**< Pulsation of the time dependant shear (only when WITH_SHEAR and TIME_DEPENDANT_SHEAR is on) */
	
	int    antialiasing;		/**< 2/3 Antialisaing rule. Could be removed if you assume is unimportant (untested feature). */

	// Output parameters
	double toutput_time;		/**< Time between two outputs in the timevar file */
	double toutput_flow;		/**< Time between two snapshot outputs */
};

#endif

