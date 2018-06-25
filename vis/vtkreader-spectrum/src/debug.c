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



#include <stdlib.h>

#include "common.h"
#include "gfft.h"

#ifdef DEBUG
	
void D_show_field(double complex * field) {
	// Print several informations about the field field
	double complex * df;
	double * dfr;
	double maxfield, minfield, avgfield, avg2field;
	int i;
	
	df = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX );
	if (df == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for df allocation");
	dfr = (double *) df;
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		df[i] =  field[i];
	}
	
	gfft_c2r(df);
	
	maxfield=dfr[0];
	minfield=dfr[0];
	avgfield=0.0;
	avg2field=0.0;
	
	for( i = 0 ; i < NTOTAL_COMPLEX * 2 ; i++) {
		if( dfr[i] > maxfield ) maxfield = dfr[i];
		if( dfr[i] < minfield ) minfield = dfr[i];
		avgfield+=dfr[i];
		avg2field+=dfr[i]*dfr[i];
	}
	
	maxfield=maxfield/ ((double) NTOTAL);
	minfield=minfield/ ((double) NTOTAL);
	avgfield=avgfield/ ((double)  NTOTAL*NTOTAL);
	avg2field=avg2field/ ((double) NTOTAL*NTOTAL*NTOTAL);
	
#ifdef MPI_SUPPORT
	reduce(&maxfield,2);
	reduce(&minfield,3);
	reduce(&avgfield,1);
	reduce(&avgfield,1);
#endif

	MPI_Printf("maxfield= %12e, minfield= %12e, avgfield= %12e, avg2field= %12e\n",maxfield, minfield, avgfield, avg2field);
	
	fftw_free(df);
}

void D_show_all(struct Field fldi) {
	MPI_Printf("   vx:");
	D_show_field(fldi.vx);
	MPI_Printf("   vy:");
	D_show_field(fldi.vy);
	MPI_Printf("   vz:");
	D_show_field(fldi.vz);
#ifdef BOUSSINESQ
	MPI_Printf("   th:");
	D_show_field(fldi.th);
#endif
#ifdef MHD
	MPI_Printf("   bx:");
	D_show_field(fldi.vx);
	MPI_Printf("   by:");
	D_show_field(fldi.vy);
	MPI_Printf("   bz:");
	D_show_field(fldi.vz);
#endif
	return;
}

void debug_start_f(const char ErrorRoutine[], const int line, const char ErrorFile[]) {
	MPI_Printf("Start of %s (Line %d of file %s)\n", ErrorRoutine, line, ErrorFile);
	return;
}

void debug_end_f(const char ErrorRoutine[], const int line, const char ErrorFile[]) {
	MPI_Printf("End of %s (Line %d of file %s)\n", ErrorRoutine, line, ErrorFile);
	return;
}

#endif
	
	