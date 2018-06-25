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
#ifdef WITH_SHEAR


double time_shift(double t) {
	double tremap;
#ifdef TIME_DEPENDANT_SHEAR
	tremap = sin(param.omega_shear * t);
#else
	tremap = fmod(t + param.ly / (2.0 * param.shear * param.lx) , param.ly / (param.shear * param.lx)) - param.ly / (2.0 * param.shear * param.lx);
#endif
	return(tremap);
}

void kvolve(const double tremap) {
	int i, j, k;
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif	
	for( i = 0; i < NX_COMPLEX; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kxt[ IDX3D ] = kx[ IDX3D ] + tremap * param.shear * ky[ IDX3D ];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
						       ky[IDX3D] * ky[IDX3D]+
							   kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}

	return;
}

#endif
