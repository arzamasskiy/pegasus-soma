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

void init_transpose();
void finish_transpose();
double read_transpose_timer();
void transpose_real(const int nxin, const int nyin, const int nzin, const int nproc, double *qin, double *qout);
void transpose_complex_XY(double complex *qin, double complex *qout);
void transpose_complex_YX(double complex *qin, double complex *qout);
void transpose_complex_ZY(double complex *qin, double complex *qout);
void transpose_complex_YZ(double complex *qin, double complex *qout);