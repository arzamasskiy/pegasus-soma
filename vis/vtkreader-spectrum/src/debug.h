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



#ifdef DEBUG

void D_show_all(struct Field fldi);
void D_show_field(double complex * field);
void debug_start_f(const char ErrorRoutine[], const int line, const char ErrorFile[]);
void debug_end_f(const char ErrorRoutine[], const int line, const char ErrorFile[]);

#define		DEBUG_START_FUNC		debug_start_f(__func__, __LINE__, __FILE__)
#define		DEBUG_END_FUNC			debug_end_f(__func__, __LINE__, __FILE__)

#else

#define		DEBUG_START_FUNC
#define		DEBUG_END_FUNC

#endif
