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



/* Error.h: Error handling, useful for debugging purposes */

#define		ERROR_HANDLER(ERROR_TYPE, ERROR_MESSAGE)		error_h(ERROR_TYPE, ERROR_MESSAGE, __func__ , __LINE__, __FILE__)
#define		ERROR_WARNING		1
#define		ERROR_NONCRITICAL	2
#define		ERROR_CRITICAL		3

void error_h (const int ErrorType,
			  const char ErrorMessage[], 
			  const char ErrorRoutine[], 
			  const int line, 
			  const char Filename[] );
