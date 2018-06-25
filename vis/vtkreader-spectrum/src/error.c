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
#include <stdio.h>
#include "common.h"
#include "error.h"

void error_h (const int ErrorType,
			  const char ErrorMessage[], 
			  const char ErrorRoutine[], 
			  const int line, 
			  const char Filename[] ) {
	MPI_Printf("*************************************************\n");
	if ( ErrorType == ERROR_WARNING )
		MPI_Printf("Warning in ");
	if ( ErrorType == ERROR_NONCRITICAL )
		MPI_Printf("Non Critical Error in ");
	if ( ErrorType == ERROR_CRITICAL )
		MPI_Printf("Critical Error in ");
	MPI_Printf("File: ");
	MPI_Printf(Filename);
	MPI_Printf(" Routine: ");
	MPI_Printf(ErrorRoutine);
		MPI_Printf(" Line: %d \n",line);
	MPI_Printf(ErrorMessage);
	
	if(ErrorType == ERROR_CRITICAL) {
		MPI_Printf("\n Terminating.\n");
		MPI_Printf("*************************************************\n");
#ifdef MPI_SUPPORT
		MPI_Abort(MPI_COMM_WORLD,1);
#endif
		exit(1);
	}
	else {
		MPI_Printf("\nResuming execution...\n");
		MPI_Printf("*************************************************\n");
	}
	return;
}
