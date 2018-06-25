
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#include "common.h"
#include "output.h"
#include "gfft.h"
#include "read_config.h"
#include "shear.h"


void vtk_func() {
	int frame_number;
	double t;
	
	// This is where the user should write his own code
	
	for(frame_number=0; frame_number<400; frame_number++) {
	
		t=read_timevar(frame_number);
		kvolve(time_shift(t));
	
		MPI_Printf("Opening t=%g\n",t);
		
		read_vtk(fld,frame_number,t);
	
		output1Dspectrum(fld, t);
	
	}
	
	return;
}
	
	
int main(int argc, char *argv[]) {

	MPI_Printf("General purpose VTK reader v1.0\n");
	MPI_Printf("Copyright (c) 2004-2009 Geoffroy Lesur (University of Cambridge, UK)\n");
	MPI_Printf("This program comes with ABSOLUTELY NO WARRANTY;\n");
	MPI_Printf("This is free software, and you are welcome to\n");
    MPI_Printf("redistribute it under certain conditions.\n");
	MPI_Printf("Compiled on %s at %s\n",__DATE__ , __TIME__);
	read_config();
	
	printf("param.toutput_flow=%g\n",param.toutput_flow);
	
	MPI_Printf("Initializing...\n");
	init_common();
	init_gfft();
	init_output();

	MPI_Printf("Running.\n");
	
	vtk_func();
	
	finish_output();
	finish_gfft();
	finish_common();
	
	MPI_Printf("Terminated.\n");
	return(0);
}
