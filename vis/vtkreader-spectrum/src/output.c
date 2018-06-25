#include "common.h"
#include "debug.h"
#include "shear.h"

#define MAX_N_BIN					10000

#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI)
#define	OUTPUT_SPECTRUM_FILENAME	"spectrum_reader.dat"

#define KS_REF			32.0
#define DKS_REF			( OUTPUT_SPECTRUM_K_BIN / 2.0 )

#ifdef WITH_SHEAR

double complex		*w1d, *w2d;						/** 1D arrays used by the remap methods */

fftw_plan	fft_1d_forward;								/**< 1D FFT transforms. Used by remap routines.*/
fftw_plan	fft_1d_backward;							/**< 1D FFT transforms. Used by remap routines.*/


/***********************************************************/
/** 
	Remap a real field from the current sheared frame to the classical 
	cartesian frame. It could be more optimized, but since it is used
	in an output routine, I don't think such an optimization is worth doing.
	
	@param wri Real array to be remapped. The remapped array will be stored here.
	@param t Current time of the simulation (used to compute what the sheared frame is).
	
*/
/***********************************************************/

void remap_output_back(	double wri[], 
					const double t) {
					
	int i,j,k;
	double tvelocity;
	double tremap;
	complex double wexp;
	complex double phase;
	
	DEBUG_START_FUNC;
	
#ifdef TIME_DEPENDANT_SHEAR
	tremap = time_shift(t);
	tvelocity = 0.0;
#else	
	tremap = time_shift(t);
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
#endif	
	
	
	for( i = 0 ; i < NX ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			for( j = 0 ; j < NY ; j++) {
				w1d[j] = wri[k + j * (NZ + 2) + NY * (NZ + 2) * i];
			}
		
			// Transform w1d, which will be stored in w2d
			fftw_execute(fft_1d_forward);
					
			for( j = 0 ; j < NY ; j++) {
			// advection phase = ky*
#ifdef TIME_DEPENDANT_SHEAR
				// There is no proper remap in this case
				phase = -(double complex) ((2.0 * M_PI) / param.ly * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i) / (double) NX - 0.5 ) * tremap ) * param.lx * param.shear );
#else
				phase = -(double complex) ((2.0 * M_PI) / param.ly * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i) / (double) NX ) * tremap - tvelocity / 2.0 ) * param.lx * param.shear);
#endif
				if(phase != phase) printf("Error\n");
				
				wexp = cexp( I * phase);
									
				w2d[ j ] = w2d[ j ] * wexp;
			}
			
			fftw_execute(fft_1d_backward);
		
			for( j = 0 ; j < NY ; j++) {
				wri[k + j * (NZ + 2) + NY * (NZ + 2) * i] = w1d[j] / NY;
			}
		}
	}
	
	DEBUG_END_FUNC;
	
	return;
}

/***********************************************************/
/** 
	Remap a real field from the current sheared frame to the classical 
	cartesian frame. It could be more optimized, but since it is used
	in an output routine, I don't think such an optimization is worth doing.
	
	@param wri Real array to be remapped. The remapped array will be stored here.
	@param t Current time of the simulation (used to compute what the sheared frame is).
	
*/
/***********************************************************/

void remap_output(	double wri[], 
					const double t) {
					
	int i,j,k;
	double tvelocity;
	double tremap;
	complex double wexp;
	complex double phase;
	
	DEBUG_START_FUNC;
	
#ifdef TIME_DEPENDANT_SHEAR
	tremap = time_shift(t);
	tvelocity = 0.0;
#else	
	tremap = time_shift(t);
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
#endif	
	
	
	for( i = 0 ; i < NX ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			for( j = 0 ; j < NY ; j++) {
				w1d[j] = wri[k + j * (NZ + 2) + NY * (NZ + 2) * i];
			}
		
			// Transform w1d, which will be stored in w2d
			fftw_execute(fft_1d_forward);
					
			for( j = 0 ; j < NY ; j++) {
			// advection phase = ky*
#ifdef TIME_DEPENDANT_SHEAR
				// There is no proper remap in this case
				phase = (double complex) ((2.0 * M_PI) / param.ly * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i) / (double) NX - 0.5 ) * tremap ) * param.lx * param.shear );
#else
				phase = (double complex) ((2.0 * M_PI) / param.ly * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i) / (double) NX ) * tremap - tvelocity / 2.0 ) * param.lx * param.shear);
#endif
				wexp = cexp( I * phase);
									
				w2d[ j ] = w2d[ j ] * wexp;
			}
			
			fftw_execute(fft_1d_backward);
		
			for( j = 0 ; j < NY ; j++) {
				wri[k + j * (NZ + 2) + NY * (NZ + 2) * i] = w1d[j] / NY;
			}
		}
	}
	
	DEBUG_END_FUNC;
	
	return;
}

#endif




/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so, 
    will convert the big endian input into a little endian output. 
	@param in_number floating point number to be converted in xxx endian */
/* *************************************************************************** */

float fix_endian(float in_number)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
	
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap)
    {
		unsigned char *bytes = (unsigned char*) &in_number;
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
	return(in_number);
}

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so, 
    it will force the data to be big-endian. 
	@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

float big_endian(float in_number)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
	
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap)
    {
		unsigned char *bytes = (unsigned char*) &in_number;
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
	return(in_number);
}

/***********************************************************/
/** 
	Output a floating point big endian array from a complex array.
	Use Forran output format (required by VTK). This routine is essentially
	useful for VTK files (hence its name...).
 	
	@param ht File handler in which the data has to be written
	@param wi Complex array containing the field to be written. This array is transformed into real space
	and remapped (if SHEAR is present) before being written.
	@param t Current time of the simulation.
	
*/
/***********************************************************/
void write_field(FILE * ht, double complex wi[], const double t) {
	// Write the data in the file handler *ht
	int i,j,k;
	float q0;

	DEBUG_START_FUNC;

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = wi[i];
	}
	
	gfft_c2r(w1);

	for( i = 0 ; i < 2 * NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
	}
	
#ifdef WITH_SHEAR
	remap_output(wr1,t);
#endif

	for( k = 0 ; k < NZ; k++) {
		for( j = 0; j < NY; j++) {		
			for( i = 0; i < NX; i++) {
				q0 = big_endian( (float) wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ] );
				fwrite(&q0, sizeof(float), 1, ht);
			}
		}
	}

	DEBUG_END_FUNC;
	
	return;
}

/***********************************************************/
/** 
	Read binary data from a VTK file
 	
	@param ht File handler in which the data has to be written
	@param wi Complex array containing the field to be written. This array is transformed into real space
	and remapped (if SHEAR is present) before being written.
	@param t Current time of the simulation.
	
*/
/***********************************************************/
void read_field(FILE * ht, double complex wi[], const double t) {
	// Write the data in the file handler *ht
	int i,j,k;
	float q0;

	DEBUG_START_FUNC;

	for( k = 0 ; k < NZ; k++) {
		for( j = 0; j < NY; j++) {
			for( i = 0; i < NX; i++) {
				fread(&q0, sizeof(float), 1, ht);
				wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ] = (double) fix_endian( q0 );
				if(fix_endian(q0) != fix_endian(q0)) ERROR_HANDLER(ERROR_CRITICAL, "Nan read from file");
			}
		}
	}
	
#ifdef WITH_SHEAR
	remap_output_back(wr1,t);
#endif
	
	chknan_r(wr1);
	gfft_r2c(wr1);
	chknan_c(w1);
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		wi[i] = w1[i];
	}
	
	DEBUG_END_FUNC;
	
	return;
}

// Geo's personnal VTK writer, using structured data points	
/***********************************************************/
/** 
	Read a legacy VTK file readable by Paraview. This routine
	will output all the variables in files data/v****.vtk.
	
	@param n Number of the file in which the output will done.
	@param t Current time of the simulation.
*/
/***********************************************************/

void read_vtk(struct Field fldo, const int n, double t) {
	FILE *ht = NULL;
	char  filename[50];
	int num_remain_field;
	char dummy_line[200];
	int nx_read, ny_read, nz_read, num_field_read;
	
	DEBUG_START_FUNC;

	sprintf(filename,"bin/fhshear.%04i.vtk",n);

	ht=fopen(filename,"r");
		
	fgets(dummy_line, 200, ht);			//vtk DataFile Version XXX
	fgets(dummy_line, 200, ht);     // CONSERVED vars at time...
	fgets(dummy_line, 200, ht);			// BINARY
	fgets(dummy_line, 200, ht);			// DATASET STRUCTURED POINTS
	fgets(dummy_line, 200, ht);			// DIMENSIONS %d %d %d\n
	
	// Check dimensions
	sscanf(dummy_line,"%*s %d %d %d", &nx_read, &ny_read, &nz_read);
	if(NX != nx_read) ERROR_HANDLER(ERROR_CRITICAL,"nx_read is different from NX");
	if(NY != ny_read) ERROR_HANDLER(ERROR_CRITICAL,"ny_read is different from NY");
	if(NZ != nz_read) ERROR_HANDLER(ERROR_CRITICAL,"nz_read is different from NZ");
	
	fgets(dummy_line, 200, ht);			// ORIGIN %g %g %g
	fgets(dummy_line, 200, ht);			// SPACING %g %g %g
	fgets(dummy_line, 200, ht);			// CELL_DATA
	fgets(dummy_line, 200, ht);			// VECTOR cell_centered_B float
//	fgets(dummy_line, 200, ht);			// LOOKUP_TABLE default
	
	read_field(ht, fldo.vx, t);
	
	num_remain_field = 2;		// Number of remaining field to be written
#ifdef MHD
	num_remain_field += 3;
#endif
#ifdef BOUSSINESQ
	num_remain_field += 1;
#endif
	
	fgets(dummy_line, 200, ht);			// FIELD FieldData %d\n
	sscanf(dummy_line,"%*s %*s %d", &num_field_read);
	if(num_field_read < num_remain_field) ERROR_HANDLER(ERROR_CRITICAL, "num_remain_field is wrong");
	
	fgets(dummy_line, 200, ht);			// xx 1 %d float
	read_field(ht, fldo.vy, t);
	
	fgets(dummy_line, 200, ht);			// xx 1 %d float
	read_field(ht, fldo.vz, t);

#ifdef MHD
	fgets(dummy_line, 200, ht);			// xx 1 %d float
	read_field(ht, fldo.bx, t);
	
	fgets(dummy_line, 200, ht);			// xx 1 %d float
	read_field(ht, fldo.by, t);
	
	fgets(dummy_line, 200, ht);			// xx 1 %d float
	read_field(ht, fldo.bz, t);
#endif
#ifdef BOUSSINESQ	
	fgets(dummy_line, 200, ht);			// xx 1 %d float
	read_field(ht, fldo.th, t);
#endif

	fclose(ht);
	
	DEBUG_END_FUNC;
	
	return;
}

// Geo's personnal VTK writer, using structured data points	
/***********************************************************/
/** 
	Output a legacy VTK file readable by Paraview. This routine
	will output all the variables in files data/v****.vtk.
	
	@param n Number of the file in which the output will done.
	@param t Current time of the simulation.
*/
/***********************************************************/

void write_vtk(double complex *wi, const char filename[], const char fieldname[], const double t) {
	FILE *ht = NULL;
	
	DEBUG_START_FUNC;

	ht=fopen(filename,"w");
	
	fprintf(ht, "# vtk DataFile Version 2.0\n");
	fprintf(ht, "Created by VTKreader\n");
	fprintf(ht, "BINARY\n");
	fprintf(ht, "DATASET STRUCTURED_POINTS\n");
	fprintf(ht, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
	fprintf(ht, "ORIGIN %g %g %g\n", -param.lx/2.0, -param.ly/2.0, -param.lz/2.0);
	fprintf(ht, "SPACING %g %g %g\n", param.lx/NX, param.ly/NY, param.lz/NZ);
	
	// Write the primary scalar (f***ing VTK legacy format...)
	fprintf(ht, "POINT_DATA %d\n",NX*NY*NZ);
	fprintf(ht, "SCALARS %s float\n", fieldname);
	fprintf(ht, "LOOKUP_TABLE default\n");

	write_field(ht,wi,t);

	fclose(ht);
	
	DEBUG_END_FUNC;
	
	return;
	
}

double read_timevar(const int n) {
	FILE *ht = NULL;
	int nloop;
	double toutput, t;
	char dummy_line[5000];
	
	toutput = ((double) n) * param.toutput_flow;		// rough estimate of toutput
	t=-1.0;									// t initial
	
	ht=fopen("timevar","r");
	nloop=0;
	while(t < toutput) {
		nloop++;
		fgets(dummy_line, 5000, ht);
		sscanf(dummy_line, "%lf",&t);
		if(nloop>32000) ERROR_HANDLER(ERROR_CRITICAL,"can't find the relevant time in timevar");
	}
	
	fclose(ht);
	
	DEBUG_END_FUNC;
	
	return(t);
}

/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_spectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
					if( k == 0) 
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + creal( 2.0 * wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
	fprintf(ht,"%08e\t", ti);
	
	for( i = 0; i < nbin; i++) 
		fprintf(ht,"%08e\t", spectrum[i]);
	
	fprintf(ht,"\n");
	
	fclose(ht);


		
	DEBUG_END_FUNC;
	return;
}

/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/

void output1Dspectrum(const struct Field fldi, const double ti) {
	int i;
	double mask_in;
	
	DEBUG_START_FUNC;
	
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	// V,B and theta spectrums
	write_spectrum(fldi.vx, fldi.vx, ti);
	write_spectrum(fldi.vy, fldi.vy, ti);
	write_spectrum(fldi.vz, fldi.vz, ti);
	
#ifdef MHD
	write_spectrum(fldi.bx, fldi.bx, ti);
	write_spectrum(fldi.by, fldi.by, ti);
	write_spectrum(fldi.bz, fldi.bz, ti);
#else
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif

#ifdef BOUSSINESQ
	write_spectrum(fldi.th, fldi.th, ti);
#else
	write_spectrum(w1, w1, ti);
#endif

	// Transport spectrums
	write_spectrum(fldi.vx,fldi.vy, ti);
#ifdef MHD
	write_spectrum(fldi.bx,fldi.by, ti);
#else
	write_spectrum(w1, w1, ti);
#endif

// **********************************************************************************	
// **********************************************************************************
// Compute what's happening in terms of flux
// **********************************************************************************	
// **********************************************************************************


	// Transfer spectrums
	// Kinetic energy transfer
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w5[i] + kz[i] * w9[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );
	}

	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);

#ifdef MHD

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the Faraday tensor components in w7-w15...

	// Let's compute th u*b tensor in wr7-wr15
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf involved in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);
	
	// The first term is the -\partial_j u_j b_i term (field transport)
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = -I * mask[i] * (kxt[i] * w7[i] + ky[i] * w10[i] + kz[i] * w13[i]);
		w2[i] = -I * mask[i] * (kxt[i] * w8[i] + ky[i] * w11[i] + kz[i] * w14[i]);
		w3[i] = -I * mask[i] * (kxt[i] * w9[i] + ky[i] * w12[i] + kz[i] * w15[i]);
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
	// Second term is \partial_j u_i b_j term (stretching term)
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * mask[i] * (kxt[i] * w7[i]  + ky[i] * w8[i]  + kz[i] * w9[i]);
		w2[i] = I * mask[i] * (kxt[i] * w10[i] + ky[i] * w11[i] + kz[i] * w12[i]);
		w3[i] = I * mask[i] * (kxt[i] * w13[i] + ky[i] * w14[i] + kz[i] * w15[i]);
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
	// Let's do the Lorentz Force
	// We already have (bx,by,bz) in w4-w6. No need to compute them again...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr2[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr3[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}

	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);


	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = I * mask[i] * (kxt[i] * w1[i] + ky[i] * w7[i] + kz[i] * w8[i]);
		w5[i] = I * mask[i] * (kxt[i] * w7[i] + ky[i] * w2[i] + kz[i] * w9[i]);
		w6[i] = I * mask[i] * (kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w3[i]);
	}
	
	write_spectrum(fldi.vx, w4, ti);
	write_spectrum(fldi.vy, w5, ti);
	write_spectrum(fldi.vz, w6, ti);
	
	// Compute the shear fluxes
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = param.shear * ky[i] * kxt[i] * pow(ik2t[i], 0.5) / OUTPUT_SPECTRUM_K_BIN * fldi.vx[i];
		w2[i] = param.shear * ky[i] * kxt[i] * pow(ik2t[i], 0.5) / OUTPUT_SPECTRUM_K_BIN * fldi.vy[i];
		w3[i] = param.shear * ky[i] * kxt[i] * pow(ik2t[i], 0.5) / OUTPUT_SPECTRUM_K_BIN * fldi.vz[i];
		w4[i] = param.shear * ky[i] * kxt[i] * pow(ik2t[i], 0.5) / OUTPUT_SPECTRUM_K_BIN * fldi.bx[i];
		w5[i] = param.shear * ky[i] * kxt[i] * pow(ik2t[i], 0.5) / OUTPUT_SPECTRUM_K_BIN * fldi.by[i];
		w6[i] = param.shear * ky[i] * kxt[i] * pow(ik2t[i], 0.5) / OUTPUT_SPECTRUM_K_BIN * fldi.bz[i];
	}
	
	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);
	write_spectrum(fldi.bx, w4, ti);
	write_spectrum(fldi.by, w5, ti);
	write_spectrum(fldi.bz, w6, ti);
	
// **********************************************************************************	
// **********************************************************************************
// Test locality
// **********************************************************************************	
// **********************************************************************************
	
	//*********************************************
	// T_uu transfert
	//*********************************************
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		if( (k2t[i] < (KS_REF + DKS_REF) * (KS_REF + DKS_REF)) && (k2t[i] > (KS_REF - DKS_REF)*(KS_REF - DKS_REF))) {
			mask_in=1.0;
		}
		else {
			mask_in=0.0;
		}
		
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
		w4[i] = mask_in * fldi.vx[i];
		w5[i] = mask_in * fldi.vy[i];
		w6[i] = mask_in * fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w10[i] + kz[i] * w13[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w11[i] + kz[i] * w14[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w9[i] + ky[i] * w12[i] + kz[i] * w15[i] );
	}

	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);

	//*********************************************
	// T_bb transfert
	//*********************************************
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		if( (k2t[i] < (KS_REF + DKS_REF) * (KS_REF + DKS_REF)) && (k2t[i] > (KS_REF - DKS_REF)*(KS_REF - DKS_REF))) {
			mask_in=1.0;
		}
		else {
			mask_in=0.0;
		}
		
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
		w4[i] = mask_in * fldi.bx[i];
		w5[i] = mask_in * fldi.by[i];
		w6[i] = mask_in * fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w10[i] + kz[i] * w13[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w11[i] + kz[i] * w14[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w9[i] + ky[i] * w12[i] + kz[i] * w15[i] );
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
	//*********************************************
	// -T_bu transfert (it's actually T_ub computed here, see Pouquet paper)
	//*********************************************
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		if( (k2t[i] < (KS_REF + DKS_REF) * (KS_REF + DKS_REF)) && (k2t[i] > (KS_REF - DKS_REF)*(KS_REF - DKS_REF))) {
			mask_in=1.0;
		}
		else {
			mask_in=0.0;
		}
		
		w1[i] = fldi.bx[i];
		w2[i] = fldi.by[i];
		w3[i] = fldi.bz[i];
		w4[i] = mask_in * fldi.vx[i];
		w5[i] = mask_in * fldi.vy[i];
		w6[i] = mask_in * fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w10[i] + kz[i] * w13[i] );
		w2[i] =  I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w11[i] + kz[i] * w14[i] );
		w3[i] =  I * mask[i] * (
					kxt[i] * w9[i] + ky[i] * w12[i] + kz[i] * w15[i] );
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);

//*********************************************
	// -T_ub transfert (actually T_bu, for same reason as above)
//*********************************************
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		if( (k2t[i] < (KS_REF + DKS_REF) * (KS_REF + DKS_REF)) && (k2t[i] > (KS_REF - DKS_REF)*(KS_REF - DKS_REF))) {
			mask_in=1.0;
		}
		else {
			mask_in=0.0;
		}
		
		w1[i] = fldi.bx[i];
		w2[i] = fldi.by[i];
		w3[i] = fldi.bz[i];
		w4[i] = mask_in * fldi.bx[i];
		w5[i] = mask_in * fldi.by[i];
		w6[i] = mask_in * fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
	gfft_r2c_t(wr12);
	gfft_r2c_t(wr13);
	gfft_r2c_t(wr14);
	gfft_r2c_t(wr15);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w10[i] + kz[i] * w13[i] );
		w2[i] =  I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w11[i] + kz[i] * w14[i] );
		w3[i] =  I * mask[i] * (
					kxt[i] * w9[i] + ky[i] * w12[i] + kz[i] * w15[i] );
	}

	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);


	// Helicity spectrums
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * ik2t[i] * (ky[i] * fldi.bz[i] - kz[i] * fldi.by[i] );
		w2[i] = I * ik2t[i] * (kz[i] * fldi.bx[i] - kxt[i]* fldi.bz[i] );
		w3[i] = I * ik2t[i] * (kxt[i]* fldi.by[i] - ky[i] * fldi.bx[i] );
	}
	
	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
#else
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	
#endif

	
	DEBUG_END_FUNC;
	
	return;
}

/**********************************************************/
/**
	Initialise the 1D spectrum output routine, used
	to output the spectrum
	This routine print the mode ks in the first line
	It also counts the number of mode in each shell and 
	output it in the second line of OUTPUT_SPECTRUM_FILENAME
*/
/*********************************************************/
void init1Dspectrum() {
	int i,j,k,m;
	int nbin;
	FILE * ht;
	double spectrum[ MAX_N_BIN ];
	
	DEBUG_START_FUNC;
	
	ht = fopen(OUTPUT_SPECTRUM_FILENAME,"w");
	
	// Count the number of bins
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( m=0; m < nbin; m++) 
		fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
	fprintf(ht,"\n");
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin)
					spectrum[ m ] = spectrum[ m ] + 1.0;
			}
		}
	}

	for( i = 0; i < nbin; i++) 
		fprintf(ht,"%08e\t", spectrum[i]);
	
	fprintf(ht,"\n");
	
	fclose(ht);
	
	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	Initialize the output variables. Should be called only in the begining .
*/
/**************************************************************************************/

void init_output() {

	DEBUG_START_FUNC;
	
#ifdef WITH_SHEAR
// Initialize 1D arrays for remaps
	w1d = (double complex *) fftw_malloc( sizeof(double complex) * NY);
	if (w1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1d allocation");
	
	w2d = (double complex *) fftw_malloc( sizeof(double complex) * NY);
	if (w2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2d allocation");

// FFT plans (we use dummy arrays since we use the "guru" interface of fft3 in the code)
// The in place/ out of place will be set automatically at this stage

#ifdef _OPENMP
	fftw_plan_with_nthreads( 1 );
#endif

	fft_1d_forward = fftw_plan_dft_1d(NY, w1d, w2d, FFTW_FORWARD, FFT_PLANNING);
	if (fft_1d_forward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D forward plan creation failed");
	
	fft_1d_backward = fftw_plan_dft_1d(NY, w2d, w1d, FFTW_BACKWARD, FFT_PLANNING);
	if (fft_1d_backward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D backward plan creation failed");
	
#endif	

	init1Dspectrum();
	
	DEBUG_END_FUNC;
	
	return;
}
/****************************************************************************/
/**
	Free the variables used by the output routines
*/
/****************************************************************************/
void finish_output() {
	DEBUG_START_FUNC;
	
#ifdef WITH_SHEAR
	fftw_free(w1d);
	fftw_free(w2d);
	
	fftw_destroy_plan(fft_1d_forward);
	fftw_destroy_plan(fft_1d_backward);	
#endif
	
	DEBUG_END_FUNC;
		
	return;
}

