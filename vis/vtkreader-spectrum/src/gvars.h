#ifndef _GVARS_
#define _GVARS_

#define		NX				64			/**< X Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NY				64			/**< Y Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NZ				64			/**< Z Dimension in real space. */

#define		MHD						/**< Uncomment to activate MHD*/

//#define		BOUSSINESQ				/**< Uncomment to activate Boussinesq */
//#define		VERTSTRAT				/**< Vertical stratification. Otherwise, Boussinesq stratification is in X */

#define		WITH_ROTATION

#define		WITH_SHEAR				/**< Uncomment to activate mean SHEAR */
//#define		TIME_DEPENDANT_SHEAR	/**< Enable Time dependant shear */

//#define		FORCING					/**< Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) */

#define		FFT_PLANNING	FFTW_MEASURE  /**< can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). Measure leads to longer initialisation of fft routines */



#endif
