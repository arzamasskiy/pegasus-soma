#include "../copyright.h"
/*============================================================================*/
/*! \file dist_func.c
 *  \brief Contains initialization for various distribution functions, as 
 *         well as routines that compute delta-f weights and set moments
 *         of the background distribution function for the delta-f method.
 *
 * PURPOSE: Contains initialization for various distribution functions, as 
 *          well as routines that compute delta-f weights and set moments
 *          of the background distribution function for the delta-f method.
 *
 * NOTE FOR BI-MAXWELLIAN FUNCTIONS: parallel direction is assumed to be x1
 *
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * 
 * - getdf_maxwell()
 * - getdf_bimaxwell()
 * - setbg_maxwell()
 * - setbg_bimaxwell()
 * - initdf_static()
 * - initdf_maxwell()
 * - initdf_bimaxwell()
 * 
 * PRIVATE FUNCTION PROTOTYPES:
 * - ran2()
 * - rangaussian()
 *
 *============================================================================*/
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../pegasus.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   ran2()        - returns a uniformly distributed random number
 *   rangaussian() - returns a normally distributed random number
 *============================================================================*/

static double ran2(long int *idum);
static double ran_gaussian(long int *idum);

#ifdef DELTA_F
/*============================================================================*/
/*---------------------------GET DELTA-F WEIGHTS-------------------------------
 *
 * getdf_maxwell()
 * getdf_bimaxwell()
 */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void getdf_maxwell(GrainS *gr, Real *df)
 *  \brief Computes F_0(t) for delta-f weights
 */
void getdf_maxwell(GrainS *gr, Real *df)
{
  Real vsq;
  
  vsq = SQR(gr->v1)+SQR(gr->v2)+SQR(gr->v3);
  *df = exp(-vsq/beta);
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real getdf_bimaxwell(GrainS *gr, Real *df)
 *  \brief Computes F_0(t) for delta-f weights w/ B along x1
 */
void getdf_bimaxwell(GrainS *gr, Real *df)
{
  Real vsq_prl,vsq_prp;
  
  vsq_prl = SQR(gr->v1);
  vsq_prp = SQR(gr->v2)+SQR(gr->v3);
  *df = exp(-vsq_prl/beta_prl)*exp(-vsq_prp/beta_prp);
  
  return;
}

/*============================================================================*/
/*--------------------------SET DELTA-F BACKGROUND-----------------------------
 *
 * setbg_maxwell()
 * setbg_bimaxwell()
 */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void setbg_maxwell(GridS *pG)
 *  \brief Sets (analytic) moments of the background Maxwellian
 */
void setbg_maxwell(GridS *pG)
{
  int i,j,k;
  GPCouple *pw;
  
  for (i=pG->is; i<=pG->ie; i++)
    for (j=pG->js; j<=pG->je; j++)
      for (k=pG->ks; k<=pG->ke; k++) {
        pw = (&(pG->Bkgrd[k][j][i]));
        pw->grid_d   = 1.0;
        pw->grid_M1  = 0.0;
        pw->grid_M2  = 0.0;
        pw->grid_M3  = 0.0;
        pw->grid_p11 = 0.5*beta;
        pw->grid_p12 = 0.0;
        pw->grid_p13 = 0.0;
        pw->grid_p22 = 0.5*beta;
        pw->grid_p23 = 0.0;
        pw->grid_p33 = 0.5*beta;
      }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void setbg_bimaxwell(GridS *pG)
 *  \brief Sets (analytic) moments of the background bi-Maxwellian w/ B along x1
 */
void setbg_bimaxwell(GridS *pG)
{
  int i,j,k;
  GPCouple *pw;
  
  for (i=pG->is; i<=pG->ie; i++)
    for (j=pG->js; j<=pG->je; j++)
      for (k=pG->ks; k<=pG->ke; k++) {
        pw = (&(pG->Bkgrd[k][j][i]));
        pw->grid_d   = 1.0;
        pw->grid_M1  = 0.0;
        pw->grid_M2  = 0.0;
        pw->grid_M3  = 0.0;
        pw->grid_p11 = 0.5*beta_prl;
        pw->grid_p12 = 0.0;
        pw->grid_p13 = 0.0;
        pw->grid_p22 = 0.5*beta_prp;
        pw->grid_p23 = 0.0;
        pw->grid_p33 = 0.5*beta_prp;
      }
  
  return;
}

#endif

/*============================================================================*/
/*------------------INITIALIZE FULL-F DISTRIBUTION FUNCTIONS-------------------
 *
 * initdf_static()
 * initdf_maxwell()
 * initdf_bimaxwell()
 */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn Real initdf_static(Real **vp, const long int Np, const long int NpTot)
 *  \brief Initializes distribution function f = \delta(v)
 */
void initdf_static(long int *idum, Real **vp, 
                   const long int Np, const long int NpTot)
{
  int i;
  long int q;
  
  for (i=0; i<3; i++) {
    for (q=0; q<Np; q++) {
      vp[i][q] = 0.0;
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real initdf_maxwell(long int *idum, Real **vp,
 *                          const long int Np, const long int NpTot)
 *  \brief Initializes Maxwellian distribution function f ~ exp(-v^2/beta)
 */

void initdf_maxwell(long int *idum, Real **vp,
                    const long int Np, const long int NpTot)
{
  int i,mpierr;
  long int q;
  double vsqsum, vcm[3];
  Real lambda;
#ifdef MPI_PARALLEL
  double gvsqsum, gvcm[3];
#endif
  
  vsqsum = 0.0;
  for (i=0; i<3; i++)
    vcm[i] = 0.0;
  
  for (q=0; q<Np; q++)
  {
    vp[0][q] = ran_gaussian(idum);
    vp[1][q] = ran_gaussian(idum);
    vp[2][q] = ran_gaussian(idum);
    vcm[0]  += vp[0][q];
    vcm[1]  += vp[1][q];
    vcm[2]  += vp[2][q];
  }

#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(vcm,gvcm,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vcm[0] = gvcm[0]; vcm[1] = gvcm[1]; vcm[2] = gvcm[2];
#endif
  
  for (i=0; i<3; i++) {
    vcm[i] /= NpTot;
    for (q=0; q<Np; q++) {
      vp[i][q] -= vcm[i];
      vsqsum   += SQR(vp[i][q]);
    }
  }
  
#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(&vsqsum,&gvsqsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vsqsum = gvsqsum;
#endif
  
  lambda = sqrt(0.5*((double)NpTot)*beta/vsqsum) * sqrt(3.0);
  
  for (i=0; i<3; i++) {
    for (q=0; q<Np; q++) {
      vp[i][q] *= lambda;
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real initdf_bimaxwell(long int *idum, Real **vp,
 *                            const long int Np, const long int NpTot)
 *  \brief Initializes bi-Maxwellian distribution function
 *         f ~ exp[-v^2_prp/beta_prp - v^2_||/beta_||] w/ B along x1
 */
void initdf_bimaxwell(long int *idum, Real **vp,
                      const long int Np, const long int NpTot)
{
  int i,mpierr;
  long int q;
  double vsqsum_prp, vsqsum_prl, vcm[3];
  Real lambda_prp, lambda_prl;
#ifdef MPI_PARALLEL
  double gvsqsum, gvcm[3];
#endif
  
  vsqsum_prp = 0.0; vsqsum_prl = 0.0;
  for (i=0; i<3; i++)
    vcm[i] = 0.0;
  
  for (q=0; q<Np; q++)
  {
    vp[0][q] = ran_gaussian(idum);
    vp[1][q] = ran_gaussian(idum);
    vp[2][q] = ran_gaussian(idum);
    vcm[0]  += vp[0][q];
    vcm[1]  += vp[1][q];
    vcm[2]  += vp[2][q];
  }
  
#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(vcm,gvcm,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vcm[0] = gvcm[0]; vcm[1] = gvcm[1]; vcm[2] = gvcm[2];
#endif /* MPI_PARALLEL */
  
  for (i=0; i<3; i++) {
    vcm[i] /= NpTot;
    for (q=0; q<Np; q++) {
      vp[i][q]   -= vcm[i];
      if (i == 0) {
        vsqsum_prl += SQR(vp[i][q]);
      } else {
        vsqsum_prp += SQR(vp[i][q]);
      }
    }
  }
  
#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(&vsqsum_prp,&gvsqsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vsqsum_prp = gvsqsum;
  mpierr = MPI_Allreduce(&vsqsum_prl,&gvsqsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vsqsum_prl = gvsqsum;
#endif /* MPI_PARALLEL */
  
  lambda_prl = sqrt(0.5*((double)NpTot)*beta_prl/vsqsum_prl);
  lambda_prp = sqrt(0.5*((double)NpTot)*beta_prp/vsqsum_prp) * sqrt(2.0);
  
  for (i=0; i<3; i++) {
    for (q=0; q<Np; q++) {
      if (i == 0) {
        vp[i][q] *= lambda_prl;
      } else{
        vp[i][q] *= lambda_prp;
      }
    }
  }

  return;
}

/*============================================================================*/
/*-----------------------------PRIVATE FUNCTIONS-------------------------------
 *
 * ran2()
 * rangaussian()
 */
/*============================================================================*/

/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. 
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

static double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;
  
  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

/*------------------------------------------------------------------------------
 * rangaussian: extracted from the Numerical Recipes in C (version 2) code.
 */

static double ran_gaussian(long int *idum)
{
  static int iset = 0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if (*idum < 0) iset = 0;
  if (iset == 0) {
    do {
      v1 = 2.0 * ran2(idum) - 1.0;
      v2 = 2.0 * ran2(idum) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1 * fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
