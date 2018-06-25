#include "copyright.h"
/*============================================================================*/
/*! \file gyrate.c
 *  \brief Problem generator for particle gyration test
 */
/*============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "pegasus.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

Real Lx,Ly; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);

static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int ixs,jxs,i,j,ip,jp;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1l,x1u,x2l,x2u;
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */
  Real L1,L2,x1min,x2min,rx1,rx2,rx3,rv1,rv2,rv3;
  long p,Npar,NparGrid,NparTot;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;

#ifdef SHEARING_BOX
  ShBoxCoord = xy;
  Omega_0 = par_getd_def("problem","omega",0.0);
  qshear  = par_getd_def("problem","qshear",1.5);
  if (Omega_0 != 0.0) {
    Shear_0 = qshear*Omega_0;
  } else {
    Shear_0 = par_getd_def("problem","shear",1.0);
  }
#endif
  
/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  iseed = -1 - (ixs + pDomain->Nx[0]*jxs);

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie+1; i++) {
    pGrid->B1i[ks][j][i] = 0.0;
    pGrid->B2i[ks][j][i] = 0.0;
    pGrid->B3i[ks][j][i] = 1.0;
  }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[ks][j][i].B1c = 0.5*(pGrid->B1i[ks][j][i]+pGrid->B1i[ks][j][i+1]);
      pGrid->U[ks][j][i].B2c = 0.5*(pGrid->B2i[ks][j][i]+pGrid->B2i[ks][j+1][i]);
      pGrid->U[ks][j][i].B3c = pGrid->B3i[ks][j][i];
    }
  }

/* insert charged particles */

  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;
  
  x2min = pGrid->MinX[1];
  L2    = pGrid->MaxX[1] - x2min;

  /* get particle number */
  Npar  = (long)(par_geti("particle","parnumgrid"));
  NparGrid = Npar;

  pGrid->nparticle = NparGrid;
  grproperty[0].num = NparGrid;
  grproperty[0].m   = 1.0/(double)Npar;

  NparTot = NparGrid;
#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(&NparGrid,&NparTot,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[problem]: MPI_Allreduce error = %d\n",mpierr);
#endif

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  /* initialize particle velocity arrays */
  if ((rv = (Real**)calloc_2d_array(3,NparGrid,sizeof(Real))) == NULL) {
    peg_error("[problem]: Error allocating memory for particle velocity\n");
  }

  /* set initial conditions for the particles */
  p = 0;

  initdf(&iseed,rv,NparGrid,NparTot);

  rx3 = 0.5*(pGrid->MinX[2]+pGrid->MaxX[2]);

  for (p=0;p<Npar;p++)
  {
    rx1 = x1min + L1*ran2(&iseed);
    rx2 = x2min + L2*ran2(&iseed);
        
    pGrid->particle[p].property = 0;

    rv1 = 1.0;
    rv2 = 0.0;
    rv3 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
    rv2 -= -Shear_0*rx1;
#endif
#endif 
          
    pGrid->particle[p].x1  = rx1;
    pGrid->particle[p].x2  = rx2;
    pGrid->particle[p].x3  = rx3;
    pGrid->particle[p].v1  = rv1;
    pGrid->particle[p].v2  = rv2;
    pGrid->particle[p].v3  = rv3;
          
    pGrid->particle[p].pos   = 1; /* grid particle */
    pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
    pGrid->particle[p].init_id = myID_Comm_world;
#endif
  }
  
/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    frst = 0;
  }

  free_2d_array(rv);
  
  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{

/* enroll new history variables */

  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");

  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

PropFun_t get_usr_par_prop(const char *name)
{
  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  return;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}

#ifdef DELTA_F
void getdf_user(GrainS *gr, Real *df)
{
  return;
}
void setbg_user(GridS *pG)
{
  return;
}
#endif
void initdf_user(long int *idum, Real **vp,
                 const long int Np, const long int NpTot)
{
  return;
}


/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
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

double ran2(long int *idum)
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
 * MHD history variables
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

/*! \fn static Real hst_Bx(const GridS *pG, const int i,const int j,const int k)
 *  \brief x-component of B-field */
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B1c;
}

/*! \fn static Real hst_By(const GridS *pG, const int i,const int j,const int k)
 *  \brief y-component of B-field */
static Real hst_By(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

static Real hst_Bz(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B3c;
}
