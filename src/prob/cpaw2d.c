#include "copyright.h"
/*============================================================================*/
/*! \file hybrid_2d_test.c
 *  \brief Problem generator for hybrid-PIC tests
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

static SolnS **RootSoln=NULL;
static Real A3(const Real x1, const Real x2);
static Real fac, sin_a, cos_a, b_par, b_perp, k_par, Lx, Ly;

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
  SolnS **Soln;
  GrainS *pq;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int ixs,jxs,i,j,ip,jp;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,x1l,x1u,x2l,x2u;
  Real amp;
  static int frst=1;  /* flag so new history variables enrolled only once */
  Real L1,L2,x1min,x2min,rx1,rx2,rx3,dv1,dv2,dv3,dvsq,cs,sn,angle;
#ifdef DELTA_F
  Real df;
#endif
  Real lambda,v_par,omega_r,omega_l,v_perp;
  long p,q,Npar,Npar2,NparGrid,NparTot;
  int  n,dir,nx1,nx2;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;

  nx1 = (ie - is + 1) + 2*nghost;
  nx2 = (je - js + 1) + 2*nghost;

  if ((Soln = (SolnS**)calloc_2d_array(nx2,nx1,sizeof(SolnS))) == NULL)
    peg_error("[problem]: Error allocating memory for Soln\n");

  if (pDomain->Level == 0){
    if ((RootSoln =(SolnS**)calloc_2d_array(nx2,nx1,sizeof(SolnS)))==NULL)
      peg_error("[problem]: Error allocating memory for RootSoln\n");
  }

/* Read problem parameters. */
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  eta_hyper = par_getd_def("problem","eta_hyper",0.0);
#endif
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

  amp = par_getd("problem","amp");
  v_perp = amp;
  b_perp = v_perp;
  v_par  = 0.0;
  b_par  = 1.0;

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  iseed = -1 - (ixs + pDomain->Nx[0]*jxs);

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

/* An angle =  0.0 is a wave aligned with the x1-direction. */
/* An angle = 90.0 is a wave aligned with the x2-direction. */

  angle = par_getd("problem","angle");

  /* Compute the sin and cos of the angle and the wavelength. */
  /* Put one wavelength in the grid */

  if (angle == 0.0) {
    sin_a = 0.0;
    cos_a = 1.0;
    lambda = Lx;
  }
  else if (angle == 90.0) {
    sin_a = 1.0;
    cos_a = 0.0;
    lambda = Ly;
  }
  else {
   
/* We put 1 wavelength in each direction.  Hence the wavelength
 *      lambda = (pDomain->RootMaxX[0] - pDomain->RootMinX[0])*cos_a
 *  AND lambda = (pDomain->RootMaxX[1] - pDomain->RootMinX[1])*sin_a;
 *  are both satisfied. */

    if(Lx == Ly){
      cos_a = sin_a = sqrt(0.5);
    }
    else{
      angle = atan((double)(Lx/Ly));
      sin_a = sin(angle);
      cos_a = cos(angle);
    }
    /* Use the larger angle to determine the wavelength */
    if (cos_a >= sin_a) {
      lambda = Lx*cos_a;
    } else {
      lambda = Ly*sin_a;
    }
  }

/* initialize wavenumbers, given input number of waves per L */

  k_par = 2.0*PI/lambda;

  dir = par_geti_def("problem","dir",1); /* right(1)/left(2) polarization */
  if (dir == 1)
    fac = 1.0;
  else
    fac = -1.0;
  
  /* for whistler wave */
  omega_r = 0.5*SQR(k_par)*(sqrt(1.0+SQR(2.0/k_par)) + 1.0);
  omega_l = 0.5*SQR(k_par)*(sqrt(1.0+SQR(2.0/k_par)) - 1.0);
  
  if (dir == 1)
    v_perp = v_perp * k_par / omega_r;
  else
    v_perp = v_perp * k_par / omega_l;


  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
      
      sn = sin((double)k_par*(x1*cos_a + x2*sin_a)) * fac;
      cs = cos((double)k_par*(x1*cos_a + x2*sin_a));
      
      Soln[j][i].d  = 1.0;
      Soln[j][i].M1 = v_par*cos_a + v_perp*sn*sin_a;
      Soln[j][i].M2 = v_par*sin_a - v_perp*sn*cos_a;
      Soln[j][i].M3 =             - v_perp*cs;

/* Initialize magnetic field */
      
      cs = cos(k_par*(x1*cos_a + x2*sin_a));

      x1 -= 0.5*pGrid->dx1;
      x2 -= 0.5*pGrid->dx2;

      pGrid->B1i[ks][j][i] = -(A3(x1,(x2+pGrid->dx2)) - A3(x1,x2))/pGrid->dx2;
      pGrid->B2i[ks][j][i] =  (A3((x1+pGrid->dx1),x2) - A3(x1,x2))/pGrid->dx1;
      pGrid->B3i[ks][j][i] =  b_perp*cs;

    }
  }
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        Soln[j][i].B1c = 0.5*(pGrid->B1i[ks][j][i]+pGrid->B1i[ks][j][i+1]);
        Soln[j][i].B2c = 0.5*(pGrid->B2i[ks][j][i]+pGrid->B2i[ks][j+1][i]);
        Soln[j][i].B3c = pGrid->B3i[ks][j][i];
        pGrid->U[ks][j][i].B1c = Soln[j][i].B1c;
        pGrid->U[ks][j][i].B2c = Soln[j][i].B2c;
        pGrid->U[ks][j][i].B3c = Soln[j][i].B3c;
      }
    }

/* insert charged particles */

  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;

  x2min = pGrid->MinX[1];
  L2    = pGrid->MaxX[1] - x2min;

  /* get particle number */
  //  NparGrid  = (long)(par_geti("particle","parnumgrid"));

  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar2 = SQR(Npar);

  NparGrid = ((long)Npar2)*pGrid->Nx[0]*pGrid->Nx[1];

  pGrid->nparticle = NparGrid;
  grproperty[0].num = NparGrid;
  grproperty[0].m   = 1.0/(double)Npar2;

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

  /* initialize distribution function */
  initdf(&iseed,rv,NparGrid,NparTot);
  
  rx3 = 0.5*(pGrid->MinX[2]+pGrid->MaxX[2]);

  for (j=js; j<=je; j++)
  {
    x2l = x2min + ((double)(j-js))*pGrid->dx2;
    x2u = x2min + ((double)(j-js+1))*pGrid->dx2;
    
    for (i=is; i<=ie; i++)
    {
      x1l = x1min + ((double)(i-is))*pGrid->dx1;
      x1u = x1min + ((double)(i-is+1))*pGrid->dx1;
      
      for (ip=0;ip<Npar;ip++)
      {
        
        rx1 = x1l+pGrid->dx1/Npar*(ip+0.5);

        for (jp=0;jp<Npar;jp++)
        {
          rx2 = x2l+pGrid->dx2/((double)Npar)*((double)jp+0.5);

          sn = sin((double)k_par*(rx1*cos_a + rx2*sin_a)) * fac;
          cs = cos((double)k_par*(rx1*cos_a + rx2*sin_a));
          
          /* compute perturbation */
          dv1  = v_par*cos_a + v_perp*sn*sin_a;
          dv2  = v_par*sin_a - v_perp*sn*cos_a;
          dv3  =             - v_perp*cs;
#ifdef SHEARING_BOX
#ifndef FARGO
          dv2 -= -Shear_0*rx1;
#endif
#endif      
          dvsq = SQR(dv1)+SQR(dv2)+SQR(dv3);
          
          pq = (&(pGrid->particle[p]));

          pq->property = 0;
          pq->x1 = rx1;
          pq->x2 = rx2;
          pq->x3 = rx3;
          
          /* draw velocities from full distribution function */
          pq->v1 = rv[0][p] + dv1;
          pq->v2 = rv[1][p] + dv2;
          pq->v3 = rv[2][p] + dv3;
          
#ifdef DELTA_F
          /* compute f(t=0) and store in structure */
          pq->v1 -= dv1;
          pq->v2 -= dv2;
          pq->v3 -= dv3;
          getdf(pq,&df);
          pq->f_0 = df;
          pq->v1 += dv1;
          pq->v2 += dv2;
          pq->v3 += dv3;
#endif
          pq->pos = 1;
          pq->my_id = p;
#ifdef MPI_PARALLEL
          pq->init_id = myID_Comm_world;
#endif
          p++;
          
        } /* jp loop */
      } /* ip loop */
    } /* i loop */
  } /* j loop */

/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    frst = 0;
  }

  /* save solution on root grid */

  if (pDomain->Level == 0) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        RootSoln[j][i].d  = Soln[j][i].d ;
        RootSoln[j][i].M1 = Soln[j][i].M1;
        RootSoln[j][i].M2 = Soln[j][i].M2;
        RootSoln[j][i].M3 = Soln[j][i].M3;
        RootSoln[j][i].B1c = Soln[j][i].B1c;
        RootSoln[j][i].B2c = Soln[j][i].B2c;
        RootSoln[j][i].B3c = Soln[j][i].B3c;
      }
    }
  }

  /* finalize */
  free_2d_array(rv);
  free_2d_array(Soln);

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
  dump_history_enroll(hst_Bz, "<Bz>");

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

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
  GridS *pGrid;
  DomainS *pD;
  GPCouple *pq;
  int i,j,is,ie,js,je,ks,Nx1,Nx2;
  Real rms_error=0.0;
  SolnS error;
  FILE *fp;
  char *fname;
  error.d = 0.0;
  error.M1 = 0.0;
  error.M2 = 0.0;
  error.M3 = 0.0;
  error.B1c = 0.0;
  error.B2c = 0.0;
  error.B3c = 0.0;

  /* Compute error only on root Grid, which is in Domain[0][0] */
  pD = (DomainS*)&(pM->Domain[0][0]);
  pGrid = pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  Nx1 = (ie-is+1);
  Nx2 = (je-js+1);
  
  deposit(pD,0);

  /* compute L1 error in each variable, and rms total error */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pq = &(pGrid->Coup[ks][j][i]);
      error.d   += fabs(pq->grid_d   - RootSoln[j][i].d  );
      error.M1  += fabs(pq->grid_M1  - RootSoln[j][i].M1 );
      error.M2  += fabs(pq->grid_M2  - RootSoln[j][i].M2 );
      error.M3  += fabs(pq->grid_M3  - RootSoln[j][i].M3 );
      error.B1c += fabs(pGrid->U[ks][j][i].B1c - RootSoln[j][i].B1c);
      error.B2c += fabs(pGrid->U[ks][j][i].B2c - RootSoln[j][i].B2c);
      error.B3c += fabs(pGrid->U[ks][j][i].B3c - RootSoln[j][i].B3c);
    }
  }

  /* Compute RMS error over all variables */

  rms_error = SQR(error.d) + SQR(error.M1) + SQR(error.M2) + SQR(error.M3);
  rms_error += SQR(error.B1c) + SQR(error.B2c) + SQR(error.B3c);
  rms_error = sqrt(rms_error)/(double)(Nx1*Nx2);

  /* Print error to file "hyb2d-errors.dat" */

  fname = "cpaw2d-errors.dat";
  /* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      peg_error("[Userwork_after_loop]: Unable to reopen file.\n");
      return;
    }
  }
  /* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      peg_error("[Userwork_after_loop]: Unable to open file.\n");
      return;
    }
    /* write out some header information */
    fprintf(fp,"# Nx1   Nx2  RMS-Error  d  M1  M2  M3");
    fprintf(fp,"  B1c  B2c  B3c");
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %e  %e  %e  %e  %e",Nx1,Nx2,rms_error,
          (error.d/(double)(Nx1*Nx2)),
          (error.M1/(double)(Nx1*Nx2)),
          (error.M2/(double)(Nx1*Nx2)),
          (error.M3/(double)(Nx1*Nx2)));
  fprintf(fp,"  %e  %e  %e",
          (error.B1c/(double)(Nx1*Nx2)),
          (error.B2c/(double)(Nx1*Nx2)),
          (error.B3c/(double)(Nx1*Nx2)));
  fprintf(fp,"\n");

  fclose(fp);
  
  free_2d_array(RootSoln);

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

static Real A3(const Real x1, const Real x2)
{
  return b_par*(x1*sin_a - x2*cos_a) 
      - (fac*b_perp/k_par)*cos(k_par*(x1*cos_a + x2*sin_a));
}
