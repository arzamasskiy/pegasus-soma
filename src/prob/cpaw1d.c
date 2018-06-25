#include "copyright.h"
/*============================================================================*/
/*! \file cpaw1d.c
 *  \brief Problem generator for cpaw in 1d
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

static SolnS *Soln=NULL;

static Real Lx; /* root grid size, global to share with output functions */

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
  GrainS *pq;
  int i,is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks;
  int ixs,ip;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,x1l,x1u;
  Real amp;
  static int frst=1;
  Real L1,x1min,rx1,rx2,rx3,dv1,dv2,dv3,dvsq,cs,sn;
#ifdef DELTA_F
  Real df;
#endif
  Real v_par,omega_r,omega_l,v_perp,x_par;
  long p,q,Npar,NparGrid,NparTot;
  int n,dir,nx1,nwx;
  Real b_par,b_perp,fac,lambda,k_par;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;

  nx1 = (ie-is+1) + 2*nghost;
  if ((Soln = (SolnS*)calloc_1d_array(nx1,sizeof(SolnS)))==NULL)
    peg_error("[problem]: Error allocating memory for Soln\n");

  v_par = par_getd_def("problem","v_par",0.0);
  amp = par_getd("problem","amp");
  v_perp = amp;
  b_perp = v_perp;
  v_par  = v_par;
  b_par  = 1.0;

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  iseed = -1 - (ixs);

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];

  /* initialize wavenumbers, given input number of waves per L */
  lambda = Lx;
  k_par  = 2.0*PI/lambda;

  nwx = par_geti_def("problem","nwx",1);
  k_par *= ((double)nwx);

  dir = par_geti_def("problem","dir",1); /* right(1)/left(2) polarization */
  if (dir == 1)
    fac = 1.0; /* right polarization */
  else
    fac = -1.0; /* left polarization */

  /* for whistler wave */
  omega_r = 0.5*SQR(k_par)*(sqrt(1.0+SQR(2.0/k_par)) + 1.0);
  omega_l = 0.5*SQR(k_par)*(sqrt(1.0+SQR(2.0/k_par)) - 1.0);
  
  if (dir == 1)
    v_perp = v_perp * k_par / omega_r;
  else
    v_perp = v_perp * k_par / omega_l;

  /* if delta-f, set (analytic) moments of the background dist. function */
#ifdef DELTA_F
  setbg(pGrid);
#endif
  
  /* initialize magnetic field */
  for (i=is; i<=ie; i++) {
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
    
    x_par = x1;
    sn = sin(k_par*x_par);
    cs = cos(k_par*x_par) * fac;
    
    Soln[i].d  = 1.0;
    Soln[i].M1 = v_par;
    Soln[i].M2 = -v_perp*sn;
    Soln[i].M3 = -v_perp*cs;
    
    pGrid->B1i[ks][js][i] = b_par;
    pGrid->B2i[ks][js][i] = b_perp*sn;
    pGrid->B3i[ks][js][i] = b_perp*cs;
    if (i==ie) pGrid->B1i[ks][js][ie+1] = b_par;
  }
  for (i=is; i<=ie; i++) {
    Soln[i].B1c = 0.5*(pGrid->B1i[ks][js][i]+pGrid->B1i[ks][js][i+1]);
    Soln[i].B2c = pGrid->B2i[ks][js][i];
    Soln[i].B3c = pGrid->B3i[ks][js][i];
    pGrid->U[ks][js][i].B1c = Soln[i].B1c;
    pGrid->U[ks][js][i].B2c = Soln[i].B2c;
    pGrid->U[ks][js][i].B3c = Soln[i].B3c;
  }

/* insert charged particles */

  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;

  /* get particle number */
  //  Npar  = (long)(par_geti("particle","parnumgrid"));

  Npar  = (int)(par_geti("particle","parnumcell"));
  NparGrid = (long)Npar*pGrid->Nx[0];
  
  pGrid->nparticle = NparGrid;
  grproperty[0].num = NparGrid;
  grproperty[0].m   = 1.0/((Real)Npar);
  
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
  rx2 = 0.5*(pGrid->MinX[1]+pGrid->MaxX[1]);
    
  /* loop over cells */
    for (i=is; i<=ie; i++)
    {
      x1l = x1min + ((Real)(i-is))*pGrid->dx1;
      x1u = x1min + ((Real)(i-is+1))*pGrid->dx1;
      
      /* loop over particles */
      for (ip=0;ip<Npar;ip++)
      {
        
        rx1 = x1l+pGrid->dx1/(Real)Npar*((Real)ip+0.5);
        x_par = rx1;
        sn = sin(k_par*x_par);
        cs = cos(k_par*x_par) * fac;

        /* compute perturbation */
        dv1    =  v_par;
        dv2    = -v_perp*sn;
        dv3    = -v_perp*cs;
        
	pq = (&(pGrid->particle[p]));

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

	pq->property = 0;
	pq->pos = 1;
	pq->my_id = p;
#ifdef MPI_PARALLEL
	pq->init_id = myID_Comm_world;
#endif 
	p++;
      } /* ip loop */
    } /* i loop */

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
  int i,is,ie,js,ks,Nx1;
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
  js = pGrid->js; 
  ks = pGrid->ks;
  Nx1 = (ie-is+1);
  
  deposit(pD,0);
  
  /* compute L1 error in each variable, and rms total error */

  for (i=is; i<=ie; i++) {
    pq = &(pGrid->Coup[ks][js][i]);
    error.d   += fabs(pq->grid_d  - Soln[i].d  );
    error.M1  += fabs(pq->grid_M1 - Soln[i].M1 );
    error.M2  += fabs(pq->grid_M2 - Soln[i].M2 );
    error.M3  += fabs(pq->grid_M3 - Soln[i].M3 );
    error.B1c += fabs(pGrid->U[ks][js][i].B1c - Soln[i].B1c);
    error.B2c += fabs(pGrid->U[ks][js][i].B2c - Soln[i].B2c);
    error.B3c += fabs(pGrid->U[ks][js][i].B3c - Soln[i].B3c);
  }
      
  /* Compute RMS error over all variables, and print out */
  
  rms_error = SQR(error.d) + SQR(error.M1) + SQR(error.M2) + SQR(error.M3);
  rms_error += SQR(error.B1c) + SQR(error.B2c) + SQR(error.B3c);
  rms_error = sqrt(rms_error)/((double)Nx1);
  
  /* Print error to file "hyb_cpaw1d-errors.dat" */
  
  fname = "cpaw1d-errors.dat";
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
      peg_perr(-1,"[Userwork_after_loop]: Unable to open file.\n");
      free(fname);
      return;
    }
    /* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
    fprintf(fp,"  B1c  B2c  B3c");
    fprintf(fp,"\n#\n");
  }
  
  fprintf(fp,"%d  %e",Nx1,rms_error);
  fprintf(fp,"  %e  %e  %e  %e",
          (error.d/(double)Nx1),
          (error.M1/(double)Nx1),
          (error.M2/(double)Nx1),
          (error.M3/(double)Nx1) );
  fprintf(fp,"  %e  %e  %e",
          (error.B1c/(double)Nx1),
          (error.B2c/(double)Nx1),
          (error.B3c/(double)Nx1));
  fprintf(fp,"\n");
  
  fclose(fp);
  
  free_1d_array(Soln);

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

