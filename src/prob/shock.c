#include "copyright.h"
/*============================================================================*/
/*! \file shock.c
 *  \brief Problem generator for shock
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

static Real Lx,Ly,Lz; /* root grid size, global to share with output functions */

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
  int j,js = pGrid->js, je = pGrid->je;
  int k,ks = pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs,ip,jp,kp;
  long int iseed;
  Real x1,x2,x3,x1l,x1u,x2l,x2u,x3l,x3u;
  static int frst=1;  /* flag so new history variables enrolled only once */
  Real L1,L2,L3,x1min,x2min,x3min,rx1,rx2,rx3,dv1,dv2,dv3,dvsq;
  long p,Npar,Npar3,NparGrid,NparTot;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;
  
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
#endif

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

  /* initialize magnetic field */
  for (k=ks; k<=ke+1; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      pGrid->B1i[k][j][i] = 1.0;
      pGrid->B2i[k][j][i] = 0.0;
      pGrid->B3i[k][j][i] = 0.0;
    }
  }
  }
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].B1c = 1.0;
      pGrid->U[k][j][i].B2c = 0.0;
      pGrid->U[k][j][i].B3c = 0.0;
    }
  }
  }
  
/* insert charged particles */

  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;

  x2min = pGrid->MinX[1];
  L2    = pGrid->MaxX[1] - x2min;
  
  x3min = pGrid->MinX[2];
  L3    = pGrid->MaxX[2] - x3min;

  /* get particle number */
  //  Npar  = (long)(par_geti("particle","parnumgrid"));

  Npar  = (int)(pow(par_geti("particle","parnumcell"),1.0/3.0));
  Npar3 = SQR(Npar)*Npar;

  NparGrid = (long)Npar3*pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2];
  
  pGrid->nparticle = NparGrid;
  grproperty[0].num = NparGrid;
  grproperty[0].m   = 1.0/((double)Npar3);

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
  
  if (par_geti("domain1","bc_ix1") == 6)
  {
    dv1= vinject;
    dv2= 0.0;
    dv3= 0.0;
  } 
  else if (par_geti("domain1","bc_ox1") == 6)
  {
    dv1=-vinject;
    dv2= 0.0;
    dv3= 0.0;
  }
  else if (par_geti("domain1","bc_ix2") == 6)
  {  
    dv1=  0.0;
    dv2=  vinject;
    dv3=  0.0;
  }
  else if (par_geti("domain1","bc_ox2") == 6)
  {  
    dv1=  0.0;
    dv2= -vinject;
    dv3=  0.0;
  }
  else if (par_geti("domain1","bc_ix3") == 6)
  { 
    dv1=  0.0;
    dv2=  0.0;
    dv3=  vinject;
  }
  else if (par_geti("domain1","bc_ox3") == 6)
  { 
    dv1=  0.0;
    dv2=  0.0;
    dv3= -vinject;
  }
  else
  {
    dv1 =  0.0;
    dv2 =  0.0;
    dv3 =  0.0;
  }

  dvsq = SQR(dv1)+SQR(dv2)+SQR(dv3);
                 
  /* loop over cells */
  for (k=ks; k<=ke; k++)
  {
    x3l = x3min + ((double)(k-ks))*pGrid->dx3;
    x3u = x3min + ((double)(k-ks+1))*pGrid->dx3;
  
  for (j=js; j<=je; j++)
  {
    x2l = x2min + ((double)(j-js))*pGrid->dx2;
    x2u = x2min + ((double)(j-js+1))*pGrid->dx2;
    
    for (i=is; i<=ie; i++)
    {
      x1l = x1min + ((double)(i-is))*pGrid->dx1;
      x1u = x1min + ((double)(i-is+1))*pGrid->dx1;
      
      /* loop over particles */
      for (ip=0;ip<Npar;ip++)
      {
        
        rx1 = x1l+pGrid->dx1/(double)Npar*((double)ip+0.5);

        for (jp=0;jp<Npar;jp++)
        {
          
          rx2 = x2l+pGrid->dx2/(double)Npar*((double)jp+0.5);

          for (kp=0;kp<Npar;kp++)
          {
            
            rx3 = x3l+pGrid->dx3/(double)Npar*((double)kp+0.5);
          
            pq = (&(pGrid->particle[p]));

            pq->x1 = rx1;
            pq->x2 = rx2;
            pq->x3 = rx3;
            pq->v1 = rv[0][p] + dv1;
            pq->v2 = rv[1][p] + dv2;
            pq->v3 = rv[2][p] + dv3;
            pq->property = 0;
            pq->pos = 1;
            pq->my_id = p;
#ifdef MPI_PARALLEL
            pq->init_id = myID_Comm_world;
#endif 
            p++;
          } /* kp loop */
        } /* jp loop */
      } /* ip loop */
    } /* i loop */
  } /* j loop */
  } /* k loop */

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


#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k, Real *eta_0)
{
  *eta_0 = 0.0;

  return;
}
#endif


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
 /*
 GridS *pG;
  Real divb;
  
  pG = (pM->Domain[0][0].Grid);
  divb = compute_div_b(pG);
  peg_pout(0,"divb = %e\n",divb);
 */ 
  return;
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

