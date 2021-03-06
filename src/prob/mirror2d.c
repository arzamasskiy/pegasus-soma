#include "copyright.h"
/*============================================================================*/
/*! \file mirror
 *  \brief Problem generator for mirror
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

static Real Lx,Ly; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);

static int property_type(const GrainS *gr, const GrainAux *grsub);
static int property_limit(const GrainS *gr, const GrainAux *grsub);

static int mytype;
static long nlis;

static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bxsqm1(const GridS *pG, const int i, const int j, const int k);
static Real InstPar(const GridS *pG, const int i, const int j, const int k);
static Real Delta(const GridS *pG, const int i, const int j, const int k);

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
  int ks = pGrid->ks;
  int ixs,jxs,ip,jp;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,x1l,x1u,x2l,x2u;
  Real amp;
  static int frst=1;
  Real L1,L2,x1min,x2min,x3min,rx1,rx2,dv1,dv2,dv3;
#ifdef DELTA_F
  Real df;
#endif
  long p,Npar,Npar2,NparGrid,NparTot;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;

#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  eta_hyper = par_getd_def("problem","eta_hyper",0.0);
#endif
  
/* Read problem parameters. */
  amp = par_getd("problem","amp");
  
/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  iseed = -1 - (ixs + pDomain->Nx[0]*jxs);

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  
  /* if delta-f, set (analytic) moments of the background dist. function */
#ifdef DELTA_F
  setbg(pGrid);
#endif
  
  /* initialize magnetic field */
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[ks][j][i].B1c = 1.0;
      pGrid->U[ks][j][i].B2c = 0.0;
      pGrid->U[ks][j][i].B3c = 0.0;
      pGrid->B1i[ks][j][i] = 1.0;
      pGrid->B2i[ks][j][i] = 0.0;
      pGrid->B3i[ks][j][i] = 0.0;
      if (i==ie) pGrid->B1i[ks][j][ie+1] = 1.0;
      if (j==je) pGrid->B2i[ks][je+1][i] = 0.0;
    }
  }

/* insert charged particles */

  /* identifiers for particle dumps */
  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]);
  mytype = par_geti_def("problem","mytype",0);
  
  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;

  x2min = pGrid->MinX[1];
  L2    = pGrid->MaxX[1] - x2min;

  x3min = pGrid->MinX[2];
  
  /* get particle number */
  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar2 = SQR(Npar);
  NparGrid = (long)Npar2*pGrid->Nx[0]*pGrid->Nx[1];

  pGrid->nparticle = NparGrid;
  grproperty[0].num = NparGrid;
  grproperty[0].m   = 1.0/((double)Npar2);

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

  for (j=js; j<=je; j++)
  {
    x2l = x2min + ((double)(j-js))*pGrid->dx2;
    x2u = x2min + ((double)(j-js+1))*pGrid->dx2;

    for (i=is; i<=ie; i++)
    {
      x1l = x1min + ((double)(i-is))*pGrid->dx1;
      x1u = x1min + ((double)(i-is+1))*pGrid->dx1;

      /* random bulk velocity fluctuations */
      dv1      =  amp*(ran2(&iseed)-1.0);
      dv2      =  amp*(ran2(&iseed)-1.0);
      dv3      =  amp*(ran2(&iseed)-1.0);
      
      for (ip=0;ip<Npar;ip++)
      {
        
        rx1 = x1l+pGrid->dx1/(double)Npar*((double)ip+0.5);

        for (jp=0;jp<Npar;jp++)
        {
          
          rx2 = x2l+pGrid->dx2/(double)Npar*((double)jp+0.5);
          
	    pq = (&(pGrid->particle[p]));

	    pq->x1 = rx1;
	    pq->x2 = rx2;
	    pq->x3 = x3min;
        
	    /* draw thermal velocities from distribution function */
	    pq->v1 = rv[0][p];
	    pq->v2 = rv[1][p];
	    pq->v3 = rv[2][p];
          
#ifdef DELTA_F
	    /* compute f(t=0) and store in structure */
	    getdf(pq,&df);
	    pq->f_0 = df;
#endif
	    /* add bulk motion */
            pq->v1 += dv1;
	    pq->v2 += dv2;
	    pq->v3 += dv3;  

	    pq->property = 0;
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
    dump_history_enroll(hst_Bxsqm1, "<Bx2m1>");
    dump_history_enroll(InstPar, "<InstPar>");
    dump_history_enroll(Delta, "<Delta>");
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
  GridS *pGrid;
  pGrid=(pM->Domain[0][0].Grid);
  
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  eta_hyper = par_getd_def("problem","eta_hyper",0.0);
#endif

  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]);
  mytype = par_geti_def("problem","mytype",0);

#ifdef DELTA_F
  setbg(pGrid);
#endif
  
/* enroll new history variables */
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_Bxsqm1, "<Bx2m1>");
  dump_history_enroll(InstPar, "<InstPar>");
  dump_history_enroll(Delta, "<Delta>");
  
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
  if (strcmp(name,"limit")==0) return property_limit;
  if (strcmp(name,"type")==0)  return property_type;
  return NULL;
}

static int property_type(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->property == mytype))
    return 1;
  else
    return 0;
}

static int property_limit(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<nlis))
    return 1;
  else
    return 0;
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

static Real hst_Bxsqm1(const GridS *pG, const int i, const int j, const int k)
{
  return SQR(pG->U[k][j][i].B1c)-1;
}

static Real Delta(const GridS *pG, const int i, const int j, const int k)
{
  GPCouple *pq;
  ConsS *pCons;
  Real p_tot,bsq,p_prl,d_tot,d_prl;

  pq = &(pG->Coup[k][j][i]);
  pCons = &(pG->U[k][j][i]);

  p_tot = ONE_3RD * ( pq->grid_p11 + pq->grid_p22 + pq->grid_p33 );
  d_tot = ONE_3RD * ( SQR(pq->grid_M1) + SQR(pq->grid_M2) + SQR(pq->grid_M3) ) / pq->grid_d;

  p_tot = p_tot - d_tot;

  bsq   = SQR(pCons->B1c) + SQR(pCons->B2c) + SQR(pCons->B3c);

  p_prl = pCons->B1c * ( pCons->B1c * pq->grid_p11 + pCons->B2c * pq->grid_p12
                         + pCons->B3c * pq->grid_p13 ) / bsq
    + pCons->B2c * ( pCons->B1c * pq->grid_p12 + pCons->B2c * pq->grid_p22
                     + pCons->B3c * pq->grid_p23 ) / bsq
    + pCons->B3c * ( pCons->B1c * pq->grid_p13 + pCons->B2c * pq->grid_p23
                     + pCons->B3c * pq->grid_p33 ) / bsq;

  d_prl = SQR( pCons->B1c * pq->grid_M1 + pCons->B2c * pq->grid_M2 + pCons->B3c * pq->grid_M3 )
    / pq->grid_d / bsq;

  p_prl = p_prl - d_prl;

  return 1.5 * ( p_tot/p_prl - 1.0);
}

static Real InstPar(const GridS *pG, const int i, const int j, const int k)
{
  GPCouple *pq;
  ConsS *pCons;
  Real p_tot,bsq,p_prl,p_prp;
  
  pq = &(pG->Coup[k][j][i]);
  pCons = &(pG->U[k][j][i]);
  
  p_tot = ONE_3RD * ( pq->grid_p11 + pq->grid_p22 + pq->grid_p33 );

  bsq   = SQR(pCons->B1c) + SQR(pCons->B2c) + SQR(pCons->B3c);

  p_prl = pCons->B1c * ( pCons->B1c * pq->grid_p11 + pCons->B2c * pq->grid_p12
                       + pCons->B3c * pq->grid_p13 ) / bsq
        + pCons->B2c * ( pCons->B1c * pq->grid_p12 + pCons->B2c * pq->grid_p22
                       + pCons->B3c * pq->grid_p23 ) / bsq
        + pCons->B3c * ( pCons->B1c * pq->grid_p13 + pCons->B2c * pq->grid_p23
                       + pCons->B3c * pq->grid_p33 ) / bsq;

  p_prp = 1.5*p_tot - 0.5*p_prl;

  return ( p_prp/p_prl - 1.0 - 0.5*bsq/p_prp );
}
