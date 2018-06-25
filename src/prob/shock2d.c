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

static Real Lx,Ly; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * hst_*            - new history variables
 *============================================================================*/

static int mytype;
static long nlis;

static double ran2(long int *idum);

static void conductwall_ox1(GridS *pGrid);

static int property_type(const GrainS *gr, const GrainAux *grsub);
static int property_limit(const GrainS *gr, const GrainAux *grsub);

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
  int ks = pGrid->ks;
  int ixs,jxs,ip,jp;
  long int iseed;
  Real x1,x2,x3,x1l,x1u,x2l,x2u,x3l;
  static int frst=1;  /* flag so new history variables enrolled only once */
  Real L1,L2,x1min,x2min,x3min,rx1,rx2,rx3,dv1,dv2,dv3;
  long p,Npar,Npar2,NparGrid,NparTot,ntrack;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;
  
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",1.0);
#endif

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  iseed = -1 - (ixs + pDomain->Nx[0]*jxs);

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  /* initialize magnetic field */
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      pGrid->B1i[ks][j][i] = 1.0;
      pGrid->B2i[ks][j][i] = 0.0;
      pGrid->B3i[ks][j][i] = 0.0;
    }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[ks][j][i].B1c = 1.0;
      pGrid->U[ks][j][i].B2c = 0.0;
      pGrid->U[ks][j][i].B3c = 0.0;
    }
  }
  
/* insert charged particles */

  /* identifiers for particle dumps */
  nlis = par_geti_def("problem","nlis",0);
  mytype = par_geti_def("problem","mytype",0);
  ntrack = 1;

  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;

  x2min = pGrid->MinX[1];
  L2    = pGrid->MaxX[1] - x2min;
  
  x3min = pGrid->MinX[2];

  /* get particle number */
  //  Npar  = (long)(par_geti("particle","parnumgrid"));

  Npar  = (int)(pow(par_geti("particle","parnumcell"),1.0/2.0));
  Npar2 = SQR(Npar);
  NparGrid = (long)Npar2*pGrid->Nx[0]*pGrid->Nx[1];
  
  pGrid->nparticle = NparGrid;
 
  grproperty[0].num = NparGrid-ntrack;
  grproperty[0].m   = 1.0/((double)Npar2);
  if (mytype != 0) {
    grproperty[1].num = ntrack;
    grproperty[1].m   = 1.0/((double)Npar2);
  }

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
  else
  {
    dv1 =  0.0;
    dv2 =  0.0;
    dv3 =  0.0;
  }

  /* loop over cells */
  x3l = x3min;
  
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

            rx3 = x3l;
          
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
        } /* jp loop */
      } /* ip loop */
    } /* i loop */
  } /* j loop */


  /* declare tracked particles -- here it's one per Grid */
  pGrid->particle[0].property = 1;

  free_2d_array(rv);

/* enroll special boundary conditions */
  
//  bvals_mhd_fun(pDomain, right_x1, conductwall_ox1);

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
  mytype = par_geti_def("problem","mytype",0);
  nlis = par_geti_def("problem","nlis",0);

#ifdef RESISTIVITY  
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
#endif

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
  if (strcmp(name,"type")==0) return property_type;
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


static void conductwall_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  int ju,ku; /* j-upper, k-upper */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i].B1c= 1.0;
        pGrid->U[k][j][ie+i].B2c= 0.0;
        pGrid->U[k][j][ie+i].B3c= 0.0;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i]= 1.0;
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i]= 0.0;
      }
    }
  }

  if (pGrid->Nx[2] >1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i]= 0.0;
      }
    }
  }

  return;
}

