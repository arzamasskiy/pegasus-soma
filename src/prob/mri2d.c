#include "copyright.h"
/*============================================================================*/
/*! \file mri2d
 *  \brief Problem generator for axisymmetric MRI
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

static Real Lx,Lz; /* root grid size, global to share with output functions */
static Real betay;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);
static double ran_gaussian(long int *idum);

static int property_type(const GrainS *gr, const GrainAux *grsub);
static int property_limit(const GrainS *gr, const GrainAux *grsub);

static int mytype;
static long nlis;

static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k);
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
  int   ks = pGrid->ks;
  int ixs,jxs,ip,jp,ifield;
  long int iseed; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,x1l,x1u,x2l,x2u;
  Real kx, amp, angle;
  static int frst=1;
  Real L1,L2,x1min,x2min,rx1,rx2,rx3,dv1,dv2,dv3;
#ifdef DELTA_F
  Real df;
#endif
  long p,q,Npar,Npar2,NparGrid,NparTot,ntrack;
  int n;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;

  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  eta_hyper = par_getd_def("problem","eta_hyper",0.0);

  /* Read problem parameters. */
  amp = par_getd("problem","amp");
  angle = par_getd_def("problem","angle",0.0);
  angle *= PI/180.0;
  ifield = par_geti_def("problem","ifield",1);
  
  /* Enroll gravitational potential function */
  ShBoxCoord = xz;
  Omega_0 = par_getd_def("problem","omega",0.01);
  qshear = par_getd_def("problem","qshear",1.5);
  if (Omega_0 == 0) {
    Shear_0 = par_getd_def("problem","shear",0.015);
  } else {
    Shear_0 = qshear*Omega_0;
  }
  
  /* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  iseed = -1 - (ixs + pDomain->Nx[0]*jxs);

  /* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Lz = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  kx = 2*PI/Lx;

  /* if delta-f, set (analytic) moments of the background dist. function */
  betay = beta * ( 2.0*Omega_0 + cos((double)angle) - Shear_0 ) / ( 2.0*Omega_0 + cos((double)angle) );
#ifdef DELTA_F
  setbg(pGrid);
#endif
  
  /* initialize magnetic field */
  /* for 2d shearing box, B1 = Bx, B2 = Bz, B3 = By                                                                         
   * ifield = 1: uniform magnetic field with By/Bz = tan(angle) [default]                                                   
   * ifield = 2: Bz = B0 sin(x1) field with zero-net-flux                                                                   
   */
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      if (ifield == 1) {
        pGrid->U[ks][j][i].B1c = 0.0;
        pGrid->U[ks][j][i].B2c = cos((double)angle);
        pGrid->U[ks][j][i].B3c = sin((double)angle);
        pGrid->B1i[ks][j][i] = 0.0;
        pGrid->B2i[ks][j][i] = cos((double)angle);
        pGrid->B3i[ks][j][i] = sin((double)angle);
        if (i==ie) pGrid->B1i[ks][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[ks][je+1][i] = cos((double)angle);
      }
      if (ifield == 2) {
        cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
        pGrid->U[ks][j][i].B1c = 0.0;
        pGrid->U[ks][j][i].B2c = sin((double)kx*x1);
        pGrid->U[ks][j][i].B3c = 0.0;
        pGrid->B1i[ks][j][i] = 0.0;
        pGrid->B2i[ks][j][i] = sin((double)kx*x1);
        pGrid->B3i[ks][j][i] = 0.0;
        if (i==ie) pGrid->B1i[ks][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[ks][je+1][i] = sin((double)kx*x1);
      }
    }
  }


/* insert charged particles */

  /* identifiers for particle dumps */
  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]);
  mytype = par_geti_def("problem","mytype",0);
  ntrack = 1;
  
  /* get grid size */
  x1min = pGrid->MinX[0];
  L1    = pGrid->MaxX[0] - x1min;

  x2min = pGrid->MinX[1];
  L2    = pGrid->MaxX[1] - x2min;

  /* get particle number */
  //  Npar  = (long)(par_geti("particle","parnumgrid"));

  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar2 = SQR(Npar);
  NparGrid = (long)Npar2*pGrid->Nx[0]*pGrid->Nx[1];

  pGrid->nparticle = NparGrid;
  grproperty[0].num = NparGrid;
  grproperty[0].m   = 1.0/((double)Npar2);
  if (mytype !=0) {
    grproperty[0].num -= ntrack;
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

  rx3 = pGrid->MinX[2];

  for (j=js; j<=je; j++)
	{
	  x2l = x2min + ((double)(j-js))*pGrid->dx2;
	  x2u = x2min + ((double)(j-js+1))*pGrid->dx2;

	  for (i=is; i<=ie; i++)
	    {
	      x1l = x1min + ((double)(i-is))*pGrid->dx1;
	      x1u = x1min + ((double)(i-is+1))*pGrid->dx1;

	      dv1  = amp*(ran2(&iseed)-1.0);
	      dv2  = amp*(ran2(&iseed)-1.0);
	      dv3  = amp*(ran2(&iseed)-1.0);
      
      for (ip=0;ip<Npar;ip++)
      {
        
        rx1 = x1l+pGrid->dx1/(double)Npar*((double)ip+0.5);

        for (jp=0;jp<Npar;jp++)
		    {
          
		      rx2 = x2l+pGrid->dx2/(double)Npar*((double)jp+0.5);

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
		    } /* jp loop */
      } /* ip loop */
    } /* i loop */
	} /* j loop */
  
  /* declare tracked particles -- here it's one per Grid */
  if (mytype != 0) pGrid->particle[0].property = 1;
  
  /* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BxBy, "<-BxBy>");
    dump_history_enroll(Delta, "<Delta>");
    frst = 0;
  }
  
  free_2d_array(rv);

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * t_usr_par_prop()      - returns a user defined particle selection function
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
  Real angle;

  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  eta_hyper = par_getd_def("problem","eta_hyper",0.0);
  
  ShBoxCoord = xz;
  Omega_0 = par_getd_def("problem","omega",0.01);
  qshear = par_getd_def("problem","qshear",1.5);
  if (Omega_0 == 0) {
    Shear_0 = par_getd_def("problem","shear",0.01);
  } else {
    Shear_0 = qshear*Omega_0;
  }

  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]);
  mytype = par_geti_def("problem","mytype",0);
  
  /* if delta-f, set (analytic) moments of the background dist. function */
#ifdef DELTA_F
  angle = par_getd_def("problem","angle",0.0);
  angle *= PI/180.0;
  betay = beta * ( 2.0*Omega_0 + cos((double)angle) - Shear_0 ) / ( 2.0*Omega_0 + cos((double)angle) );
  setbg(pGrid);
#endif

  /* enroll new history variables */
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-BxBy>");
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
  return pG->U[k][j][i].B3c;
}

static Real hst_Bz(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k)
{
  return -pG->U[k][j][i].B1c * pG->U[k][j][i].B3c;
}

static Real Delta(const GridS *pG, const int i, const int j, const int k)
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
  
  p_prp = 1.5 * p_tot - 0.5 * p_prl;

  return (p_prp/p_tot - p_prl/p_tot);
}



#ifdef DELTA_F
void getdf_user(GrainS *gr, Real *df)
{
  Real vsq,vysq;

  vsq  = SQR(gr->v1)+SQR(gr->v2);
  vysq = SQR(gr->v3);
                    
  *df = exp(-vsq/beta)*exp(-vysq/betay);
                    
  return;
}

void setbg_user(GridS *pG)
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
        pw->grid_p33 = 0.5*betay;
      }
  
  return;
}
#endif

void initdf_user(long int *idum, Real **vp,
                                   const long int Np, const long int NpTot)
{
  int i,mpierr;
  long int q;
  Real vsqsum_xz, vsqsum_y, vcm[3], lambda_xz, lambda_y;
#ifdef MPI_PARALLEL
  Real gvsqsum, gvcm[3];
#endif
  
  vsqsum_xz = 0.0; vsqsum_y = 0.0;
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
      if (i == 2) {
        vsqsum_y += SQR(vp[i][q]);
      } else {
        vsqsum_xz += SQR(vp[i][q]);
      }
    }
  }
  
#ifdef MPI_PARALLEL
  mpierr = MPI_Allreduce(&vsqsum_xz,&gvsqsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vsqsum_xz = gvsqsum;
  mpierr = MPI_Allreduce(&vsqsum_y,&gvsqsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpierr) peg_error("[dist_func]: MPI_Allreduce error = %d\n", mpierr);
  vsqsum_y = gvsqsum;
#endif /* MPI_PARALLEL */
  
  lambda_xz = sqrt(0.5*((double)NpTot)*beta /vsqsum_xz) * sqrt(2.0);
  lambda_y  = sqrt(0.5*((double)NpTot)*betay/vsqsum_y);
  
  for (i=0; i<3; i++) {
    for (q=0; q<Np; q++) {
      if (i == 2) {
        vp[i][q] *= lambda_y;
      } else{
        vp[i][q] *= lambda_xz;
      }
    }
  }
  
  return;
}

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
