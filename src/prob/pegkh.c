#include "copyright.h"
/*============================================================================*/
/*! \file pegkh.c
 *  \brief Problem generator for KH instability using Pegasus. 
 *
 * PURPOSE: Problem generator for KH instability using Pegasus in analogy to
 * the Athena KH problem generator.
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

static Real Lx, Ly; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * hst_*            - new history variables
 *============================================================================*/

static int mytype;
static long nlis;

static int property_type(const GrainS *gr, const GrainAux *grsub);
static int property_limit(const GrainS *gr, const GrainAux *grsub);

static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real Delta(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  GrainS *pq;
  
  int i,j,n,ip,jp;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  Real x1u, x1l, x2u, x2l;
  Real amp,vflow,b0x,b0z,lshear,lmode,sigma,x1,x2,x3;
  long int iseed = -1;
  static int frst=1;  /* flag so new history variables enrolled only once */
  Real x1min,x2min,rx1,rx2,rx3,dv1,dv2,dv3;
  long p,q,Npar,Npar2,NparGrid,NparTot,ntrack;
  int my_iproc,my_jproc,my_kproc;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;
  
  /* Read problem parameters */
  vflow = par_getd("problem","vflow");
  amp   = par_getd("problem","amp");             // in units of d_i\omega_i
  b0x   = par_getd_def("problem", "b0x", 1.0);
  b0z   = par_getd_def("problem", "b0z", 0.0);
  lshear= par_getd_def("problem", "lshear", 3.0);     // in units of d_i
  lmode = par_getd_def("problem", "lmode", 100.0);
  sigma = par_getd_def("problem", "sigma", 0.2);
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem", "eta_O", 0.001);
  eta_hyper = par_getd_def("problem","eta_hyper",0.0);
#endif
  
  /* Ensure a different initial random seed for each process in an MPI calc. */
  int ixs = pGrid->Disp[0];
  int jxs = pGrid->Disp[1];
  iseed = -1 - (ixs + pDomain->Nx[0]*jxs);
  
  /* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  
  /* Initialize velocity field and magnetic field */
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->B1i[ks][j][i] = b0x;
      pGrid->B2i[ks][j][i] = 0.0;
      pGrid->B3i[ks][j][i] = b0z;
      pGrid->U[ks][j][i].B1c = b0x;
      pGrid->U[ks][j][i].B2c = 0.0;
      pGrid->U[ks][j][i].B3c = b0z;
    }
    pGrid->B1i[ks][j][ie+1] = b0x;
  }
  for (i=is; i<=ie; i++) {
    pGrid->B2i[ks][je+1][i] = 0.0;
  }
  
  /* insert charged particles */
  

  /* identifiers for particle dumps */
  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]);
  mytype = par_geti_def("problem","mytype",0);
  ntrack = 1;

  /* get grid size */
  x1min = pGrid->MinX[0];
  x2min = pGrid->MinX[1];
  
  /* get particle number */
  //  Npar  = (long)(par_geti("particle","parnumgrid"));
  
  Npar  = (long)(sqrt(par_geti("particle","parnumcell")));
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
      
      for (ip=0;ip<Npar;ip++)
	    {        
	      rx1 = x1l+pGrid->dx1/((double)Npar)*((double)ip+0.5);
        
	      for (jp=0;jp<Npar;jp++)
        {         
          rx2 = x2l+pGrid->dx2/((double)Npar)*((double)jp+0.5);
          
          pq = (&(pGrid->particle[p]));
          
          pq->x1 = rx1;
          pq->x2 = rx2;
          pq->x3 = rx3;
          
          dv1  = 0.5*vflow*(tanh((rx2-0.25*Ly)/lshear) - tanh((rx2-0.75*Ly)/lshear) - 1.0);
          dv2  = amp*sin(2.0*PI*rx1/lmode)*exp(-SQR(rx2-0.25*Ly)/(Ly*Ly*sigma*sigma)) + amp*sin(2.0*PI*rx1/lmode)*exp(-SQR(rx2-0.75*Ly)/(Ly*Ly*sigma*sigma));
          dv3  = 0.0;
          
          /* draw velocities from full distribution function */
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
  if (mytype != 0) pGrid->particle[0].property = 1;
      
  /* enroll new history variables, only once  */
  
  if (frst == 1) {
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(Delta, "<Delta>");
    frst = 0;
  }
  
  free_2d_array(rv);
  
  return;
}

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
  
  /* enroll new history variables */
  
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
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
