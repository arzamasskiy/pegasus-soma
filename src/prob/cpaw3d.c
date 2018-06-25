#include "copyright.h"
/*============================================================================*/
/*! \file hyb_cpaw3d.c
 *  \brief Problem generator for circularly polarized Alfven wave in 3D with
 *         hybrid integrator
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

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
static SolnS ***RootSoln=NULL;

/* Parameters which define initial solution -- made global so that they can be
 * shared with functions A1,2,3 which compute vector potentials */
static Real b_par, b_perp;
static Real ang_2, ang_3;
static Real fac, sin_a2, cos_a2, sin_a3, cos_a3;
static Real lambda, k_par;
static Real Lx,Ly,Lz;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * A1() - 1-component of vector potential for initial conditions
 * A2() - 2-component of vector potential for initial conditions
 * A3() - 3-component of vector potential for initial conditions
 *============================================================================*/

static Real A1(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A3(const Real x1, const Real x2, const Real x3);

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
  SolnS ***Soln;
  GrainS *pq;
  int i,is = pGrid->is, ie = pGrid->ie;
  int j,js = pGrid->js, je = pGrid->je;
  int k,ks = pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs,ip,jp,kp;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real dx1 = pGrid->dx1, dx2 = pGrid->dx2, dx3 = pGrid->dx3;
  Real hdx1 = 0.5*pGrid->dx1, hdx2 = 0.5*pGrid->dx2, hdx3 = 0.5*pGrid->dx3;
  Real x1,x2,x3,x1l,x1u,x2l,x2u,x3l,x3u;
  Real amp;
#ifdef DELTA_F
  Real df;
#endif
  static int frst=1;  /* flag so new history variables enrolled only once */
  Real L1,L2,L3,x1min,x2min,x3min,rx1,rx2,rx3,dv1,dv2,dv3,cs,sn;
  Real v_par,omega_r,omega_l,v_perp,x_par;
  long p,q,Npar3,NparGrid,NparTot;
  int  Npar,n,dir,nx1,nx2,nx3;
#ifdef MPI_PARALLEL
  int mpierr;
#endif
  Real **rv;
  
  

  nx1 = (ie - is + 1) + 2*nghost;
  nx2 = (je - js + 1) + 2*nghost;
  nx3 = (ke - ks + 1) + 2*nghost;

  if(pGrid->Nx[2] <= 1)
    peg_error("[problem]: hyb_cpaw3d assumes a 3D grid\n");
  
  if ((Soln = (SolnS***)calloc_3d_array(nx3,nx2,nx1,sizeof(SolnS))) == NULL)
    peg_error("[problem]: Error allocating memory for Soln\n");

  if (pDomain->Level == 0){
    if ((RootSoln =(SolnS***)calloc_3d_array(nx3,nx2,nx1,sizeof(SolnS)))==NULL)
      peg_error("[problem]: Error allocating memory for RootSoln\n");
  }

/* Read problem parameters.  Note Omega_0 set to 1.0 by default */
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
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
  
#ifdef DELTA_F
  setbg(pGrid);
#endif

/* Imposing periodicity and one wavelength along each grid direction */
  
  ang_3  = atan(Lx/Ly);
  //ang_3  = 0.0;
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);
  
  ang_2  = atan(0.5*(Lx*cos_a3 + Ly*sin_a3)/Lz);
  //ang_2  = 0.0;
  sin_a2 = sin(ang_2);
  cos_a2 = cos(ang_2);
  
  x1 = Lx*cos_a2*cos_a3;
  x2 = Ly*cos_a2*sin_a3;
  x3 = Lz*sin_a2;
  

/* For lambda choose the smaller of the 3 */
  
  lambda = x1;
  lambda = MIN(lambda,x2);
  lambda = MIN(lambda,x3);

/* Initialize k_parallel */
  
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


  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        
        x_par = cos_a2*(x1*cos_a3 + x2*sin_a3) + x3*sin_a2;
        
        sn = sin(k_par*x_par);
        cs = cos(k_par*x_par) * fac;
        
        Soln[k][j][i].d  = 1.0;
        Soln[k][j][i].M1 = v_par*cos_a2*cos_a3 + v_perp*sn*sin_a3 +
                           v_perp*cs*sin_a2*cos_a3;
        Soln[k][j][i].M2 = v_par*cos_a2*sin_a3 - v_perp*sn*cos_a3 +
                           v_perp*cs*sin_a2*sin_a3;
        Soln[k][j][i].M3 = v_par*sin_a2 - v_perp*cs*cos_a2;

/* Initialize magnetic field */

        x1 -= 0.5*pGrid->dx1;
        x2 -= 0.5*pGrid->dx2;
        x3 -= 0.5*pGrid->dx3;

        pGrid->B1i[k][j][i] = (A3(x1,x2+dx2 ,x3+hdx3) - A3(x1,x2,x3+hdx3))/dx2 -
                              (A2(x1,x2+hdx2,x3+dx3 ) - A2(x1,x2+hdx2,x3))/dx3;

        pGrid->B2i[k][j][i] = (A1(x1+hdx1,x2,x3+dx3 ) - A1(x1+hdx1,x2,x3))/dx3 -
                              (A3(x1+dx1 ,x2,x3+hdx3) - A3(x1,x2,x3+hdx3))/dx1;

        pGrid->B3i[k][j][i] = (A2(x1+dx1,x2+hdx2,x3) - A2(x1,x2+hdx2,x3))/dx1 -
                              (A1(x1+hdx1,x2+dx2,x3) - A1(x1+hdx1,x2,x3))/dx2;

      }
    }
  }
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        Soln[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
        Soln[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
        Soln[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
        pGrid->U[k][j][i].B1c = Soln[k][j][i].B1c;
        pGrid->U[k][j][i].B2c = Soln[k][j][i].B2c;
        pGrid->U[k][j][i].B3c = Soln[k][j][i].B3c;
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
  //  NparGrid  = (long)(par_geti("particle","parnumgrid"));

  Npar = (int)(round(pow(par_geti("particle","parnumcell"),1.0/3.0)));
  Npar3 = Npar*SQR(Npar);

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

  initdf(&iseed,rv,NparGrid,NparTot);
  
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
      
        for (ip=0;ip<Npar;ip++)
        {
          rx1 = x1l+pGrid->dx1/Npar*(ip+0.5);

          for (jp=0;jp<Npar;jp++)
          {
            rx2 = x2l+pGrid->dx2/Npar*(jp+0.5);

            for (kp=0;kp<Npar;kp++)
            {
              rx3 = x3l+pGrid->dx3/Npar*(kp+0.5);

              x_par = cos_a2*(rx1*cos_a3 + rx2*sin_a3) + rx3*sin_a2;

              sn = sin(k_par*x_par);
              cs = cos(k_par*x_par) * fac;
              
              /* compute perturbation */
              dv1 = v_par*cos_a2*cos_a3 + v_perp*sn*sin_a3 + v_perp*cs*sin_a2*cos_a3;
              dv2 = v_par*cos_a2*sin_a3 - v_perp*sn*cos_a3 + v_perp*cs*sin_a2*sin_a3;
              dv3 = v_par*sin_a2 - v_perp*cs*cos_a2;
#ifdef SHEARING_BOX
#ifndef FARGO
              dv2 -= Shear_0*rx1;
#endif
#endif
              
              pq = &(pGrid->particle[p]);
              
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

  /* save solution on root grid */

  if (pDomain->Level == 0) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          RootSoln[k][j][i].d  = Soln[k][j][i].d ;
          RootSoln[k][j][i].M1 = Soln[k][j][i].M1;
          RootSoln[k][j][i].M2 = Soln[k][j][i].M2;
          RootSoln[k][j][i].M3 = Soln[k][j][i].M3;
          RootSoln[k][j][i].B1c = Soln[k][j][i].B1c;
          RootSoln[k][j][i].B2c = Soln[k][j][i].B2c;
          RootSoln[k][j][i].B3c = Soln[k][j][i].B3c;
        }}}
  }

  /* finalize */
  free_2d_array(rv);
  free_3d_array(Soln);

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
/* enroll new history variables */

  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");

/* enroll gravitational potential function */

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

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k, Real *eta_0)
{
  *eta_0 = 0.0;

  return;
}
#endif

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
  GridS *pGrid;
  DomainS *pD;
  GPCouple *pq;
  int i,j,k,is,ie,js,je,ks,ke,Nx1,Nx2,Nx3,count;
  Real rms_error=0.0;
  SolnS error,total_error;
  FILE *fp;
  char *fname;
#if defined MPI_PARALLEL
  double err[8], tot_err[8];
  int ierr,myID;
#endif
  total_error.d = 0.0;
  total_error.M1 = 0.0;
  total_error.M2 = 0.0;
  total_error.M3 = 0.0;
  total_error.B1c = 0.0;
  total_error.B2c = 0.0;
  total_error.B3c = 0.0;

  /* Compute error only on root Grid, which is in Domain[0][0] */
  pD = (DomainS*)&(pM->Domain[0][0]);
  pGrid = pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  
  deposit(pD,0);

  /* compute L1 error in each variable, and rms total error */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
    error.B1c = 0.0;
    error.B2c = 0.0;
    error.B3c = 0.0;
    for (i=is; i<=ie; i++) {
      pq = &(pGrid->Coup[k][j][i]);
      error.d   += fabs(pq->grid_d  - RootSoln[k][j][i].d  );
      error.M1  += fabs(pq->grid_M1 - RootSoln[k][j][i].M1 );
      error.M2  += fabs(pq->grid_M2 - RootSoln[k][j][i].M2 );
      error.M3  += fabs(pq->grid_M3 - RootSoln[k][j][i].M3 );
      error.B1c += fabs(pGrid->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
      error.B2c += fabs(pGrid->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
      error.B3c += fabs(pGrid->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
    }
    
    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
    total_error.B1c += error.B1c;
    total_error.B2c += error.B2c;
    total_error.B3c += error.B3c;
  }}
  
#if defined MPI_PARALLEL
  Nx1 = pM->Domain[0][0].Nx[0];
  Nx2 = pM->Domain[0][0].Nx[1];
  Nx3 = pM->Domain[0][0].Nx[2];
#else
  Nx1 = (ie-is+1);
  Nx2 = (je-js+1);
  Nx3 = (ke-ks+1);
#endif
  count = Nx1*Nx2*Nx3;

#ifdef MPI_PARALLEL 
  /* Now we have to use an All_Reduce to get the total error over all the MPI
   * grids.  Begin by copying the error into the err[] array */
  
  err[0] = total_error.d;
  err[1] = total_error.M1;
  err[2] = total_error.M2;
  err[3] = total_error.M3;
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
  
  ierr = MPI_Reduce(err,tot_err,7,MPI_DOUBLE,MPI_SUM,0,
                    pM->Domain[0][0].Comm_Domain);
  
  /* If I'm the parent, copy the sum back to the total_error variable */
  
  ierr = MPI_Comm_rank(pM->Domain[0][0].Comm_Domain, &myID);
  if(myID == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
  }
  else return; /* The child grids do not do any of the following code */
  
#endif /* MPI_PARALLEL */
  
  /* Compute RMS error over all variables, and print out */
  
  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
            + SQR(total_error.M3);
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
             + SQR(total_error.B3c);
  rms_error = sqrt(rms_error)/(double)count;
  
  /* Print error to file "hyb_cpaw3d-errors.dat" */
  
#ifdef MPI_PARALLEL
  fname = "../cpaw3d-errors.dat";
#else
  fname = "cpaw3d-errors.dat";
#endif
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
  
  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);
  fprintf(fp,"  %e  %e  %e  %e",
          (total_error.d/(double)count),
          (total_error.M1/(double)count),
          (total_error.M2/(double)count),
          (total_error.M3/(double)count) );
  fprintf(fp,"  %e  %e  %e",
          (total_error.B1c/(double)count),
          (total_error.B2c/(double)count),
          (total_error.B3c/(double)count));
  fprintf(fp,"\n");
  
  fclose(fp);
  
  free_3d_array(RootSoln);

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

/*----------------------------------------------------------------------------*/
/*! \fn static Real A1(const Real x1, const Real x2, const Real x3)
 *  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0,
 * and Ay, Az are functions of x and y alone.
 */

static Real A1(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Ay, Az;
  
  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;
  
  Ay = fac*(b_perp/k_par)*sin(k_par*x);
  Az = (b_perp/k_par)*cos(k_par*x) + b_par*y;
  
  return -Ay*sin_a3 - Az*sin_a2*cos_a3;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real A2(const Real x1, const Real x2, const Real x3)
 *  \brief A2: 2-component of vector potential
 */

static Real A2(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Ay, Az;
  
  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;
  
  Ay = fac*(b_perp/k_par)*sin(k_par*x);
  Az = (b_perp/k_par)*cos(k_par*x) + b_par*y;
  
  return Ay*cos_a3 - Az*sin_a2*sin_a3;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real A3(const Real x1, const Real x2, const Real x3)
 *  \brief A3: 3-component of vector potential
 */

static Real A3(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Az;
  
  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;
  
  Az = (b_perp/k_par)*cos(k_par*x) + b_par*y;
  
  return Az*cos_a2;
}
