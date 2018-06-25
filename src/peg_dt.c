#include "copyright.h"
/*============================================================================*/
/*! \file peg_dt.c
 *  \brief Computes timestep for simulation.
 *
 * PURPOSE: Computes timestep for simulations
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - peg_dt() - computes dt						      */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "pegasus.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void peg_dt(MeshS *pM)
 *  \brief Computes timestep */

void peg_dt(MeshS *pM)
{
  GridS *pG;
  GrainS *gr;
  int nl,nd,i,j,k;
  long q;
  Real di,v1,v2,v3,qsq,b1,b2,b3,bsq,vh1,vh2,vh3;
  Real max_v1=0.0,max_v2=0.0,max_v3=0.0,max_dti=0.0,max_b=0.0,max_di=0.0;
  Real new_dt;
#ifdef MPI_PARALLEL
  double my_dt,dt;
  int ierr;
#endif
/*
#ifdef DELTA_F
  Real f0_t,particle_df,max_df=0.0;
#endif
*/

/* Check if timestep is small enough */
  
  for (nl=0; nl<(pM->NLevels); nl++) {
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
      
  if (pM->Domain[nl][nd].Grid != NULL) {
    pG=(pM->Domain[nl][nd].Grid);
        
    for (k=pG->ks; k<=pG->ke; k++) {
      for (j=pG->js; j<=pG->je; j++) {
        for (i=pG->is; i<=pG->ie; i++) {
          di = 1.0/(pG->Coup[k][j][i].grid_d);
          v1 = pG->Coup[k][j][i].grid_M1*di;
          v2 = pG->Coup[k][j][i].grid_M2*di;
          v3 = pG->Coup[k][j][i].grid_M3*di;
          qsq = v1*v1 + v2*v2 + v3*v3;

          max_di = MAX(di,max_di);
          
          b1 = pG->U[k][j][i].B1c
             + fabs((double)(pG->B1i[k][j][i] - pG->U[k][j][i].B1c));
          b2 = pG->U[k][j][i].B2c
             + fabs((double)(pG->B2i[k][j][i] - pG->U[k][j][i].B2c));
          b3 = pG->U[k][j][i].B3c
             + fabs((double)(pG->B3i[k][j][i] - pG->U[k][j][i].B3c));
          bsq = b1*b1 + b2*b2 + b3*b3;
          max_b = MAX(max_b,sqrt((double)bsq));
           
          if (pG->Nx[0] > 1) {
            vh1 = sqrt((double)bsq)*di/pG->dx1;
            max_v1 = MAX(max_v1,fabs(v1)+sqrt((double)bsq)+fabs(vh1));
          }
          if (pG->Nx[1] > 1) {
            vh2 = sqrt((double)bsq)*di/pG->dx2;
            max_v2 = MAX(max_v2,fabs(v2)+sqrt((double)bsq)+fabs(vh2));
          }
          if (pG->Nx[2] > 1) {
            vh3 = sqrt((double)bsq)*di/pG->dx3;
            max_v3 = MAX(max_v3,fabs(v3)+sqrt((double)bsq)+fabs(vh3));
          }

        }
      }
    }

    for (q=0; q<pG->nparticle; q++) {
      gr = &(pG->particle[q]);
      if (pG->Nx[0] > 1)
        max_v1 = MAX(max_v1, CourNo*gr->v1);
      if (pG->Nx[1] > 1)
        max_v2 = MAX(max_v2, CourNo*gr->v2);
      if (pG->Nx[2] > 1)
        max_v3 = MAX(max_v3, CourNo*gr->v3);
      /*
#ifdef DELTA_F
      getdf(gr, &f0_t);
      particle_df = fabs(1.0-(f0_t/gr->f_0));
      max_df = MAX(max_df, particle_df);
#endif
      */
    }

    if (pG->Nx[0] > 1)
      max_dti = MAX(max_dti, max_v1/pG->dx1);
    if (pG->Nx[1] > 1)
      max_dti = MAX(max_dti, max_v2/pG->dx2);
    if (pG->Nx[2] > 1)
      max_dti = MAX(max_dti, max_v3/pG->dx3);
    
    max_dti = MAX(max_dti, 0.5*max_b/PI);

  }}} /*--- End loop over Domains --------------------------------------------*/

 new_dt = CourNo/max_dti;

#ifdef MPI_PARALLEL
  my_dt = new_dt;
  ierr = MPI_Allreduce(&my_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  new_dt = dt;
#endif

/*  
#ifdef DELTA_F
  peg_pout(0,"maximum particle df/f=%e\n",max_df);
#endif
*/

  if (new_dt >= pM->dt) {
    peg_pout(0,"Timestep can be up to %f.\n",new_dt);
  } else {
    peg_pout(0,"Timestep needs to be reduced below %f.\n",new_dt);
    peg_perr(-1,"max_b = %f, max_di = %f, max_v1 = %f, max_v2 = %f, max_v3 = %f.\n",
    max_b,max_di,max_v1,max_v2,max_v3);
    exit(EXIT_FAILURE);
  }

  return;
}
