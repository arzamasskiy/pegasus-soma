#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_3d_hybrid.c
 *  \brief Integrate coupled induction equation and particles
 *
 * PURPOSE: 
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_3d_hybrid()
 * - integrate_init_3d()
 * - integrate_destruct_3d()
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "../defs.h"
#include "../pegasus.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "../particles/particle.h"

/* The interface electric fields; cell-centered electric fields */
static Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
static Real ***B1i=NULL, ***B2i=NULL, ***B3i=NULL;
static Real ***B1c=NULL, ***B2c=NULL, ***B3c=NULL;
static Real ***J1i=NULL, ***J2i=NULL, ***J3i=NULL;


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   calculate_emfs() - calculates electric field
 *   putmhdvalues() - puts E & B info in particle coupling structure
 *============================================================================*/

static void putmhdvalues(GridS *pG);
static void calculate_emfs(const GridS *pG, short lab);
static void calculate_currents(const GridS *pG, short lab);
static void clear_emfs(void);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_3d_hybrid(DomainS *pD)
 *  \brief  */

void integrate_3d_hybrid(DomainS *pD)
{
  GridS *pG = pD->Grid;
  Real q1 = pG->dt/pG->dx1;
  Real q2 = pG->dt/pG->dx2;
  Real q3 = pG->dt/pG->dx3;
  int i,is = pG->is, ie = pG->ie;
  int j,js = pG->js, je = pG->je;
  int k,ks = pG->ks, ke = pG->ke;
  int my_iproc, my_jproc, my_kproc;
  struct timeval tv_init, tvs;
  double step_time;  
/*=== Step 0: Preparations ===================================================*/

  gettimeofday(&tvs,NULL);
  tv_init=tvs;

/*--- Step 0a ------------------------------------------------------------------
 * Copy face- and cell-centered magnetic fields at t^{n} into local arrays
 */
 
  for (k=klp; k<=kup; k++) {
    for (j=jlp; j<=jup; j++) {
      for (i=ilp; i<=iup; i++) {
        B1i[k][j][i] = pG->B1i[k][j][i];
        B2i[k][j][i] = pG->B2i[k][j][i];
        B3i[k][j][i] = pG->B3i[k][j][i];
        B1c[k][j][i] = pG->U[k][j][i].B1c;
        B2c[k][j][i] = pG->U[k][j][i].B2c;
        B3c[k][j][i] = pG->U[k][j][i].B3c;
      }
    }
  }
  
/*=== Step 1: Compute values at t^{n+1/2} ====================================*/
  
/*--- Step 1a ------------------------------------------------------------------
 * Deposit particles on grid, correct bndry values, and set bndry conditions.
 *
 * gives grid[ks-5:ke+5][js-5:je+5][is-5:ie+5]
 * for Steps 1b, 1e
 */

// moved to main.c
//  deposit(pD,1);

/*--- Step 1b ------------------------------------------------------------------
 * Calculate cell- and edge-centered electric fields at t^{n}. Store them
 * as grid quantities for use later in Step 1j.
 *
 * requires grid[ks-5:ke+5][js-5:je+5][is-5:ie+5]
 * from Step 1a
 *
 * requires B1c [ks-5:ke+5][js-5:je+5][is-4:ie+4]
 *          B2c [ks-5:ke+5][js-4:je+4][is-5:ie+5]
 *          B3c [ks-4:ke+4][js-5:je+5][is-5:ie+5]
 *          B1i [ks-5:ke+5][js-5:je+5][is-4:ie+5]
 *          B2i [ks-5:ke+5][js-4:je+5][is-5:ie+5]
 *          B3i [ks-4:ke+5][js-5:je+5][is-5:ie+5]
 * from Step 0a
 *
 * gives emf1[ks-4:ke+5][js-4:je+5][is-4:ie+4]
 *       emf2[ks-4:ke+5][js-4:je+4][is-4:ie+5]
 *       emf3[ks-4:ke+4][js-4:je+5][is-4:ie+5]
 * for Step 1c
 *
 * gives emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for storage
 *
 * stores emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *        emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *        emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1k
 *
 * stores emf1[ks:ke+1][js:je+1][is:ie  ]
 *        emf2[ks:ke+1][js:je  ][is:ie+1]
 *        emf3[ks:ke  ][js:je+1][is:ie+1]
 * for Step 1k
 */
  
  calculate_currents(pG,0);

  calculate_emfs(pG,0);


/*--- Step 1c ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 *
 * requires emf1[ks-4:ke+5][js-4:je+5][is-4:ie+4]
 *          emf2[ks-4:ke+5][js-4:je+4][is-4:ie+5]
 *          emf3[ks-4:ke+4][js-4:je+5][is-4:ie+5]
 * from Step 1b
 *
 * gives B1i[ks-4:ke+4][js-4:je+4][is-3:ie+4]
 *       B2i[ks-4:ke+4][js-3:je+4][is-4:ie+4]
 *       B3i[ks-3:ke+4][js-4:je+4][is-4:ie+4]
 * for Steps 1d and 1e
 */

  for (k=ks-4; k<=ke+4; k++) {
    for (j=js-4; j<=je+4; j++) {
      for (i=is-4; i<=ie+4; i++) {
        B1i[k][j][i]  += q3 * (emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                         q2 * (emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        B2i[k][j][i]  += q1 * (emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                         q3 * (emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        B3i[k][j][i]  += q2 * (emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                         q1 * (emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
    }
  }

/*--- Step 1d ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 *
 * requires B1i[ks-4:ke+4][js-4:je+4][is-3:ie+4]
 *          B2i[ks-4:ke+4][js-3:je+4][is-4:ie+4]
 *          B3i[ks-3:ke+4][js-4:je+4][is-4:ie+4]
 * from Step 1c
 *
 * gives B1c[ks-4:ke+4][js-4:je+4][is-3:ie+3]
 *       B2c[ks-4:ke+4][js-3:je+3][is-4:ie+4]
 *       B3c[ks-3:ke+3][js-4:je+4][is-4:ie+4]
 * for Step 1e
 *
 * gives B1c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B2c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B3c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1f
 */
  
  for (k=ks-4; k<=ke+4; k++) {
    for (j=js-4; j<=je+4; j++) {
      for (i=is-4; i<=ie+4; i++) {
        B1c[k][j][i] = 0.5 * ( B1i[k][j][i] + B1i[k][j][i+1] );
        B2c[k][j][i] = 0.5 * ( B2i[k][j][i] + B2i[k][j+1][i] );
        B3c[k][j][i] = 0.5 * ( B3i[k][j][i] + B3i[k+1][j][i] );
      }
    }
  }
  
/*--- Step 1e ------------------------------------------------------------------
 * Re-calculate cell- and edge-centered electric fields using B^{n+1}. Then 
 * average to get electric fields at t^{n+1/2}.
 *
 * requires grid[ks-4:ke+4][js-4:je+4][is-4:ie+4]
 * from Step 1a
 *
 * requires B1c[ks-4:ke+4][js-4:je+4][is-3:ie+3]
 *          B2c[ks-4:ke+4][js-3:je+3][is-4:ie+4]
 *          B3c[ks-3:ke+3][js-4:je+4][is-4:ie+4]
 * from Step 1d
 *
 * requires B1i[ks-4:ke+4][js-4:je+4][is-3:ie+4]
 *          B2i[ks-4:ke+4][js-3:je+4][is-4:ie+4]
 *          B3i[ks-3:ke+4][js-4:je+4][is-4:ie+4]
 * from Step 1c
 *
 * gives emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1g
 *
 * gives emf1[ks-3:ke+4][js-3:je+4][js-3:je+3]
 *       emf2[ks-3:ke+4][js-3:je+3][js-3:je+4]
 *       emf3[ks-3:ke+3][js-3:je+4][js-3:je+4]
 * for Step 1i
 */
  
  calculate_currents(pG,1);

  calculate_emfs(pG,1);


/*--- Step 1f ------------------------------------------------------------------
 * Time-average to get magnetic fields at t^{n+1/2}.
 *
 * requires B1c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B2c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B3c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * from Step 1d
 *
 * gives B1c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B2c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B3c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1f
 */  
  
  for (k=ks-4; k<=ke+4; k++) {
    for (j=js-4; j<=je+4; j++) {
      for (i=is-4; i<=ie+4; i++) {
        B1c[k][j][i] = 0.5 * ( pG->U[k][j][i].B1c + B1c[k][j][i] );
        B2c[k][j][i] = 0.5 * ( pG->U[k][j][i].B2c + B2c[k][j][i] );
        B3c[k][j][i] = 0.5 * ( pG->U[k][j][i].B3c + B3c[k][j][i] );
      }
    }
  }
  
/*--- Step 1g ------------------------------------------------------------------
 * Put cell-centered electric and magnetic fields in particle coupling array.
 *
 * requires emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B1c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B2c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B3c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * from Steps 1e and 1f
 *
 * gives emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B1c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B2c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B3c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1h
 */
  
  putmhdvalues(pG);
  
/*--- Step 1h ------------------------------------------------------------------
 * Integrate temporary copy of particles to t^{n+1} and deposit on grid.
 *
 * requires emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B1c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B2c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B3c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * from Step 1g
 *
 * gives grid[ks-3:ke+3][js-3:je+3][is-3:ie+3]
 * for Step 1k
 */

  gettimeofday(&tvs,NULL);
  step_time = (double)(tvs.tv_sec - tv_init.tv_sec) +
      1.0e-6*(double)(tvs.tv_usec - tv_init.tv_usec);
  peg_pout(0,"time for Maxwell solver before first push       --->%e s\n",step_time);
  tv_init=tvs;  

  Integrate_Particles_1(pD);

  gettimeofday(&tvs,NULL);
  step_time = (double)(tvs.tv_sec - tv_init.tv_sec) +
      1.0e-6*(double)(tvs.tv_usec - tv_init.tv_usec);
  peg_pout(0,"time for predictor particle push                --->%e s\n",step_time);
  tv_init=tvs;  

/*--- Step 1i ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 *
 * requires emf1[ks-3:ke+4][js-3:je+4][js-3:je+3]
 *          emf2[ks-3:ke+4][js-3:je+3][js-3:je+4]
 *          emf3[ks-3:ke+3][js-3:je+4][js-3:je+4]
 * from Step 1e
 *
 * gives B1i[ks-3:ke+3][js-3:je+3][is-2:ie+3]
 *       B2i[ks-3:ke+3][js-2:je+3][is-3:ie+3]
 *       B3i[ks-2:ke+3][js-3:je+3][is-3:ie+3]
 * for Step 1j
 *
 * gives B1i[ks-1:ke+1][js-1:je+1][is  :ie+1]
 *       B2i[ks-1:ke+1][js  :je+1][is-1:ie+1]
 *       B3i[ks  :ke+1][js-1:je+1][is-1:ie+1]
 * for Step 1k
 */

  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
        B1i[k][j][i] = pG->B1i[k][j][i] +
                       q3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                       q2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        B2i[k][j][i] = pG->B2i[k][j][i] +
                       q1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                       q3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        B3i[k][j][i] = pG->B3i[k][j][i] + 
                       q2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                       q1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
    }
  }
  
/*--- Step 1j ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 *
 * requires B1i[ks-3:ke+3][js-3:je+3][is-2:ie+3]
 *          B2i[ks-3:ke+3][js-2:je+3][is-3:ie+3]
 *          B3i[ks-2:ke+3][js-3:je+3][is-3:ie+3]
 * from Step 1i
 *
 * gives B1c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B2c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B3c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1l
 *
 * gives B1c[ks-3:ke+3][js-3:je+3][is-2:ie+2]
 *       B2c[ks-3:ke+3][js-2:je+2][is-3:ie+3]
 *       B3c[ks-2:ke+2][js-3:je+3][is-3:ie+3]
 * for Step 1k
 */
  
  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
        B1c[k][j][i] = 0.5 * ( B1i[k][j][i] + B1i[k][j][i+1] );
        B2c[k][j][i] = 0.5 * ( B2i[k][j][i] + B2i[k][j+1][i] );
        B3c[k][j][i] = 0.5 * ( B3i[k][j][i] + B3i[k+1][j][i] );
      }
    }
  }
  
/*--- Step 1k ------------------------------------------------------------------
 * Calculate cell- and edge-centered electric field at t^{n+1}. Then average
 * to get electric fields at t^{n+1/2}.
 *
 * requires grid[ks-3:ke+3][js-3:je+3][is-3:ie+3]
 * from Step 1h
 *
 * requires B1i[ks-1:ke+1][js-1:je+1][is  :ie+1]
 *          B2i[ks-1:ke+1][js  :je+1][is-1:ie+1]
 *          B3i[ks  :ke+1][js-1:je+1][is-1:ie+1]
 * from Step 1i
 *
 * requires B1c[ks-3:ke+3][js-3:je+3][is-2:ie+2]
 *          B2c[ks-3:ke+3][js-2:je+2][is-3:ie+3]
 *          B3c[ks-2:ke+2][js-3:je+3][is-3:ie+3]
 * from Step 1j
 *
 * gives emf1[ks:ke+1][js:je+1][is:ie  ]
 *       emf2[ks:ke+1][js:je  ][is:ie+1]
 *       emf3[ks:ke  ][js:je+1][is:ie+1]
 * for Step 2b
 *
 * gives emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1m
 */
  
  calculate_currents(pG,2);
  
  calculate_emfs(pG,2);

/*--- Step 1l ------------------------------------------------------------------
 * Time-average to get magnetic fields at t^{n+1/2}.
 *
 * requires B1c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B2c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B3c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * from Step 1j
 *
 * gives B1c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B2c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B3c[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 1m
 */  

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        B1c[k][j][i] = 0.5 * ( pG->U[k][j][i].B1c + B1c[k][j][i] );
        B2c[k][j][i] = 0.5 * ( pG->U[k][j][i].B2c + B2c[k][j][i] );
        B3c[k][j][i] = 0.5 * ( pG->U[k][j][i].B3c + B3c[k][j][i] );
      }
    }
  }
  
/*--- Step 1m ------------------------------------------------------------------
 * Put cell-centered electric and magnetic fields in particle coupling array.
 *
 * requires emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B1c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B2c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B3c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * from Steps 1k and 1l
 *
 * gives emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B1c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B2c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *       B3c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * for Step 2a
 */
  
  putmhdvalues(pG);
  
/*=== Step 2: Compute values at t^{n+1} ======================================*/
  
/*--- Step 2a ------------------------------------------------------------------
 * Integrate particles to t^{n+1}.
 *
 * requires emf1_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf2_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          emf3_cc[ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B1c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B2c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 *          B3c    [ks-2:ke+2][js-2:je+2][is-2:ie+2]
 * -3:+3 if injector used
 *
 * from Step 1m
 */
  
  gettimeofday(&tvs,NULL);
  step_time = (double)(tvs.tv_sec - tv_init.tv_sec) +
      1.0e-6*(double)(tvs.tv_usec - tv_init.tv_usec);
  peg_pout(0,"time for Maxwell solver before second push      --->%e s\n",step_time);
  tv_init=tvs;

  Integrate_Particles_2(pD);
  
  gettimeofday(&tvs,NULL);
  step_time = (double)(tvs.tv_sec - tv_init.tv_sec) +
      1.0e-6*(double)(tvs.tv_usec - tv_init.tv_usec);
  peg_pout(0,"time for corrector particle push                --->%e s\n",step_time);
  tv_init=tvs;

/*--- Step 2b ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 *
 * NOTE: In order to conserve Bz, the edge-centered emf_y must be 
 *       averaged with a remapped value from the other side of the domain.
 *
 * requires emf1[ks:ke+1][js:je+1][is:ie  ]
 *          emf2[ks:ke+1][js:je  ][is:ie+1]
 *          emf3[ks:ke  ][js:je+1][is:ie+1]
 * from Step 1l
 *
 * gives pG->B1i[ks:ke  ][js:je  ][is:ie+1]
 *       pG->B2i[ks:ke  ][js:je+1][is:ie  ]
 *       pG->B3i[ks:ke+1][js:je  ][is:ie  ]
 * for Step 2c
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->B1i[k][j][i] += 
                       q3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                       q2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        pG->B2i[k][j][i] +=
                       q1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                       q3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        pG->B3i[k][j][i] +=
                       q2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                       q1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
      pG->B1i[k][j][ie+1] +=
                        q3*(emf2[k+1][j  ][ie+1] - emf2[k][j][ie+1]) -
                        q2*(emf3[k  ][j+1][ie+1] - emf3[k][j][ie+1]);
    }
    for (i=is; i<=ie; i++) {
      pG->B2i[k][je+1][i] +=
                        q1*(emf3[k  ][je+1][i+1] - emf3[k][je+1][i]) -
                        q3*(emf1[k+1][je+1][i  ] - emf1[k][je+1][i]);
    }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B3i[ke+1][j][i] += 
                        q2*(emf1[ke+1][j+1][i  ] - emf1[ke+1][j][i]) -
                        q1*(emf2[ke+1][j  ][i+1] - emf2[ke+1][j][i]);
    }
  }
  
/*--- Step 2c ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 *
 * requires pG->B1i[ks:ke  ][js:je  ][is:ie+1]
 *          pG->B2i[ks:ke  ][js:je+1][is:ie  ]
 *          pG->B3i[ks:ke+1][js:je  ][is:ie  ]
 * from Step 2b
 */
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
      }
    }
  }
  
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

static void clear_emfs(void)
{
  int i,j,k;
  
  for (k=klp; k<=kup; k++) {
    for (j=jlp; j<=jup; j++) {
      for (i=ilp; i<=iup; i++) {
        emf1[k][j][i] = 0.0;
        emf2[k][j][i] = 0.0;
        emf3[k][j][i] = 0.0;
        emf1_cc[k][j][i] = 0.0;
        emf2_cc[k][j][i] = 0.0;
        emf3_cc[k][j][i] = 0.0;
      }
    }
  }
  
  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void calculate_current(const GridS *pG, char *stage)
 *  \brief Calculate edge-centered current densities
 */
static void calculate_currents(const GridS *pG, short lab)
{
  int i,j,k,kb,kt,jb,jt,ib,it;
  Real dx1i = 1.0/pG->dx1, dx2i = 1.0/pG->dx2, dx3i = 1.0/pG->dx3;
  
  switch (lab) {
    case 0:
     /* need to calculate
      *   J1i[ks-4:ke+5][js-4:je+5][is-5:ie+5]
      *   J2i[ks-4:ke+5][js-5:je+5][is-4:ie+5]
      *   J3i[ks-5:ke+5][js-4:je+5][is-4:ie+5]
      */
      kb=pG->ks-4; kt=pG->ke+4;
      jb=pG->js-4; jt=pG->je+4;
      ib=pG->is-4; it=pG->ie+4;
      break;
      
    case 1:
     /* need to calculate
      *  J1i[ks-3:ke+4][js-3:je+4][is-4:ie+4]
      *  J2i[ks-3:ke+4][js-4:je+4][is-3:ie+4]
      *  J3i[ks-4:ke+4][js-3:je+4][is-3:ie+4]
      */
      kb=pG->ks-3; kt=pG->ke+3;
      jb=pG->js-3; jt=pG->je+3;
      ib=pG->is-3; it=pG->ie+3;
      break;
      
    case 2:
     /* need to calculate
      *  J1i[ks  :ke+1][js  :je+1][is-1:ie+1]
      *  J2i[ks  :ke+1][js-1:je+1][is  :ie+1]
      *  J3i[ks-1:ke+1][js  :je+1][is  :ie+1]
      */
      kb=pG->ks-2; kt=pG->ke+2;
      jb=pG->js-2; jt=pG->je+2;
      ib=pG->is-2; it=pG->ie+2;
      break;
      
    default:
      peg_error("[integrate_2d_hybrid]: stage undefined!\n");
  }
  
  /* calculate_currents gives:
   *   J1i[kb  :kt+1][jb  :jt+1][ib-1:it+1]
   *   J2i[kb  :kt+1][jb-1:jt+1][ib  :it+1]
   *   J3i[kb-1:kt+1][jb  :jt+1][ib  :it+1]
   * these require:
   *   B1i[kb-1:kt+1][jb-1:jt+1][ib  :it+1]
   *   B2i[kb-1:kt+1][jb  :jt+1][ib-1:it+1]
   *   B3i[kb  :kt+1][jb-1:jt+1][ib-1:it+1]
   */
  
  for (k=kb; k<=kt+1; k++) {
    for (j=jb; j<=jt+1; j++) {
      for (i=ib; i<=it+1; i++) {
        J1i[k][j][i] =  dx2i * ( B3i[k][j][i] - B3i[k][j-1][i] )
                       -dx3i * ( B2i[k][j][i] - B2i[k-1][j][i] );
        J2i[k][j][i] =  dx3i * ( B1i[k][j][i] - B1i[k-1][j][i] )
                       -dx1i * ( B3i[k][j][i] - B3i[k][j][i-1] );
        J3i[k][j][i] =  dx1i * ( B2i[k][j][i] - B2i[k][j][i-1] )
                       -dx2i * ( B1i[k][j][i] - B1i[k][j-1][i] );
      }
      i=ib-1;
      J1i[k][j][i] =  dx2i * ( B3i[k][j][i] - B3i[k][j-1][i] )
                     -dx3i * ( B2i[k][j][i] - B2i[k-1][j][i] );
    }
    j=jb-1;
    for (i=ib; i<=it+1; i++) {
      J2i[k][j][i] =  dx3i * ( B1i[k][j][i] - B1i[k-1][j][i] )
                     -dx1i * ( B3i[k][j][i] - B3i[k][j][i-1] );
    }
  }
  k=kb-1;
  for (j=jb; j<=jt+1; j++) {
    for (i=ib; i<=it+1; i++) {
      J3i[k][j][i] =  dx1i * ( B2i[k][j][i] - B2i[k][j][i-1] )
                     -dx2i * ( B1i[k][j][i] - B1i[k][j-1][i] );
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void calculate_emfs(GridS *pG, char *stage)
 *  \brief Calculate cell- and edge-centered electric fields
 * 
 * NOTE: ions are not pushed by E = eta J
 *
 */

static void calculate_emfs(const GridS *pG, short lab)
{
  int k,ks=pG->ks,ke=pG->ke;
  int j,js=pG->js,je=pG->je;
  int i,is=pG->is,ie=pG->ie;
  int kb,kt,jb,jt,ib,it,il,iu,jl,ju,kl,ku;
  Real den,dinv,M1,M2,M3,vf1,vf2,vf3,j1,j2,j3,e1,e2,e3;
  Real B1u,B2u,B3u,pemf1,pemf2,pemf3;
  Real dx1i = 1.0/pG->dx1, dx2i = 1.0/pG->dx2, dx3i = 1.0/pG->dx3;
#ifdef ADIABATIC
  Real adiab;
#endif
  
  switch (lab) {
    case 0:
      kb=ks-4; kt=ke+4;
      jb=js-4; jt=je+4;
      ib=is-4; it=ie+4;
      kl=ks-2; ku=ke+2;
      jl=js-2; ju=je+2;
      il=is-2; iu=ie+2;
      break;

    case 1:
      kb=ks-3; kt=ke+3;
      jb=js-3; jt=je+3;
      ib=is-3; it=ie+3;
      if (vinject == 0) {
        il=is-2; iu=ie+2;
        jl=js-2; ju=je+2;
        kl=ks-2; ku=ke+2;
      } else {
        il=is-3; iu=ie+3;
        jl=js-3; ju=je+3;
        kl=ks-3; ku=ke+3;
      }
      break;
      
    case 2:
      kb=ks; kt=ke;
      jb=js; jt=je;
      ib=is; it=ie;
      kl=ks-2; ku=ke+2;
      jl=js-2; ju=je+2;
      il=is-2; iu=ie+2;
      break;
      
    default:
      peg_error("[integrate_3d_hybrid]: stage undefined!\n");
  }

  
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        
        j1   = 0.25 * ( J1i[k][j][i] + J1i[k][j+1][i] + J1i[k+1][j][i] + J1i[k+1][j+1][i]);
        j2   = 0.25 * ( J2i[k][j][i] + J2i[k][j][i+1] + J2i[k+1][j][i] + J2i[k+1][j][i+1]);
        j3   = 0.25 * ( J3i[k][j][i] + J3i[k][j][i+1] + J3i[k][j+1][i] + J3i[k][j+1][i+1]);
        
        dinv = 1.0/pG->Coup[k][j][i].grid_d;
        vf1  = (pG->Coup[k][j][i].grid_M1 - j1)*dinv;
        vf2  = (pG->Coup[k][j][i].grid_M2 - j2)*dinv;
        vf3  = (pG->Coup[k][j][i].grid_M3 - j3)*dinv;
        
        pemf1 = 0.5 * dx1i * dinv * ( ZTeTi * 0.5 * beta ) *
              ( pG->Coup[k][j][i+1].grid_d - pG->Coup[k][j][i-1].grid_d );
        pemf2 = 0.5 * dx2i * dinv * ( ZTeTi * 0.5 * beta ) *
              ( pG->Coup[k][j+1][i].grid_d - pG->Coup[k][j-1][i].grid_d );
        pemf3 = 0.5 * dx3i * dinv * ( ZTeTi * 0.5 * beta ) *
              ( pG->Coup[k+1][j][i].grid_d - pG->Coup[k-1][j][i].grid_d );

#ifdef ADIABATIC
        adiab = Gamma * pow(pG->Coup[k][j][i].grid_d,Gamma_1);
        pemf1 *= adiab;
        pemf2 *= adiab;
        pemf3 *= adiab;
#endif
        
        e1 = B2c[k][j][i] * vf3 - B3c[k][j][i] * vf2 - pemf1;
        e2 = B3c[k][j][i] * vf1 - B1c[k][j][i] * vf3 - pemf2;
        e3 = B1c[k][j][i] * vf2 - B2c[k][j][i] * vf1 - pemf3;
        
        switch (lab) {
          case 0:
            pG->U[k][j][i].E1c = e1;
            pG->U[k][j][i].E2c = e2;
            pG->U[k][j][i].E3c = e3;
            break;
            
          case 1:
            emf1_cc[k][j][i] = 0.5 * ( pG->U[k][j][i].E1c + e1 );
            emf2_cc[k][j][i] = 0.5 * ( pG->U[k][j][i].E2c + e2 );
            emf3_cc[k][j][i] = 0.5 * ( pG->U[k][j][i].E3c + e3 );
            break;
            
          case 2:
            emf1_cc[k][j][i] = 0.5 * ( pG->U[k][j][i].E1c + e1 );
            emf2_cc[k][j][i] = 0.5 * ( pG->U[k][j][i].E2c + e2 );
            emf3_cc[k][j][i] = 0.5 * ( pG->U[k][j][i].E3c + e3 );
            break;
            
          default:
            peg_error("[integrate_3d_hybrid]: stage undefined!\n");
        }
        
      }
    }
  }
  
  /* requires:
   *  J1i[kb  :kt+1][jb  :jt+1][ib  :it  ]
   *  J2i[kb  :kt+1][jb-1:jt+1][ib  :it+1]
   *  J3i[kb-1:kt+1][jb  :jt+1][ib  :it+1]
   */
  for (k=kb; k<=kt+1; k++) {
    for (j=jb; j<=jt+1; j++) {
      for (i=ib; i<=it; i++) {
        
        den = 0.25 * ( pG->Coup[k  ][j][i].grid_d  + pG->Coup[k  ][j-1][i].grid_d  
                     + pG->Coup[k-1][j][i].grid_d  + pG->Coup[k-1][j-1][i].grid_d  );
        M2  = 0.25 * ( pG->Coup[k  ][j][i].grid_M2 + pG->Coup[k  ][j-1][i].grid_M2
                     + pG->Coup[k-1][j][i].grid_M2 + pG->Coup[k-1][j-1][i].grid_M2 );
        M3  = 0.25 * ( pG->Coup[k  ][j][i].grid_M3 + pG->Coup[k  ][j-1][i].grid_M3 
                     + pG->Coup[k-1][j][i].grid_M3 + pG->Coup[k-1][j-1][i].grid_M3 );
      
        j2 = 0.25 * ( J2i[k][j][i] + J2i[k][j][i+1] + J2i[k][j-1][i] + J2i[k][j-1][i+1] );
        j3 = 0.25 * ( J3i[k][j][i] + J3i[k][j][i+1] + J3i[k-1][j][i] + J3i[k-1][j][i+1] );

        vf2 = (M2 - j2)/den;
        vf3 = (M3 - j3)/den;
        
        B3u = 0.5 * ( B3i[k][j][i] + B3i[k][j-1][i] );
        B2u = 0.5 * ( B2i[k][j][i] + B2i[k-1][j][i] );
        
        e1 = B2u * vf3 - B3u * vf2;
        
        switch (lab) {
          case 0:
            emf1[k][j][i] = e1;
            pG->E1i[k][j][i] = e1;
            break;
            
          case 1:
            emf1[k][j][i] = 0.5 * ( pG->E1i[k][j][i] + e1 );
            break;
            
          case 2:
            emf1[k][j][i] = 0.5 * ( pG->E1i[k][j][i] + e1 );

            break;
            
          default:
            peg_error("[integrate_3d_hybrid]: stage undefined!\n");
        }
        
      }
    }
  }
  
  /* requires:
   *  J1i[kb  :kt+1][jb  :jt+1][ib-1:it+1]
   *  J2i[kb  :kt+1][jb  :jt  ][ib  :it+1]
   *  J3i[kb-1:kt+1][jb  :jt+1][ib  :it+1]
   */
  for (k=kb; k<=kt+1; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=ib; i<=it+1; i++) {
        
        den = 0.25 * ( pG->Coup[k  ][j][i].grid_d  + pG->Coup[k  ][j][i-1].grid_d  
                     + pG->Coup[k-1][j][i].grid_d  + pG->Coup[k-1][j][i-1].grid_d  );
        M1  = 0.25 * ( pG->Coup[k  ][j][i].grid_M1 + pG->Coup[k  ][j][i-1].grid_M1 
                     + pG->Coup[k-1][j][i].grid_M1 + pG->Coup[k-1][j][i-1].grid_M1 );
        M3  = 0.25 * ( pG->Coup[k  ][j][i].grid_M3 + pG->Coup[k  ][j][i-1].grid_M3 
                     + pG->Coup[k-1][j][i].grid_M3 + pG->Coup[k-1][j][i-1].grid_M3 );
      
        j1 = 0.25 * ( J1i[k][j][i] + J1i[k][j+1][i] + J1i[k][j][i-1] + J1i[k][j+1][i-1] );
        j3 = 0.25 * ( J3i[k][j][i] + J3i[k][j+1][i] + J3i[k-1][j][i] + J3i[k-1][j+1][i] );

        vf1 = (M1 - j1)/den;
        vf3 = (M3 - j3)/den;
        
        B3u = 0.5 * ( B3i[k][j][i] + B3i[k][j][i-1] );
        B1u = 0.5 * ( B1i[k][j][i] + B1i[k-1][j][i] );
        
        e2 = B3u * vf1 - B1u * vf3;

        switch (lab) {
          case 0:
            emf2[k][j][i] = e2;
            pG->E2i[k][j][i] = e2;
            break;
            
          case 1:
            emf2[k][j][i] = 0.5 * ( pG->E2i[k][j][i] + e2 );
            break;
            
          case 2:
            emf2[k][j][i] = 0.5 * ( pG->E2i[k][j][i] + e2 );

            break;
            
          default:
            peg_error("[integrate_3d_hybrid]: stage undefined!\n");
        }
        
      }
    }
  }
  
  /* requires:
   *  J1i[kb  :kt+1][jb  :jt+1][ib-1:it+1]
   *  J2i[kb  :kt+1][jb-1:jt+1][ib  :it+1]
   *  J3i[kb  :kt  ][jb  :jt+1][ib  :it+1]
   */
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt+1; j++) {
      for (i=ib; i<=it+1; i++) {
        
        den = 0.25 * ( pG->Coup[k][j  ][i].grid_d  + pG->Coup[k][j  ][i-1].grid_d
                     + pG->Coup[k][j-1][i].grid_d  + pG->Coup[k][j-1][i-1].grid_d  );
        M1  = 0.25 * ( pG->Coup[k][j  ][i].grid_M1 + pG->Coup[k][j  ][i-1].grid_M1
                     + pG->Coup[k][j-1][i].grid_M1 + pG->Coup[k][j-1][i-1].grid_M1 );
        M2  = 0.25 * ( pG->Coup[k][j  ][i].grid_M2 + pG->Coup[k][j  ][i-1].grid_M2
                     + pG->Coup[k][j-1][i].grid_M2 + pG->Coup[k][j-1][i-1].grid_M2 );
      
        j1 = 0.25 * ( J1i[k][j][i] + J1i[k][j][i-1] + J1i[k+1][j][i] + J1i[k+1][j][i-1] );
        j2 = 0.25 * ( J2i[k][j][i] + J2i[k][j-1][i] + J2i[k+1][j][i] + J2i[k+1][j-1][i] );

        vf1 = (M1 - j1)/den;
        vf2 = (M2 - j2)/den;
        
        B1u = 0.5 * ( B1i[k][j][i] + B1i[k][j-1][i] );
        B2u = 0.5 * ( B2i[k][j][i] + B2i[k][j][i-1] );
        
        e3 = B1u * vf2 - B2u * vf1;
        
        switch (lab) {
          case 0:
            emf3[k][j][i] = e3;
            pG->E3i[k][j][i] = e3;
            break;
            
          case 1:
            emf3[k][j][i] = 0.5 * ( pG->E3i[k][j][i] + e3 );
            break;
            
          case 2:
            emf3[k][j][i] = 0.5 * ( pG->E3i[k][j][i] + e3 );

            break;
            
          default:
            peg_error("[integrate_3d_hybrid]: stage undefined!\n");
        }
        
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void putmhdvalues(GridS *pG)
 *  \brief Put magnetic and electric fields in particle coupling structure
 *
 * il+3 through iu-3 needed when injector not used
 * il+2 through iu-2 needed for predictor step when injector is used
 * nothing is harmed by generally doing il+2 through iu-2
 */

static void putmhdvalues(GridS *pG)
{
  int i,j,k;
  GPCouple *pq;
  
  for (k=pG->ks-2; k<=pG->ke+2; k++) {
    for (j=pG->js-2; j<=pG->je+2; j++) {
      for (i=pG->is-2; i<=pG->ie+2; i++) {
        pq = &(pG->Coup[k][j][i]);
        pq->grid_b1 = B1c[k][j][i];
        pq->grid_b2 = B2c[k][j][i];
        pq->grid_b3 = B3c[k][j][i];
        pq->grid_e1 = emf1_cc[k][j][i];
        pq->grid_e2 = emf2_cc[k][j][i];
        pq->grid_e3 = emf3_cc[k][j][i];
        pq->grid_eb = B1c[k][j][i] * emf1_cc[k][j][i]
                    + B2c[k][j][i] * emf2_cc[k][j][i]
                    + B3c[k][j][i] * emf3_cc[k][j][i];
      }
    }
  } 
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void upwind(const Real Ul, const Real Uc, const Real Ur, const Real xi)
 *  \brief Return upwinded value
 */

static Real upwind(const Real Ul, const Real Uc, const Real Ur, const Real xi)
{
  Real dUl,dUr,dU;
  
  dUl = Uc-Ul;
  dUr = Ur-Uc;
  if (dUl*dUr > 0.0) {
    dU = (dUl*dUr)*SIGN(dUl+dUr)/MAX(fabs(dUl+dUr),TINY_NUMBER);
  } else {
    dU = 0.0;
  }
  return (Uc + SIGN(xi)*(1.0-fabs(xi))*dU);
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_3d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
 */
void integrate_init_3d(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;
  
  /* Cycle over all Grids on this processor to find maximum Nx1, Nx2 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }
  
  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX(MAX(size1,size2),size3);
  
  if ((emf1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((emf1_cc = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((B1i = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B2i = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B3i = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((B1c = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B2c = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B3c = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((J1i = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((J2i = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((J3i = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  return;
  
on_error:
  integrate_destruct();
  peg_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_2d(void)
 *  \brief Free temporary integration arrays 
 */
void integrate_destruct_3d(void)
{
  
  if (emf1    != NULL) free_3d_array(emf1);
  if (emf2    != NULL) free_3d_array(emf2);
  if (emf3    != NULL) free_3d_array(emf3);
  if (emf1_cc != NULL) free_3d_array(emf1_cc);
  if (emf2_cc != NULL) free_3d_array(emf2_cc);
  if (emf3_cc != NULL) free_3d_array(emf3_cc);
  if (B1i     != NULL) free_3d_array(B1i);
  if (B2i     != NULL) free_3d_array(B2i);
  if (B3i     != NULL) free_3d_array(B3i);
  if (B1c     != NULL) free_3d_array(B1c);
  if (B2c     != NULL) free_3d_array(B2c);
  if (B3c     != NULL) free_3d_array(B3c);
  if (J1i     != NULL) free_3d_array(J1i);
  if (J2i     != NULL) free_3d_array(J2i);
  if (J3i     != NULL) free_3d_array(J3i);
  
  return;
}
