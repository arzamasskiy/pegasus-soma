 #include "../copyright.h"
/*============================================================================*/
/*! \file integrate_2d_hybrid.c
 *  \brief Integrate coupled induction equation and particles
 *
 * PURPOSE: 
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_2d_hybrid()
 * - integrate_init_2d()
 * - integrate_destruct_2d()
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../pegasus.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "../particles/particle.h"

/* The interface electric fields; cell-centered electric fields */
static Real **emf1=NULL, **emf2=NULL, **emf3=NULL;
static Real **emf1_cc=NULL, **emf2_cc=NULL, **emf3_cc=NULL;
static Real **B1i=NULL, **B2i=NULL, **B3i=NULL;
static Real **B1c=NULL, **B2c=NULL, **B3c=NULL;
static Real **J1i=NULL, **J2i=NULL, **J3i=NULL;

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
/*! \fn void integrate_2d_hybrid(DomainS *pD)
 *  \brief 2d hybrid-PIC integrator */

void integrate_2d_hybrid(DomainS *pD)
{
  GridS *pG = pD->Grid;
  Real q1 = pG->dt/pG->dx1;
  Real q2 = pG->dt/pG->dx2;
  int i,is = pG->is, ie = pG->ie;
  int j,js = pG->js, je = pG->je;
  int ks = pG->ks;
  int my_iproc, my_jproc, my_kproc;
#ifdef SHEARING_BOX
q1*=-1;q2*=-1;
#endif 
/*=== Step 0: Preparations ===================================================*/
  
/*--- Step 0a ------------------------------------------------------------------
 * Copy face- and cell-centered magnetic fields at t^{n} into local arrays
 */

    for (j=jlp; j<=jup; j++) {
      for (i=ilp; i<=iup; i++) {
        B1i[j][i] = pG->B1i[ks][j][i];
        B2i[j][i] = pG->B2i[ks][j][i];
        B3i[j][i] = pG->B3i[ks][j][i];
        B1c[j][i] = pG->U[ks][j][i].B1c;
        B2c[j][i] = pG->U[ks][j][i].B2c;
        B3c[j][i] = pG->U[ks][j][i].B3c;
      }
    }
  
/*=== Step 1: Compute values at t^{n+1/2} ====================================*/
  
/*--- Step 1a ------------------------------------------------------------------
 * Deposit particles on grid, correct bndry values, and set bndry conditions
 */

  deposit(pD,1);
  
/*--- Step 1b ------------------------------------------------------------------
 * Calculate cell- and edge-centered electric fields at t^{n}. Store them
 * as grid quantities for use later.  Note that remapping and symmetrization 
 * must be done on J at shearing periodic boundaries.
 */
  
  calculate_currents(pG,0);

  calculate_emfs(pG,0);

/*--- Step 1c ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 */
  
  for (j=js-4; j<=je+4; j++) {
    for (i=is-4; i<=ie+4; i++) {
      B1i[j][i]  -= q2 * ( emf3[j+1][i  ] - emf3[j][i] );
      B2i[j][i]  += q1 * ( emf3[j  ][i+1] - emf3[j][i] );
      B3i[j][i]  += q2 * ( emf1[j+1][i  ] - emf1[j][i] )
                  - q1 * ( emf2[j  ][i+1] - emf2[j][i] );
    }
  }
  
/*--- Step 1d ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 */

  for (j=js-4; j<=je+4; j++) {
    for (i=is-4; i<=ie+4; i++) {
      B1c[j][i] = 0.5 * ( B1i[j][i] + B1i[j][i+1] );
      B2c[j][i] = 0.5 * ( B2i[j][i] + B2i[j+1][i] );
      B3c[j][i] =         B3i[j][i];
    }
  }

/*--- Step 1e ------------------------------------------------------------------
 * Re-calculate cell- and edge-centered electric fields using B^{n+1}. Then
 * average to get electric fields at t^{n+1/2}. Note that remapping and 
 * symmetrization must be done on J at shearing periodic boundaries.
 */
  
  calculate_currents(pG,1);

  calculate_emfs(pG,1);

/*--- Step 1f ------------------------------------------------------------------
 * Time-average to get magnetic fields at t^{n+1/2}.
 */
  
  for (j=js-4; j<=je+4; j++) {
    for (i=is-4; i<=ie+4; i++) {
      B1c[j][i] = 0.5 * ( pG->U[ks][j][i].B1c + B1c[j][i] );
      B2c[j][i] = 0.5 * ( pG->U[ks][j][i].B2c + B2c[j][i] );
      B3c[j][i] = 0.5 * ( pG->U[ks][j][i].B3c + B3c[j][i] );
    }
  }
  
/*--- Step 1g ------------------------------------------------------------------
 * Put cell-centered electric and magnetic fields in particle coupling array.
 */
  
  putmhdvalues(pG);
  
/*--- Step 1h ------------------------------------------------------------------
 * Integrate temporary copy of particles to t^{n+1} and deposit on grid.
 */
  
  Integrate_Particles_1(pD);

/*--- Step 1i ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 */
  
  for (j=js-3; j<=je+3; j++) {
    for (i=is-3; i<=ie+3; i++) {
      B1i[j][i]  = pG->B1i[ks][j][i] -
                   q2 * ( emf3[j+1][i  ] - emf3[j][i] );
      B2i[j][i]  = pG->B2i[ks][j][i] +
                   q1 * ( emf3[j  ][i+1] - emf3[j][i] );
      B3i[j][i]  = pG->B3i[ks][j][i] +
                   q2 * ( emf1[j+1][i  ] - emf1[j][i] ) -
                   q1 * ( emf2[j  ][i+1] - emf2[j][i] );
    }
  }
  
/*--- Step 1j ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 */
  
  for (j=js-3; j<=je+3; j++) {
    for (i=is-3; i<=ie+3; i++) {
      B1c[j][i] = 0.5 * ( B1i[j][i] + B1i[j][i+1] );
      B2c[j][i] = 0.5 * ( B2i[j][i] + B2i[j+1][i] );
      B3c[j][i] =         B3i[j][i];
    }
  }
  
/*--- Step 1k ------------------------------------------------------------------
 * Calculate cell- and edge-centered electric field at t^{n+1}. Then average
 * to get electric fields at t^{n+1/2}. Note that remapping and symmetrization
 * must be done on J at shearing periodic boundaries.
 */
  
  calculate_currents(pG,2);

  calculate_emfs(pG,2);

/*--- Step 1l ------------------------------------------------------------------
 * Average to get magnetic fields at t^{n+1/2}.
 */  
  
  for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+2; i++) {
      B1c[j][i] = 0.5 * ( pG->U[ks][j][i].B1c + B1c[j][i] );
      B2c[j][i] = 0.5 * ( pG->U[ks][j][i].B2c + B2c[j][i] );
      B3c[j][i] = 0.5 * ( pG->U[ks][j][i].B3c + B3c[j][i] );
    }
  }
  
/*--- Step 1m ------------------------------------------------------------------
 * Put cell-centered electric and magnetic fields in particle coupling array.
 */
  
  putmhdvalues(pG);

/*=== Step 2: Compute values at t^{n+1} ======================================*/
  
/*--- Step 2a ------------------------------------------------------------------
 * Integrate particles to t^{n+1}.
 */
  
  Integrate_Particles_2(pD);
   
/*--- Step 2b ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 */
  
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B1i[ks][j][i]  -= q2 * ( emf3[j+1][i  ] - emf3[j][i] );
      pG->B2i[ks][j][i]  += q1 * ( emf3[j  ][i+1] - emf3[j][i] );
      pG->B3i[ks][j][i]  += q2 * ( emf1[j+1][i  ] - emf1[j][i] ) -
                            q1 * ( emf2[j  ][i+1] - emf2[j][i] );
    }
    pG->B1i[ks][j][ie+1] -= q2 * ( emf3[j+1][ie+1] - emf3[j][ie+1] );
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += q1 * ( emf3[je+1][i+1] - emf3[je+1][i] );
  }

/*--- Step 2c ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 */
  
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i]+pG->B1i[ks][j][i+1]);
        pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i]+pG->B2i[ks][j+1][i]);
        pG->U[ks][j][i].B3c = pG->B3i[ks][j][i];
      }
    }
  
  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
static void clear_emfs(void)
{
  int i,j;

  for (j=jlp; j<=jup; j++) {
    for (i=ilp; i<=iup; i++) {
      emf1[j][i] = 0.0;
      emf2[j][i] = 0.0;
      emf3[j][i] = 0.0;
      emf1_cc[j][i] = 0.0;
      emf2_cc[j][i] = 0.0;
      emf3_cc[j][i] = 0.0;
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
  int i,j,jb,jt,ib,it, ks = pG->ks;
  Real dx1i = 1.0/pG->dx1, dx2i = 1.0/pG->dx2;
#ifdef SHEARING_BOX
dx1i*=-1;dx2i*=-1;
#endif 
  switch (lab) {
    case 0:
      jb=pG->js-4; jt=pG->je+4;
      ib=pG->is-4; it=pG->ie+4;
      break;
      
    case 1:
      jb=pG->js-3; jt=pG->je+3;
      ib=pG->is-3; it=pG->ie+3;
      break;
      
    case 2:
      jb=pG->js-2; jt=pG->je+2;
      ib=pG->is-2; it=pG->ie+2;
      break;
      
    default:
      peg_error("[integrate_2d_hybrid]: stage undefined!\n");
  }
  
  /* calculate_emfs requires:
   *   J1i[jb  :jt+1][ib-1:it+1]
   *   J2i[jb-1:jt+1][ib  :it+1]
   *   J3i[jb  :jt+1][ib  :it+1]
   * these require:
   *   B1i[jb-1:jt+1][ib  :it+1]
   *   B2i[jb  :jt+1][ib-1:it+1]
   *   B3i[jb-1:jt+1][ib-1:it+1]
   */
  
  for (j=jb; j<=jt+1; j++) {
    for (i=ib; i<=it+1; i++) {
      J1i[j][i] =  dx2i * ( B3i[j][i] - B3i[j-1][i] );
      J2i[j][i] = -dx1i * ( B3i[j][i] - B3i[j][i-1] );
      J3i[j][i] =  dx1i * ( B2i[j][i] - B2i[j][i-1] )
                  -dx2i * ( B1i[j][i] - B1i[j-1][i] );
    }
    J1i[j][ib-1] =  dx2i * ( B3i[j][ib-1] - B3i[j-1][ib-1] );
  }
  for (i=ib; i<=it+1; i++) {
    J2i[jb-1][i] = -dx1i * ( B3i[jb-1][i] - B3i[jb-1][i-1] );
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void calculate_emfs(GridS *pG, char *stage)
 *  \brief Calculate cell- and edge-centered electric fields
 *
 * NOTE: ions are not pushed by E = eta J
 *
 * For cell-centered:
 * il+3 through iu-3 needed when injector not used
 * il+2 through iu-2 needed for predictor step when injector is used
 * nothing is harmed by generally doing il+2 through iu-2
 *
 */

static void calculate_emfs(const GridS *pG, short lab)
{
  int ks=pG->ks;
  int j,js=pG->js,je=pG->je;
  int i,is=pG->is,ie=pG->ie;
  int jb,jt,ib,it,il,iu,jl,ju;
  Real den,dinv,M1,M2,M3,vf1,vf2,vf3,j1,j2,j3,e1,e2,e3;
  Real B1u,B2u,B3u,pemf1,pemf2,adiab;
  Real dx1i = 1.0/pG->dx1, dx2i = 1.0/pG->dx2;
  Real x1,x2,x3;  
  Real dx1isq = SQR(dx1i), dx2isq = SQR(dx2i);
  switch (lab) {
    case 0:
      jb=js-4; jt=je+4;
      ib=is-4; it=ie+4;
      jl=js-2; ju=je+2;
      il=is-2; iu=ie+2;
      break;
      
    case 1:
      jb=js-3; jt=je+3;
      ib=is-3; it=ie+3;
      if (vinject == 0) {
        il=is-2; iu=ie+2;
        jl=js-2; ju=je+2;
      } else {
        il=is-3; iu=ie+3;
        jl=js-3; ju=je+3;
      }
      break;
      
    case 2:
      jb=js; jt=je;
      ib=is; it=ie;
      jl=js-2; ju=je+2;
      il=is-2; iu=ie+2;
      break;
      
    default:
      peg_error("[integrate_2d_hybrid]: stage undefined!\n");
  }
    
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      
      j1    = 0.5  * ( J1i[j][i] + J1i[j+1][i] );
      j2    = 0.5  * ( J2i[j][i] + J2i[j][i+1] );
      j3    = 0.25 * ( J3i[j][i] + J3i[j][i+1] + J3i[j+1][i] + J3i[j+1][i+1] );
      
      dinv  = 1.0/pG->Coup[ks][j][i].grid_d;
      vf1   = (pG->Coup[ks][j][i].grid_M1 - j1)*dinv;
      vf2   = (pG->Coup[ks][j][i].grid_M2 - j2)*dinv;
      vf3   = (pG->Coup[ks][j][i].grid_M3 - j3)*dinv;
      
      pemf1 = 0.5 * dinv * dx1i * ( ZTeTi * 0.5 * beta ) *
            ( pG->Coup[ks][j][i+1].grid_d - pG->Coup[ks][j][i-1].grid_d );
      pemf2 = 0.5 * dinv * dx2i * ( ZTeTi * 0.5 * beta ) *
            ( pG->Coup[ks][j+1][i].grid_d - pG->Coup[ks][j-1][i].grid_d );

#ifdef ADIABATIC
      adiab = Gamma * pow(pG->Coup[ks][j][i].grid_d,Gamma_1);
      pemf1 *= adiab;
      pemf2 *= adiab;
#endif

#ifdef SHEARING_BOX
      e1 = B3c[j][i] * vf2 - B2c[j][i] * vf3 - pemf1;
      e2 = B1c[j][i] * vf3 - B3c[j][i] * vf1 - pemf2;
      e3 = B2c[j][i] * vf1 - B1c[j][i] * vf2;
#else
      e1 = B2c[j][i] * vf3 - B3c[j][i] * vf2 - pemf1;
      e2 = B3c[j][i] * vf1 - B1c[j][i] * vf3 - pemf2;
      e3 = B1c[j][i] * vf2 - B2c[j][i] * vf1;
#endif 
      switch (lab) {
        case 0:
          pG->U[ks][j][i].E1c = e1;
          pG->U[ks][j][i].E2c = e2;
          pG->U[ks][j][i].E3c = e3;
          break;
          
        case 1:
          emf1_cc[j][i] = 0.5 * ( pG->U[ks][j][i].E1c + e1 );
          emf2_cc[j][i] = 0.5 * ( pG->U[ks][j][i].E2c + e2 );
          emf3_cc[j][i] = 0.5 * ( pG->U[ks][j][i].E3c + e3 );
          break;
          
        case 2:
          emf1_cc[j][i] = 0.5 * ( pG->U[ks][j][i].E1c + e1 );
          emf2_cc[j][i] = 0.5 * ( pG->U[ks][j][i].E2c + e2 );
          emf3_cc[j][i] = 0.5 * ( pG->U[ks][j][i].E3c + e3 );
          break;
          
        default:
          peg_error("[integrate_2d_hybrid]: stage undefined!\n");
      }
      
    }
  }
  
  for (j=jb; j<=jt+1; j++) {
    for (i=ib; i<=it; i++) {
      
      den = 0.5 * ( pG->Coup[ks][j][i].grid_d  + pG->Coup[ks][j-1][i].grid_d  );
      M2  = 0.5 * ( pG->Coup[ks][j][i].grid_M2 + pG->Coup[ks][j-1][i].grid_M2 );
      M3  = 0.5 * ( pG->Coup[ks][j][i].grid_M3 + pG->Coup[ks][j-1][i].grid_M3 );
        
      j2 = 0.25 * ( J2i[j][i] + J2i[j][i+1] + J2i[j-1][i] + J2i[j-1][i+1] );
      j3 = 0.5  * ( J3i[j][i] + J3i[j][i+1] );
        
      vf2 = (M2 - j2)/den;
      vf3 = (M3 - j3)/den;

      B3u = 0.5 * ( B3i[j][i] + B3i[j-1][i] );
      B2u = B2i[j][i];
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
#ifdef SHEARING_BOX
      e1 = B3u * vf2 - B2u * (vf3 - Shear_0 * x1);
      if (eta_Ohm > 0) e1 += eta_Ohm * J1i[j][i];
#else
      e1 = B2u * vf3 - B3u * vf2;
#endif     
 
      switch (lab) {
        case 0:
          emf1[j][i] = e1;
          pG->E1i[ks][j][i] = e1;
          break;
          
        case 1:
          emf1[j][i] = 0.5 * ( pG->E1i[ks][j][i] + e1 );
          break;
          
        case 2:
          emf1[j][i] = 0.5 * ( pG->E1i[ks][j][i] + e1 );
#ifdef SHEARING_BOX
          if (eta_hyper > 0) emf1[j][i] -= eta_hyper * ((J1i[j][i-1] - 2.0*J1i[j][i] + J1i[j][i+1]) * dx1isq + (J1i[j-1][i] - 2.0*J1i[j][i] + J1i[j+1][i]) * dx2isq);
#endif
	  break;
          
        default:
          peg_error("[integrate_2d_hybrid]: stage undefined!\n");
      }
      
    }
  }
  
  for (j=jb; j<=jt; j++) {
    for (i=ib; i<=it+1; i++) {
      
      den = 0.5 * ( pG->Coup[ks][j][i].grid_d  + pG->Coup[ks][j][i-1].grid_d  );
      M1  = 0.5 * ( pG->Coup[ks][j][i].grid_M1 + pG->Coup[ks][j][i-1].grid_M1 );
      M3  = 0.5 * ( pG->Coup[ks][j][i].grid_M3 + pG->Coup[ks][j][i-1].grid_M3 );
        
      j1 = 0.25 * ( J1i[j][i] + J1i[j+1][i] + J1i[j][i-1] + J1i[j+1][i-1] );
      j3 = 0.5  * ( J3i[j][i] + J3i[j+1][i] );
        
      vf1 = (M1 - j1)/den;
      vf3 = (M3 - j3)/den;

      B3u = 0.5 * ( B3i[j][i] + B3i[j][i-1] );
      B1u = B1i[j][i];
      fc_pos(pG,i,js,ks,&x1,&x2,&x3);
#ifdef SHEARING_BOX
      e2 = B1u * (vf3 - Shear_0*x1) - B3u * vf1;
      if (eta_Ohm > 0) e2 += eta_Ohm * J2i[j][i];
#else      
      e2 = B3u * vf1 - B1u * vf3;
#endif  
    
     switch (lab) {
        case 0:
          emf2[j][i] = e2;
          pG->E2i[ks][j][i] = e2;
          break;
          
        case 1:
          emf2[j][i] = 0.5 * ( pG->E2i[ks][j][i] + e2 );
          break;
          
        case 2:
          emf2[j][i] = 0.5 * ( pG->E2i[ks][j][i] + e2 );
#ifdef SHEARING_BOX
          if (eta_hyper > 0) emf2[j][i] -= eta_hyper * ((J2i[j][i-1] - 2.0*J2i[j][i] + J2i[j][i+1]) * dx1isq + (J2i[j-1][i] - 2.0*J2i[j][i] + J2i[j+1][i]) * dx2isq);
#endif
	  break;
          
        default:
          peg_error("[integrate_2d_hybrid]: stage undefined!\n");
      }
      
    }
  }

  for (j=jb; j<=jt+1; j++) {
    for (i=ib; i<=it+1; i++) {
      
      den = 0.25 * ( pG->Coup[ks][j  ][i].grid_d  + pG->Coup[ks][j  ][i-1].grid_d
                   + pG->Coup[ks][j-1][i].grid_d  + pG->Coup[ks][j-1][i-1].grid_d  );
      M1  = 0.25 * ( pG->Coup[ks][j  ][i].grid_M1 + pG->Coup[ks][j  ][i-1].grid_M1
                   + pG->Coup[ks][j-1][i].grid_M1 + pG->Coup[ks][j-1][i-1].grid_M1 );
      M2  = 0.25 * ( pG->Coup[ks][j  ][i].grid_M2 + pG->Coup[ks][j  ][i-1].grid_M2
                   + pG->Coup[ks][j-1][i].grid_M2 + pG->Coup[ks][j-1][i-1].grid_M2 );
        
      j1 = 0.5 * ( J1i[j][i] + J1i[j][i-1] );
      j2 = 0.5 * ( J2i[j][i] + J2i[j-1][i] );
      
      vf1 = (M1 - j1)/den;
      vf2 = (M2 - j2)/den;

      B1u = 0.5 * ( B1i[j][i] + B1i[j-1][i] );
      B2u = 0.5 * ( B2i[j][i] + B2i[j][i-1] );
#ifdef SHEARING_BOX
      e3 = B2u * vf1 - B1u * vf2;
      if (eta_Ohm > 0) e3 += eta_Ohm * J3i[j][i];
#else     
      e3 = B1u * vf2 - B2u * vf1;
#endif
      
      switch (lab) {
        case 0:
          emf3[j][i] = e3;
          pG->E3i[ks][j][i] = e3;
          break;
          
        case 1:
          emf3[j][i] = 0.5 * ( pG->E3i[ks][j][i] + e3 );
          break;
          
        case 2:
          emf3[j][i] = 0.5 * ( pG->E3i[ks][j][i] + e3 );
#ifdef SHEARING_BOX
          if (eta_hyper > 0) emf3[j][i] -= eta_hyper * ((J3i[j][i-1] - 2.0*J3i[j][i] + J3i[j][i+1]) * dx1isq + (J3i[j-1][i] - 2.0*J3i[j][i] + J3i[j+1][i]) * dx2isq);
#endif
          break;
          
        default:
          peg_error("[integrate_2d_hybrid]: stage undefined!\n");
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
  int i,j,ks = pG->ks;
  GPCouple *pq;
  
    for (j=pG->js-2; j<=pG->je+2; j++) {
      for (i=pG->is-2; i<=pG->ie+2; i++) {
        pq = &(pG->Coup[ks][j][i]);
        pq->grid_b1 = B1c[j][i];
        pq->grid_b2 = B2c[j][i];
        pq->grid_b3 = B3c[j][i];
        pq->grid_e1 = emf1_cc[j][i];
        pq->grid_e2 = emf2_cc[j][i];
        pq->grid_e3 = emf3_cc[j][i];
        pq->grid_eb = B1c[j][i] * emf1_cc[j][i]
                    + B2c[j][i] * emf2_cc[j][i]
                    + B3c[j][i] * emf3_cc[j][i];
      }
    } 
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_3d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
 */
void integrate_init_2d(MeshS *pM)
{
  int nmax,size1=0,size2=0,nl,nd;
  
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
      }
    }
  }
  
  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  nmax = MAX(size1,size2);
  
  if ((emf1 = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((emf1_cc = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((B1i = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B2i = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B3i = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((B1c = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B2c = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B3c = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((J1i = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((J2i = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((J3i = (Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
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
void integrate_destruct_2d(void)
{
  
  if (emf1    != NULL) free_2d_array(emf1);
  if (emf2    != NULL) free_2d_array(emf2);
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf1_cc != NULL) free_2d_array(emf1_cc);
  if (emf2_cc != NULL) free_2d_array(emf2_cc);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
  if (B1i     != NULL) free_2d_array(B1i);
  if (B2i     != NULL) free_2d_array(B2i);
  if (B3i     != NULL) free_2d_array(B3i);
  if (B1c     != NULL) free_2d_array(B1c);
  if (B2c     != NULL) free_2d_array(B2c);
  if (B3c     != NULL) free_2d_array(B3c);
  if (J1i     != NULL) free_2d_array(J1i);
  if (J2i     != NULL) free_2d_array(J2i);
  if (J3i     != NULL) free_2d_array(J3i);
  
  return;
}
