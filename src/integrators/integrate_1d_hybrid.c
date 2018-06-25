#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_1d_hybrid.c
 *  \brief Integrate coupled induction equation and particles
 *
 * PURPOSE: 
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_1d_hybrid()
 * - integrate_init_1d()
 * - integrate_destruct_1d()
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
static Real *emf1=NULL, *emf2=NULL, *emf3=NULL;
static Real *emf1_cc=NULL, *emf2_cc=NULL, *emf3_cc=NULL;
static Real *B1i=NULL, *B2i=NULL, *B3i=NULL;
static Real *B1c=NULL, *B2c=NULL, *B3c=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   calculate_emfs() - calculates electric field
 *   putmhdvalues() - puts E & B info in particle coupling structure
 *============================================================================*/

static void putmhdvalues(GridS *pG);
static void calculate_emfs(const GridS *pG, short lab);
static void clear_emfs(void);


static Real upwind(const Real Ul, const Real Uc, const Real Ur, const Real xi);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_1d_rlf(DomainS *pD)
 *  \brief 1D RLF integrator */

void integrate_1d_hybrid(DomainS *pD)
{
  GridS *pG = pD->Grid;
  Real q1 = pG->dt/pG->dx1;
  int i,is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  
/*=== Step 0: Preparations ===================================================*/
  
/*--- Step 0a ------------------------------------------------------------------
 * Copy face- and cell-centered magnetic fields at t^{n} into local arrays
 */
  
  for (i=ilp; i<=iup; i++) {
    B1i[i] = pG->B1i[ks][js][i];
    B2i[i] = pG->B2i[ks][js][i];
    B3i[i] = pG->B3i[ks][js][i];
    B1c[i] = pG->U[ks][js][i].B1c;
    B2c[i] = pG->U[ks][js][i].B2c;
    B3c[i] = pG->U[ks][js][i].B3c;
  }

/*=== Step 1: Compute values at t^{n+1/2} ====================================*/
  
/*--- Step 1a ------------------------------------------------------------------
 * Deposit particles on grid, correct bndry values, and set bndry conditions.
 */
  
  deposit(pD,1);
  
/*--- Step 1b ------------------------------------------------------------------
 * Calculate cell- and edge-centered electric fields at t^{n}. Store them
 * as grid quantities for use later.
 */
  
  calculate_emfs(pG,0);

/*--- Step 1c ------------------------------------------------------------------
 * Calculate face-centered magnetic fields at t^{n+1} using CT.
 */

  for (i=is-4; i<=ie+4; i++) {
    B2i[i] += q1 * ( emf3[i+1] - emf3[i] );
    B3i[i] -= q1 * ( emf2[i+1] - emf2[i] );
  }
  
 
/*--- Step 1d ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 */
  
  for (i=is-4; i<=ie+4; i++) {
    B1c[i] = 0.5 * ( B1i[i] + B1i[i+1] );
    B2c[i] =         B2i[i];
    B3c[i] =         B3i[i];
  }
  
/*--- Step 1e ------------------------------------------------------------------
 * Re-calculate cell- and edge-centered electric fields using B^{n+1}. Then
 * average to get electric fields at t^{n+1/2}.
 */
  
  calculate_emfs(pG,1);
  
/*--- Step 1f ------------------------------------------------------------------
 * Time-average to get magnetic fields at t^{n+1/2}.
 */
  
  for (i=is-4; i<=ie+4; i++) {
    B1c[i] = 0.5 * ( pG->U[ks][js][i].B1c + B1c[i] );
    B2c[i] = 0.5 * ( pG->U[ks][js][i].B2c + B2c[i] );
    B3c[i] = 0.5 * ( pG->U[ks][js][i].B3c + B3c[i] );
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
  
  for (i=is-3; i<=ie+3; i++) {
    B2i[i] = pG->B2i[ks][js][i] + q1 * ( emf3[i+1] - emf3[i] );
    B3i[i] = pG->B3i[ks][js][i] - q1 * ( emf2[i+1] - emf2[i] );
  }
  
/*--- Step 1j ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields
 */
  
  for (i=is-3; i<=ie+3; i++) {
    B1c[i] = 0.5 * ( B1i[i] + B1i[i+1] );
    B2c[i] =         B2i[i];
    B3c[i] =         B3i[i];
  }
  
/*--- Step 1k ------------------------------------------------------------------
 * Calculate cell- and edge-centered electric field at t^{n+1}. Then average
 * to get electric fields at t^{n+1/2}.
 */
  
  calculate_emfs(pG,2);
  
/*--- Step 1l ------------------------------------------------------------------
 * Average to get magnetic fields at t^{n+1/2}.
 */  

  for (i=is-2; i<=ie+2; i++) {
    B1c[i] = 0.5 * ( pG->U[ks][js][i].B1c + B1c[i] );
    B2c[i] = 0.5 * ( pG->U[ks][js][i].B2c + B2c[i] );
    B3c[i] = 0.5 * ( pG->U[ks][js][i].B3c + B3c[i] );
  }
  
/*--- Step 1m ------------------------------------------------------------------
 * Put cell-centered electric and magnetic fields in particle coupling array
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
  
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][js][i] += q1 * ( emf3[i+1] - emf3[i] );
    pG->B3i[ks][js][i] -= q1 * ( emf2[i+1] - emf2[i] );
  }
  
/*--- Step 2c ------------------------------------------------------------------
 * Set cell-centered magnetic fields to average of face-centered fields.
 */
  
  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].B1c = 0.5*(pG->B1i[ks][js][i]+pG->B1i[ks][js][i+1]);
    pG->U[ks][js][i].B2c = pG->B2i[ks][js][i];
    pG->U[ks][js][i].B3c = pG->B3i[ks][js][i];
  }
  
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

static void clear_emfs(void)
{
  int i;
  
  for (i=ilp; i<=iup; i++) {
    emf1[i] = 0.0;
    emf2[i] = 0.0;
    emf3[i] = 0.0;
    emf1_cc[i] = 0.0;
    emf2_cc[i] = 0.0;
    emf3_cc[i] = 0.0;
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
 * is-2 through ie+2 needed when injector not used
 * is-3 through ie+3 needed for predictor step when injector is used
 */

static void calculate_emfs(const GridS *pG, short lab)
{
  int ks=pG->ks;
  int js=pG->js;
  int i,is=pG->is,ie=pG->ie;
  int ib,it,il,iu,adiab;
  Real den,dinv,M1,M2,M3,vf1,vf2,vf3,j2,j3,e1,e2,e3;
  Real B1u,B2u,B3u,pemf1;
  Real dx1i = 1.0/pG->dx1;

  /*
  Real dtodx1=dx1i*pG->dt,Ul,Ur,Uc;
  int uw;
  */
 
  switch (lab) {
    case 0:
      ib=is-4; it=ie+4;
      il=is-2; iu=ie+2;
      break;
      
    case 1:
      ib=is-3; it=ie+3;
      if (vinject == 0) {
        il=is-2; iu=ie+2;
      } else { 
        il=is-3; iu=ie+3;
      }
      break;
      
    case 2:
      ib=is; it=ie;
      il=is-2; iu=ie+2;
      break;
      
    default:
      peg_error("[integrate_1d_hybrid]: stage undefined!\n");
  }
  
  
  for (i=il; i<=iu; i++) {
    
    j2    = 0.5 * dx1i * ( B3c[i-1] - B3c[i+1] );
    j3    = 0.5 * dx1i * ( B2c[i+1] - B2c[i-1] );
    
    dinv  =  1.0/pG->Coup[ks][js][i].grid_d;
    vf1   = (pG->Coup[ks][js][i].grid_M1     )*dinv;
    vf2   = (pG->Coup[ks][js][i].grid_M2 - j2)*dinv;
    vf3   = (pG->Coup[ks][js][i].grid_M3 - j3)*dinv;
    
    pemf1 = 0.5 * dx1i * dinv * ( ZTeTi * 0.5 * beta ) *
          ( pG->Coup[ks][js][i+1].grid_d - pG->Coup[ks][js][i-1].grid_d );
    
#ifdef ADIABATIC
    adiab = Gamma * pow(pG->Coup[ks][js][i].grid_d,Gamma_1);
    pemf1 *= adiab;
#endif
    
    e1 = B2c[i] * vf3 - B3c[i] * vf2 - pemf1;
    e2 = B3c[i] * vf1 - B1c[i] * vf3;
    e3 = B1c[i] * vf2 - B2c[i] * vf1;
    
    switch (lab) {
      case 0:
        pG->U[ks][js][i].E1c = e1;
        pG->U[ks][js][i].E2c = e2;
        pG->U[ks][js][i].E3c = e3;
        break;
        
      case 1:
        emf1_cc[i] = 0.5 * ( pG->U[ks][js][i].E1c + e1 );
        emf2_cc[i] = 0.5 * ( pG->U[ks][js][i].E2c + e2 );
        emf3_cc[i] = 0.5 * ( pG->U[ks][js][i].E3c + e3 );
        break;
        
      case 2:
        emf1_cc[i] = 0.5 * ( pG->U[ks][js][i].E1c + e1 );
        emf2_cc[i] = 0.5 * ( pG->U[ks][js][i].E2c + e2 );
        emf3_cc[i] = 0.5 * ( pG->U[ks][js][i].E3c + e3 );
        break;
        
      default:
        peg_error("[integrate_1d_hybrid]: stage undefined!\n");
    }
    
  }

  for (i=ib; i<=it; i++) {
    
    den = pG->Coup[ks][js][i].grid_d ;
    M2  = pG->Coup[ks][js][i].grid_M2;
    M3  = pG->Coup[ks][js][i].grid_M3;
        
    j2  = dx1i * 0.5 * ( B3i[i-1] - B3i[i+1] );
    j3  = dx1i * 0.5 * ( B2i[i+1] - B2i[i-1] );
        
    vf2 = (M2 - j2)/den;
    vf3 = (M3 - j3)/den;
        
    B3u = B3i[i];
    B2u = B2i[i];
    
    e1 = B2u * vf3 - B3u * vf2;
        
    switch (lab) {
      case 0:
        emf1[i] = e1;
        pG->E1i[ks][js][i] = e1;
        break;
        
      case 1:
        emf1[i] = 0.5 * ( pG->E1i[ks][js][i] + e1 );
        break;
        
      case 2:
        emf1[i] = 0.5 * ( pG->E1i[ks][js][i] + e1 );
        break;
        
      default:
        peg_error("[integrate_1d_hybrid]: stage undefined!\n");
    }

  }

  for (i=ib; i<=it+1; i++) {
    
    den = 0.5 * ( pG->Coup[ks][js][i].grid_d  + pG->Coup[ks][js][i-1].grid_d  );
    M1  = 0.5 * ( pG->Coup[ks][js][i].grid_M1 + pG->Coup[ks][js][i-1].grid_M1 );
    M3  = 0.5 * ( pG->Coup[ks][js][i].grid_M3 + pG->Coup[ks][js][i-1].grid_M3 );
     
    j2 = dx1i * ( B3i[i-1] - B3i[i] );
    j3 = dx1i * ( B2c[i] - B2c[i-1] );
        
    vf1 = M1/den;
    vf3 = (M3 - j3)/den;
    
    /*
    if (vf1 >= 0.0)
      uw = i-1;
    else
      uw = i;
    Ul = B3i[uw-1]; Uc = B3i[uw]; Ur = B3i[uw+1];
    B3u = upwind(Ul,Uc,Ur,vf1*dtodx1);
    */
    B3u = 0.5 * ( B3i[i] + B3i[i-1] );   

    B1u = B1i[i];
      
    e2 = B3u * vf1 - B1u * vf3;
    
    switch (lab) {
      case 0:
        emf2[i] = e2;
        pG->E2i[ks][js][i] = e2;
        break;
        
      case 1:
        emf2[i] = 0.5 * ( pG->E2i[ks][js][i] + e2 );
        break;
        
      case 2:
        emf2[i] = 0.5 * ( pG->E2i[ks][js][i] + e2 );
        break;
        
      default:
        peg_error("[integrate_1d_hybrid]: stage undefined!\n");
    }
    
  }
  
  for (i=ib; i<=it+1; i++) {
    
    den = 0.5 * ( pG->Coup[ks][js][i].grid_d  + pG->Coup[ks][js][i-1].grid_d  );
    M1  = 0.5 * ( pG->Coup[ks][js][i].grid_M1 + pG->Coup[ks][js][i-1].grid_M1 );
    M2  = 0.5 * ( pG->Coup[ks][js][i].grid_M2 + pG->Coup[ks][js][i-1].grid_M2 );
        
    j2  = dx1i * ( B3c[i-1] - B3c[i] );
    j3  = dx1i * ( B2i[i] - B2i[i-1] );
    
    vf1 = M1/den;
    vf2 = (M2 - j2)/den;

    /*
    if (vf1 >= 0.0)
      uw = i-1;
    else
      uw = i;
    Ul = B2i[uw-1]; Uc = B2i[uw]; Ur = B2i[uw+1];
    B2u = upwind(Ul,Uc,Ur,vf1*dtodx1);
    */

    B2u = 0.5 * ( B2i[i] + B2i[i-1] );       


    B1u = B1i[i];
    
    e3 = B1u * vf2 - B2u * vf1;
    
    switch (lab) {
      case 0:
        emf3[i] = e3;
        pG->E3i[ks][js][i] = e3;
        break;
        
      case 1:
        emf3[i] = 0.5 * ( pG->E3i[ks][js][i] + e3 );
        break;
        
      case 2:
        emf3[i] = 0.5 * ( pG->E3i[ks][js][i] + e3 );
        break;
        
      default:
        peg_error("[integrate_1d_hybrid]: stage undefined!\n");
    }
    
  }
  
  return;
}


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
/*! \fn void putmhdvalues(GridS *pG)
 *  \brief Put magnetic and electric fields in particle coupling structure
 *
 * il+3 through iu-3 needed when injector not used
 * il+2 through iu-2 needed for predictor step when injector is used
 * nothing is harmed by generally doing il+2 through iu-2
 */

static void putmhdvalues(GridS *pG)
{
  int i;
  GPCouple *pq;
  
  for (i=pG->is-2; i<=pG->ie+2; i++) {
    pq = &(pG->Coup[pG->ks][pG->js][i]);
    pq->grid_b1 = B1c[i];
    pq->grid_b2 = B2c[i];
    pq->grid_b3 = B3c[i];
    pq->grid_e1 = emf1_cc[i];
    pq->grid_e2 = emf2_cc[i];
    pq->grid_e3 = emf3_cc[i];
    pq->grid_eb = B1c[i] * emf1_cc[i]
                + B2c[i] * emf2_cc[i]
                + B3c[i] * emf3_cc[i];
  }  
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_1d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
 */
void integrate_init_1d(MeshS *pM)
{
  int size1=0,nl,nd;
  
  /* Cycle over all Grids on this processor to find maximum Nx1 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if ((pM->Domain[nl][nd].Grid->Nx[0]) > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
      }
    }
  }
  
  size1 = size1 + 2*nghost;
  
  if ((emf1 = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((emf1_cc = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((B1i = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B2i = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B3i = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  
  if ((B1c = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B2c = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((B3c = (Real*)calloc_1d_array(size1,sizeof(Real)))==NULL)
    goto on_error;
  
  return;
  
on_error:
  integrate_destruct();
  peg_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_1d(void)
 *  \brief Free temporary integration arrays 
 */
void integrate_destruct_1d(void)
{
  
  if (emf1    != NULL) free_1d_array(emf1);
  if (emf2    != NULL) free_1d_array(emf2);
  if (emf3    != NULL) free_1d_array(emf3);
  if (emf1_cc != NULL) free_1d_array(emf1_cc);
  if (emf2_cc != NULL) free_1d_array(emf2_cc);
  if (emf3_cc != NULL) free_1d_array(emf3_cc);
  if (B1i     != NULL) free_1d_array(B1i);
  if (B2i     != NULL) free_1d_array(B2i);
  if (B3i     != NULL) free_1d_array(B3i);
  if (B1c     != NULL) free_1d_array(B1c);
  if (B2c     != NULL) free_1d_array(B2c);
  if (B3c     != NULL) free_1d_array(B3c);

  return;
}
