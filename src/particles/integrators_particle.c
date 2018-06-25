#include "../copyright.h"
/*===========================================================================*/
/*! \file integrators_particle.c
 *  \brief Provide four kinds of particle integrators.
 *
 * PURPOSE: provide four kinds of particle integrators, namely, 2nd order
 *   explicit, 2nd order semi-implicit, 2nd order fully implicit, and
 *   2nd order semi-implicit with modified Boris push .
 * 
 * CONTAINS PUBLIC FUNCTIONS:
 * - Integrate_Particles();
 * - int_par_exp   ()
 * - int_par_semimp()
 * - int_par_fulimp()
 * - int_par_boris ()
 * - feedback_predictor()
 * - feedback_corrector()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - JudgeCrossing()  - judge if the particle cross the grid boundary
 * - Get_Drag()       - calculate the drag force
 * - Get_Force()      - calculate forces other than the drag
 * - Get_MHDForce() - calculate forces for modified Boris push
 *
 * REFERENCE:
 *   X.-N. Bai & J.M. Stone, 2010, ApJS, 190, 297 									      
 *   MWKunz, Aug 21 2012 for modified Boris push
  *  SOMA2018 version modified by L Arzamasskiy                              */
/*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../pegasus.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   JudgeCrossing()  - judge if the particle cross the grid boundary
  *   Get_Force()     - calculate forces other than electric-field force
  *============================================================================*/
static void JudgeCrossing(GridS *pG, Real x1, Real x2, Real x3, GrainS *gr);
static Real3Vect Get_Force(GridS *pG, Real x1, Real x2, Real x3,
                           Real v1, Real v2, Real v3, Real B1, Real B2, Real B3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*---------------------------- Main Integrator -------------------------------*/
/*! \fn void Integrate_Particles(DomainS *pD)
 *  \brief Main hybrid particle integrator.
 *
 * Input: 
 * Output: 
 * Note: 
 */
void Integrate_Particles_1(DomainS *pD)
{
  int i,j,k,is,js,ks,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  int n0 = ncell-1;
  long p;
#ifdef DELTA_F
  Real f0_t;
#endif
  Real dv1, dv2, dv3;       /* amount of velocity update */
  Real weight[3][3][3],drho;
  Real3Vect cell1;
  GPCouple *pq;
  GrainS *curG, *curP, mygr;    /* pointer of the current working position */
  
  GridS *pG = pD->Grid;         /* set ptr to Grid */
  
  curP = &(mygr);               /* temporary particle */
  
  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;

  /* set density and velocity in coupling array */
  set_gpcouple(pG,1);
    
  for(p=0; p<pG->nparticle; p++)
  {/* loop over all particles */
    
    /* calculate change in velocity */
    
    curG = &(pG->particle[p]);
    int_par_boris(pG, curG, cell1, &dv1, &dv2, &dv3);

    /* particle update to curP */

    /* velocity update */
    curP->v1 = curG->v1 + dv1;
    curP->v2 = curG->v2 + dv2;
    curP->v3 = curG->v3 + dv3;

    /* position update */
    if (pG->Nx[0] > 1)
      curP->x1 = curG->x1 + 0.5*pG->dt*(curG->v1 + curP->v1);
    else /* do not move if this dimension collapses */
      curP->x1 = curG->x1;
    
    if (pG->Nx[1] > 1)
      curP->x2 = curG->x2 + 0.5*pG->dt*(curG->v2 + curP->v2);
    else /* do not move if this dimension collapses */
      curP->x2 = curG->x2;
    
    if (pG->Nx[2] > 1)
      curP->x3 = curG->x3 + 0.5*pG->dt*(curG->v3 + curP->v3);
    else /* do not move if this dimension collapses */
      curP->x3 = curG->x3;
    /* if prediction stage, deposit particles onto grid */
    drho = grproperty[curG->property].m;
#ifdef DELTA_F
  getdf(curP, &f0_t);
  drho *= (1.0 = (f0_t/curG->f_0))
#endif
    getweight(pG, curP->x1, curP->x2, curP->x3, cell1, weight, &is, &js, &ks);
      
    k1 = MAX(ks, klp);  k2 = MIN(ks+n0, kup);
    j1 = MAX(js, jlp);  j2 = MIN(js+n0, jup);
    i1 = MAX(is, ilp);  i2 = MIN(is+n0, iup);
      
    for (k=k1; k<=k2; k++) {
      k0 = k-k1;
      for (j=j1; j<=j2; j++) {
        j0 = j-j1;
        for (i=i1; i<=i2; i++) {
          i0 = i-i1;
          pq = &(pG->Coup[k][j][i]);
            
          /* interpolate the temporary particles to the grid */
          pq->grid_d  += weight[k0][j0][i0]*drho;
          pq->grid_M1 += weight[k0][j0][i0]*drho*curP->v1;
          pq->grid_M2 += weight[k0][j0][i0]*drho*curP->v2;
          pq->grid_M3 += weight[k0][j0][i0]*drho*curP->v3;
        }
      }
    }  
  } /* end of the for loop */

  /* if end of prediction stage, exchange with boundaries and smooth */
  exchange_gpcouple(pD,2);
  if (nfpass > 0) smooth_gpcouple(pD,2);

  
  return;
} 

void Integrate_Particles_2(DomainS *pD)
{
  int i,j,k,is,js,ks,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  int n0 = ncell-1;
  long p;
#ifdef DELTA_F
  Real f0_t;
#endif
  Real dv1, dv2, dv3;       /* amount of velocity update */
  Real weight[3][3][3],drho;
  Real3Vect cell1;
  GPCouple *pq;
  GrainS *curG, *curP, mygr;    /* pointer of the current working position */
  
  GridS *pG = pD->Grid;         /* set ptr to Grid */
  
  curP = &(mygr);               /* temporary particle */
  
  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;

  for(p=0; p<pG->nparticle; p++)
  {/* loop over all particles */
    
    /* calculate change in velocity */
    
    curG = &(pG->particle[p]);
    int_par_boris(pG, curG, cell1, &dv1, &dv2, &dv3);

    /* particle update to curP */

    /* velocity update */
    curP->v1 = curG->v1 + dv1;
    curP->v2 = curG->v2 + dv2;
    curP->v3 = curG->v3 + dv3;
    /* position update */
    if (pG->Nx[0] > 1)
      curP->x1 = curG->x1 + 0.5*pG->dt*(curG->v1 + curP->v1);
    else /* do not move if this dimension collapses */
      curP->x1 = curG->x1;
    
    if (pG->Nx[1] > 1)
      curP->x2 = curG->x2 + 0.5*pG->dt*(curG->v2 + curP->v2);
    else /* do not move if this dimension collapses */
      curP->x2 = curG->x2;
    
    if (pG->Nx[2] > 1)
      curP->x3 = curG->x3 + 0.5*pG->dt*(curG->v3 + curP->v3);
    else /* do not move if this dimension collapses */
      curP->x3 = curG->x3;
    
    /* if correction stage, update particle status */
    JudgeCrossing(pG, curP->x1, curP->x2, curP->x3, curG);

    /* if correction stage, update the particle */
    curG->x1 = curP->x1;
    curG->x2 = curP->x2;
    curG->x3 = curP->x3;
    curG->v1 = curP->v1;
    curG->v2 = curP->v2;
    curG->v3 = curP->v3;
        
    //p++; OLD_VERSION
    
  } /* end of the for loop */
  return;
} 

/* ------- 2nd order semi implicit particle integrator w/ Boris push ---------*/
/*! \fn void int_par_boris(Grid *pG, Grain *curG, Real3Vect cell1, 
 *                            Real *dv1, Real *dv2, Real *dv3)
 *  \brief 2nd order semi implicit particle integrator with Boris push
 *         This is essentially the Boris push, modified to take into 
 *         consideration Coriolis and shear terms while preserving 2nd
 *         order accuracy.
 *
 * Input: 
 *   grid pointer (pG), grain pointer (curG), cell size indicator (cell1)
 * Output:
 *   dv1,dv2,dv3: velocity update
 */

void int_par_boris(GridS *pG, GrainS *curG, Real3Vect cell1,
                   Real *dv1, Real *dv2, Real *dv3)
{
    int is, js, ks;
    Real3Vect ft;             /* total force */
    Real qom, b, b2;          /* other shortcut expressions */
    Real x1n, x2n, x3n;       /* first order new position at half a time step */
    Real v1n, v2n, v3n;       /* velocities after half-push from electric field */
    Real M11,M12,M13,M21,M22,M23,M31,M32,M33; /* inverse matrix elements */
#ifdef SHEARING_BOX
  Real b1, oh, sh;
#endif
    Real w1h, w2h, w3h, wh;   /* (qB/mc)*h */
    Real weight[3][3][3];     /* weight function */
    Real bfld1, bfld2, bfld3, efld1, efld2, efld3; 
                              /* fields interpolated to particle position */

    /* step 1: prediction of the particle position after half a time step */
    if (pG->Nx[0] > 1)  x1n = curG->x1+0.5*curG->v1*pG->dt;
    else x1n = curG->x1;
    if (pG->Nx[1] > 1)  x2n = curG->x2+0.5*curG->v2*pG->dt;
    else x2n = curG->x2;
    if (pG->Nx[2] > 1)  x3n = curG->x3+0.5*curG->v3*pG->dt;
    else x3n = curG->x3;

    
    /* Step 2: interpolation to get magnetic and electric fields
     *  at predicted position.
     */
    getweight(pG, x1n, x2n, x3n, cell1, weight, &is, &js, &ks);
    if (getvalues(pG, weight, is, js, ks
                ,&bfld1, &bfld2, &bfld3, &efld1, &efld2, &efld3
                ) == 0)
    {/* particle is on the grid */
      qom = grproperty[curG->property].qomc;
      bfld1 *= qom;  bfld2 *= qom;  bfld3 *= qom;
      efld1 *= qom;  efld2 *= qom;  efld3 *= qom;
    }
    else
    { /* particle out of the grid, free motion, with warning sign */
      peg_perr(0, "Particle moved out of grid %d with position (%f,%f,%f)!\n",
               myID_Comm_world,x1n,x2n,x3n); /* warning! */
    }

    /* Add forcing to electric field */

  
    /* Add user-defined force to electric field */
    /*
     Userforce_particle(&ft, x1, x2, x3);
     efld1 += ft.x1;
     efld2 += ft.x2;
     efld3 += ft.x3;
    */
  
    /* Do a half-push due to electric field */
    v1n = curG->v1 + 0.5*pG->dt*efld1;
    v2n = curG->v2 + 0.5*pG->dt*efld2;
    v3n = curG->v3 + 0.5*pG->dt*efld3;
  
    /* Get Lorentz, Coriolis, and shear forces for modified Boris rotation */
    ft = Get_Force(pG,x1n,x2n,x3n,v1n,v2n,v3n,bfld1,bfld2,bfld3);
  
    /* step 3: calculate velocity update */
  
    /* shortcut expressions */
    w1h = bfld1*pG->dt; w2h = bfld2*pG->dt; w3h = bfld3*pG->dt;
    wh  = sqrt(SQR(w1h)+SQR(w2h)+SQR(w3h));
    b   = 2.0;

    /* velocity evolution */
#ifdef SHEARING_BOX
  oh = Omega_0 * pG->dt;
  sh = Shear_0 * pG->dt;
#endif

#ifdef SHEARING_BOX
    b1 = 1.0/(SQR(b) +2.0*oh*(2.0*oh-sh)
       + SQR(wh)+w2h*(4.0*oh-sh)-w1h*w3h*sh/b);
    b2  = b*b1;
    M11 = 1.0+SQR(w1h/b);
    M12 = w1h*(w2h+2.0*oh)/SQR(b)-w3h/b;
    M13 = w1h*w3h/SQR(b)+(w2h+2*oh)/b;
    M21 = w1h*(w2h-sh+2.0*oh)/SQR(b)+w3h/b;
    M22 = 1.0+(w2h+2.0*oh)*(w2h-sh+2.0*oh)/SQR(b);
    M23 = (w2h+2.0*oh)*w3h/SQR(b)-w1h/b;
    M31 = w1h*w3h/SQR(b)-(w2h-sh+2.0*oh)/b;
    M32 = (w2h-sh+2.0*oh)*w3h/SQR(b)+w1h/b;
    M33 = 1.0+SQR(w3h/b);
#else
    b2  = b/(SQR(b)+SQR(wh));
    M11 = 1.0+SQR(w1h/b);
    M12 = w1h*w2h/SQR(b)+w3h/b;
    M13 = w1h*w3h/SQR(b)-w2h/b;
    M21 = w1h*w2h/SQR(b)-w3h/b;
    M22 = 1.0+SQR(w2h/b);
    M23 = w2h*w3h/SQR(b)+w1h/b;
    M31 = w1h*w3h/SQR(b)+w2h/b;
    M32 = w2h*w3h/SQR(b)-w1h/b;
    M33 = 1.0+SQR(w3h/b);
#endif
    *dv1 = pG->dt*efld1 + pG->dt*2.0*b2*(M11*ft.x1+M12*ft.x2+M13*ft.x3);
    *dv2 = pG->dt*efld2 + pG->dt*2.0*b2*(M21*ft.x1+M22*ft.x2+M23*ft.x3);
    *dv3 = pG->dt*efld3 + pG->dt*2.0*b2*(M31*ft.x1+M32*ft.x2+M33*ft.x3);
  
    return;
}
  
/*=========================== PRIVATE FUNCTIONS ==============================*/

/*--------------------------------------------------------------------------- */
/*! \fn void JudgeCrossing(GridS *pG, Real x1, Real x2, Real x3, GrainS *gr)
 *  \brief Judge if the particle is a crossing particle */
void JudgeCrossing(GridS *pG, Real x1, Real x2, Real x3, GrainS *gr)
{
    /* if it crosses the grid boundary, mark it as a crossing out particle */
    if ((x1>=x1upar) || (x1<x1lpar) || (x2>=x2upar) || (x2<x2lpar) ||
                                       (x3>=x3upar) || (x3<x3lpar) )
      gr->pos = 10;


    return;
}

/*--------------------------------------------------------------------------- */
/*! \fn Real3Vect Get_Force(GridS *pG, Real v1, Real v2, Real v3,
 *                                     Real B1, Real B2, Real B3)
 *  \brief Calculate the forces to the particle for Boris pusher
 *
 * Input:
 *   pG: grid;
 *   x1,x2,x3: particle position;
 *   v1,v2,v3: particle velocity;
 *   B1,B2,B3: magnetic field at particle position;
 * Return:
 *   forces;
 */
Real3Vect Get_Force(GridS *pG, Real x1, Real x2, Real x3,
                               Real v1, Real v2, Real v3,
                               Real B1, Real B2, Real B3)
{
  Real3Vect ft;
  
  ft.x1 = v2*B3 - v3*B2;
  ft.x2 = v3*B1 - v1*B3;
  ft.x3 = v1*B2 - v2*B1;
   
  #ifdef SHEARING_BOX
  ft.x1 = v3*B2 - v2*B3 + 2.0*v3*Omega_0;
  ft.x2 = v1*B3 - v3*B1;
  ft.x3 = v2*B1 - v1*B2 + (Shear_0 - 2.0*Omega_0)*v1;

  #endif
 
  return ft;
}
