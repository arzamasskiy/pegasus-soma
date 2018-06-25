#include "../copyright.h"
/*============================================================================*/
/*! \file utils_particle.c
 *  \brief Contains most of the utilities for the particle code.
 *
 * PURPOSE: Contains most of the utilities for the particle code: all the
 *   interpolation functions, stopping time calculation functions, shuffle
 *   algorithms. Also contained are the default (and trivial) gas velocity
 *   shift function. The get_gasinfo(Grid *pG) routine is used for test
 *   purposes only.
 * 
 * CONTAINS PUBLIC FUNCTIONS:
 * 
 * - getwei_linear()
 * - getwei_TSC   ()
 * - getwei_QP    ()
 * - getvalues()
 * - get_ts_epstein()
 * - get_ts_general()
 * - get_ts_fixed  ()
 * - get_gasinfo()
 * - feedback_clear()
 * - distrFB      ()
 * - void shuffle()
 * - void gasvshift_zero()
 * 
 * PRIVATE FUNCTION PROTOTYPES:
 * - compare_gr()         - compare the location of the two particles
 * - quicksort_particle() - sort the particles using the quicksort
 *
 *============================================================================*/
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
 *   compare_gr()         - compare the location of the two particles
 *   quicksort_particle() - sort the particles using the quicksort
 *============================================================================*/
static int compare_gr(GridS *pG, Real3Vect cell1, GrainS gr1, GrainS gr2);
static void quicksort_particle(GridS *pG, Real3Vect cell1, long start, long end);
static void countingsort_particle(GridS *pG, Real3Vect cell1);
static void inline filter_gpcouple(GPCouple *A, GPCouple *B, GPCouple *C,
                                   GPCouple *tmp, short lab);
static void inline clear_gpcouple(GPCouple *A, short lab);

/*============================== ALL FUNCTIONS ===============================*/
/*------------------------------------------------------------------------------
 * interpolation functions;
 * deposit related functions;
 * shuffle related functions;
 * filter related functions;
 *----------------------------------------------------------------------------*/

/*============================================================================*/
/*-----------------------------INTERPOLATION----------------------------------
 *
 * getwei_linear()
 * getwei_TSC()
 * getwei_QP ()
 * getvalues();
 */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void getwei_linear(GridS *pG, Real x1, Real x2, Real x3,
 *           Real3Vect cell1, Real weight[3][3][3], int *is, int *js, int *ks)
 *  \brief Get weight using linear interpolation
 *
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */

void getwei_linear(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c;		/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[2], wei2[2], wei3[2];/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (cell1.x1 > 0.0) {
    i = celli(pG, x1, cell1.x1, &i1, &a);	/* x1 index */
    i1 = i+i1-1;	*is = i1;		/* starting x1 index */
    wei1[1] = a - i1 - 0.5;			/* one direction weight */
    wei1[0] = 1.0 - wei1[1];			/* 0: left; 1: right */
  }
  else { /* x1 dimension collapses */
    *is = pG->is;
    wei1[1] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (cell1.x2 > 0.0) {
    j = cellj(pG, x2, cell1.x2, &j1, &b);	/* x2 index */
    j1 = j+j1-1;	*js = j1;		/* starting x2 index */
    wei2[1] = b - j1 - 0.5;			/* one direction weight */
    wei2[0] = 1.0 - wei2[1];			/* 0: left; 1: right */
  }
  else { /* x2 dimension collapses */
    *js = pG->js;
    wei2[1] = 0.0;
    wei2[0] = 1.0;
  }

  /* x3 direction */
  if (cell1.x3 > 0.0) {
    k = cellk(pG, x3, cell1.x3, &k1, &c);	/* x3 index */
    k1 = k+k1-1;	*ks = k1;		/* starting x3 index */
    wei3[1] = c - k1 - 0.5;			/* one direction weight */
    wei3[0] = 1.0 - wei3[1];			/* 0: left; 1: right */
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;
    wei3[1] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<2; k++)
    for (j=0; j<2; j++)
      for (i=0; i<2; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void getwei_TSC(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
 *                      Real weight[3][3][3], int *is, int *js, int *ks)
 *
 *  \brief Get weight using Triangular Shaped Cloud (TSC) interpolation 
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */
void getwei_TSC(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                           Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c, d;	/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[3], wei2[3], wei3[3];/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (cell1.x1 > 0.0) {
    celli(pG, x1, cell1.x1, &i, &a);		/* x1 index */
    *is = i - 1;				/* starting x1 index, wei[0] */
    d = a - i;
    wei1[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
    wei1[1] = 0.75-SQR(d-0.5);			/* one direction weight */
    wei1[2] = 0.5*SQR(d);
  }
  else { /* x1 dimension collapses */
    *is = pG->is;
    wei1[1] = 0.0;	wei1[2] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (cell1.x2 > 0.0) {
    cellj(pG, x2, cell1.x2, &j, &b);		/* x2 index */
    *js = j - 1;				/* starting x2 index */
    d = b - j;
    wei2[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
    wei2[1] = 0.75-SQR(d-0.5);			/* one direction weight */
    wei2[2] = 0.5*SQR(d);
  }
  else { /* x2 dimension collapses */
    *js = pG->js;
    wei2[1] = 0.0;	wei2[2] = 0.0;
    wei2[0] = 1.0;
  }

  /* x3 direction */
  if (cell1.x3 > 0.0) {
    cellk(pG, x3, cell1.x3, &k, &c);		/* x3 index */
    *ks = k - 1;				/* starting x3 index */
    d = c - k;
    wei3[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
    wei3[1] = 0.75-SQR(d-0.5);			/* one direction weight */
    wei3[2] = 0.5*SQR(d);
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;
    wei3[1] = 0.0;	wei3[2] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void getwei_QP(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
 *                       Real weight[3][3][3], int *is, int *js, int *ks)
 *  \brief Get weight using quadratic polynomial interpolation 
 *
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */
void getwei_QP(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                          Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c, d;	/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[3], wei2[3], wei3[3];/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (cell1.x1 > 0.0) {
    celli(pG, x1, cell1.x1, &i, &a);		/* x1 index */
    *is = i - 1;				/* starting x1 index, wei[0] */
    d = a - i;
    wei1[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
    wei1[1] = 1.0-SQR(d-0.5);			/* one direction weight */
    wei1[2] = 0.5*(d-0.5)*(d+0.5);
  }
  else { /* x1 dimension collapses */
    *is = pG->is;
    wei1[1] = 0.0;	wei1[2] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (cell1.x2 > 0.0) {
    cellj(pG, x2, cell1.x2, &j, &b);		/* x2 index */
    *js = j - 1;				/* starting x2 index */
    d = b - j;
    wei2[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
    wei2[1] = 1.0-SQR(d-0.5);			/* one direction weight */
    wei2[2] = 0.5*(d-0.5)*(d+0.5);
  }
  else { /* x2 dimension collapses */
    *js = pG->js;
    wei2[1] = 0.0;	wei2[2] = 0.0;
    wei2[0] = 1.0;
  }

  /* x3 direction */
  if (cell1.x3 > 0.0) {
    cellk(pG, x3, cell1.x3, &k, &c);		/* x3 index */
    *ks = k - 1;				/* starting x3 index */
    d = c - k;
    wei3[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
    wei3[1] = 1.0-SQR(d-0.5);			/* one direction weight */
    wei3[2] = 0.5*(d-0.5)*(d+0.5);
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;
    wei3[1] = 0.0;	wei3[2] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn int getvalues()
 *  \brief Get interpolated value using the weight
 *
 * Input:
 * - pG: grid; weight: weight function;
 * - is,js,ks: starting cell indices in the grid.
 * Output:
 * - interpolated values of B-field, and E-field
 *
 * Return: 0: normal exit; -1: particle lie out of the grid, cannot interpolate
 * Note: this interpolation works in any 1-3 dimensions.
 */

int getvalues(GridS *pG, Real weight[3][3][3], int is, int js, int ks
              ,Real *bfld1, Real *bfld2, Real *bfld3
              ,Real *efld1, Real *efld2, Real *efld3
#ifdef DRIVING
              ,Real *ffld1, Real *ffld2, Real *ffld3
#endif
              ){
  int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  Real b1, b2, b3, e1, e2, e3, eb;		/* E & B-fields */
#ifdef DRIVING
  Real f1, f2, f3;
#endif
  GPCouple *pq;
  Real totwei, totwei1;		/* total weight (in case of edge cells) */
  Real bsq, edotb;
  
  /* linear interpolation */
  b1 = 0.0; b2 = 0.0; b3 = 0.0; e1 = 0.0; e2 = 0.0; e3 = 0.0; eb = 0.0;
#ifdef DRIVING
  f1 = 0.0; f2 = 0.0; f3 = 0.0;
#endif
  totwei = 0.0;		totwei1 = 1.0;
  
  /* Interpolate E & B fields */
  /* Note: in lower dimensions only wei[0] is non-zero */
  n0 = ncell-1;
  k1 = MAX(ks, klp);	k2 = MIN(ks+n0, kup);
  j1 = MAX(js, jlp);	j2 = MIN(js+n0, jup);
  i1 = MAX(is, ilp);	i2 = MIN(is+n0, iup);
  
  for (k=k1; k<=k2; k++) {
    k0=k-k1;
    for (j=j1; j<=j2; j++) {
      j0=j-j1;
      for (i=i1; i<=i2; i++) {
        i0=i-i1;
        
        pq = &(pG->Coup[k][j][i]);
        b1 += weight[k0][j0][i0] * pq->grid_b1;
        b2 += weight[k0][j0][i0] * pq->grid_b2;
        b3 += weight[k0][j0][i0] * pq->grid_b3;
        e1 += weight[k0][j0][i0] * pq->grid_e1;
        e2 += weight[k0][j0][i0] * pq->grid_e2;
        e3 += weight[k0][j0][i0] * pq->grid_e3;
        eb += weight[k0][j0][i0] * pq->grid_eb;
#ifdef DRIVING
        f1 += weight[k0][j0][i0] * pG->force[k][j][i].x1;
        f2 += weight[k0][j0][i0] * pG->force[k][j][i].x2;
        f3 += weight[k0][j0][i0] * pG->force[k][j][i].x3;
#endif
        totwei += weight[k0][j0][i0];
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;
  
  /* ensure interpolated E dot interpolated B = interpolated E dot B */
  bsq   = SQR(b1)+SQR(b2)+SQR(b3);
  edotb = e1*b1+e2*b2+e3*b3;
  e1 += (eb*totwei-edotb)*b1/MAX(bsq,TINY_NUMBER);
  e2 += (eb*totwei-edotb)*b2/MAX(bsq,TINY_NUMBER);
  e3 += (eb*totwei-edotb)*b3/MAX(bsq,TINY_NUMBER);
  
  totwei1 = 1.0/totwei;
  *bfld1 = b1*totwei1; *bfld2 = b2*totwei1; *bfld3 = b3*totwei1;
  *efld1 = e1*totwei1; *efld2 = e2*totwei1; *efld3 = e3*totwei1;
#ifdef DRIVING
  *ffld1 = f1*totwei1; *ffld2 = f2*totwei1; *ffld3 = f3*totwei1;
#endif
  
  return 0;
}

/*----------------------------------------------------------------------------*/
/*! \fn void set_gpcouple(GridS *pG, short lab)
 *  \brief clear the gridded density and velocity; set to background if delta-f
 */
void set_gpcouple(GridS *pG, short lab)
{
  int i,j,k;
  GPCouple *pq, *pw;
    
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        clear_gpcouple(&(pG->Coup[k][j][i]),lab);
      }
  
#ifdef DELTA_F
  for (k=pG->ks; k<=pG->ke; k++)
    for (j=pG->js; j<=pG->je; j++)
      for (i=pG->is; i<=pG->ie; i++) {
        pq = (&(pG->Coup[k][j][i]));
        pw = (&(pG->Bkgrd[k][j][i]));
        
        pq->grid_d  += pw->grid_d;
        pq->grid_M1 += pw->grid_M1;
        pq->grid_M2 += pw->grid_M2;
        pq->grid_M3 += pw->grid_M3;
        if (lab == 0) {
          pq->grid_p11 += pw->grid_p11;
          pq->grid_p12 += pw->grid_p12;
          pq->grid_p13 += pw->grid_p13;
          pq->grid_p22 += pw->grid_p22;
          pq->grid_p23 += pw->grid_p23;
          pq->grid_p33 += pw->grid_p33;
        }
      }
#endif
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void deposit(GridS *pG)
 *  \brief Deposit primary particles on grid and store in Coup structure
 */
void deposit(DomainS *pD, short lab)
{
  GridS *pG=pD->Grid;
  int i,j,k,is,js,ks,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  int n0 = ncell-1;
  long p;
  Real weight[3][3][3],drho;
  Real3Vect cell1;
  GrainS *gr;
  GPCouple *pq;
#ifdef DELTA_F
  Real f0_t;
#endif
  
  if (pG->Nx[0] > 1) cell1.x1 = 1.0/pG->dx1;  else  cell1.x1 = 0.0;
  if (pG->Nx[1] > 1) cell1.x2 = 1.0/pG->dx2;  else  cell1.x2 = 0.0;
  if (pG->Nx[2] > 1) cell1.x3 = 1.0/pG->dx3;  else  cell1.x3 = 0.0;
  
  /* initialization */
  set_gpcouple(pG,lab);

  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);
    
    drho = grproperty[gr->property].m;
    
#ifdef DELTA_F
    getdf(gr, &f0_t);
    drho *= (1.0 - (f0_t/gr->f_0));
#endif

    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);
    
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

          /* interpolate the particles to the grid */
          pq->grid_d  += weight[k0][j0][i0]*drho;
          pq->grid_M1 += weight[k0][j0][i0]*drho*gr->v1;
          pq->grid_M2 += weight[k0][j0][i0]*drho*gr->v2;
          pq->grid_M3 += weight[k0][j0][i0]*drho*gr->v3;
          if (lab == 0) {
            pq->grid_p11 += weight[k0][j0][i0]*drho*gr->v1*gr->v1;
            pq->grid_p12 += weight[k0][j0][i0]*drho*gr->v1*gr->v2;
            pq->grid_p13 += weight[k0][j0][i0]*drho*gr->v1*gr->v3;
            pq->grid_p22 += weight[k0][j0][i0]*drho*gr->v2*gr->v2;
            pq->grid_p23 += weight[k0][j0][i0]*drho*gr->v2*gr->v3;
            pq->grid_p33 += weight[k0][j0][i0]*drho*gr->v3*gr->v3;
          }
        }
      }
    }
        
  }
  
  exchange_gpcouple(pD,lab);
  if (nfpass > 0) smooth_gpcouple(pD,lab);
  
  return;
}

/*============================================================================*/
/*---------------------------------SHUFFLE------------------------------------
 *
 * shuffle()
 * compare_gr()
 * quicksort_particle()
 */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void Delete_Ghost(GridS *pG)
 *  \brief Delete ghost particles */
void Delete_Ghost(GridS *pG)
{
  long p;
  GrainS *gr;
  
  p = 0;
  while (p<pG->nparticle)
  {/* loop over all particles */
    gr = &(pG->particle[p]);
    
    if (gr->pos == 0)
    {/* gr is a ghost particle */
      pG->nparticle -= 1;
      grproperty[gr->property].num -= 1;
      pG->particle[p] = pG->particle[pG->nparticle];
    }
    else {
      p++;
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shuffle(GridS *pG)
 *  \brief Shuffle the particles
 *
 * Input: pG: grid with particles;
 * Output: pG: particles in the linked list are rearranged by the order of their
 *         locations that are consistent with grid cell storage.
 */
void shuffle(GridS *pG)
{
  Real3Vect cell1;

  if (pG->Nx[0] > 1) cell1.x1 = 1.0/pG->dx1;  else  cell1.x1 = 0.0;
  if (pG->Nx[1] > 1) cell1.x2 = 1.0/pG->dx2;  else  cell1.x2 = 0.0;
  if (pG->Nx[2] > 1) cell1.x3 = 1.0/pG->dx3;  else  cell1.x3 = 0.0;

  /* output status */
  peg_pout(0, "Resorting particles...\n");

  /* sort the particles according to their positions */
  //quicksort_particle(pG, cell1, 0, pG->nparticle-1);
  countingsort_particle(pG, cell1);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn int compare_gr(GridS *pG, Real3Vect cell1, GrainS gr1, GrainS gr2)
 *  \brief Compare the order of two particles according to their positions in 
 *  the grid
 *
 * Input: pG: grid; 
 * -      cell1: 1/dx1,1/dx2,1/dx3, or 0 if that dimension collapses.
 * -      gr1,gr2: pointers of the two particles to be compared.
 * Output: pointer of the particle that should be put in front of the other.
 */
static int compare_gr(GridS *pG, Real3Vect cell1, GrainS gr1, GrainS gr2)
{
  int i1,j1,k1, i2,j2,k2;

  k1 = (int)((gr1.x3 - pG->MinX[2]) * cell1.x3);	/* x3 index of gr1 */
  k2 = (int)((gr2.x3 - pG->MinX[2]) * cell1.x3);	/* x3 index of gr2 */
  if (k1 < k2) return 1;
  if (k1 > k2) return 2;

  j1 = (int)((gr1.x2 - pG->MinX[1]) * cell1.x2);	/* x2 index of gr1 */
  j2 = (int)((gr2.x2 - pG->MinX[1]) * cell1.x2);	/* x2 index of gr2 */
  if (j1 < j2) return 1;
  if (j1 > j2) return 2;

  i1 = (int)((gr1.x1 - pG->MinX[0]) * cell1.x1);	/* x1 index of gr1 */
  i2 = (int)((gr2.x1 - pG->MinX[0]) * cell1.x1);	/* x1 index of gr2 */
  if (i1 < i2) return 1;
  if (i1 > i2) return 2;

  /* if they have equal indices, arbitrarily choose gr1 */
  return 1;
}


/*----------------------------------------------------------------------------*/
/*! \fn int getCellID(double *x, double *xMin, int *n, Real3Vect cell1)
 *  \brief Calculate the cell ID (position in memory) of a single particle on a 
 *  grid 
 *  
 *
 * Input: x: particle positions 
 * -      xMin: left-most edge locations of grid
 * -      n: number of grid points (without ghosts)
 *        cell1: Vector of inverse grid sizes (1/dx)
 * Output: cellID (position in particle array)
 */
static inline int getCellID(double *x, double *xMin, int *n, Real3Vect cell1){
  int iz = (int) ((x[2] - xMin[2]) * cell1.x3);
  int iy = (int) ((x[1] - xMin[1]) * cell1.x2);
  int ix = (int) ((x[0] - xMin[0]) * cell1.x1);

  return (n[1]*n[0])*iz + n[0]*iy + ix;
}
/*----------------------------------------------------------------------------*/
/*! \fn void countingsort_particle(Grid *pG, Vector cell1)
 *  \brief Counting sort algorithm to shuffle the particles
 *
 * Input: pG, cell1: for CellID subroutine only. See above.
 */


static void countingsort_particle(GridS *pG, Real3Vect cell1)
{
  static int *pa = NULL, *pa_save = NULL;
  int i,j,k,se, cellid;
  double pos[3];

  int ncell = pG->Nx[0]*pG->Nx[1]*pG->Nx[2];

  GrainS *in, *out, *stop, *tmp;
  GrainS tmp_s;
  tmp = &tmp_s;

  if(pa == NULL || pa_save == NULL){
    pa = malloc(sizeof(int)*(ncell + 1));
    pa_save = malloc(sizeof(int)*(ncell + 1));
  }
  
  for(i = 0; i < ncell; i++) pa[i] = 0;

  for(i = 0; i < pG->nparticle; i++){
      in = &pG->particle[i];
      pos[0] = in->x1;
      pos[1] = in->x2;
      pos[2] = in->x3;
      cellid = getCellID(pos, pG->MinX,pG->Nx,cell1);
      pa[cellid]++;
  }

  for(i=k=0; i <= ncell; i++){
    j = pa[i];
    pa_save[i] = pa[i] = k;
    k += j;
  }

  i = 0;
  while( i < ncell){
    if(pa[i] >= pa_save[i+1]){
      i++;
    } else {
      in = stop = &pG->particle[pa[i]];
      do{
        pos[0] = in->x1;
        pos[1] = in->x2;
        pos[2] = in->x3;
        cellid = getCellID(pos, pG->MinX,pG->Nx,cell1);
        out = &pG->particle[pa[cellid]];
        pa[cellid]++;
        if(out != stop){
          tmp_s = *out;
          *out = *in;
          *in = tmp_s;
        } else if(out != in){
          *out = *in;
        }
      } while( out != stop);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void quicksort_particle(Grid *pG, Vector cell1, long start, long end)
 *  \brief Quick sort algorithm to shuffle the particles
 *
 * Input: pG, cell1: for compare_gr subroutine only. See above.
 *        head, rear: head and rear of the linked list.
 *          They do not contain data, or equal the pivot in the recursion.
 *        length: length of the linked list (does not contain head or rear).
 * Output: *head: linked list with shuffling finished.
 */
static void quicksort_particle(GridS *pG, Real3Vect cell1, long start, long end)
{
  long i, pivot;
  GrainS gr;
  if (end <= start) return;	/* automatically sorted already */

  /* location of the pivot at half chain length */
  pivot = (long)((start+end+1)/2);

  /* move the pivot to the start */
  gr = pG->particle[pivot];
  pG->particle[pivot] = pG->particle[start];
  pG->particle[start] = gr;

  /* initial configuration */
  pivot = start;
  i = start + 1;

  /* move the particles that are "smaller" than the pivot before it */
  while (i <= end) {
    if (compare_gr(pG, cell1, pG->particle[pivot], pG->particle[i]) == 2)
    {/* the ith particle is smaller, move it before the pivot */
      gr = pG->particle[pivot];
      pG->particle[pivot] = pG->particle[i];
      pG->particle[i] = pG->particle[pivot+1];
      pG->particle[pivot+1] = gr;
      pivot += 1;
    }
    i += 1;
  }

  /* recursively call this routine to complete sorting */
  quicksort_particle(pG, cell1, start, pivot-1);
  quicksort_particle(pG, cell1, pivot+1, end);

  return;
}

/*============================================================================*/
/*---------------------------------FILTER--------------------------------------
 *
 * smooth_gpcouple()
 * filter_gpcouple()
 * clear_gpcouple()
 */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void smooth_gpcouple(DomainS *pD, short lab)
 *  \brief handles smoothing and (possible) re-filling of ghost zones
 *
 *  Notes: Each filter pass ruins one layer of ghost zones. Five layers of 
 *         ghost zones holding deposited quantities are required for the hybrid 
 *         integrator to work at the prediction step; three layers of ghost zones
 *         are required at the correction step. There are 5 layers of ghost zones 
 *         in total. In prediction step, we filter requested number of layers, ruin-
 *         ing them as we go with the knowledge that they'll be refilled. In the 
 *         correction step, if =< 2 filter passes are requested, there are still 
 *         enough accurate ghost zones after filtering. If nfpass > 2, then we have 
 *         to clear the ghost zones and re-fill them with accurate values after the 
 *         filter passes are complete.
 *         If smoothing for output purposes (lab = 0), then there is no need to
 *         re-fill ghost zones, but we must include the pressure tensor.
 */
void smooth_gpcouple(DomainS *pD, short lab)
{ 
  GridS *pG=pD->Grid;
  int i,j,k,ib,it,jb,jt,kb,kt,n,di,dj,dk;
  GPCouple tmp1,tmp2,tmp3;
  
  if (nfpass == 0) return;
  
  switch (lab) {
      
    case 0:
      ib = pG->is-nfpass;
      it = pG->ie+nfpass;
      if (pG->Nx[1] > 1) {
        jb = pG->js-nfpass;
        jt = pG->je+nfpass;
      } else {
        jb = pG->js;
        jt = pG->je;
      }
      if (pG->Nx[2] > 1) {
        kb = pG->ks-nfpass;
        kt = pG->ke+nfpass;
      } else {
        kb = pG->ks;
        kt = pG->ke;
      }
      break;

    case 1:
      ib = pG->is-nfpass;
      it = pG->ie+nfpass;
      if (pG->Nx[1] > 1) {
        jb = pG->js-nfpass;
        jt = pG->je+nfpass;
      } else {
        jb = pG->js;
        jt = pG->je;
      }
      if (pG->Nx[2] > 1) {
        kb = pG->ks-nfpass;
        kt = pG->ke+nfpass;
      } else {
        kb = pG->ks;
        kt = pG->ke;
      }
      break;
      
    case 2:
      ib = nfpass > 2 ? pG->is-nfpass : pG->is-nfpass-3;
      it = nfpass > 2 ? pG->ie+nfpass : pG->ie+nfpass+3;
      if (pG->Nx[1] > 1) {
        jb =  nfpass > 2 ? pG->js-nfpass : pG->js-nfpass-3;
        jt =  nfpass > 2 ? pG->je+nfpass : pG->je+nfpass+3;
      } else {
        jb = pG->js;
        jt = pG->je;
      }
      if (pG->Nx[2] > 1) {
        kb =  nfpass > 2 ? pG->ks-nfpass: pG->ks-nfpass-3;
        kt =  nfpass > 2 ? pG->ke+nfpass: pG->ke+nfpass+3;
      } else {
        kb = pG->ks;
        kt = pG->ke;
      }
      break;

    default:
      peg_error("[smooth_gpcouple]: stage undefined!\n");

  }
    
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (n=1; n<=nfpass; n++) {
        // filter x1-pencil of constant x2 and x3
        tmp2 = pG->Coup[k][j][ib];
        filter_gpcouple(&(pG->Coup[k][j][it-2]),
                        &(pG->Coup[k][j][it-1]),
                        &(pG->Coup[k][j][it]),&(tmp3),lab);
        i = ib+1;
        while (i<=it-2) {
          filter_gpcouple(&(pG->Coup[k][j][i-1]),
                          &(pG->Coup[k][j][i  ]),
                          &(pG->Coup[k][j][i+1]),&(tmp1),lab);
          pG->Coup[k][j][i-1] = tmp2;
          i++;
          filter_gpcouple(&(pG->Coup[k][j][i-1]),
                          &(pG->Coup[k][j][i  ]),
                          &(pG->Coup[k][j][i+1]),&(tmp2),lab);
          pG->Coup[k][j][i-1] = tmp1;
          i++;
        }
        pG->Coup[k][j][i-1]  = tmp2;
        pG->Coup[k][j][it-1] = tmp3; //Slight redundancy for even grids to fix odd grids
      }
    }
  }
    
  if (pG->Nx[1] > 1) {
    
    for (k=kb; k<=kt; k++) {
      for (i=ib+nfpass; i<=it-nfpass; i++) {
        for (n=1; n<=nfpass; n++) {
          // filter x2-pencil of constant x1 and x3
          tmp2 = pG->Coup[k][jb][i];
          filter_gpcouple(&(pG->Coup[k][jt-2][i]),
                          &(pG->Coup[k][jt-1][i]),
                          &(pG->Coup[k][jt][i]),&(tmp3),lab);
          j = jb+1;
          while (j<=jt-2) {
            filter_gpcouple(&(pG->Coup[k][j-1][i]),
                            &(pG->Coup[k][j  ][i]),
                            &(pG->Coup[k][j+1][i]),&(tmp1),lab);
            pG->Coup[k][j-1][i] = tmp2;
            j++;
            filter_gpcouple(&(pG->Coup[k][j-1][i]),
                            &(pG->Coup[k][j  ][i]),
                            &(pG->Coup[k][j+1][i]),&(tmp2),lab);
            pG->Coup[k][j-1][i] = tmp1;
            j++;
          }
          pG->Coup[k][j-1][i]  = tmp2;
          pG->Coup[k][jt-1][i] = tmp3; //Slight redundancy for even grids to fix odd grids
        }
      }
    }
    
  }
  
  if (pG->Nx[2] > 1) {
    
    for (j=jb+nfpass; j<=jt-nfpass; j++) {
      for (i=ib+nfpass; i<=it-nfpass; i++) {
        for (n=1; n<=nfpass; n++) {
          // filter x3-pencil of constant x1 and x2
          tmp2 = pG->Coup[kb][j][i];
          filter_gpcouple(&(pG->Coup[kt-2][j][i]),
                          &(pG->Coup[kt-1][j][i]),
                          &(pG->Coup[kt][j][i]),&(tmp3),lab);
          k = kb+1;
          while (k<=kt-2) {
            filter_gpcouple(&(pG->Coup[k-1][j][i]),
                            &(pG->Coup[k  ][j][i]),
                            &(pG->Coup[k+1][j][i]),&(tmp1),lab);
            pG->Coup[k-1][j][i] = tmp2;
            k++;
            filter_gpcouple(&(pG->Coup[k-1][j][i]),
                            &(pG->Coup[k  ][j][i]),
                            &(pG->Coup[k+1][j][i]),&(tmp2),lab);
            pG->Coup[k-1][j][i] = tmp1;
            k++;
          }
          pG->Coup[k-1][j][i]  = tmp2;
          pG->Coup[kt-1][j][i] = tmp3; //Slight redundancy for even grids to fix odd grids
        }
      }
    }
    
  }

  if ((lab == 1) || (lab == 2 && nfpass > 2)) {
    
    /* clear ghost zones of density and momentum */
    di = pG->ie+1-ilp;
    for (k=klp; k<=kup; k++) {
      for (j=jlp; j<=jup; j++) {
        for (i=ilp; i<pG->is; i++) {
          clear_gpcouple(&(pG->Coup[k][j][i  ]),lab);
          clear_gpcouple(&(pG->Coup[k][j][i+di]),lab);
        }
      }
    }
   
    if (pG->Nx[1] > 1) {
      dj = pG->je+1-jlp;
      for (k=klp; k<=kup; k++) {
        for (j=jlp; j<pG->js; j++) {
          for (i=pG->is; i<=pG->ie; i++) {
            clear_gpcouple(&(pG->Coup[k][j   ][i]),lab);
            clear_gpcouple(&(pG->Coup[k][j+dj][i]),lab);
          }
        }
      }    
    }
  
    if (pG->Nx[2] > 1) {  
      dk = pG->ke+1-klp;
      for (k=klp; k<pG->ks; k++) {
        for (j=pG->js; j<=pG->je; j++) {
          for (i=pG->is; i<=pG->ie; i++) {
            clear_gpcouple(&(pG->Coup[k   ][j][i]),lab);
            clear_gpcouple(&(pG->Coup[k+dk][j][i]),lab);
          }
        }
      }
    }
  
    /* re-fill ghost zones with density and momentum */
    exchange_gpcouple(pD,2+lab);
        
  } /* end of if statement for lab != 0 && nfpass !=1 */
  
  return;
}

static void inline filter_gpcouple(GPCouple *A, GPCouple *B, GPCouple *C, 
                                   GPCouple *tmp, short lab) 
{
  tmp->grid_d  = 0.25 * A->grid_d  + 0.5 * B->grid_d  + 0.25 * C->grid_d;
  tmp->grid_M1 = 0.25 * A->grid_M1 + 0.5 * B->grid_M1 + 0.25 * C->grid_M1;
  tmp->grid_M2 = 0.25 * A->grid_M2 + 0.5 * B->grid_M2 + 0.25 * C->grid_M2;
  tmp->grid_M3 = 0.25 * A->grid_M3 + 0.5 * B->grid_M3 + 0.25 * C->grid_M3;
  if (lab == 0) {
    tmp->grid_p11 = 0.25 * A->grid_p11 + 0.5 * B->grid_p11 + 0.25 * C->grid_p11;
    tmp->grid_p12 = 0.25 * A->grid_p12 + 0.5 * B->grid_p12 + 0.25 * C->grid_p12;
    tmp->grid_p13 = 0.25 * A->grid_p13 + 0.5 * B->grid_p13 + 0.25 * C->grid_p13;
    tmp->grid_p22 = 0.25 * A->grid_p22 + 0.5 * B->grid_p22 + 0.25 * C->grid_p22;
    tmp->grid_p23 = 0.25 * A->grid_p23 + 0.5 * B->grid_p23 + 0.25 * C->grid_p23;
    tmp->grid_p33 = 0.25 * A->grid_p33 + 0.5 * B->grid_p33 + 0.25 * C->grid_p33;
  }
  
  return;
}

static void inline clear_gpcouple(GPCouple *A, short lab)
{
  A->grid_d  = 0.0;
  A->grid_M1 = 0.0;
  A->grid_M2 = 0.0;
  A->grid_M3 = 0.0;
  if (lab == 0) {
    A->grid_p11 = 0.0;
    A->grid_p12 = 0.0;
    A->grid_p13 = 0.0;
    A->grid_p22 = 0.0;
    A->grid_p23 = 0.0;
    A->grid_p33 = 0.0;
  }
  
  return;
}
