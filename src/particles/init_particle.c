#include "../copyright.h"
/*===========================================================================*/
/*! \file init_particle.c
 *  \brief Initialize particle related structures and functions.
 *
 * PURPOSE: Initialize particle related structures and functions. Particle
 *   integrator is enrolled by calling integrate_particle_init(int type). In
 *   init_particle(Grid *pG, Domain *pD), all the particle related arrays are
 *   allocated, interpolation and stopping time calculation functions are
 *   enrolled, most global variables for particles are evaluated. Also contains
 *   functions for particle array reallocation and destruction.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - integrate_particle_init();
 * - init_particle();
 * - particle_destruct();
 * - particle_realloc();
 *                                                                            */
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
 *   grid_limit()   - get the limit coordinates of the grid
 *============================================================================*/
void grid_limit(MeshS *pM);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*! \fn void init_particle(MeshS *pM)
 *  \brief Initialization for particles
 *
 * Allocate memory for the gas velocity/sound speed array, feedback array.
 * Note we enforce that each type has equal number of particles to ensure equal
 * resolution.
 */
void init_particle(MeshS *pM)
{
  int i, N1T, N2T, N3T, integratortype, interp, tsmode, distfunc;
  GrainS *GrArray;
  long size = 1000, size1 = 1, size2 = 1;
  DomainS *pD;
  GridS   *pG;

  pD = (DomainS*)&(pM->Domain[0][0]);  /* set ptr to Domain */
  pG = pD->Grid;          /* set ptr to Grid */

  /* get coordinate limit */
  grid_limit(pM);
  N1T = iup-ilp+1;
  N2T = jup-jlp+1;
  N3T = kup-klp+1;

  /* check particle types */
  npartypes = par_geti_def("particle","partypes",1);

  if (npartypes < 0)
    peg_error("[init_particle]: Particle types must not be negative!\n");

  /* initialize the particle array */
  if(par_exist("particle","parnumcell"))
  {
    /* if we consider number of particles per cell */
    size1 = N1T*N2T*N3T*(long)(par_geti("particle","parnumcell"));
    if (size1 < 0)
      peg_error("[init_particle]: Particle number must be non-negative!\n");
  }

  if(par_exist("particle","parnumgrid"))
  {
    /* if we consider number of particles per grid */
    size2 = (long)(npartypes*par_geti("particle","parnumgrid"));
    if (size2 < 0)
      peg_error("[init_particle]: Particle number must be non-negative!\n");
    /* account for the ghost cells and different types */
    size2 = (long)(size2/((double)(pG->Nx[0]*pG->Nx[1]*pG->Nx[2]))*N1T*N2T*N3T);

  }

  size = MAX(size, MAX(size1, size2));
  pG->arrsize = (long)(1.2*size);   /* account for number fluctuations */

  pG->particle = (GrainS*)calloc_1d_array(pG->arrsize, sizeof(GrainS));
  if (pG->particle == NULL) goto on_error;

  pG->parsub   = (GrainAux*)calloc_1d_array(pG->arrsize,sizeof(GrainAux));
  if (pG->parsub == NULL) goto on_error;

  /* allocate memory for particle properties */
  grproperty = (Grain_Property*)calloc_1d_array(npartypes,
                                              sizeof(Grain_Property));
  if (grproperty == NULL) goto on_error;

  /* Set particle integrator to according to particle types */
  for (i=0; i<npartypes; i++)
    grproperty[i].integrator =par_geti_def("particle","integrator",1);

  /* set the interpolation function pointer */
  interp = par_geti_def("particle","interp",2);
  if (interp == 1)
  { /* linear interpolation */
    getweight = getwei_linear;
    ncell = 2;
  }
  else if (interp == 2)
  { /* TSC interpolation */
    getweight = getwei_TSC;
    ncell = 3;
  }
  else if (interp == 3)
  { /* Quadratic polynomial interpolation */
    getweight = getwei_QP;
    ncell = 3;
  }
  else
    peg_error("[init_particle]: Value of interp must be 1, 2 or 3!\n");

  /* set the distribution function pointers */
  distfunc = par_geti_def("problem","distfunc",2);
  switch (distfunc)
  {
  case 1: /* static distribution */
    initdf = initdf_static;
    break;

  case 2: /* Maxwellian distribution */
    initdf = initdf_maxwell;
#ifdef DELTA_F
    getdf = getdf_maxwell;
    setbg = setbg_maxwell;
#endif
    break;

  case 3: /* Bi-Maxwellian distribution */
    initdf = initdf_bimaxwell;
#ifdef DELTA_F
    getdf = getdf_bimaxwell;
    setbg = setbg_bimaxwell;
#endif
    break;

  case 4: /* User-defined distrbution */
    initdf = initdf_user;
#ifdef DELTA_F
    getdf = getdf_user;
    setbg = setbg_user;
#endif
    break;

  default: peg_error("[init_particle]: Value of distfunc must be 1, 2, 3, or 4!\n");
    return;
  }

  /* set number of filter passes */
  nfpass = par_geti_def("particle","nfpass",0);
  if (nfpass < 0)
    peg_error("[init_particle]: Number of filter passes must not be negative!\n");

  /* limit number of passes to number of ghost zones */
  nfpass = MIN(nfpass,nghost);
  
  /* set charge-to-mass ratio for particles */
  for (i=0; i<npartypes; i++)
    grproperty[i].qomc = par_getd_def("particle","qomc",1.0);
  
  /* allocate the memory for gas-particle coupling array */
  pG->Coup = (GPCouple***)calloc_3d_array(N3T,N2T,N1T, sizeof(GPCouple));
  if (pG->Coup == NULL) goto on_error;
  
#ifdef DELTA_F
  /* allocate the memory for the coupling array of background quanities */
  pG->Bkgrd = (GPCouple***) calloc_3d_array(N3T,N2T,N1T, sizeof(GPCouple));
  if (pG->Bkgrd == NULL) goto on_error;
#endif
  
  return;

  on_error:
    peg_error("[init_particle]: Error allocating memory.\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void particle_destruct(MeshS *pM)
 *  \brief Finalization for particles
 */
void particle_destruct(MeshS *pM)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]); 
  GridS *pG = pD->Grid; 

  free_1d_array(pG->particle);
  free_1d_array(pG->parsub);
  free_1d_array(grproperty);

  /* free memory for feedback and delta-f background arrays */
  if (pG->Coup != NULL) free_3d_array(pG->Coup);
#ifdef DELTA_F
  if (pG->Bkgrd != NULL) free_3d_array(pG->Bkgrd);
#endif

  initdf = NULL;
#ifdef DELTA_F
  getdf = NULL;
  setbg = NULL;
#endif
  getweight = NULL;

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void particle_realloc(Grid *pG, long n)
 *  \brief Enlarge the particle array
 */
void particle_realloc(GridS *pG, long n)
{
  pG->arrsize = MAX((long)(1.2*pG->arrsize), n);

  /* for the main particle array */
  if ((pG->particle = (GrainS*)realloc(pG->particle,
                                      pG->arrsize*sizeof(GrainS))) == NULL)
  {
    peg_error("[init_particle]: Error re-allocating memory with array size\
 %ld.\n", n);
  }

  /* for the auxilary array */
  if ((pG->parsub = (GrainAux*)realloc(pG->parsub,
                                      pG->arrsize*sizeof(GrainAux))) == NULL)
  {
    peg_error("[init_particle]: Error re-allocating memory with array size\
 %ld.\n", n);
  }

  return;
}

/*============================================================================*/
/*----------------------------- Private Functions ----------------------------*/

/*----------------------------------------------------------------------------*/
/*! \fn void grid_limit(MeshS *pM)
 *  \brief Calculate the left and right grid limit
 *
 * Input: pG: grid;
 * Output: ilp,iup,jlp,jup,klp,kup: grid limit indices;
 *         x1lpar,x1upar,x2lpar,x2upar,x3lpar,x3upar: grid boundary coordinates
 */
void grid_limit(MeshS *pM)
{
  DomainS *pD;
  GridS   *pG;
  pD = (DomainS*)&(pM->Domain[0][0]);
  pG = pD->Grid; 

  int m1, m2, m3;	/* dimension flags */
  int my_iproc, my_jproc, my_kproc;

  if (pG->Nx[0] > 1) m1 = 1;
  else m1 = 0;

  if (pG->Nx[1] > 1) m2 = 1;
  else m2 = 0;

  if (pG->Nx[2] > 1) m3 = 1;
  else m3 = 0;

/* set left and right grid indices */
  ilp = pG->is - m1*nghost;
  iup = pG->ie + m1*nghost;

  jlp = pG->js - m2*nghost;
  jup = pG->je + m2*nghost;

  klp = pG->ks - m3*nghost;
  kup = pG->ke + m3*nghost;

/* set left and right boundary for removing particles
 * Note: for outflow B.C. (ibc=2), we only remove the particles in
 * the outermost layer of the ghost cells
 */
  get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

  if ((pM->BCFlag_ix1 == 2) && (my_iproc == 0))
    x1lpar = pG->MinX[0] - (nghost-1)*m1*pG->dx1;
  else
    x1lpar = pG->MinX[0];

  if ((pM->BCFlag_ox1 == 2) && (my_iproc == pD->NGrid[0]-1))
    x1upar = pG->MaxX[0] + (nghost-1)*m1*pG->dx1;
  else
    x1upar = pG->MaxX[0];

  if ((pM->BCFlag_ix2 == 2) && (my_jproc == 0))
    x2lpar = pG->MinX[1] - (nghost-1)*m2*pG->dx2;
  else
    x2lpar = pG->MinX[1];

  if ((pM->BCFlag_ox2 == 2) && (my_jproc == pD->NGrid[1]-1))
    x2upar = pG->MaxX[1] + (nghost-1)*m2*pG->dx2;
  else
    x2upar = pG->MaxX[1];

  if ((pM->BCFlag_ix3 == 2) && (my_kproc == 0))
    x3lpar = pG->MinX[2] - (nghost-1)*m3*pG->dx3;
  else
    x3lpar = pG->MinX[2];

  if ((pM->BCFlag_ox3 == 2) && (my_kproc == pD->NGrid[2]-1))
    x3upar = pG->MaxX[2] + (nghost-1)*m3*pG->dx3;
  else
    x3upar = pG->MaxX[2];

  /* increase x?upar in non-existing dimensions  *
   * to avoid incorrectly labeling the particles */
  if (pG->Nx[0] == 1) x1upar += 1.0;
  if (pG->Nx[1] == 1) x2upar += 1.0;
  if (pG->Nx[2] == 1) x3upar += 1.0;

  return;
}
