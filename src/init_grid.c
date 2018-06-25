#include "copyright.h"
/*============================================================================*/
/*! \file init_grid.c 
 *  \brief Initializes most variables in the Grid structure.
 *
 * PURPOSE: Initializes most variables in the Grid structure.  Allocates memory
 *   for 3D arrays of Cons, interface B, etc.  With SMR, finds all overlaps
 *   between child and parent Grids, and initializes data needed for restriction
 *   flux-correction, and prolongation steps.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - init_grid()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - checkOverlap() - checks for overlap of cubes, and returns overlap coords
 * - checkOverlapTouch() - same as above, but checks for overlap and/or touch */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "pegasus.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void init_grid(MeshS *pM)
 *  \brief Initializes most variables in the Grid structure.
 */

void init_grid(MeshS *pM)
{
  DomainS *pD;
  GridS *pG;
  int nDim,nl,nd,myL,myM,myN;
  int i,l,m,n,n1z,n2z,n3z,n1p,n2p,n1r,n2r;

/* number of dimensions in Grid. */
  nDim=1;
  for (i=1; i<3; i++) if (pM->Nx[i]>1) nDim++;

/* Loop over all levels and domains per level */

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */
      pG = pM->Domain[nl][nd].Grid;          /* set ptr to Grid */

      pG->time = pM->time;

/* get (l,m,n) coordinates of Grid being updated on this processor */

      get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);

/* ---------------------  Intialize grid in 1-direction --------------------- */
/* Initialize is,ie,dx1
 * Compute Disp, MinX[0], and MaxX[0] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[0] = pD->GData[myN][myM][myL].Nx[0];

      if(pG->Nx[0] > 1) {
        pG->is = nghost;
        pG->ie = pG->Nx[0] + nghost - 1;
      }
      else
        pG->is = pG->ie = 0;
    
      pG->dx1 = pD->dx[0];
    
      pG->Disp[0] = pD->Disp[0];
      pG->MinX[0] = pD->MinX[0];
      for (l=1; l<=myL; l++) {
        pG->Disp[0] +=        pD->GData[myN][myM][l-1].Nx[0];
        pG->MinX[0] += (Real)(pD->GData[myN][myM][l-1].Nx[0])*pG->dx1;
      }
      pG->MaxX[0] = pG->MinX[0] + (Real)(pG->Nx[0])*pG->dx1;
    
/* ---------------------  Intialize grid in 2-direction --------------------- */
/* Initialize js,je,dx2
 * Compute Disp, MinX[1], and MaxX[1] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[1] = pD->GData[myN][myM][myL].Nx[1];
    
      if(pG->Nx[1] > 1) {
        pG->js = nghost;
        pG->je = pG->Nx[1] + nghost - 1;
      }
      else
        pG->js = pG->je = 0;
    
      pG->dx2 = pD->dx[1];

      pG->Disp[1] = pD->Disp[1];
      pG->MinX[1] = pD->MinX[1];
      for (m=1; m<=myM; m++) {
        pG->Disp[1] +=        pD->GData[myN][m-1][myL].Nx[1];
        pG->MinX[1] += (Real)(pD->GData[myN][m-1][myL].Nx[1])*pG->dx2;
      }
      pG->MaxX[1] = pG->MinX[1] + (Real)(pG->Nx[1])*pG->dx2;

/* ---------------------  Intialize grid in 3-direction --------------------- */
/* Initialize ks,ke,dx3
 * Compute Disp, MinX[2], and MaxX[2] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[2] = pD->GData[myN][myM][myL].Nx[2];

      if(pG->Nx[2] > 1) {
        pG->ks = nghost;
        pG->ke = pG->Nx[2] + nghost - 1;
      }
      else
        pG->ks = pG->ke = 0;

      pG->dx3 = pD->dx[2];

      pG->Disp[2] = pD->Disp[2];
      pG->MinX[2] = pD->MinX[2];
      for (n=1; n<=myN; n++) {
        pG->Disp[2] +=        pD->GData[n-1][myM][myL].Nx[2];
        pG->MinX[2] += (Real)(pD->GData[n-1][myM][myL].Nx[2])*pG->dx3;
      }
      pG->MaxX[2] = pG->MinX[2] + (Real)(pG->Nx[2])*pG->dx3;

/* ---------  Allocate 3D arrays to hold Cons based on size of grid --------- */

      if (pG->Nx[0] > 1)
        n1z = pG->Nx[0] + 2*nghost;
      else
        n1z = 1;

      if (pG->Nx[1] > 1)
        n2z = pG->Nx[1] + 2*nghost;
      else
        n2z = 1;

      if (pG->Nx[2] > 1)
        n3z = pG->Nx[2] + 2*nghost;
      else
        n3z = 1;

/* Build a 3D array of type ConsS */

      pG->U = (ConsS***)calloc_3d_array(n3z, n2z, n1z, sizeof(ConsS));
      if (pG->U == NULL) goto on_error1;

    
/* Build 3D arrays to hold interface fields */

      pG->B1i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->B1i == NULL) goto on_error2;

      pG->B2i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->B2i == NULL) goto on_error3;

      pG->B3i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->B3i == NULL) goto on_error4;
      
      pG->E1i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->E1i == NULL) goto on_error5;
      
      pG->E2i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->E2i == NULL) goto on_error6;
      
      pG->E3i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->E3i == NULL) goto on_error7;

/*-- Get IDs of neighboring Grids in Domain communicator ---------------------*/
/* If Grid is at the edge of the Domain (so it is either a physical boundary,
 * or an internal boundary between fine/coarse grids), then ID is set to -1
 */

/* Left-x1 */
      if(myL > 0) pG->lx1_id = pD->GData[myN][myM][myL-1].ID_Comm_Domain;
      else pG->lx1_id = -1;

/* Right-x1 */
      if(myL <(pD->NGrid[0])-1)
        pG->rx1_id = pD->GData[myN][myM][myL+1].ID_Comm_Domain;
      else pG->rx1_id = -1;

/* Left-x2 */
      if(myM > 0) pG->lx2_id = pD->GData[myN][myM-1][myL].ID_Comm_Domain;
      else pG->lx2_id = -1;

/* Right-x2 */
      if(myM <(pD->NGrid[1])-1)
        pG->rx2_id = pD->GData[myN][myM+1][myL].ID_Comm_Domain;
      else pG->rx2_id = -1;

/* Left-x3 */
      if(myN > 0) pG->lx3_id = pD->GData[myN-1][myM][myL].ID_Comm_Domain;
      else pG->lx3_id = -1;

/* Right-x3 */
      if(myN <(pD->NGrid[2])-1)
        pG->rx3_id = pD->GData[myN+1][myM][myL].ID_Comm_Domain;
      else pG->rx3_id = -1;

    }
  }}

  return;

/*--- Error messages ---------------------------------------------------------*/

  on_error7:
    free_3d_array(pG->E3i);
  on_error6:
    free_3d_array(pG->E2i);
  on_error5:
    free_3d_array(pG->E1i);
  on_error4:
    free_3d_array(pG->B3i);
  on_error3:
    free_3d_array(pG->B2i);
  on_error2:
    free_3d_array(pG->B1i);
  on_error1:
    free_3d_array(pG->U);
    peg_error("[init_grid]: Error allocating memory\n");
}
