#include "../copyright.h"
/*============================================================================*/
/*! \file exchange.c
 *  \brief Exchange the gas-particle coupling array between ghost and
 *         boundary zones.
 *
 * PURPOSE: Particles near grid boundaries deposit their physical properties
 *   partially to the ghost zones. This part of the deposit is to be mapped
 *   to the grid zone. The procedure is opposite to setting boundary conditions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - exchange_gpcouple()
 * - exchange_gpcouple_init()
 * - exchange_gpcouple_fun()
 * - exchange_gpcouple_destruct()
 *
 * PRIVATE FUNCTIONS:
 * - reflecting_???
 * - outflow_???
 * - periodic_???
 * - send_???
 * - receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
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

/* Define the maximum number of variables in the gas-particle coupling array
 * to be exchanged at any step (0-3) */
#define NVar_Max 10
/* Define the maximum number of zones to be exchanged (NExc = 2)
 *      + the number of ghost zones to be filled (nghost = 5) */
#define NLayer_Max 7

/*! \struct GPExc
 *  \brief Define structure which holds variables for gas-particle exchange */
typedef struct GPExc_s{
  Real U[NVar_Max];
}GPExc;

/* number of variables to be exchanged for each cell in a specific step */
static int NVar;
/* number of ghost layers to be exchanged */
static int NExc;
/* number of grid layers (offset from boundary) to be copied
 * For pure exchange, set to 0;
 * To fill (copy) particle deposits to ghost zones, set to nghost
 * Total number of layers is thus NLayer = NExc + NOfst.  */
static int NOfst;

/* grid index limit for the exchange */
static int il,iu, jl,ju, kl,ku;
static int ib,it, jb,jt, kb,kt;

/* Temporary array where the exchange operation is executed */
static GPExc ***myCoup=NULL;

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

/*====================== PROTOTYPE OF PRIVATE FUNCTIONS ======================*/
/*----------------------------------------------------------------------------*/

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - apply reflecting BCs at boundary ???
 *   outflow_???()  - apply outflow BCs at boundary ???
 *   periodic_???() - apply periodic BCs at boundary ???
 *   pack_???()     - pack data at ??? boundary
 *   unpack_???()   - unpack data at ??? boundary
 *============================================================================*/

static void reflect_ix1_exchange(GridS *pG);
static void reflect_ox1_exchange(GridS *pG);
static void reflect_ix2_exchange(GridS *pG);
static void reflect_ox2_exchange(GridS *pG);
static void reflect_ix3_exchange(GridS *pG);
static void reflect_ox3_exchange(GridS *pG);

static void outflow_ix1_exchange(GridS *pG);
static void outflow_ox1_exchange(GridS *pG);
static void outflow_ix2_exchange(GridS *pG);
static void outflow_ox2_exchange(GridS *pG);
static void outflow_ix3_exchange(GridS *pG);
static void outflow_ox3_exchange(GridS *pG);

static void inflow_ix1_exchange(GridS *pG);
static void inflow_ox1_exchange(GridS *pG);
static void inflow_ix2_exchange(GridS *pG);
static void inflow_ox2_exchange(GridS *pG);
static void inflow_ix3_exchange(GridS *pG);
static void inflow_ox3_exchange(GridS *pG);

static void periodic_ix1_exchange(GridS *pG);
static void periodic_ox1_exchange(GridS *pG);
static void periodic_ix2_exchange(GridS *pG);
static void periodic_ox2_exchange(GridS *pG);
static void periodic_ix3_exchange(GridS *pG);
static void periodic_ox3_exchange(GridS *pG);


#ifdef MPI_PARALLEL
static void pack_ix1_exchange(GridS *pG);
static void pack_ox1_exchange(GridS *pG);
static void pack_ix2_exchange(GridS *pG);
static void pack_ox2_exchange(GridS *pG);
static void pack_ix3_exchange(GridS *pG);
static void pack_ox3_exchange(GridS *pG);

static void unpack_ix1_exchange(GridS *pG);
static void unpack_ox1_exchange(GridS *pG);
static void unpack_ix2_exchange(GridS *pG);
static void unpack_ox2_exchange(GridS *pG);
static void unpack_ix3_exchange(GridS *pG);
static void unpack_ox3_exchange(GridS *pG);

#endif /* MPI_PARALLEL */


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void exchange_feedback(DomainS *pD, short step)
 *  \brief Calls appropriate functions to copy feedback in the ghost
 *    zones back to the grid.
 *
 *    The function pointers (*pD->???_EBCFun) are set during
 *    initialization by exchange_feedback_init() to be either a user-defined
 *    function, or one of the functions corresponding to reflecting, periodic,
 *    or outflow.  If the left- or right-Grid ID numbers are >= 1 (neighboring
 *    grids exist), then MPI calls are used.
 *
 *  Order for updating boundary conditions must always be x3-x2-x1 in order to
 *  fill the corner cells properly (opposite to setting MHD B.C.!)
 */

void exchange_gpcouple(DomainS *pD, short lab)
{
  GridS *pG = pD->Grid;
  int i,j,k;

#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, ierr, mIndex;
#endif /* MPI_PARALLEL */
	
  /*--- Step 1. ------------------------------------------------------------------	
   * Copy the information in the Gas-Particle coupling array into temporary array
   * This step depends on the parameter "lab", where for
   * lab = 0: particle binning for output purpose
   * lab = 1: predictor step of feedback exchange
   * lab = 2: corrector step of feedback exchange
   * All the operations in this routine are performed on the temporary array,
   * which will be copied back to the main array GPCoup at the end.
   *----------------------------------------------------------------------------*/
  
  switch (lab) {
    case 0: /* particle binning for output purpose */
		  
      NVar = 10; NExc = 1; NOfst = nfpass;
      
      for (k=klp; k<=kup; k++) {
        for (j=jlp; j<=jup; j++) {
          for (i=ilp; i<=iup; i++) {
            myCoup[k][j][i].U[0]=pG->Coup[k][j][i].grid_M1;
            myCoup[k][j][i].U[1]=pG->Coup[k][j][i].grid_M2;
            myCoup[k][j][i].U[2]=pG->Coup[k][j][i].grid_M3;
            myCoup[k][j][i].U[3]=pG->Coup[k][j][i].grid_d;
            myCoup[k][j][i].U[4]=pG->Coup[k][j][i].grid_p11;
            myCoup[k][j][i].U[5]=pG->Coup[k][j][i].grid_p12;
            myCoup[k][j][i].U[6]=pG->Coup[k][j][i].grid_p13;
            myCoup[k][j][i].U[7]=pG->Coup[k][j][i].grid_p22;
            myCoup[k][j][i].U[8]=pG->Coup[k][j][i].grid_p23;
            myCoup[k][j][i].U[9]=pG->Coup[k][j][i].grid_p33;
          }}}
      break;
		  
    case 1: /* for predictor step of hybrid exchange */
      /* 5 filled ghost layers needed (if no filtering)
       or nfpass ghost layers needed (if filtering) */
      
      NVar = 4; NExc = 1; NOfst = nfpass > 0 ? nfpass : 5;
      
      for (k=klp; k<=kup; k++) {
        for (j=jlp; j<=jup; j++) {
          for (i=ilp; i<=iup; i++) {
            myCoup[k][j][i].U[0]=pG->Coup[k][j][i].grid_M1;
            myCoup[k][j][i].U[1]=pG->Coup[k][j][i].grid_M2;
            myCoup[k][j][i].U[2]=pG->Coup[k][j][i].grid_M3;
            myCoup[k][j][i].U[3]=pG->Coup[k][j][i].grid_d;
          }}}
      break;
      
    case 2: /* for corrector step of hybrid exchange */
      /* 3+nfpass filled ghost layers needed or nfpass 
       ghost layers needed (if nfpass > 2) */
      
      NVar = 4; NExc = 2; NOfst = nfpass > 2 ? nfpass : nfpass+3;

      for (k=klp; k<=kup; k++) {
        for (j=jlp; j<=jup; j++) {
          for (i=ilp; i<=iup; i++) {
            myCoup[k][j][i].U[0]=pG->Coup[k][j][i].grid_M1;
            myCoup[k][j][i].U[1]=pG->Coup[k][j][i].grid_M2;
            myCoup[k][j][i].U[2]=pG->Coup[k][j][i].grid_M3;
            myCoup[k][j][i].U[3]=pG->Coup[k][j][i].grid_d;
          }}}
      break;
      
    case 3: /* re-set ghost zones after filter */
      /* only called if predictor step */
      
      NVar = 4; NExc = 0; NOfst = 5;
      
      for (k=klp; k<=kup; k++) {
        for (j=jlp; j<=jup; j++) {
          for (i=ilp; i<=iup; i++) {
            myCoup[k][j][i].U[0]=pG->Coup[k][j][i].grid_M1;
            myCoup[k][j][i].U[1]=pG->Coup[k][j][i].grid_M2;
            myCoup[k][j][i].U[2]=pG->Coup[k][j][i].grid_M3;
            myCoup[k][j][i].U[3]=pG->Coup[k][j][i].grid_d;
          }}}
      break;
      
    case 4: /* re-set ghost zones after filter */
      /* only called if corrector step and nfpass > 2 */
      
      NVar = 4; NExc = 0; NOfst = 3;
      
      for (k=klp; k<=kup; k++) {
        for (j=jlp; j<=jup; j++) {
          for (i=ilp; i<=iup; i++) {
            myCoup[k][j][i].U[0]=pG->Coup[k][j][i].grid_M1;
            myCoup[k][j][i].U[1]=pG->Coup[k][j][i].grid_M2;
            myCoup[k][j][i].U[2]=pG->Coup[k][j][i].grid_M3;
            myCoup[k][j][i].U[3]=pG->Coup[k][j][i].grid_d;
          }}}
      break;
      
    default:
      peg_perr(-1,"[exchange_GPCouple]: lab must be equal to 0, 1, 2, 3, or 4!\n");
      
  }
  
  /* set left and right grid indices */
  if (pG->Nx[0] > 1) {
    il = pG->is - NExc;         iu = pG->ie + NExc;
    ib = pG->is - NOfst;        it = pG->ie + NOfst;
  } else {
    il = ib = pG->is;           iu = it = pG->ie;
  }   
  
  if (pG->Nx[1] > 1) {
    jl = pG->js - NExc;         ju = pG->je + NExc;
    jb = pG->js - NOfst;        jt = pG->je + NOfst;
  } else {
    jl = jb = pG->js;           ju = jt = pG->je;
  }   
  
  if (pG->Nx[2] > 1) {
    kl = pG->ks - NExc;         ku = pG->ke + NExc;
    kb = pG->ks - NOfst;        kt = pG->ke + NOfst;
  } else {
    kl = kb = pG->ks;           ku = kt = pG->ke;
  }
  
  /*--- Step 2. ------------------------------------------------------------------
   * Feedback exchange in x3-direction */
  
  if (pG->Nx[2] > 1){
    
#ifdef MPI_PARALLEL
    cnt1 = pG->Nx[0] > 1 ? pG->Nx[0] + 2*NExc : 1;
    cnt2 = pG->Nx[1] > 1 ? pG->Nx[1] + 2*NExc : 1;
    cnt = (NExc+NOfst)*cnt1*cnt2*NVar;
    
    /* MPI blocks to both left and right */
    if (pG->rx3_id >= 0 && pG->lx3_id >= 0) {
      
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,LtoR_tag,
                       pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,RtoL_tag,
                       pD->Comm_Domain, &(recv_rq[1]));
      
      /* pack and send data L and R */
      pack_ix3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,RtoL_tag,
                       pD->Comm_Domain, &(send_rq[0]));
      
      pack_ox3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,LtoR_tag,
                       pD->Comm_Domain, &(send_rq[1]));
      
      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
      
      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_exchange(pG);
      if (mIndex == 1) unpack_ox3_exchange(pG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_exchange(pG);
      if (mIndex == 1) unpack_ox3_exchange(pG);
      
    }
    
    /* Physical boundary on left, MPI block on right */
    if (pG->rx3_id >= 0 && pG->lx3_id < 0) {
      
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,RtoL_tag,
                       pD->Comm_Domain, &(recv_rq[1]));
      
      /* pack and send data R */
      pack_ox3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,LtoR_tag,
                       pD->Comm_Domain, &(send_rq[1]));
      
      /* set physical boundary */
      (*(pD->ix3_EBCFun))(pG);
      
      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3_exchange(pG);
      
    }
    
    /* MPI block on left, Physical boundary on right */
    if (pG->rx3_id < 0 && pG->lx3_id >= 0) {
      
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,LtoR_tag,
                       pD->Comm_Domain, &(recv_rq[0]));
      
      /* pack and send data L */
      pack_ix3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,RtoL_tag,
                       pD->Comm_Domain, &(send_rq[0]));
      
      /* set physical boundary */
      (*(pD->ox3_EBCFun))(pG);
      
      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3_exchange(pG);
      
    }
#endif /* MPI_PARALLEL */
    
    /* Physical boundaries on both left and right */
    if (pG->rx3_id < 0 && pG->lx3_id < 0) {
      (*(pD->ix3_EBCFun))(pG);
      (*(pD->ox3_EBCFun))(pG);
    }
    
  }
  
  /*--- Step 3. ------------------------------------------------------------------
   * Feedback exchange in x2-direction */
  
  if (pG->Nx[1] > 1){
    
#ifdef MPI_PARALLEL
    cnt1 = pG->Nx[0] > 1 ? pG->Nx[0] + 2*NExc : 1;
    cnt3 = pG->Nx[2] > 1 ? pG->Nx[2] + 2*NOfst: 1;
    cnt = (NExc+NOfst)*cnt1*cnt3*NVar;
    
    /* MPI blocks to both left and right */
    if (pG->rx2_id >= 0 && pG->lx2_id >= 0) {
      
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
                       pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
                       pD->Comm_Domain, &(recv_rq[1]));
      
      /* pack and send data L and R */
      pack_ix2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
                       pD->Comm_Domain, &(send_rq[0]));
      
      pack_ox2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
                       pD->Comm_Domain, &(send_rq[1]));
      
      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
      
      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_exchange(pG);
      if (mIndex == 1) unpack_ox2_exchange(pG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_exchange(pG);
      if (mIndex == 1) unpack_ox2_exchange(pG);
      
    }
    
    /* Physical boundary on left, MPI block on right */
    if (pG->rx2_id >= 0 && pG->lx2_id < 0) {
      
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
                       pD->Comm_Domain, &(recv_rq[1]));
      
      /* pack and send data R */
      pack_ox2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
                       pD->Comm_Domain, &(send_rq[1]));
      
      /* set physical boundary */
      (*(pD->ix2_EBCFun))(pG);
      
      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2_exchange(pG);
      
    }
    
    /* MPI block on left, Physical boundary on right */
    if (pG->rx2_id < 0 && pG->lx2_id >= 0) {
      
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
                       pD->Comm_Domain, &(recv_rq[0]));
      
      /* pack and send data L */
      pack_ix2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
                       pD->Comm_Domain, &(send_rq[0]));
      
      /* set physical boundary */
      (*(pD->ox2_EBCFun))(pG);
      
      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2_exchange(pG);
      
    }
#endif /* MPI_PARALLEL */
    
    /* Physical boundaries on both left and right */
    if (pG->rx2_id < 0 && pG->lx2_id < 0) {
      (*(pD->ix2_EBCFun))(pG);
      (*(pD->ox2_EBCFun))(pG);
    }
    
   }
  
  /*--- Step 4. ------------------------------------------------------------------
   * Feedback exchange in x1-direction */
  
  if (pG->Nx[0] > 1){
    
#ifdef MPI_PARALLEL
    cnt2 = pG->Nx[1] > 1 ? pG->Nx[1] + 2*NOfst : 1;
    cnt3 = pG->Nx[2] > 1 ? pG->Nx[2] + 2*NOfst : 1;
    cnt = (NExc+NOfst)*cnt2*cnt3*NVar;
    
    /* MPI blocks to both left and right */
    if (pG->rx1_id >= 0 && pG->lx1_id >= 0) {
      
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,LtoR_tag,
                       pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,RtoL_tag,
                       pD->Comm_Domain, &(recv_rq[1]));
      
      /* pack and send data L and R */
      pack_ix1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,RtoL_tag,
                       pD->Comm_Domain, &(send_rq[0]));
      
      pack_ox1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,LtoR_tag,
                       pD->Comm_Domain, &(send_rq[1]));
      
      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);
      
      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_exchange(pG);
      if (mIndex == 1) unpack_ox1_exchange(pG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_exchange(pG);
      if (mIndex == 1) unpack_ox1_exchange(pG);
    }
    
    /* Physical boundary on left, MPI block on right */
    if (pG->rx1_id >= 0 && pG->lx1_id < 0) {
      
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,RtoL_tag,
                       pD->Comm_Domain, &(recv_rq[1]));
      
      /* pack and send data R */
      pack_ox1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,LtoR_tag,
                       pD->Comm_Domain, &(send_rq[1]));
      
      /* set physical boundary */
      (*(pD->ix1_EBCFun))(pG);
      
      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1_exchange(pG);
      
    }
    
    /* MPI block on left, Physical boundary on right */
    if (pG->rx1_id < 0 && pG->lx1_id >= 0) {
      
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,LtoR_tag,
                       pD->Comm_Domain, &(recv_rq[0]));
      
      /* pack and send data L */
      pack_ix1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,RtoL_tag,
                       pD->Comm_Domain, &(send_rq[0]));
      
      /* set physical boundary */
      (*(pD->ox1_EBCFun))(pG);
      
      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1_exchange(pG);
      
    }
#endif /* MPI_PARALLEL */
    
    /* Physical boundaries on both left and right */
    if (pG->rx1_id < 0 && pG->lx1_id < 0) {
      (*(pD->ix1_EBCFun))(pG);
      (*(pD->ox1_EBCFun))(pG);
    } 
  }
  
  /*--- Step 5. ------------------------------------------------------------------	
   * Copy the variables from the temporary array where exchange has finished back
   * to the Gas-Particle coupling array. Again, for
   * lab = 0: particle binning for output purpose
   * lab = 1: predictor step of feedback exchange
   * lab = 2: corrector step of feedback exchange
   *----------------------------------------------------------------------------*/
	
  switch (lab) {
    case 0:	/* particle binning for output purpose */
      for (k=kb; k<=kt; k++) {
        for (j=jb; j<=jt; j++) {
          for (i=ib; i<=it; i++) {
            pG->Coup[k][j][i].grid_M1 = myCoup[k][j][i].U[0];
            pG->Coup[k][j][i].grid_M2 = myCoup[k][j][i].U[1];
            pG->Coup[k][j][i].grid_M3 = myCoup[k][j][i].U[2];
            pG->Coup[k][j][i].grid_d  = myCoup[k][j][i].U[3];
            pG->Coup[k][j][i].grid_p11= myCoup[k][j][i].U[4];
            pG->Coup[k][j][i].grid_p12= myCoup[k][j][i].U[5];
            pG->Coup[k][j][i].grid_p13= myCoup[k][j][i].U[6];
            pG->Coup[k][j][i].grid_p22= myCoup[k][j][i].U[7];
            pG->Coup[k][j][i].grid_p23= myCoup[k][j][i].U[8];
            pG->Coup[k][j][i].grid_p33= myCoup[k][j][i].U[9];
          }}}
      break;
			
    case 1: /* predictor step of hybrid exchange */
      for (k=kb; k<=kt; k++) {
        for (j=jb; j<=jt; j++) {
          for (i=ib; i<=it; i++) {
            pG->Coup[k][j][i].grid_M1= myCoup[k][j][i].U[0];
            pG->Coup[k][j][i].grid_M2= myCoup[k][j][i].U[1];
            pG->Coup[k][j][i].grid_M3= myCoup[k][j][i].U[2];
            pG->Coup[k][j][i].grid_d = myCoup[k][j][i].U[3];
          }}}
      break;
      
    case 2: /* corrector step of hybrid exchange */
      for (k=kb; k<=kt; k++) {
        for (j=jb; j<=jt; j++) {
          for (i=ib; i<=it; i++) {
            pG->Coup[k][j][i].grid_M1= myCoup[k][j][i].U[0];
            pG->Coup[k][j][i].grid_M2= myCoup[k][j][i].U[1];
            pG->Coup[k][j][i].grid_M3= myCoup[k][j][i].U[2];
            pG->Coup[k][j][i].grid_d = myCoup[k][j][i].U[3];
          }}}
      break;
      
    case 3: /* predictor step after filter */
      for (k=kb; k<=kt; k++) {
        for (j=jb; j<=jt; j++) {
          for (i=ib; i<=it; i++) {
            pG->Coup[k][j][i].grid_M1= myCoup[k][j][i].U[0];
            pG->Coup[k][j][i].grid_M2= myCoup[k][j][i].U[1];
            pG->Coup[k][j][i].grid_M3= myCoup[k][j][i].U[2];
            pG->Coup[k][j][i].grid_d = myCoup[k][j][i].U[3];
          }}}
      break;
      
    case 4: /* corrector step after filter */
      for (k=kb; k<=kt; k++) {
        for (j=jb; j<=jt; j++) {
          for (i=ib; i<=it; i++) {
            pG->Coup[k][j][i].grid_M1= myCoup[k][j][i].U[0];
            pG->Coup[k][j][i].grid_M2= myCoup[k][j][i].U[1];
            pG->Coup[k][j][i].grid_M3= myCoup[k][j][i].U[2];
            pG->Coup[k][j][i].grid_d = myCoup[k][j][i].U[3];
          }}}
      break;
      
    default:
      peg_perr(-1,"[exchange_GPCouple]: lab must be equal to 0, 1, 2, 3, or 4!\n");
      
  }
	
  return;
  
}

/*----------------------------------------------------------------------------*/
/*! \fn void exchange_gpcouple_init(MeshS *pM)
 *  \brief Sets function pointers for the exchange of gas-particle coupling
 *   array, allocates memory, etc. for the initialization
 */
void exchange_gpcouple_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int ibc_x1, obc_x1; /* x1 inner and outer boundary condition flag */
  int ibc_x2, obc_x2; /* x2 inner and outer boundary condition flag */
  int ibc_x3, obc_x3; /* x3 inner and outer boundary condition flag */
  int N1T,N2T,N3T;
#ifdef MPI_PARALLEL
  int i,j,k;
  int x1cnt, x2cnt, x3cnt; /* Number of Gas passed in x1-, x2-, x3-dir. */
  int nx1t, nx2t, nx3t, size;
#endif /* MPI_PARALLEL */
  
  if (pM->NLevels > 1)
    peg_error("[exchange_init]: particle module does not suport SMR\n");
  
  pD = &(pM->Domain[0][0]);
  pG = pD->Grid;
  
  /* Set function pointers for physical boundaries in x1-direction */
  
  if(pG->Nx[0] > 1) {
    if(pD->ix1_EBCFun == NULL){  /* EBCFun ptr was not set in prob gen */
      
      switch(pM->BCFlag_ix1){
          
        case 1: /* Reflecting */
          pD->ix1_EBCFun = reflect_ix1_exchange;
          break;
          
        case 5: /* Reflecting */
          pD->ix1_EBCFun = reflect_ix1_exchange;
          break;
          
        case 2: /* Outflow */
          pD->ix1_EBCFun = outflow_ix1_exchange;
          break;
          
        case 4: /* Periodic */
          pD->ix1_EBCFun = periodic_ix1_exchange;
          break;
          
        case 6: /* Inflow */
          pD->ix1_EBCFun = inflow_ix1_exchange;
          break;
          
        default:
          peg_perr(-1,"[exchange_init]: bc_ix1 = %d unknown\n",
                   pM->BCFlag_ix1);
          exit(EXIT_FAILURE);
      }
    }
    
    if(pD->ox1_EBCFun == NULL){
      
      switch(pM->BCFlag_ox1){
          
        case 1: /* Reflecting */
          pD->ox1_EBCFun = reflect_ox1_exchange;
          break;
          
        case 5: /* Reflecting */
          pD->ox1_EBCFun = reflect_ox1_exchange;
          break;
          
        case 2: /* Outflow */
          pD->ox1_EBCFun = outflow_ox1_exchange;
          break;
          
        case 4: /* Periodic */
          pD->ox1_EBCFun = periodic_ox1_exchange;
          break;
          
        case 6: /* Inflow */
          pD->ox1_EBCFun = inflow_ox1_exchange;
          break;
          
        default:
          peg_perr(-1,"[exchange_init]: bc_ox1 = %d unknown\n",
                   pM->BCFlag_ox1);
          exit(EXIT_FAILURE);
      }
    }
  }
  
  /* Set function pointers for physical boundaries in x2-direction */
  
  if(pG->Nx[1] > 1) {
    if(pD->ix2_EBCFun == NULL){
      
      switch(pM->BCFlag_ix2){
          
        case 1: /* Reflecting */
          pD->ix2_EBCFun = reflect_ix2_exchange;
          break;
          
        case 5: /* Reflecting */
          pD->ix2_EBCFun = reflect_ix2_exchange;
          break;
          
        case 2: /* Outflow */
          pD->ix2_EBCFun = outflow_ix2_exchange;
          break;
          
        case 4: /* Periodic */
          pD->ix2_EBCFun = periodic_ix2_exchange;
          break;
          
        case 6: /* Inflow */
          pD->ix3_EBCFun = inflow_ix2_exchange;
          break;
          
        default:
          peg_perr(-1,"[exchange_init]: bc_ix2 = %d unknown\n",
                   pM->BCFlag_ix2);
          exit(EXIT_FAILURE);
      }
    }
    
    if(pD->ox2_EBCFun == NULL){
      
      switch(pM->BCFlag_ox2){
          
        case 1: /* Reflecting */
          pD->ox2_EBCFun = reflect_ox2_exchange;
          break;
          
        case 5: /* Reflecting */
          pD->ox2_EBCFun = reflect_ox2_exchange;
          break;
          
        case 2: /* Outflow */
          pD->ox2_EBCFun = outflow_ox2_exchange;
          break;
          
        case 4: /* Periodic */
          pD->ox2_EBCFun = periodic_ox2_exchange;
          break;
          
        case 6: /* Inflow */
          pD->ox2_EBCFun = inflow_ox2_exchange;
          break;
          
        default:
          peg_perr(-1,"[exchange_init]: bc_ox2 = %d unknown\n",
                   pM->BCFlag_ox2);
          exit(EXIT_FAILURE);
      }
    }
  }
  
  /* Set function pointers for physical boundaries in x3-direction */
  
  if(pG->Nx[2] > 1) {
    if(pD->ix3_EBCFun == NULL){
      
      switch(pM->BCFlag_ix3){
          
        case 1: /* Reflecting */
          pD->ix3_EBCFun = reflect_ix3_exchange;
          break;
          
        case 5: /* Reflecting */
          pD->ix3_EBCFun = reflect_ix3_exchange;
          break;
          
        case 2: /* Outflow */
          pD->ix3_EBCFun = outflow_ix3_exchange;
          break;
          
        case 4: /* Periodic */
          pD->ix3_EBCFun = periodic_ix3_exchange;
          break;
          
        case 6: /* Inflow */
          pD->ix3_EBCFun = inflow_ix3_exchange;
          break;
          
        default:
          peg_perr(-1,"[exchange_init]: bc_ix3 = %d unknown\n",
                   pM->BCFlag_ix3);
          exit(EXIT_FAILURE);
      }
    }
    
    if(pD->ox3_EBCFun == NULL){
      
      switch(pM->BCFlag_ox3){
          
        case 1: /* Reflecting */
          pD->ox3_EBCFun = reflect_ox3_exchange;
          break;
          
        case 5: /* Reflecting */
          pD->ox3_EBCFun = reflect_ox3_exchange;
          break;
          
        case 2: /* Outflow */
          pD->ox3_EBCFun = outflow_ox3_exchange;
          break;
          
        case 4: /* Periodic */
          pD->ox3_EBCFun = periodic_ox3_exchange;
          break;
          
        case 6: /* Inflow */
          pD->ox3_EBCFun = inflow_ox3_exchange;
          break;
          
        default:
          peg_perr(-1,"[exchange_init]: bc_ox3 = %d unknown\n",
                   pM->BCFlag_ox3);
          exit(EXIT_FAILURE);
      }
    }
  }
	
  if (pG->Nx[0] > 1) {
    N1T = pG->Nx[0]+2*nghost;
  } else {
    N1T=1;
  }
  if (pG->Nx[1] > 1) {
    N2T = pG->Nx[1]+2*nghost;
  } else {
    N2T=1;
  }
  if (pG->Nx[2] > 1) {
    N3T = pG->Nx[2]+2*nghost;
  } else {
    N3T=1;
  }
  
  if ((myCoup = (GPExc***)calloc_3d_array(N3T,N2T,N1T, sizeof(GPExc))) == NULL)
    peg_error("[exchange_init]: Failed to allocate the myCoup array.\n");
	
#ifdef MPI_PARALLEL
  x1cnt = x2cnt = x3cnt = 0;
  
  for (k=0; k<(pD->NGrid[2]); k++){
    for (j=0; j<(pD->NGrid[1]); j++){
      for (i=0; i<(pD->NGrid[0]); i++){
        if(pD->NGrid[2] > 1){
          nx1t = pD->GData[k][j][i].Nx[0];
          if(nx1t > 1) nx1t += 2*NLayer_Max;
          
          nx2t = pD->GData[k][j][i].Nx[1];
          if(nx2t > 1) nx2t += 2*NLayer_Max;
          
          x3cnt = nx1t*nx2t > x3cnt ? nx1t*nx2t : x3cnt;
        }
        
        if(pD->NGrid[1] > 1){
          nx1t = pD->GData[k][j][i].Nx[0];
          if(nx1t > 1) nx1t += 2*NLayer_Max;
          
          nx3t = pD->GData[k][j][i].Nx[2];
          if(nx3t > 1) nx3t += 2*NLayer_Max;
          
          x2cnt = nx1t*nx3t > x2cnt ? nx1t*nx3t : x2cnt;
        }
        
        if(pD->NGrid[0] > 1){
          nx2t = pD->GData[k][j][i].Nx[1];
          if(nx2t > 1) nx2t += 2*NLayer_Max;
          
          nx3t = pD->GData[k][j][i].Nx[2];
          if(nx3t > 1) nx3t += 2*NLayer_Max;
          
          x1cnt = nx2t*nx3t > x1cnt ? nx2t*nx3t : x1cnt;
        }
        
      }
    }
  }
  
  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;
  
  size *= NLayer_Max*NVar_Max; /* Multiply by the third dimension */
  
  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      peg_error("[exchange_init]: Failed to allocate send buffer\n");
    
    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      peg_error("[exchange_init]: Failed to allocate recv buffer\n");
  }
  
  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    peg_error("[exchange_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    peg_error("[exchange_init]: Failed to allocate send MPI_Request array\n");
  
#endif /* MPI_PARALLEL */
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void exchange_gpcouple_fun(enum Direction dir, VGFun_t prob_bc)
 *  \brief Sets function pointers for user-defined exchange gas-particle
 *   coupling in the problem generator
 */

void exchange_gpcouple_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
    case left_x1:
      pD->ix1_EBCFun = prob_bc;
      break;
    case right_x1:
      pD->ox1_EBCFun = prob_bc;
      break;
    case left_x2:
      pD->ix2_EBCFun = prob_bc;
      break;
    case right_x2:
      pD->ox2_EBCFun = prob_bc;
      break;
    case left_x3:
      pD->ix3_EBCFun = prob_bc;
      break;
    case right_x3:
      pD->ox3_EBCFun = prob_bc;
      break;
    default:
      peg_perr(-1,"[bvals_exchange_fun]: Unknown direction = %d\n",dir);
      exit(EXIT_FAILURE);
  }
  return;
}

/*! \fn void exchange_gpcouple_destruct(GridS *pG, Domain *pD)
 *  \brief Finalize the exchange of gas-particle coupling */
void exchange_gpcouple_destruct(MeshS *pM)
{
  free_3d_array(myCoup);

#ifdef MPI_PARALLEL
  if(send_buf != NULL) free_2d_array(send_buf);
  if(recv_buf != NULL) free_2d_array(recv_buf);
#endif
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???
 *   outflow_???
 *   periodic_???
 *   send_???
 *   receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1,5)
 */

static void reflect_ix3_exchange(GridS *pG)
{
  int kr;
  int i,j,k,n;
  
  for (k=kl; k<pG->ks; k++) {
    kr = 2*pG->ks-k-1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[kr][j][i].U[2] -= myCoup[k][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
  
  for (k=kb; k<pG->ks; k++) {
    kr = 2*pG->ks-k-1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[k][j][i].U[2] = -myCoup[kr][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
      } 
    }   
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1,5)
 */

static void reflect_ox3_exchange(GridS *pG)
{
  int kr;
  int i,j,k,n;
  
  for (k=pG->ke+1; k<=ku; k++) {
    kr = 2*pG->ke-k+1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[kr][j][i].U[2] -= myCoup[k][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
  
  for (k=pG->ke+1; k<=kt; k++) {
    kr = 2*pG->ke-k+1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[k][j][i].U[2] = -myCoup[kr][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1,5)
 */
static void reflect_ix2_exchange(GridS *pG)
{
  int jr;
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=jl; j<pG->js; j++) {
      jr = 2*pG->js-j-1;
      for (i=il; i<=iu; i++) {
        myCoup[k][jr][i].U[0] += myCoup[k][j][i].U[0];
        myCoup[k][jr][i].U[1] -= myCoup[k][j][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][jr][i].U[n] += myCoup[k][j][i].U[n];
      }
    } 
  } 
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<pG->js; j++) {
      jr = 2*pG->js-j-1;
      for (i=il; i<=iu; i++) { 
        myCoup[k][j][i].U[0] = myCoup[k][jr][i].U[0];
        myCoup[k][j][i].U[1] = -myCoup[k][jr][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][jr][i].U[n];
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1,5)
 */

static void reflect_ox2_exchange(GridS *pG)
{
  int jr;
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=ju; j++) {
      jr = 2*pG->je-j+1;
      for (i=il; i<=iu; i++) {
        myCoup[k][jr][i].U[0] += myCoup[k][j][i].U[0];
        myCoup[k][jr][i].U[1] -= myCoup[k][j][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][jr][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=jt; j++) {
      jr = 2*pG->je-j+1;
      for (i=il; i<=iu; i++) {
        myCoup[k][j][i].U[0] = myCoup[k][jr][i].U[0];
        myCoup[k][j][i].U[1] = -myCoup[k][jr][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][jr][i].U[n];
      }
    } 
  } 
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix1_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1,5)
 */

static void reflect_ix1_exchange(GridS *pG)
{
  int ir;
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=il; i<pG->is; i++) {
        ir = 2*pG->is-i-1;
        myCoup[k][j][ir].U[0] -= myCoup[k][j][i].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][ir].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=ib; i<pG->is; i++) {
        ir = 2*pG->is-i-1;
        myCoup[k][j][i].U[0] = -myCoup[k][j][ir].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][j][ir].U[n];
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox1_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1,5)
 */

static void reflect_ox1_exchange(GridS *pG)
{
  int ir;
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) { 
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=iu; i++) {
        ir = 2*pG->ie-i+1;
        myCoup[k][j][ir].U[0] -= myCoup[k][j][i].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][ir].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
  
  for (k=kb; k<=kt; k++) { 
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=it; i++) {
        ir = 2*pG->ie-i+1;
        myCoup[k][j][i].U[0] = -myCoup[k][j][ir].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][j][ir].U[n];
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix3_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Inner x3 boundary (ibc=2)
 */

static void outflow_ix3_exchange(GridS *pG)
{
  int kr;
  int i,j,k,n;
  
  /*
  for (k=kl; k<pG->ks; k++) {
    kr = 2*pG->ks-k-1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
   */
  
  for (k=kb; k<pG->ks; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[pG->ks][j][i].U[n];
      } 
    }   
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox3_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Outer x3 boundary (obc_x3=2)
 */

static void outflow_ox3_exchange(GridS *pG)
{
  int kr;
  int i,j,k,n;
  
  /*
  for (k=pG->ke+1; k<=ku; k++) {
    kr = 2*pG->ke-k+1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
   */
  
  for (k=pG->ke+1; k<=kt; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[pG->ke][j][i].U[n];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix2_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Inner x2 boundary (ibc_x2=2)
 */
static void outflow_ix2_exchange(GridS *pG)
{
  int jr;
  int i,j,k,n;
  
  /*
  for (k=kb; k<=kt; k++) {
    for (j=jl; j<pG->js; j++) {
      jr = 2*pG->js-j-1;
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][jr][i].U[n] += myCoup[k][j][i].U[n];
      }
    } 
  } 
   */
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<pG->js; j++) {
      for (i=il; i<=iu; i++) { 
        for (n=0;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][pG->js][i].U[n];
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox2_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Outer x2 boundary (obc_x2=2)
 */

static void outflow_ox2_exchange(GridS *pG)
{
  int jr;
  int i,j,k,n;
  
  /*
  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=ju; j++) {
      jr = 2*pG->je-j+1;
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][jr][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
   */
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=jt; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][pG->je][i].U[n];
      }
    } 
  } 
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix1_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Inner x1 boundary (ibc_x1=2)
 */

static void outflow_ix1_exchange(GridS *pG)
{
  int ir;
  int i,j,k,n;
  
  /*
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=il; i<pG->is; i++) {
        ir = 2*pG->is-i-1;
        for (n=0;n<NVar; n++)
          myCoup[k][j][ir].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
   */
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=ib; i<pG->is; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][j][pG->is].U[n];
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox1_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions, Outer x1 boundary (obc_x1=2)
 */

static void outflow_ox1_exchange(GridS *pG)
{
  int ir;
  int i,j,k,n;
  
  /*
  for (k=kb; k<=kt; k++) { 
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=iu; i++) {
        ir = 2*pG->ie-i+1;
        for (n=0;n<NVar; n++)
          myCoup[k][j][ir].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
   */
  
  for (k=kb; k<=kt; k++) { 
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=it; i++) {
        for (n=0;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][j][pG->ie].U[n];
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void inflow_ix3_exchange(GridS *pG)
 *  \brief INFLOW boundary conditions, Inner x3 boundary (ibc_x3=6)
 */

static void inflow_ix3_exchange(GridS *pG)
{
  int i,j,k,n;
  
  for (k=kb; k<pG->ks; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=ib; i<=it; i++) {
        myCoup[k][j][i].U[0] = vinject;
        myCoup[k][j][i].U[1] = 0.0;
        myCoup[k][j][i].U[2] = 0.0;
        myCoup[k][j][i].U[3] = 1.0;
        if (NVar > 4) {
          myCoup[k][j][i].U[4] = 0.5*beta;
          myCoup[k][j][i].U[5] = 0.0;
          myCoup[k][j][i].U[6] = 0.0;
          myCoup[k][j][i].U[7] = 0.5*beta;
          myCoup[k][j][i].U[8] = 0.0;
          myCoup[k][j][i].U[9] = 0.5*beta;
        }
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void inflow_ox3_exchange(GridS *pG)
 *  \brief INFLOW boundary conditions, Outer x3 boundary (obc_x3=6)
 */

static void inflow_ox3_exchange(GridS *pG)
{
  int i,j,k,n;
  
  for (k=pG->ke+1; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=ib; i<=it; i++) {
        myCoup[k][j][i].U[0] = -vinject;
        myCoup[k][j][i].U[1] =  0.0;
        myCoup[k][j][i].U[2] =  0.0;
        myCoup[k][j][i].U[3] =  1.0;
        if (NVar > 4) {
          myCoup[k][j][i].U[4] = 0.5*beta;
          myCoup[k][j][i].U[5] = 0.0;
          myCoup[k][j][i].U[6] = 0.0;
          myCoup[k][j][i].U[7] = 0.5*beta;
          myCoup[k][j][i].U[8] = 0.0;
          myCoup[k][j][i].U[9] = 0.5*beta;
        }
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void inflow_ix2_exchange(GridS *pG)
 *  \brief INFLOW boundary conditions, Inner x2 boundary (ibc_x2=6)
 */

static void inflow_ix2_exchange(GridS *pG)
{
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<pG->js; j++) {
      for (i=ib; i<=it; i++) {
        myCoup[k][j][i].U[0] = vinject;
        myCoup[k][j][i].U[1] = 0.0;
        myCoup[k][j][i].U[2] = 0.0;
        myCoup[k][j][i].U[3] = 1.0;
        if (NVar > 4) {
          myCoup[k][j][i].U[4] = 0.5*beta;
          myCoup[k][j][i].U[5] = 0.0;
          myCoup[k][j][i].U[6] = 0.0;
          myCoup[k][j][i].U[7] = 0.5*beta;
          myCoup[k][j][i].U[8] = 0.0;
          myCoup[k][j][i].U[9] = 0.5*beta;
        }
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void inflow_ox2_exchange(GridS *pG)
 *  \brief INFLOW boundary conditions, Outer x2 boundary (obc_x2=6)
 */

static void inflow_ox2_exchange(GridS *pG)
{
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=jt; j++) {
      for (i=ib; i<=it; i++) {
        myCoup[k][j][i].U[0] = -vinject;
        myCoup[k][j][i].U[1] =  0.0;
        myCoup[k][j][i].U[2] =  0.0;
        myCoup[k][j][i].U[3] =  1.0;
        if (NVar > 4) {
          myCoup[k][j][i].U[4] = 0.5*beta;
          myCoup[k][j][i].U[5] = 0.0;
          myCoup[k][j][i].U[6] = 0.0;
          myCoup[k][j][i].U[7] = 0.5*beta;
          myCoup[k][j][i].U[8] = 0.0;
          myCoup[k][j][i].U[9] = 0.5*beta;
        }
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void inflow_ix1_exchange(GridS *pG)
 *  \brief INFLOW boundary conditions, Inner x1 boundary (ibc_x1=6)
 */

static void inflow_ix1_exchange(GridS *pG)
{
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=ib; i<pG->is; i++) {
        myCoup[k][j][i].U[0] = vinject;
        myCoup[k][j][i].U[1] = 0.0;
        myCoup[k][j][i].U[2] = 0.0;
        myCoup[k][j][i].U[3] = 1.0;
        if (NVar > 4) {
          myCoup[k][j][i].U[4] = 0.5*beta;
          myCoup[k][j][i].U[5] = 0.0;
          myCoup[k][j][i].U[6] = 0.0;
          myCoup[k][j][i].U[7] = 0.5*beta;
          myCoup[k][j][i].U[8] = 0.0;
          myCoup[k][j][i].U[9] = 0.5*beta;
        }
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void inflow_ox1_exchange(GridS *pG)
 *  \brief INFLOW boundary conditions, Outer x1 boundary (obc_x1=6)
 */

static void inflow_ox1_exchange(GridS *pG)
{
  int i,j,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=it; i++) {
        myCoup[k][j][i].U[0] = -vinject;
        myCoup[k][j][i].U[1] =  0.0;
        myCoup[k][j][i].U[2] =  0.0;
        myCoup[k][j][i].U[3] =  1.0;
        if (NVar > 4) {
          myCoup[k][j][i].U[4] = 0.5*beta;
          myCoup[k][j][i].U[5] = 0.0;
          myCoup[k][j][i].U[6] = 0.0;
          myCoup[k][j][i].U[7] = 0.5*beta;
          myCoup[k][j][i].U[8] = 0.0;
          myCoup[k][j][i].U[9] = 0.5*beta;
        }
      }
    }
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix3_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x3 boundary (ibc_x3=4)
 */
static void periodic_ix3_exchange(GridS *pG)
{
  int dk = pG->Nx[2];
  int i,j,k,k1,n;
  
  for (k=pG->ks; k<pG->ks+NExc; k++) {
    k1 = k+dk-NExc;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          myCoup[k ][j][i].U[n] += myCoup[k +dk][j][i].U[n];
          myCoup[k1][j][i].U[n] += myCoup[k1-dk][j][i].U[n];
        }}}}
  
  for (k=pG->ks-NOfst; k<pG->ks; k++) {
    k1 = k+dk+NOfst;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          myCoup[k ][j][i].U[n] = myCoup[k +dk][j][i].U[n];
          myCoup[k1][j][i].U[n] = myCoup[k1-dk][j][i].U[n];
        }}}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox3_exchange(GridS *pG) 
 *  \brief PERIODIC boundary conditions, Outer x3 boundary (obc_x3=4)
 */

static void periodic_ox3_exchange(GridS *pG)
{
  /* Do nothing. All are handled by periodic_ix3_exchange */
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix2_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x2 boundary (ibc_x2=4)
 */

static void periodic_ix2_exchange(GridS *pG)
{
  int dj = pG->Nx[1];
  int i,j,j1,k,n;
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->js; j<pG->js+NExc; j++) {
      j1 = j+dj-NExc;
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          myCoup[k][j ][i].U[n] += myCoup[k][j +dj][i].U[n];
          myCoup[k][j1][i].U[n] += myCoup[k][j1-dj][i].U[n];
        }}}}
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->js-NOfst; j<pG->js; j++) {
      j1 = j+dj+NOfst;
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          myCoup[k][j ][i].U[n] = myCoup[k][j +dj][i].U[n];
          myCoup[k][j1][i].U[n] = myCoup[k][j1-dj][i].U[n];
        }}}   }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox2_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions,Outer x2 boundary (obc_x2=4)
 */

static void periodic_ox2_exchange(GridS *pG)
{
  /* Do nothing. All are handled by periodic_ix2_exchange */
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix1_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_ix1_exchange(GridS *pG)
{
  int di = pG->Nx[0];
  int i,i1,j,k,n;
  
    for (k=kb; k<=kt; k++) {
      for (j=jb; j<=jt; j++) {
        for (i=pG->is; i<pG->is+NExc; i++) {
          i1 = i+di-NExc;
          for (n=0;n<NVar; n++) {
            myCoup[k][j][i ].U[n] += myCoup[k][j][i +di].U[n];
            myCoup[k][j][i1].U[n] += myCoup[k][j][i1-di].U[n];
          }}}}
    
    for (k=kb; k<=kt; k++) {
      for (j=jb; j<=jt; j++) {
        for (i=pG->is-NOfst; i<pG->is; i++) {
          i1 = i+di+NOfst;
          for (n=0;n<NVar; n++) {
            myCoup[k][j][i ].U[n] = myCoup[k][j][i +di].U[n];
            myCoup[k][j][i1].U[n] = myCoup[k][j][i1-di].U[n];
          }}}}
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox1_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Outer x1 boundary (obc_x1=4)
 */

static void periodic_ox1_exchange(GridS *pG)
{
  /* Do nothing. All are handled by periodic_ix1_exchange */
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 16 funs; ~400 lines */

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix3_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at inner x3 boundary
 */

static void pack_ix3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[0][0]);
  
  for (k=kl; k<=pG->ks+NOfst-1; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          *(pd++) = myCoup[k][j][i].U[n];
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox3_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x3 boundary
 */

static void pack_ox3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[1][0]);
  
  for (k=pG->ke-NOfst+1; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          *(pd++) = myCoup[k][j][i].U[n];
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix2_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at inner x2 boundary
 */

static void pack_ix2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[0][0]);
  
  for (k=kb; k<=kt; k++) {
    for (j=jl; j<=pG->js+NOfst-1; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          *(pd++) = myCoup[k][j][i].U[n];
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox2_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x2 boundary
 */

static void pack_ox2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[1][0]);
  
  for (k=kb; k<=kt; k++) {
    for (j=pG->je-NOfst+1; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        for (n=0;n<NVar; n++) {
          *(pd++) = myCoup[k][j][i].U[n];
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix1_exchange(GridS *pG)
 *  \brief pack the coupling array to bufer at inner x1 boundary
 */

static void pack_ix1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(send_buf[0][0]);

    for (k=kb; k<=kt; k++) {
      for (j=jb; j<=jt; j++) {
        for (i=il; i<=pG->is+NOfst-1; i++) {
          for (n=0;n<NVar; n++) {
            *(pd++) = myCoup[k][j][i].U[n];
          }}}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox1_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x1 boundary
 */

static void pack_ox1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(send_buf[1][0]);
  
    for (k=kb; k<=kt; k++) {
      for (j=jb; j<=jt; j++) {
        for (i=pG->ie-NOfst+1; i<=iu; i++) {
          for (n=0;n<NVar; n++) {
            *(pd++) = myCoup[k][j][i].U[n];
          }}}}
  
  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix3_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at inner x3 boundary
 */

static void unpack_ix3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[0][0]);
  
  for (k=pG->ks-NOfst; k<pG->ks+NExc; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        for (n=0;n<NVar; n++) {
          myCoup[k][j][i].U[n] += *(pd++);
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox3_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at outer x3 boundary
 */

static void unpack_ox3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[1][0]);
  
  for (k=pG->ke-NExc+1; k<=pG->ke+NOfst; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        for (n=0;n<NVar; n++) {
          myCoup[k][j][i].U[n] += *(pd++);
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix2_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at inner x2 boundary
 */

static void unpack_ix2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[0][0]);
  
  for (k=kb; k<=kt; k++){
    for (j=pG->js-NOfst; j<pG->js+NExc; j++){
      for (i=il; i<=iu; i++){
        for (n=0;n<NVar; n++) {
          myCoup[k][j][i].U[n] += *(pd++);
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox2_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at outer x2 boundary
 */

static void unpack_ox2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[1][0]);
  
  for (k=kb; k<=kt; k++){
    for (j=pG->je-NExc+1; j<=pG->je+NOfst; j++){
      for (i=il; i<=iu; i++){
        for (n=0;n<NVar; n++) {
          myCoup[k][j][i].U[n] += *(pd++);
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix1_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at inner x1 boundary
 */

static void unpack_ix1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[0][0]);
  
  for (k=kb; k<=kt; k++){
    for (j=jb; j<=jt; j++){
      for (i=pG->is-NOfst; i<pG->is+NExc; i++){
        for (n=0;n<NVar; n++) {
          myCoup[k][j][i].U[n] += *(pd++);
        }
      }}}
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox1_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at outer x1 boundary
 */

static void unpack_ox1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[1][0]);
  
  for (k=kb; k<=kt; k++){
    for (j=jb; j<=jt; j++){
      for (i=pG->ie-NExc+1; i<=pG->ie+NOfst; i++){
        for (n=0;n<NVar; n++) {
          myCoup[k][j][i].U[n] += *(pd++);
        }
      }}}
  
  return;
}

#endif /* MPI_PARALLEL */
