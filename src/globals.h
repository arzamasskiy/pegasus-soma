#ifndef GLOBALS_H
#define GLOBALS_H
/*============================================================================*/
/*! \file globals.h
 *  \brief Contains global variables.
 *
 * PURPOSE: Contains global variables:
 *   The first occurence in this file is included in main.c and defines the
 *   variables.  The second is included everywhere else.		      */
/*============================================================================*/

#ifdef MAIN_C

Real CourNo;                 /*!< Courant, Friedrichs, & Lewy (CFL) number */
#ifdef VARIABLE_DT
Real Safety;		     /*!< Safety factor for Variable dt cadence */
#endif
Real beta,beta_prp,beta_prl; /*!< ion plasma beta parameter: 2 n_i T_i / B^2 */
Real ZTeTi;                  /*!< Z T_e / T_i */
Real vinject;                /*!< Injection speed */
#if defined ADIABATIC
Real Gamma;                  /*!< adiabatic index (ratio of specific heats) */
Real Gamma_1, Gamma_2;       /*!< (Gamma)-1 and (Gamma)-2 */
#endif
int myID_Comm_world; /*!< Rank (proc ID) in MPI_COMM_WORLD, 0 for single proc */
Real d_MIN = TINY_NUMBER;    /*!< density floor */

#ifdef SHEARING_BOX
Real Omega_0, qshear;
Real Shear_0;
enum SS2DCoord ShBoxCoord;
Real eta_Ohm, eta_hyper;
#endif

WeightFun_t getweight = NULL;     /*!< get weight function */

#ifdef DELTA_F
GetDistFuncFun_t getdf = NULL;
SetBackgroundFun_t setbg = NULL;
#endif
InitDistFuncFun_t initdf = NULL;


/*----------------------------------------------------------------------------*/
/* definitions included everywhere except main.c  */

#else /* MAIN_C */

extern Real CourNo;
#ifdef VARIABLE_DT
extern Real Safety;
#endif
Real beta,beta_prp,beta_prl,ZTeTi;
Real vinject;
#if defined ADIABATIC
Real Gamma,Gamma_1, Gamma_2;
#endif
extern int myID_Comm_world;
extern Real d_MIN;

extern WeightFun_t getweight;

#ifdef DELTA_F
extern GetDistFuncFun_t getdf;
extern SetBackgroundFun_t setbg;
#endif
extern InitDistFuncFun_t initdf;

#ifdef SHEARING_BOX
extern Real Omega_0, qshear;
extern Real Shear_0;
extern enum SS2DCoord ShBoxCoord;
extern Real eta_Ohm, eta_hyper;
#endif


#endif /* MAIN_C */
#endif /* GLOBALS_H */
