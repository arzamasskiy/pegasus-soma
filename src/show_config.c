#include "copyright.h"
/*============================================================================*/
/*! \file show_config.c 
 *  \brief Outputs information on configuration of Pegasus.
 *
 * PURPOSE: Outputs information on configuration of Pegasus.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - show_config()     - prints diagnostic message showinf code configuration 
 * - show_config_par() - adds configuration information to database used by par
 *============================================================================*/

#include <stdio.h>
#include "defs.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void show_config(void)
 *  \brief The packages and features reported on by this function should
 *   be kept consistent with the optional packages and features added by the
 *   file configure.ac in the top-level directory   */

void show_config(void)
{
  int nscal;
  peg_pout(0,"\nConfiguration details:\n\n");
  peg_pout(0," Problem:                 %s\n",A_PROBLEM);

#if defined(ADIABATIC)
  peg_pout(0," Electron Eqn of State:   ADIABATIC\n");
#elif defined(ISOTHERMAL)
  peg_pout(0," Electron Eqn of State:   ISOTHERMAL\n");
#else
  peg_pout(0," Electron Eqn of State:   " EOS_STR "\n");
#endif


#if defined(SINGLE_PREC)
  peg_pout(0," Precision:               SINGLE_PREC\n");
#elif defined(DOUBLE_PREC)
  peg_pout(0," Precision:               DOUBLE_PREC\n");
#endif

#ifdef WRITE_GHOST_CELLS
  peg_pout(0," Ghost cell Output:       ON\n");
#else
  peg_pout(0," Ghost cell Output:       OFF\n");
#endif

#if defined(MPI_PARALLEL)
  peg_pout(0," Parallel Modes: MPI:     ON\n");
#else
  peg_pout(0," Parallel Modes: MPI:     OFF\n");
#endif

#ifdef DELTA_F
  peg_pout(0," Delta-f Method:          ON\n");
#else
  peg_pout(0," Delta-f Method:          OFF\n");
#endif


#ifdef VARIABLE_DT
  peg_pout(0," Variable dt:             ON\n");
#else
  peg_pout(0," Variable dt:             OFF\n");
#endif
}

/*----------------------------------------------------------------------------*/
/*! \fn void show_config_par(void)
 *  \brief Add the configure block to the parameter database used
 *    by the functions in par.c.  */

void show_config_par(void)
{
  par_sets("configure","problem",A_PROBLEM,"Name of the problem file");

#if defined(ADIABATIC)
  par_sets("configure","eq_state","adiabatic","electron eqn of state");
#elif defined(ISOTHERMAL)
  par_sets("configure","eq_state","isothermal","electron eqn of state");
#else
  par_sets("configure","eq_state",EOS_STR,"electron eqn of state");
#endif


#if defined(SINGLE_PREC)
  par_sets("configure","precision","single","Type of Real variables");
#elif defined(DOUBLE_PREC)
  par_sets("configure","precision","double","Type of Real variables");
#endif

#ifdef WRITE_GHOST_CELLS
  par_sets("configure","write_ghost","yes","Ghost cells included in output?");
#else
  par_sets("configure","write_ghost","no","Ghost cells included in output?");
#endif

#if defined(MPI_PARALLEL)
  par_sets("configure","mpi","yes","Is code MPI parallel enabled?");
#else
  par_sets("configure","mpi","no","Is code MPI parallel enabled?");
#endif

#ifdef DELTA_F
  par_sets("configure","Delta-f","yes","Delta-f method enabled?");
#else
  par_sets("configure","Delta-f","no","Delta-f method enabled?");
#endif

#ifdef VARIABLE_DT
  par_sets("configure","dt","variable","Variable or fixed dt?");
#else
  par_sets("configure","dt","fixed","Variable or fixed dt?");
#endif
  return;
}
