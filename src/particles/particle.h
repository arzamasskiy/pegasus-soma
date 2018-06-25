#ifndef PARTICLE_H
#define PARTICLE_H 
#include "../copyright.h"
/*============================================================================*/
/*! \file particle.h
 *  \brief Global variables for all functions in in the src/particles
 *   directory.								      */
/*============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../pegasus.h"
#include "../defs.h"

#include "../config.h"

/*----------------------------------------------------------------------------*/

/*---------------------------- Particle Properties ---------------------------*/
int npartypes;               /*!< number of particle types */
Grain_Property *grproperty;  /*!< array of particle properties of all types */

/*--------------------------- grid limit quantities --------------------------*/
/* left and right limit of grid indices */
int ilp,iup, jlp,jup, klp,kup;
/* left and right limit of grid boundary */
Real x1lpar, x1upar, x2lpar, x2upar, x3lpar, x3upar;

/*! \var int ncell
 *  \brief number of neighbouring cells involved in 1D interpolation */
int ncell;

/*-------------------------- Number of Filter Passes -------------------------*/
/*! \var int nfpass
 *  \brief number of passes through the filter */
int nfpass;

#endif /* PARTICLE_H */
