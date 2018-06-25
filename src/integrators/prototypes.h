#ifndef INTEGRATORS_PROTOTYPES_H
#define INTEGRATORS_PROTOTYPES_H 
#include "../copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions in the /src/integrators dir */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../pegasus.h"
#include "../defs.h"

#include "../config.h"

/* integrate.c */
VDFun_t integrate_init(MeshS *pM);
void integrate_destruct(void);

/* integrate_1d_ctu.c and integrate_1d_vl.c and integrate_1d_hybrid.c */
void integrate_destruct_1d(void);
void integrate_init_1d(MeshS *pM);
void integrate_1d_hybrid(DomainS *pD);

/* integrate_2d_ctu.c and integrate_2d_vl.c and integrate_2d_hybrid.c */
void integrate_destruct_2d(void);
void integrate_init_2d(MeshS *pM);
void integrate_2d_hybrid(DomainS *pD);

/* integrate_3d_ctu.c and integrate_3d_vl.c and integrate_3d_hybrid.c */
void integrate_destruct_3d(void);
void integrate_init_3d(MeshS *pM);
void integrate_3d_hybrid(DomainS *pD);
#ifdef SHEARING_BOX
void integrate_2d_hybrid_xz(DomainS *pD);
#endif
#endif /* INTEGRATORS_PROTOTYPES_H */
