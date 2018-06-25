#ifndef PARTICLES_PROTOTYPES_H
#define PARTICLES_PROTOTYPES_H 
#include "../copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions in the /src/particles dir      */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../pegasus.h"
#include "../defs.h"
#include "../config.h"

/* bvals_particle.c */
void bvals_particle(DomainS *pD);
#ifdef FARGO
void advect_particles(DomainS *pD);
#endif
void bvals_particle_init(MeshS *pM);
void bvals_particle_fun(enum BCDirection dir, VGFun_t prob_bc);
void bvals_final_particle(MeshS *pM);
void bvals_particle_destruct(MeshS *pM);

/* exchange.c */
void exchange_gpcouple(DomainS *pD, short lab);
void exchange_gpcouple_init(MeshS *pM);
void exchange_gpcouple_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc);
void exchange_gpcouple_destruct(MeshS *pM);

/* init_particle.c */
void init_particle(MeshS *pM);
void particle_destruct(MeshS *pM);
void particle_realloc(GridS *pG, long n);

/* injectors.c */
void inject_particles(DomainS *pD);
void inject_init(MeshS *pM);

/* integrators_particle.c */
void Integrate_Particles_1(DomainS *pD);
void Integrate_Particles_2(DomainS *pD);
void int_par_boris(GridS *pG, GrainS *curG, Real3Vect cell1,
                   Real *dv1, Real *dv2, Real *dv3);

/* output_particle.c */
void particle_to_grid(DomainS *pD, PropFun_t par_prop);
void dump_particle_track(MeshS *pM, OutputS *pOut);
int  property_all(const GrainS *gr, const GrainAux *grsub);
void dump_particle_history(MeshS *pM, OutputS *pOut);
void dump_particle_spectrum(MeshS *pM, OutputS *pOut);

/* utils_particle.c */
void getwei_linear(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);
void getwei_TSC   (GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);
void getwei_QP    (GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);

int getvalues(GridS *pG, Real weight[3][3][3], int is, int js, int ks
              ,Real *bfld1, Real *bfld2, Real *bfld3
              ,Real *efld1, Real *efld2, Real *efld3
              );

void shuffle(GridS *pG);
void Delete_Ghost(GridS *pG);
void deposit(DomainS *pD, short lab);
void set_gpcouple(GridS *pG, short lab);
void smooth_gpcouple(DomainS *pD, short lab);

/* distfunc.c */
#ifdef DELTA_F
void getdf_maxwell(GrainS *gr, Real *df);
void getdf_bimaxwell(GrainS *gr, Real *df);
void setbg_maxwell(GridS *pG);
void setbg_bimaxwell(GridS *pG);
#endif /* DELTA_F */
void initdf_static(long int *idum, Real **vp, const long int Np, const long int NpTot);
void initdf_maxwell(long int *idum, Real **vp, const long int Np, const long int NpTot);
void initdf_bimaxwell(long int *idum, Real **vp, const long int Np, const long int NpTot);

#endif /* PARTICLES_PROTOTYPES_H */
