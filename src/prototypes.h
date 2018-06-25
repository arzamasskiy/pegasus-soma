#ifndef PROTOTYPES_H
#define PROTOTYPES_H 
#include "copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions from the /src directory.      */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "pegasus.h"
#include "defs.h"
#include "config.h"

#include "integrators/prototypes.h"
#include "particles/prototypes.h"

/*----------------------------------------------------------------------------*/
/* main.c */
int pegasus_main(int argc, char *argv[]);

/*----------------------------------------------------------------------------*/
/* peg_array.c */
void*   calloc_1d_array(                      size_t nc, size_t size);
void**  calloc_2d_array(           size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);
void free_1d_array(void *array);
void free_2d_array(void *array);
void free_3d_array(void *array);

/*----------------------------------------------------------------------------*/
/* peg_log.c */
void peg_log_set_level(const int out, const int err);
void peg_log_open(const char *basename, const int lazy, const char *mode);
void peg_log_close(void);
FILE *pegout_fp(void);
FILE *pegerr_fp(void);
void peg_flush_out(void);
void peg_flush_err(void);
int peg_perr(const int level, const char *fmt, ...);
int peg_pout(const int level, const char *fmt, ...);

/*----------------------------------------------------------------------------*/
/* peg_files.c */
char *peg_fname(const char *path, const char *basename,
                const char *levstr, const char *domstr,
                const int dlen, const int idump, 
                const char *id, const char *ext);

/*----------------------------------------------------------------------------*/
/* peg_signal.c */
void peg_sig_init(void);
int  peg_sig_act(int *piquit);

/*----------------------------------------------------------------------------*/
/* baton.c */
void baton_start(const int Nb, const int tag);
void baton_stop(const int Nb, const int tag);

/*----------------------------------------------------------------------------*/
/* bvals_mhd.c  */
void bvals_mhd_init(MeshS *pM);
void bvals_mhd_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc);
void bvals_mhd(DomainS *pDomain);

/*----------------------------------------------------------------------------*/
/* bvals_heat.c  */
void bvals_heat_init(MeshS *pM);
void bvals_heat_fun(DomainS *pD, enum BCDirection dir, HVGFun_t prob_bc);
void bvals_heat(DomainS *pD, Real3Vect ***Efiltered);

/*----------------------------------------------------------------------------*/
/* cc_pos.c */
void cc_pos(const GridS *pG, const int i, const int j,const int k,
            Real *px1, Real *px2, Real *px3);
void fc_pos(const GridS *pG, const int i, const int j,const int k,
            Real *px1, Real *px2, Real *px3);
int celli(const GridS *pG, const Real x, const Real dx1_1, int *i, Real *a);
Real x1cc(const GridS *pG, const int i);
int cellj(const GridS *pG, const Real y, const Real dx2_1, int *j, Real *b);
Real x2cc(const GridS *pG, const int j);
int cellk(const GridS *pG, const Real z, const Real dx3_1, int *k, Real *c);
Real x3cc(const GridS *pG, const int k);

/*----------------------------------------------------------------------------*/
/* init_grid.c */
void init_grid(MeshS *pM);

/*----------------------------------------------------------------------------*/
/* init_mesh.c */
void init_mesh(MeshS *pM);
void get_myGridIndex(DomainS *pD, const int my_id, int *pi, int *pj, int *pk);

/*----------------------------------------------------------------------------*/
/* peg_dt.c */
void peg_dt(MeshS *pM);

/*----------------------------------------------------------------------------*/
/* new_dt.c */
void new_dt(MeshS *pM, int dtncheck);

/*----------------------------------------------------------------------------*/
/* output.c - and related files */
void init_output(MeshS *pM);
void data_output(MeshS *pM, const int flag);
void add_rst_out(OutputS *new_out);
void data_output_destruct(void);
void dump_history_enroll(const ConsFun_t pfun, const char *label);
Real ***OutData3(GridS *pGrid, OutputS *pOut, int *Nx1, int *Nx2, int *Nx3);
Real  **OutData2(GridS *pGrid, OutputS *pOut, int *Nx1, int *Nx2);
Real   *OutData1(GridS *pGrid, OutputS *pOut, int *Nx1);

void output_pdf  (MeshS *pM, OutputS *pOut);
void output_pgm  (MeshS *pM, OutputS *pOut);
void output_ppm  (MeshS *pM, OutputS *pOut);
void output_vtk  (MeshS *pM, OutputS *pOut);
void output_tab  (MeshS *pM, OutputS *pOut);

void dump_binary (MeshS *pM, OutputS *pOut);
void dump_history(MeshS *pM, OutputS *pOut);
void dump_tab    (MeshS *pM, OutputS *pOut);
void dump_vtk    (MeshS *pM, OutputS *pOut);

/*----------------------------------------------------------------------------*/
/* par.c */
void   par_open(char *filename);
void   par_cmdline(int argc, char *argv[]);
int    par_exist(char *block, char *name);

char  *par_gets(char *block, char *name);
int    par_geti(char *block, char *name);
double par_getd(char *block, char *name);

char  *par_gets_def(char *block, char *name, char   *def);
int    par_geti_def(char *block, char *name, int    def);
double par_getd_def(char *block, char *name, double def);

void   par_sets(char *block, char *name, char *sval, char *comment);
void   par_seti(char *block, char *name, char *fmt, int ival, char *comment);
void   par_setd(char *block, char *name, char *fmt, double dval, char *comment);

void   par_dump(int mode, FILE *fp);
void   par_close(void);

#ifdef MPI_PARALLEL
void par_dist_mpi(const int mytid, MPI_Comm comm);
#endif

/*----------------------------------------------------------------------------*/
/* prob/PROBLEM.c ; linked to problem.c */
void problem(DomainS *pD);
void Userwork_in_loop(MeshS *pM);
void Userwork_after_loop(MeshS *pM);
void problem_read_restart(MeshS *pM, FILE *fp);
void problem_write_restart(MeshS *pM, FILE *fp);
ConsFun_t get_usr_expr(const char *expr);
VOutFun_t get_usr_out_fun(const char *name);

PropFun_t get_usr_par_prop(const char *name);
void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2,
                  const Real x3, const Real v1, const Real v2, const Real v3);

#ifdef DELTA_F
void getdf_user(GrainS *gr, Real *df);
void setbg_user(GridS *pG);
#endif
void initdf_user(long int *idum, Real **vp, const long int Np, 
		 const long int NpTot);

/*----------------------------------------------------------------------------*/
/* restart.c  */
void dump_restart(MeshS *pM, OutputS *pout);
void restart_grids(char *res_file, MeshS *pM);

/*----------------------------------------------------------------------------*/
/* show_config.c */
void show_config(void);
void show_config_par(void);

/*----------------------------------------------------------------------------*/
/* utils.c */
char *peg_strdup(const char *in);
int peg_gcd(int a, int b);
int peg_big_endian(void);
void peg_bswap(void *vdat, int sizeof_len, int cnt);
void peg_error(char *fmt, ...);
void minmax1(Real   *data, int nx1,                   Real *dmin, Real *dmax);
void minmax2(Real  **data, int nx2, int nx1,          Real *dmin, Real *dmax);
void minmax3(Real ***data, int nx3, int nx2, int nx1, Real *dmin, Real *dmax);
void do_nothing_bc(GridS *pG);
Real compute_div_b(GridS *pG);
int sign_change(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *a, Real *b);
int bisection(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *root);
Real trapzd(Real (*func)(Real), const Real a, const Real b, const int n, const Real s);
Real qsimp(Real (*func)(Real), const Real a, const Real b);
Real avg1d(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real avg2d(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real avg3d(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real avgXZ(Real (*func)(Real, Real, Real), const GridS *pG, const int i, const int j, const int k);
Real vecpot2b1i(Real (*A2)(Real,Real,Real), Real (*A3)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k);
Real vecpot2b2i(Real (*A1)(Real,Real,Real), Real (*A3)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k);
Real vecpot2b3i(Real (*A1)(Real,Real,Real), Real (*A2)(Real,Real,Real),
                const GridS *pG, const int i, const int j, const int k);
void InverseMatrix(Real **a, int n, Real **b);
void MatrixMult(Real **a, Real **b, int m, int n, int l, Real **c);
#endif /* PROTOTYPES_H */
