#ifndef PEGASUS_H
#define PEGASUS_H 
/*============================================================================*/
/*! \file pegasus.h
 *  \brief Contains definitions of many data types and structures.
 *
 * PURPOSE: Contains definitions of the following data types and structures:
 * - Real    - either float or double, depending on configure option
 * - ConsS   - cell-centered variables
 * - GrainS  - basic properties of particles
 * - GridS   - everything in a single Grid
 * - DomainS - everything in a single Domain (potentially many Grids)
 * - MeshS   - everything across whole Mesh (potentially many Domains)
 * - OutputS - everything associated with an individual output		      */
/*============================================================================*/
#include "defs.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

/*! \typedef Real
 *  \brief Variable precision float, depends on macro set by configure.
 */ 
#if defined(SINGLE_PREC)
#ifdef MPI_PARALLEL
#error: MPI requires double precision
#endif /*MPI_PARALLEL */
typedef float  Real;
#elif defined(DOUBLE_PREC)
typedef double Real;
#else
# error "Not a valid precision flag"
#endif

/*! \struct Real3Vect
 *  \brief General 3-vectors of Reals.
 */
typedef struct Real3Vect_s{
  Real x1, x2, x3;
}Real3Vect;
/*! \struct Int3Vect
 *  \brief General 3-vectors of ints.
 */
typedef struct Int3Vect_s{
  int i, j, k;
}Int3Vect;

/*! \struct SideS
 *  \brief Sides of a cube, used to find overlaps between Grids 
 *   at different levels.
 */
typedef struct Side_s{
  int ijkl[3];    /*!< indices of left-sides  in each dir [0,1,2]=[i,j,k] */ 
  int ijkr[3];    /*!< indices of right-sides in each dir [0,1,2]=[i,j,k] */ 
}SideS;

/*! \struct GridsDataS
 *  \brief Number of zones in, and identifying information about, a Grid.
 */
typedef struct GridsData_s{
  int Nx[3];       /*!< number of zones in each dir [0,1,2]=[x1,x2,x3] */
  int Disp[3];   /*!< i,j,k displacements from origin of root [0,1,2]=[i,j,k] */
  int ID_Comm_world;      /*!< ID of process for this Grid in MPI_COMM_WORLD */
  int ID_Comm_Domain;     /*!< ID of process for this Grid in Comm_Domain    */
}GridsDataS;

/*----------------------------------------------------------------------------*/
/*! \struct ConsS
 *  \brief Field Variables.
 *  IMPORTANT!! The order of the elements in ConsS CANNOT be changed.
 */
typedef struct Cons_s{
  Real B1c;			/*!< cell centered magnetic field in 1-dir*/
  Real B2c;			/*!< cell centered magnetic field in 2-dir*/
  Real B3c;			/*!< cell centered magnetic field in 3-dir*/
  Real E1c;                     /*!< cell centered electric field in 1-dir*/
  Real E2c;                     /*!< cell centered electric field in 2-dir*/
  Real E3c;                     /*!< cell centered electric field in 3-dir*/
}ConsS;

typedef struct Soln_s{
  Real d;
  Real M1;
  Real M2;
  Real M3;
  Real B1c;
  Real B2c;
  Real B3c;
}SolnS;

/*----------------------------------------------------------------------------*/
/*! \struct GrainS
 *  \brief Basic quantities for one pseudo-particle.
 */

/* Physical quantities of a dust particle */
typedef struct Grain_s{
  Real x1,x2,x3;	/*!< coordinate in X,Y,Z */
  Real v1,v2,v3;	/*!< velocity in X,Y,Z */
#ifdef DELTA_F
  Real f_0;        /*!< full dist. func. evaluated at particle's (x(0),v(0)) */
#endif
  int property;		/*!< index of grain properties */
  short pos;      /*!< position: 0: ghost; 1: grid; >=10: cross out/in; */
  long my_id;     /*!< particle id */
#ifdef MPI_PARALLEL
  int init_id;    /*!< particle's initial host processor id */
#endif
}GrainS;

/*! \struct GrainAux
 *  \brief Auxilary quantities for a dust particle. */
typedef struct GrainAux_s{
  Real dpar;            /*!< local particle density */

}GrainAux;

/*! \struct Grain_Property
 *  \brief List of physical grain properties. */
typedef struct Grain_Property_s{
  Real m;		  /*!< mass of this type of particle */
  Real qomc;  /*!< charge-to-mass ratio of this type of particle */
  long num;		/*!< number of particles with this property */
  short integrator;	/*!< integrator type: semi Boris (1) */
}Grain_Property;

/*! \struct GPCouple
 *  \brief Grid elements for gas-particle coupling. */
typedef struct GPCouple_s{
  Real grid_d;          /*!< gas density */
  Real grid_M1;         /*!< gas 1-velocity */
  Real grid_M2;         /*!< gas 2-velocity */
  Real grid_M3;         /*!< gas 3-velocity */
  Real grid_p11;        /*!< 11-component of pressure tensor */
  Real grid_p12;        /*!< 12-component of pressure tensor */
  Real grid_p13;        /*!< 13-component of pressure tensor */
  Real grid_p22;        /*!< 22-component of pressure tensor */
  Real grid_p23;        /*!< 23-component of pressure tensor */
  Real grid_p33;        /*!< 33-component of pressure tensor */
  Real grid_b1;         /*!< 1-component of magnetic field */
  Real grid_b2;         /*!< 2-component of magnetic field */
  Real grid_b3;         /*!< 3-component of magnetic field */
  Real grid_e1;         /*!< 1-component of electric field */
  Real grid_e2;         /*!< 2-component of electric field */
  Real grid_e3;         /*!< 3-component of electric field */
  Real grid_eb;         /*!< electric field . magnetic field */
}GPCouple;

/*----------------------------------------------------------------------------*/
/*! \struct GridS
 *  \brief 3D arrays of dependent variables, plus grid data, plus particle data,
 *   plus data about child and parent Grids, plus MPI rank information for a
 *   Grid.
 *
 *   Remember a Grid is defined as the region of a Domain at some
 *   refinement level being updated by a single processor.  Uses an array of
 *   ConsS, rather than arrays of each variable, to increase locality of data
 *   for a given cell in memory.  */

typedef struct Grid_s{
  ConsS ***U;                   /*!< field variables */
  Real ***B1i,***B2i,***B3i;    /*!< interface magnetic fields */
  Real ***E1i,***E2i,***E3i;    /*!< cell-centered electric fields */

  Real MinX[3];       /*!< min(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
  Real MaxX[3];       /*!< max(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
  Real dx1,dx2,dx3;   /*!< cell size on this Grid */
  Real time, dt;           /*!< current time and timestep  */
  int nstep;               /*!< number of integration steps taken */
  int is,ie;		   /*!< start/end cell index in x1 direction */
  int js,je;		   /*!< start/end cell index in x2 direction */
  int ks,ke;		   /*!< start/end cell index in x3 direction */
  int Nx[3];     /*!< # of zones in each dir on Grid [0,1,2]=[x1,x2,x3] */
  int Disp[3];   /*!< i,j,k displacements of Grid from origin [0,1,2]=[i,j,k] */

  int rx1_id, lx1_id;  /*!< ID of Grid to R/L in x1-dir (default=-1; no Grid) */
  int rx2_id, lx2_id;  /*!< ID of Grid to R/L in x2-dir (default=-1; no Grid) */
  int rx3_id, lx3_id;  /*!< ID of Grid to R/L in x3-dir (default=-1; no Grid) */

  long nparticle;            /*!< number of particles */
  long arrsize;              /*!< size of the particle array */
  GrainS *particle;          /*!< array of all particles */
  GrainAux *parsub;          /*!< supplemental particle information */
  GPCouple ***Coup;          /*!< array of gas-particle coupling */
#ifdef DELTA_F
  GPCouple ***Bkgrd;         /*!< array of background quantities: moments of
                              * background distribution function f0 at t = 0) */
#endif
}GridS;

/*! \fn void (*VGFun_t)(GridS *pG)
 *  \brief Generic void function of Grid. */
typedef void (*VGFun_t)(GridS *pG);    /* generic void function of Grid */
typedef void (*HVGFun_t)(GridS *pG, Real3Vect ***Efiltered);    /* generic void function of Grid */

/*----------------------------------------------------------------------------*/
/*! \struct DomainS
 *  \brief Information about one region of Mesh at some particular level.
 *
 * Contains pointer to a single Grid, even though the Domain may contain many
 * Grids, because for any general parallelization mode, no more than one Grid
 * can exist per Domain per processor.
 *
 * The i,j,k displacements are measured in units of grid cells on this Domain
 */
typedef struct Domain_s{
  Real RootMinX[3]; /*!< min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real RootMaxX[3]; /*!< max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real MinX[3];     /*!< min(x) in each dir on this Domain [0,1,2]=[x1,x2,x3] */
  Real MaxX[3];     /*!< max(x) in each dir on this Domain [0,1,2]=[x1,x2,x3] */
  Real dx[3];                /*!< cell size in this Domain [0,1,2]=[x1,x2,x3] */
  int Nx[3];    /*!< # of zones in each dir in this Domain [0,1,2]=[x1,x2,x3] */
  int NGrid[3]; /*!< # of Grids in each dir in this Domain [0,1,2]=[x1,x2,x3] */
  int Disp[3]; /*!< i,j,k displacements of Domain from origin [0,1,2]=[i,j,k] */
  int Level,DomNumber;   /*!< level and ID number of this Domain */
  int InputBlock;      /*!< # of <domain> block in input file for this Domain */
  GridS *Grid;     /*!< pointer to Grid in this Dom updated on this processor */

  GridsDataS ***GData;/*!< size,location,& processor IDs of Grids in this Dom */

  VGFun_t ix1_BCFun, ox1_BCFun;/*!< ix1/ox1 BC function pointers for this Dom */
  VGFun_t ix2_BCFun, ox2_BCFun;/*!< ix2/ox2 BC function pointers for this Dom */
  VGFun_t ix3_BCFun, ox3_BCFun;/*!< ix3/ox3 BC function pointers for this Dom */
  
  VGFun_t ix1_EBCFun, ox1_EBCFun;/*!< ix1/ox1 BC function pointers for particle exchange */
  VGFun_t ix2_EBCFun, ox2_EBCFun;/*!< ix2/ox2 BC function pointers for particle exchange */
  VGFun_t ix3_EBCFun, ox3_EBCFun;/*!< ix3/ox3 BC function pointers for particle exchange */

  HVGFun_t ix1_HBCFun, ox1_HBCFun;/*!< ix1/ox1 BC function pointers for particle exchange */
  HVGFun_t ix2_HBCFun, ox2_HBCFun;/*!< ix2/ox2 BC function pointers for particle exchange */
  HVGFun_t ix3_HBCFun, ox3_HBCFun;/*!< ix3/ox3 BC function pointers for particle exchange */
  
  
#ifdef MPI_PARALLEL
  MPI_Comm Comm_Domain;      /*!< MPI communicator between Grids on this Dom */
  MPI_Group Group_Domain;    /*!< MPI group for Domain communicator */
#endif /* MPI_PARALLEL */
}DomainS;

/*! \fn void (*VDFun_t)(DomainS *pD)
 *  \brief Generic void function of Domain. */
typedef void (*VDFun_t)(DomainS *pD);  /* generic void function of Domain */

/*----------------------------------------------------------------------------*/
/*! \struct MeshS
 *  \brief Information about entire mesh hierarchy, including array of Domains.
 */

typedef struct Mesh_s{
  Real RootMinX[3]; /*!< min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real RootMaxX[3]; /*!< max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real dx[3];     /*!< cell size on root Domain [0,1,2]=[x1,x2,x3] */
  Real time, dt;  /*!< current time and timestep for entire Mesh */
  int Nx[3];    /*!< # of zones in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  int nstep;                 /*!< number of integration steps taken */
  int BCFlag_ix1, BCFlag_ox1;  /*!< BC flag on root domain for inner/outer x1 */
  int BCFlag_ix2, BCFlag_ox2;  /*!< BC flag on root domain for inner/outer x2 */
  int BCFlag_ix3, BCFlag_ox3;  /*!< BC flag on root domain for inner/outer x3 */
  int NLevels;               /*!< overall number of refinement levels in mesh */
  int *DomainsPerLevel;      /*!< number of Domains per level (DPL) */
  DomainS **Domain;        /*!< array of Domains, indexed over levels and DPL */
  char *outfilename;         /*!< basename for output files containing -id#  */
}MeshS;

/*----------------------------------------------------------------------------*/
/* OutputS: Everything for outputs. */

struct Output_s;
/*! \fn void (*VOutFun_t)(MeshS *pM, struct Output_s *pout)
 *  \brief Output function pointer. */
typedef void (*VOutFun_t)(MeshS *pM, struct Output_s *pout);
/*! \fn void (*VResFun_t)(MeshS *pM, struct Output_s *pout)
 *  \brief Restart function pointer. */
typedef void (*VResFun_t)(MeshS *pM, struct Output_s *pout);
/*! \fn Real (*ConsFun_t)(const GridS *pG, const int i,const int j,const int k) 
 *  \brief Pointer to expression that computes quant for output.*/
typedef Real (*ConsFun_t)(const GridS *pG, const int i,const int j,const int k);
/*! \fn int (*PropFun_t)(const GrainS *gr, const GrainAux *grsub)
 *  \brief Particle property selection function */
typedef int (*PropFun_t)(const GrainS *gr, const GrainAux *grsub);
typedef Real (*Parfun_t)(const GridS *pG, const GrainS *gr);

/*! \struct OutputS
 *  \brief Everything for outputs. */
typedef struct Output_s{
  int n;          /*!< the N from the <outputN> block of this output */
  Real dt;        /*!< time interval between outputs  */
  Real t;         /*!< next time to output */
  int dcycle;     /*!< cycle interval between outputs as an alternative to time */
  int cycle;      /*!< next cycle to output */
  int num;        /*!< dump number (0=first) */
  char *out;      /*!< variable (or user fun) to be output */
  char *id;       /*!< filename is of the form <basename>[.idump][.id].<ext> */
  int out_pargrid;    /*!< bin particles to grid (=1) or not (=0) */
  PropFun_t par_prop; /*!< particle property selection function */

/* level and domain number of output (default = [-1,-1] = output all levels) */

  int nlevel, ndomain;

/* variables which describe data min/max */
  Real dmin,dmax;   /*!< user defined min/max for scaling data */
  Real gmin,gmax;   /*!< computed global min/max (over all output data) */
  int sdmin,sdmax;  /*!< 0 = auto scale, otherwise use dmin/dmax */

/* variables which describe coordinates of output data volume */
  int ndim;       /*!< 3=cube 2=slice 1=vector 0=scalar */
  int reduce_x1;  /*!< flag to denote reduction in x1 (0=no reduction) */
  int reduce_x2;  /*!< flag to denote reduction in x2 (0=no reduction) */
  int reduce_x3;  /*!< flag to denote reduction in x3 (0=no reduction) */
  Real x1l, x1u;  /*!< lower/upper x1 range for data slice  */
  Real x2l, x2u;  /*!< lower/upper x2 range for data slice  */
  Real x3l, x3u;  /*!< lower/upper x3 range for data slice  */

/* variables which describe output format */
  char *out_fmt;  /*!< output format = {bin, tab, hdf, hst, pgm, ppm, ...} */
  char *dat_fmt;  /*!< format string for tabular type output, e.g. "%10.5e" */
  char *palette;  /*!< name of palette for RGB conversions */
  float *rgb;     /*!< array of RGB[256*3] values derived from palette */
  float *der;     /*!< helper array of derivatives for interpolation into RGB */

/* pointers to output functions; data expressions */
  VOutFun_t out_fun; /*!< output function pointer */
  VResFun_t res_fun; /*!< restart function pointer */
  ConsFun_t expr;   /*!< pointer to expression that computes quant for output */

}OutputS;

/*----------------------------------------------------------------------------*/
/* typedefs for functions:
 */


/* function type for interpolation schemes */
/*! \fn void (*WeightFun_t)(GridS *pG, Real x1, Real x2, Real x3,
  Real3Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
 *  \brief Interpolation scheme for particles. */
typedef void (*WeightFun_t)(GridS *pG, Real x1, Real x2, Real x3,
  Real3Vect cell1, Real weight[3][3][3], int *is, int *js, int *ks);

#ifdef DELTA_F
/* function types for delta-f scheme */
/*! \fn void (*GetDistFuncFun_t)(GrainS *gr, Real *df);
 *  \brief Computes delta-f weight for particles. */
typedef void (*GetDistFuncFun_t)(GrainS *gr, Real *df); 
/*! \fn void (*SetBackgroundFun_t)(GridS *pG)
 *  \brief Sets (analytic) moments of the background distribution function. */
typedef void (*SetBackgroundFun_t)(GridS *pG);
#endif
/* function type for initializing distribution functions */
/*! \fn void (*InitDistFuncFun_t)(long int *idum, Real **vp,
 *                                const long int Np, const long int NpTot);
 *  \brief Initializes distribution function for particles. */
typedef void (*InitDistFuncFun_t)(long int *idum, Real **vp, 
                                  const long int Np, const long int NpTot);

/*----------------------------------------------------------------------------*/
/*! \enum BCDirection
 *  \brief Directions for the set_bvals_fun() function */
enum BCDirection {left_x1, right_x1, left_x2, right_x2, left_x3, right_x3};

#endif /* PEGASUS_H */
