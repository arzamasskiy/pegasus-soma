#include "../copyright.h"
/*===========================================================================*/
/*! \file output_particle.c
 *  \brief Contains routines necessary for outputting particles.
 *
 * PURPOSE: contains routines necessary for outputting particles.
 *   There are two basic formats:
 *  - 1. Bin particles to the grid, then output particles as a grid.
 *  - 2. Dump the particle list directly.
 *
 *   For particle binning, there can be many choices since particles may have
 *   different properties. We provide a default (and trivial) particle selection
 *   function, which select all the particles with different properties. The 
 *   user can define their own particle selection functions in the problem
 *   generator and pass them to the main code.
 *
 *   The output quantities include, density, momentum density and velocity of
 *   the selected particles averaged in one grid cell. The binned data are
 *   saved in arrays dpar, and grid_M. The latter is borrowed from particle.c
 *   to save memory. The expression functions expr_??? are used to pick the 
 *   relevant quantities, which is part of the output data structure. The way
 *   to output these binned particle quantities are then exactly the same as
 *   other gas quantities.
 *
 *   Dumping particle list has not been developed yet since we need to figure
 *   out how to do visualization.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - particle_to_grid();
 * - dump_particle_spectrum();
 * - dump_particle_binary();
 * - dump_particle_heating();
 * - dump_spacecraft()
 * - dump_particle_track();
 * - dump_particle_history();
 * - property_all();
 * 
 *============================================================================*/
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../pegasus.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

/* Maximum Number of default history dump columns. */
#ifdef DELTA_F
#define NSCAL 11
#else
#define NSCAL 9
#endif
#define NBINS 200
#define NSPECBINS 80000

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   expr_*
 *============================================================================*/
static int getbfield(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
              Real *bfld1, Real *bfld2, Real *bfld3);

static int getefield(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
                     Real *efld1, Real *efld2, Real *efld3);

static int getdens(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
                     Real *dens);

Real expr_dpar (const GridS *pG, const int i, const int j, const int k);
Real expr_M1par(const GridS *pG, const int i, const int j, const int k);
Real expr_M2par(const GridS *pG, const int i, const int j, const int k);
Real expr_M3par(const GridS *pG, const int i, const int j, const int k);
Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
Real expr_V3par(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*! \fn void particle_to_grid(Grid *pG, PropFun_t par_prop)
 *  \brief Bin the particles to grid cells
 */
void particle_to_grid(DomainS *pD, PropFun_t par_prop)
{

  deposit(pD,0);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_particle_spectrum(MeshS *pM, OutputS *pOut)
 *  \brief Dump particle spectrum
 */
void dump_particle_spectrum(MeshS *pM, OutputS *pOut)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS   *pG = pD->Grid;
  FILE *pfile;
  char *fname,fmt[80];
  long p;
  int i,mhst,ind;
  //Real vsq,vloc;
  int pis,pjs,pks,iloc,jloc,kloc;
  Real rho,bmag,wx,wy,wz,wprl,wprp,wsq;
  Real bfld1,bfld2,bfld3,bsq;
  Real weight[3][3][3];
  int nbr_prp_bins,nbr_prl_bins,prlloc,prploc,nbr_bins;
  Real3Vect cell1;
  Real prlsize,prpsize,prpmax,prpmin,prlmax,prlmin;
  Real hist[NSPECBINS];
  GrainS *gr;
  //Real efld1,efld2,efld3,edotv;
  Real edotvperp, edotvprl, edotvtot, efld1, efld2, efld3;

  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;

  /* Loop over all Domains in Mesh, and output spectrum data */

  nbr_prp_bins = 200; //100;
  nbr_prl_bins = 400; //200;
  nbr_bins = nbr_prp_bins*nbr_prl_bins;
 
  /* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%14.6e"); /* Use a default format */
  }
  else {
    sprintf(fmt," %s",pOut->dat_fmt);
  }
  /* zero out histogram bins before summation */
  for (i=0; i<nbr_bins; i++) {
    hist[i] = 0.0;
  }
  
  prpmax = 4.0*sqrt(beta);
  prpmin = 0.0;
  prlmax = 4.0*sqrt(beta);
  prlmin = -prlmax;
  
  prpsize = (prpmax-prpmin)/((double)nbr_prp_bins);
  prlsize = (prlmax-prlmin)/((double)nbr_prl_bins);

  /* sort particles into bins */
  for (p=0; p<pG->nparticle; p++)
  {   
    gr = &(pG->particle[p]);

      /*
  vsq = SQR(gr->v1) + SQR(gr->v2) + SQR(gr->v3);
  vloc = (int)floor((log10(vsq)-binmin)/binsize);
  if ((vloc >= 0) && (vloc < nbr_bins))
  hist[vloc+1] += 1;

      */

    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &pis, &pjs, &pks);
    if(pG->Nx[2]>1) kloc = pks+1; else kloc = pks;
    if(pG->Nx[1]>1) jloc = pjs+1; else jloc = pjs;
    if(pG->Nx[0]>1) iloc = pis+1; else iloc = pis;

    rho  = pG->Coup[kloc][jloc][iloc].grid_d;
    bmag = SQR(pG->U[kloc][jloc][iloc].B1c) + SQR(pG->U[kloc][jloc][iloc].B2c) + SQR(pG->U[kloc][jloc][iloc].B3c);
    bmag = sqrt(bmag);
      
    wx = gr->v1 - pG->Coup[kloc][jloc][iloc].grid_M1/rho;
    wy = gr->v2 - pG->Coup[kloc][jloc][iloc].grid_M2/rho;
    wz = gr->v3 - pG->Coup[kloc][jloc][iloc].grid_M3/rho;
    wsq = wx*wx + wy*wy + wz*wz;
      
    wprl = wx*pG->U[kloc][jloc][iloc].B1c + wy*pG->U[kloc][jloc][iloc].B2c + wz*pG->U[kloc][jloc][iloc].B3c;
    wprl /= bmag;
      
    wprp = sqrt(wsq-wprl*wprl);

    if (getefield(pG, weight, pis, pjs, pks, &efld1, &efld2, &efld3) == 0)
    {
    }
    else 
    {
        peg_perr(0,"Requested particle is off grid %d with position (%f,%f,%f)!\n",
                 myID_Comm_world,gr->x1,gr->x2,gr->x3); 
    }
    //edotv = gr->v1*efld1 + gr->v2*efld2 + gr->v3*efld3;
    prlloc = (int)floor((wprl-prlmin)/prlsize);
    prploc = (int)floor((wprp-prpmin)/prpsize);
    
    if ((prlloc >= 0) && (prlloc < nbr_prl_bins) && (prploc >=0) && (prploc < nbr_prp_bins)) {
      ind = prploc*nbr_prl_bins + prlloc;
      hist[ind] += 1;
  /*
    hist[ind] += edotv;
   */
    }
  }
  /* Create filename and open file. */
  fname = peg_fname(NULL,pM->outfilename,NULL,NULL,0,0,NULL,"spec");
  if(fname == NULL){
    peg_perr(-1,"[dump_pspectrum]: Unable to create spec filename\n");
  }
  pfile = fopen(fname,"a");
  if(pfile == NULL){
    peg_perr(-1,"[dump_pspectrum]: Unable to open the spec file\n");
  }
  
  /* Write out column headers, but only for first dump */
  
  if(pOut->num == 0){
    fprintf(pfile,
      "# Pegasus particle distribution: %i vprl-bins from %.3f to %.3f and %i vprp-bins from %.3f to %.3f\n",
      nbr_prl_bins,prlmin,prlmax,nbr_prp_bins,prpmin,prpmax);
    mhst=1;
    fprintf(pfile,"#   [%i]=time   ",mhst);
    for (mhst=2; mhst<=nbr_bins+1; mhst++) {
      fprintf(pfile,"   [%i]=bin%i   ",mhst,mhst-1);
    }
    fprintf(pfile,"\n#\n");
  }
  
  /* Write out data and close file */
  fprintf(pfile,fmt,pM->time);
  for (i=0; i<nbr_bins; i++) {
    fprintf(pfile,fmt,hist[i]);
  }
  fprintf(pfile,"\n");
  fclose(pfile);

  return;
}  

/*----------------------------------------------------------------------------*/
/*! \fn void dump_particle_track(MeshS *pM, OutputS *pOut)
 *  \brief Dump tracked particles in ascii format
 */
void dump_particle_track(MeshS *pM, OutputS *pOut)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS   *pG = pD->Grid;
  FILE *pfile;
  char *fname;
  long p;
  int pid,init_id = 0;
  int i,is,js,ks, myID;
  
  static int init = 1;

  Real mu,qom,bfld1,bfld2,bfld3,bsq,vsq,vprl,vprp,vprpsq;
  Real3Vect cell1;
  Real weight[3][3][3];         /* weight function */
  GrainS *gr;

  /* Get grid limit related quantities */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else               cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else               cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else               cell1.x3 = 0.0;

/* Create track Directory */
  if(init){
    init=0;

#ifdef MPI_PARALLEL
    if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &myID))
      peg_error("[main]: Error on calling MPI_Comm_rank\n");

/* Only rank=0 processor creates track directory */
    if(myID == 0){
      mkdir("../track", 0775); /* Create track directory. May return error (directory already exists) */
    }
    MPI_Barrier(MPI_COMM_WORLD); /* Wait for rank 0 to mkdir() */
#else
    mkdir("../track", 0775); /* Create track directory. May return error (directory already exists) */
#endif
  }

/* Write all the selected particles */

  /* Write particle information */
  for (p=0; p<pG->nparticle; p++)
  {
    gr = &(pG->particle[p]);
    if ((*(pOut->par_prop))(gr,&(pG->parsub[p]))) { /* 1: true; 0: false */

      /* compute first adiabatic invariant */
      getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);
      if (getbfield(pG, weight, is, js, ks, &bfld1, &bfld2, &bfld3) == 0)
      {/* particle is on the grid */
        qom = grproperty[gr->property].qomc;
        bfld1 *= qom;  bfld2 *= qom;  bfld3 *= qom;
      }
      else
      {/* particle is off grid, warning sign */
        peg_perr(0,"Requested particle is off grid %d with position (%f,%f,%f)!\n",
                 myID_Comm_world,gr->x1,gr->x2,gr->x3); /* warning! */
      }
      bsq    = bfld1 *bfld1  + bfld2 *bfld2  + bfld3 *bfld3 ;
      vsq    = gr->v1*gr->v1 + gr->v2*gr->v2 + gr->v3*gr->v3;
      vprl   = gr->v1*bfld1  + gr->v2*bfld2  + gr->v3*bfld3 ;
      vprl   = vprl/sqrt(bsq);
      vprpsq = vsq - vprl*vprl;
      vprp   = sqrt(vprpsq);
      mu     = vprpsq/sqrt(bsq);

      pid      = gr->my_id;
#ifdef MPI_PARALLEL
      init_id  = gr->init_id;
#endif

      if((fname = peg_fname("../track","turb",NULL,NULL,6,100000*pid+init_id,NULL,"lis"))==NULL){
        peg_error("[dump_particle_binary]: Error constructing filename\n");
      }
      if((pfile = fopen(fname,"a"))==NULL){
        peg_error("[dump_particle_binary]: Unable to open particle file %s\n",fname);
      }

      if(pOut->num == 0){
        fprintf(pfile,"# Pegasus tracked particle information for particle number %i starting on grid %i\n",pid,init_id);
        fprintf(pfile,"#   [%i]=time   ",1);
        fprintf(pfile,"#   [%i]=x1     ",2);
        fprintf(pfile,"#   [%i]=x2     ",3);
        fprintf(pfile,"#   [%i]=x3     ",4);
        fprintf(pfile,"#   [%i]=v1     ",5);
        fprintf(pfile,"#   [%i]=v2     ",6);
        fprintf(pfile,"#   [%i]=v3     ",7);
        fprintf(pfile,"#   [%i]=mu     ",8);
        fprintf(pfile,"#   [%i]=vprl   ",9);
        fprintf(pfile,"#   [%i]=vprp   ",10);
        fprintf(pfile,"\n#\n");
      }

      fprintf(pfile,"%f\t",pG->time);
      fprintf(pfile,"%f\t",(float)(gr->x1));
      fprintf(pfile,"%f\t",(float)(gr->x2));
      fprintf(pfile,"%f\t",(float)(gr->x3));
      fprintf(pfile,"%f\t",(float)(gr->v1));
      fprintf(pfile,"%f\t",(float)(gr->v2));
      fprintf(pfile,"%f\t",(float)(gr->v3));
      fprintf(pfile,"%f\t",(float)(mu));
      fprintf(pfile,"%f\t",(float)(vprl));
      fprintf(pfile,"%f\t",(float)(vprp));
      fprintf(pfile,"\n");
      fclose(pfile);

    }
  }

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void dump_particle_history(MeshS *pM, OutputS *pOut)
 *  \brief Dump average particle properties to history file
 */
void dump_particle_history(MeshS *pM, OutputS *pOut)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS   *pG = pD->Grid;
  FILE *pfile;
  char *fname,fmt[80];
  long p;
  int i,is,js,ks;
  int mhst,total_phst_cnt,myID_Comm_Domain=1;
  Real pcount,scal[NSCAL+1];     /* +1 for particle count */
  Real bfld1,bfld2,bfld3,bsq,bmag,vprl,vprlsq,vprpsq;
  Real3Vect cell1;
  Real weight[3][3][3];
  GrainS *gr;
#ifdef MPI_PARALLEL
  Real my_scal[NSCAL+1];  /* +1 for particle count */
  int ierr;
#endif

#ifdef DELTA_F
  Real f0_t;
  Real dfw,dfw_sum, dfw_sq;
  total_phst_cnt = 11;
#else
  total_phst_cnt = 9;
#endif
  
  /* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%14.6e"); /* Use a default format */
  }
  else {
    sprintf(fmt," %s",pOut->dat_fmt);
  }
  
  /* store time and dt in first two elements of output vector */
  
  scal[0] = pM->time;
  scal[1] = pM->dt;
  
  /* Get grid limit related quantities */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else               cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else               cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else               cell1.x3 = 0.0;
  
  /* zero out average quantities before summation */
  for (i=2; i<=total_phst_cnt; i++) {
    scal[i] = 0.0;
  }
  
  /* sum properties over all particles on grid */
  for (p=0; p<pG->nparticle; p++)
  {
    
    gr = &(pG->particle[p]);
    
    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);
    if (getbfield(pG, weight, is, js, ks, &bfld1, &bfld2, &bfld3) == 0)
    {/* particle is on the grid */

      bsq    = bfld1 *bfld1  + bfld2 *bfld2  + bfld3 *bfld3 ;
      bmag   = sqrt(bsq);
      vprl   = gr->v1*bfld1  + gr->v2*bfld2  + gr->v3*bfld3 ;
      vprlsq = SQR(vprl)/bsq;
      vprpsq = SQR(gr->v1) + SQR(gr->v2) + SQR(gr->v3) - vprlsq;


      mhst = 2;
      scal[mhst] += SQR(gr->v1);
      mhst++;
      scal[mhst] += SQR(gr->v2);
      mhst++;
      scal[mhst] += SQR(gr->v3);
      mhst++;
      scal[mhst] += vprlsq;
      mhst++;
      scal[mhst] += vprpsq;
      mhst++;
      scal[mhst] += vprpsq/bmag;
      mhst++;
      scal[mhst] += sqrt(vprpsq)/bmag;
      
#ifdef DELTA_F
      getdf(gr, &f0_t);
      dfw = (1.0 - (f0_t/gr->f_0));
      mhst++;
      scal[mhst] += dfw;
      mhst++;
      scal[mhst] += dfw*dfw;
#endif
    }
    else
    {/* particle is off grid, warning sign */
      peg_perr(0,"Requested particle is off grid %d with position (%f,%f,%f)!\n",
               myID_Comm_world,gr->x1,gr->x2,gr->x3); /* warning! */
    }
    
  }
  
  /* Compute the sum over all Grids in Domain. While doing this, also get total number
   * of particles on Mesh and collect in root process. */
  
  pcount = (double)pG->nparticle;
#ifdef MPI_PARALLEL
  for(i=2; i<total_phst_cnt; i++){
    my_scal[i] = scal[i];
  }
  my_scal[total_phst_cnt] = pcount;
  ierr = MPI_Reduce(&(my_scal[2]),&(scal[2]),(total_phst_cnt-1),
                    MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
  pcount = scal[total_phst_cnt];
#endif /* MPI_PARALLEL */
  
  /* Only the parent (rank=0) process computes the average and writes output.
   * For single-processor jobs, myID_Comm_world is always zero */
  
#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
#endif
  if((myID_Comm_Domain==0) || (myID_Comm_world==0)){ /* I'm the parent */
    
    /* Compute particle averages */
    for(i=2; i<total_phst_cnt; i++){
      scal[i] /= pcount;
    }
    
    /* Create filename and open file. History files are always written in directory
     * of root process (rank=0 in MPI_COMM_WORLD) */
    
    fname = peg_fname(NULL,pM->outfilename,NULL,NULL,0,0,NULL,"phst");
    if(fname == NULL){
      peg_perr(-1,"[dump_phistory]: Unable to create phistory filename\n");
    }
    pfile = fopen(fname,"a");
    if(pfile == NULL){
      peg_perr(-1,"[dump_phistory]: Unable to open the phistory file\n");
    }
    
    /* Write out column headers, but only for first dump */
    
    mhst = 0;
    if(pOut->num == 0){
      fprintf(pfile,
              "# Pegasus particle history dump for %e particles\n",pcount);
      mhst++;
      fprintf(pfile,"#   [%i]=time   ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=dt      ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=v1sq    ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=v2sq    ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=v3sq    ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=vprlsq  ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=vprpsq  ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=mu      ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=rho     ",mhst);
#ifdef DELTA_F
      mhst++;
      fprintf(pfile,"   [%i]=<w_p>   ",mhst);
      mhst++;
      fprintf(pfile,"   [%i]=<w_p^2> ",mhst);
#endif
      
      fprintf(pfile,"\n#\n");
    }
    
    /* Write out data and close file */
    
    for (i=0; i<total_phst_cnt; i++) {
      fprintf(pfile,fmt,scal[i]);
    }
    fprintf(pfile,"\n");
    fclose(pfile);
  }
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn int property_all(const Grain *gr, const GrainAux *grsub)
 *  \brief Default choice for binning particles to the grid: 
 * All the particles are binned, return true for any value.
 */
int property_all(const GrainS *gr, const GrainAux *grsub)
{
    return 1;  /* always true */
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/*! \fn int getbfield()
 *  \brief Get interpolated value of B-field using the weight
 *
 * Input:
 * - pG: grid; weight: weight function;
 * - is,js,ks: starting cell indices in the grid.
 * Output:
 * - interpolated values of B-field
 *
 * Return: 0: normal exit; -1: particle lie out of the grid, cannot interpolate
 * Note: this interpolation works in any 1-3 dimensions.
 */

int getbfield(GridS *pG, Real weight[3][3][3], int is, int js, int ks, 
              Real *bfld1, Real *bfld2, Real *bfld3
              ){
  int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  Real b1, b2, b3, bsq;    /* B-field */
  Real totwei, totwei1;    /* total weight (in case of edge cells) */
  ConsS *pCons;
  
  /* linear interpolation */
  b1 = 0.0; b2 = 0.0; b3 = 0.0;
  totwei = 0.0;    totwei1 = 1.0;
  
  /* Interpolate B field */
  /* Note: in lower dimensions only wei[0] is non-zero */
  n0 = ncell-1;
  k1 = MAX(ks, klp);  k2 = MIN(ks+n0, kup);
  j1 = MAX(js, jlp);  j2 = MIN(js+n0, jup);
  i1 = MAX(is, ilp);  i2 = MIN(is+n0, iup);
  
  for (k=k1; k<=k2; k++) {
    k0=k-k1;
    for (j=j1; j<=j2; j++) {
      j0=j-j1;
      for (i=i1; i<=i2; i++) {
        i0=i-i1;
        
        pCons = &(pG->U[k][j][i]);
        b1 += weight[k0][j0][i0] * pCons->B1c;
        b2 += weight[k0][j0][i0] * pCons->B2c;
        b3 += weight[k0][j0][i0] * pCons->B3c;
        totwei += weight[k0][j0][i0];
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;
    
  totwei1 = 1.0/totwei;
  *bfld1 = b1*totwei1; *bfld2 = b2*totwei1; *bfld3 = b3*totwei1;
  
  return 0;
}

/*! \fn int getefield()
 *  \brief Get interpolated value of E-field using the weight
 *
 * Input:
 * - pG: grid; weight: weight function;
 * - is,js,ks: starting cell indices in the grid.
 * Output:
 * - interpolated values of E-field
 *
 * Return: 0: normal exit; -1: particle lie out of the grid, cannot interpolate
 * Note: this interpolation works in any 1-3 dimensions.
 */

int getefield(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
              Real *efld1, Real *efld2, Real *efld3
              ){
  int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  Real e1, e2, e3, esq;    /* E-field */
  Real totwei, totwei1;    /* total weight (in case of edge cells) */
  ConsS *pCons;
  
  /* linear interpolation */
  e1 = 0.0; e2 = 0.0; e3 = 0.0;
  totwei = 0.0;    totwei1 = 1.0;
  
  /* Interpolate E field */
  /* Note: in lower dimensions only wei[0] is non-zero */
  n0 = ncell-1;
  k1 = MAX(ks, klp);  k2 = MIN(ks+n0, kup);
  j1 = MAX(js, jlp);  j2 = MIN(js+n0, jup);
  i1 = MAX(is, ilp);  i2 = MIN(is+n0, iup);
  
  for (k=k1; k<=k2; k++) {
    k0=k-k1;
    for (j=j1; j<=j2; j++) {
      j0=j-j1;
      for (i=i1; i<=i2; i++) {
        i0=i-i1;
        
        pCons = &(pG->U[k][j][i]);
        e1 += weight[k0][j0][i0] * pCons->E1c;
        e2 += weight[k0][j0][i0] * pCons->E2c;
        e3 += weight[k0][j0][i0] * pCons->E3c;
        totwei += weight[k0][j0][i0];
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;
  
  totwei1 = 1.0/totwei;
  *efld1 = e1*totwei1; *efld2 = e2*totwei1; *efld3 = e3*totwei1;
  
  return 0;
}

/*! \fn int getdens()
 *  \brief Get interpolated value of density using the weight
 *
 * Input:
 * - pG: grid; weight: weight function;
 * - is,js,ks: starting cell indices in the grid.
 * Output:
 * - interpolated values of density
 *
 * Return: 0: normal exit; -1: particle lie out of the grid, cannot interpolate
 * Note: this interpolation works in any 1-3 dimensions.
 */

int getdens(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
              Real *dens
              ){
  int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  Real d;    /* density */
  Real totwei, totwei1;    /* total weight (in case of edge cells) */
  /* linear interpolation */
  d = 0.0;
  totwei = 0.0;    totwei1 = 1.0;
  
  /* Interpolate density */
  /* Note: in lower dimensions only wei[0] is non-zero */
  n0 = ncell-1;
  k1 = MAX(ks, klp);  k2 = MIN(ks+n0, kup);
  j1 = MAX(js, jlp);  j2 = MIN(js+n0, jup);
  i1 = MAX(is, ilp);  i2 = MIN(is+n0, iup);
  
  for (k=k1; k<=k2; k++) {
    k0=k-k1;
    for (j=j1; j<=j2; j++) {
      j0=j-j1;
      for (i=i1; i<=i2; i++) {
        i0=i-i1;
        if(pG->Coup[k][j][i].grid_d != 0.0){
          d += weight[k0][j0][i0] * pG->Coup[k][j][i].grid_d;
          totwei += weight[k0][j0][i0];
        }
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;
  
  totwei1 = 1.0/totwei;
  *dens = d*totwei1;
  
  return 0;
}

/* expr_*: where * are variables d,M1,M2,M3,V1,V2,V3 for particles */

/*! \fn Real expr_dpar(const Grid *pG, const int i, const int j, const int k) 
 *  \brief Wrapper for particle density */
Real expr_dpar(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_d;
}
/*! \fn Real expr_M1par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 1-momentum */
Real expr_M1par(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_M1;
}

/*! \fn Real expr_M2par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 2-momentum */
Real expr_M2par(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_M2;
}
/*! \fn Real expr_M3par(const Grid *pG, const int i, const int j, const int k) 
 *  \brief Wrapper for particle 3-momentum */
Real expr_M3par(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_M3;
}
/*! \fn Real expr_V1par(const Grid *pG, const int i, const int j, const int k) 
 *  \brief Wrapper for particle 1-velocity */
Real expr_V1par(const GridS *pG, const int i, const int j, const int k) {
  if (pG->Coup[k][j][i].grid_d>0.0)
    return pG->Coup[k][j][i].grid_M1/pG->Coup[k][j][i].grid_d;
  else return 0.0;
}
/*! \fn Real expr_V2par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 2-velocity */
Real expr_V2par(const GridS *pG, const int i, const int j, const int k) {
  if (pG->Coup[k][j][i].grid_d>0.0)
    return pG->Coup[k][j][i].grid_M2/pG->Coup[k][j][i].grid_d;
  else return 0.0;
}
/*! \fn Real expr_V3par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 3-velocity */
Real expr_V3par(const GridS *pG, const int i, const int j, const int k) {
  if (pG->Coup[k][j][i].grid_d>0.0)
    return pG->Coup[k][j][i].grid_M3/pG->Coup[k][j][i].grid_d;
  else return 0.0;
}

#undef NSCAL
#undef NBINS

