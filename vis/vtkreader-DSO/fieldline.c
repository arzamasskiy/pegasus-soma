#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "read_vtk.h"
#include "math.h"
#include "utils.h"
#include "time.h"

#include <mpi.h>
#include <fftw3.h>


#define TINY_NUMBER 1.0e-20

#define NCELL 3
//#define COUNT 500
#define WEIGHTS getwei_TSC
//#define WEIGHTS getwei_linear

#define Real double
#define MPIREAL MPI_DOUBLE
#define SQR(x) ((x)*(x))

struct Triplet{
  short indx;
  short indy;
  short indz;
};

typedef struct Real3Vect_s{
  Real x1,x2,x3;
}Real3Vect;


static const Real steps[] = {0.0, 0.5, 0.5, 1.0};
static const Real coeff[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

int rank = 0, numtasks;

void shiftArray(float ***in,float ***out, int nx,int ny,int nz, int lx,int ly,int lz);
void shuffle_triplets(struct Triplet *trips, long size, long int idum);

void trilinear_int(struct Parameters p, float*** arr, Real3Vect idx,
                   Real x0,Real y0, Real z0,Real *b1);
int getvalues(struct Fields* fld, struct Moments *mom, struct Parameters* par
              ,Real weight[3][3][3], int is, int js, int ks
              ,Real *bfld1, Real *bfld2, Real *bfld3
              ,Real *ufld1, Real *ufld2, Real *ufld3);
void getwei_linear(const struct Parameters *par, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);
void getwei_TSC(const struct Parameters *par, Real x1, Real x2, Real x3, Real3Vect cell1,
                           Real weight[3][3][3], int *is, int *js, int *ks);

int celli(const struct Parameters* par, const Real x, const Real dx1_1, int *i, Real *a);
int cellj(const struct Parameters* par, const Real x, const Real dx1_1, int *i, Real *a);
int cellk(const struct Parameters* par, const Real x, const Real dx1_1, int *i, Real *a);

static double ran2(long int *idum);

/*------------------------------------------------------------------------------
 *  This file takes B and u on a Eulerian grid and traces a fieldline, 
 *  calculating d u_parallel / d l_\parallel as it goes.  
 * 
 *  It should then compute the parallel velocity gradient power spectrum,
 *  hopefully demonstrating sharp cut-off.
 *
 *  It should also hopefully ensemble average all of this.
 *
 */
int main(int argc, char* argv[]){
  struct Fields fld;
  struct Moments mom;
  struct Parameters p;
  struct Triplet *trips;

  double *upar, *par_ps_ave, *par_ps_root;
  fftw_complex *par_ps;
  fftw_plan par_plan;
  
  double *uperp, *perp_ps_ave, *perp_ps_root;
  fftw_complex *perp_ps;
  fftw_plan perp_plan;
  
  FILE* fp;
  char fname[50];

  long num, ncell;
  long int idum;

  int dlog, nlog;

  int i,j,k, nsteps, count, ave_c;
  int nx,ny,nz,indx,indy,indz;

  Real lmax, dk; 
  Real Lx, Ly, Lz;
  Real x0,y0,z0;
  Real gdx, gdy, gdz;
  Real sum1, sum2;
  Real isum1, isum2;

  //RK4 temp
  Real xi, yi, zi, xj, yj,zj;
  Real kxi, kyi, kzi;


  //Real weights[3][3][3];
  Real3Vect idx;

  Real dl;

  Real bx,by,bz,ux,uy,uz;
  Real Bmag, iBmag;

  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == 0){
    if(argc != 6){
      printf("Usage: %s [fld_vtk] [mom_vtk] [dl] [lmax] [count]\n",argv[0]);
      exit(0);
    }

    read_fld(argv[1],0,&p,&fld);
    read_mom(argv[2],0,&p,&mom);

    dl = atof(argv[3]);
    lmax = atof(argv[4]);
    count = atoi(argv[5]);

    free_contiguous_3dArray(mom.n);
    mom.n = NULL;

    idum = (long int) time(NULL);
    
    sprintf(fname,"single_field_%.0f.dat",p.time);
    fp = fopen(fname,"w");
  } 

  MPI_Bcast(&p.nx,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.ny,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.nz,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.dx,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.dy,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.dz,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&dl,1,MPIREAL,0,MPI_COMM_WORLD);
  MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&lmax,1,MPIREAL,0,MPI_COMM_WORLD);
  MPI_Bcast(&idum,1,MPI_LONG,0,MPI_COMM_WORLD); //make sure use the same random number

  nx = p.nx; ny = p.ny; nz = p.nz;
  gdx = p.dx; gdy = p.dy; gdz = p.dz;

  Lx = p.nx*gdx;
  Ly = p.ny*gdy;
  Lz = p.nz*gdz;
  

  dk = 2.0*M_PI/lmax;

  idx.x1 = 1.0/gdx;
  idx.x2 = 1.0/gdy;
  idx.x3 = 1.0/gdz;

  if(rank != 0){
    fld.b1   = new_contiguous_3dArray(nz,ny,nx);
    fld.b2   = new_contiguous_3dArray(nz,ny,nx);
    fld.b3   = new_contiguous_3dArray(nz,ny,nx);
    fld.e1 = NULL; fld.e2 = NULL; fld.e3 = NULL;

    mom.n = NULL;
    mom.u1   = new_contiguous_3dArray(nz,ny,nx);
    mom.u2   = new_contiguous_3dArray(nz,ny,nx);
    mom.u3   = new_contiguous_3dArray(nz,ny,nx);
    mom.p11 = NULL; mom.p12 = NULL; mom.p13 = NULL;
    mom.p22 = NULL; mom.p23 = NULL; mom.p33 = NULL;
  }
  

  MPI_Bcast(&(fld.b1[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(fld.b2[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(fld.b3[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(mom.u1[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(mom.u2[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(mom.u3[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);

  nsteps = lmax/dl;

  dlog = (int)(gdx/(3.0*dl));
  nlog = (nsteps/dlog) + 1;
  
  if(dlog == 0) dlog = 1;

  upar = fftw_malloc(sizeof(double)*2*(nlog+1));
  par_ps  = fftw_malloc(sizeof(fftw_complex)*(nlog+2));
  par_ps_ave = fftw_malloc(sizeof(double)*(nlog+2));
  par_ps_root = fftw_malloc(sizeof(double)*(nlog+2));

  uperp = fftw_malloc(sizeof(double)*2*(nlog+1));
  perp_ps  = fftw_malloc(sizeof(fftw_complex)*(nlog+2));
  perp_ps_ave = fftw_malloc(sizeof(double)*(nlog+2));
  perp_ps_root = fftw_malloc(sizeof(double)*(nlog+2));

  par_plan  = fftw_plan_dft_r2c_1d(2*(nlog+1), upar,par_ps,FFTW_ESTIMATE);
  perp_plan = fftw_plan_dft_r2c_1d(2*(nlog+1), uperp,perp_ps,FFTW_ESTIMATE);

  ncell = nx*ny*nz;
  if(ncell % numtasks != 0){
    printf("Please choose a divisor of %ld for the number of processes.\n",ncell);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(1);
  }

  
  for(i = 0; i < nlog+2; i++){
    par_ps_ave[i] = 0.0;
    perp_ps_ave[i] = 0.0;
  }

  trips = malloc(sizeof(struct Triplet)*nx*ny*nz);

  // contruct triplets for starting points
  for(k = 0; k < nz; k++){
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        struct Triplet t;
        t.indx = i;
        t.indy = j;
        t.indz = k;
        trips[k*(ny*nx)+j*nx+i]=t;
      }
    }
  }

  //random sort here
  shuffle_triplets(trips,ncell,idum);


  //main loop
  long offset = (ncell/numtasks)*rank;

  //for(num = 0; num <  ncell/numtasks; num++){
  for(num = 0; num < count; num++){

    for(i = 0; i < 2*(nlog+1); i++) upar[i] = 0.0;
    for(i = 0; i < 2*(nlog+1); i++) uperp[i] = 0.0;

    sum1 = 0.0; sum2 =0.0;
    ave_c = 0;

    indx = trips[num + offset].indx;
    indy = trips[num + offset].indy;
    indz = trips[num + offset].indz;

    x0 = gdx*indx + p.x0;
    y0 = gdy*indy + p.y0;
    z0 = gdz*indz + p.z0;

    for(i = 0; i < nsteps; i++){
      xi = x0;
      yi = y0;
      zi = z0;
      xj = x0;
      yj = y0;
      zj = z0;
      kxi = 0.0;
      kyi = 0.0;
      kzi = 0.0;

      for ( j = 0; j < 4; j++){

        xi = xj + dl*steps[j]*kxi;
        yi = yj + dl*steps[j]*kyi;
        zi = zj + dl*steps[j]*kzi;

        if(xi < p.x0)   xi += Lx;
        if(yi < p.y0)   yi += Ly;
        if(zi < p.z0)   zi += Lz;

        if(xi >= Lx + p.x0) xi -= Lx;
        if(yi >= Ly + p.y0) yi -= Ly;
        if(zi >= Lz + p.z0) zi -= Lz;

        trilinear_int(p,fld.b1,idx,xi,yi,zi,&bx);
        trilinear_int(p,fld.b2,idx,xi,yi,zi,&by);
        trilinear_int(p,fld.b3,idx,xi,yi,zi,&bz);

        trilinear_int(p,mom.u1,idx,xi,yi,zi,&ux);
        trilinear_int(p,mom.u2,idx,xi,yi,zi,&uy);
        trilinear_int(p,mom.u3,idx,xi,yi,zi,&uz);

        //get B and u by interpolation
       // WEIGHTS(&p, xi, yi, zi, idx, weights, &is, &js, &ks);
       // getvalues(&fld, &mom, &p, weights, is, js, ks, &bx, &by, &bz, &ux, &uy, &uz);

        Bmag = sqrt(bx*bx + by*by + bz*bz);
        iBmag = 1.0 / Bmag;

        if(j==0 ){
          if(i % dlog == 0){
            k = i / dlog; 
            upar[k] = iBmag*(ux*bx + uy*by + uz*bz);
            uperp[k] = sqrt(ux*ux + uy*uy + uz*uz - upar[k]*upar[k]);
            sum1 += upar[k];
            sum2 += uperp[k];
            ave_c++;
            if(num == 0 && rank == 0) 
              fprintf(fp,"%e %e %e %e %e %e %e %e %e %e\n",dl*i,x0,y0,z0,Bmag,bx,by,bz,upar[k], uperp[k]);
          }
        }

        kxi = bx*iBmag;
        kyi = by*iBmag;
        kzi = bz*iBmag;

        x0 += dl*coeff[j]*kxi;
        y0 += dl*coeff[j]*kyi;
        z0 += dl*coeff[j]*kzi;
      }
      //dx = dl*bx*iBmag;
      //dy = dl*by*iBmag;
      //dz = dl*bz*iBmag;
      //x0 += dx;
      //y0 += dy;
      //z0 += dz;
    
      if(x0 < p.x0)   x0 += Lx;
      if(y0 < p.y0)   y0 += Ly;
      if(z0 < p.z0)   z0 += Lz;

      if(x0 >= Lx + p.x0) x0 -= Lx;
      if(y0 >= Ly + p.y0) y0 -= Ly;
      if(z0 >= Lz + p.z0) z0 -= Lz;
    }

   
    isum1 = sum1/((double)ave_c);
    isum2 = sum2/((double)ave_c);

    for(i = 0; i < ave_c; i++) upar[i]  -= isum1;
    for(i = 0; i < ave_c; i++) uperp[i] -= isum2;

    fftw_execute(par_plan);
    fftw_execute(perp_plan);
    for(i = 0; i < nlog+2; i++){
        par_ps_ave[i] +=  par_ps[i][0]*par_ps[i][0] + par_ps[i][1]*par_ps[i][1];
        perp_ps_ave[i] +=  perp_ps[i][0]*perp_ps[i][0] + perp_ps[i][1]*perp_ps[i][1];
    }
  } 

  MPI_Reduce(par_ps_ave,par_ps_root,nlog+2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(perp_ps_ave,perp_ps_root,nlog+2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank == 0){ 
    fclose(fp);
    sprintf(fname,"spectra_%.0f.dat",p.time);
    fp = fopen(fname,"w");
    for(i = 0; i < nlog+2; i++){
      par_ps_root[i] /=  (4.0*numtasks*count*(nlog+1.0)*(nlog+1.0));
      perp_ps_root[i] /=  (4.0*numtasks*count*(nlog+1.0)*(nlog+1.0));
      fprintf(fp,"%e %e %e\n",dk*i,par_ps_root[i], perp_ps_root[i]);  
    }
    fclose(fp);
  }

  freeFieldStruct(&fld);
  freeMomentStruct(&mom);
  free(trips);
  free(upar);
  free(uperp);
  MPI_Finalize();
  return 0;
}

/*------------------------------------------------------------------------------
 *  Randomly shuffle array of triplets
 */
void shuffle_triplets(struct Triplet *trips, long size, long int idum){
  idum= idum > 0 ? -idum : idum;

  long i;
  for(i = 0; i < size; i++){
    struct Triplet tmp;
    double rand = ran2(&idum);
    long swap = i +(long)(rand*(size-i));
    tmp = trips[swap];
    trips[swap] = trips[i];
    trips[i] = tmp;
  }
}

void shiftArray(float ***in,float ***out,int nx,int ny,int nz, int lx,int ly,int lz){
  int i,j,k;
  int indx,indy,indz;
  for(i = 0; i < nx; i++){
    indx = (i + lx) % nx;
    for(j = 0; j < ny; j++){
      indy = (j + ly) % ny;
      for(k = 0; k < nz; k++){
        indz = (k + lz) % nz;
        out[i][j][k] =in[indx][indy][indz];
      } 
    } 
  }
}

/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

void trilinear_int(struct Parameters p, float*** arr, Real3Vect idx,
                   Real x0,Real y0, Real z0,Real *b1){
    int i0,j0,k0,i1,j1,k1;
    Real xl,yl,zl,xd,yd,zd;
    Real c00,c01,c10,c11,c0,c1;

    i0 = (int) (idx.x1*(x0 - p.x0));
    j0 = (int) (idx.x2*(y0 - p.y0));
    k0 = (int) (idx.x3*(z0 - p.z0));

    i1 = i0 + 1;
    j1 = j0 + 1;
    k1 = k0 + 1;
  
    if(i1 >= p.nx) i1 -= p.nx;
    if(j1 >= p.ny) j1 -= p.ny;
    if(k1 >= p.nz) k1 -= p.nz;

    xl = i0 * p.dx + p.x0;
    yl = j0 * p.dy + p.y0;
    zl = k0 * p.dz + p.z0;

    xd = (x0 - xl)*idx.x1;
    yd = (y0 - yl)*idx.x2;
    zd = (z0 - zl)*idx.x3;

    c00 = arr[k0][j0][i0]*(1.0-xd) + arr[k0][j0][i1]*xd;
    c01 = arr[k1][j0][i0]*(1.0-xd) + arr[k1][j0][i1]*xd;
    c10 = arr[k0][j1][i0]*(1.0-xd) + arr[k0][j1][i1]*xd;
    c11 = arr[k1][j1][i0]*(1.0-xd) + arr[k1][j1][i1]*xd;

    c0 = c00*(1.0-yd) + c10*yd;
    c1 = c01*(1.0-yd) + c11*yd;

    *b1 = c0*(1.0-zd) + c1*zd;
}

int getvalues(struct Fields* fld, struct Moments *mom, struct Parameters* par
              ,Real weight[3][3][3], int is, int js, int ks
              ,Real *bfld1, Real *bfld2, Real *bfld3
              ,Real *ufld1, Real *ufld2, Real *ufld3
              ){
  int n0,i,j,k;
  Real b1, b2, b3;
  Real u1, u2, u3;
  Real totwei, totwei1;    /* total weight (in case of edge cells) */

  int iind, jind, kind;
  
  /* linear interpolation */
  b1 = 0.0; b2 = 0.0; b3 = 0.0; 
  u1 = 0.0; u2 = 0.0; u3 = 0.0; 
  totwei = 0.0;    totwei1 = 1.0;
  
  /* Interpolate E & B fields */
  /* Note: in lower dimensions only wei[0] is non-zero */
  n0 = NCELL-1;
  
  for (i=0; i < NCELL; i++) {
    for (j=0; j < NCELL; j++) {
      for (k=0; k < NCELL; k++) {
        
        iind = i + is; 
        jind = j + js; 
        kind = k + ks; 

        if(iind < 0) iind += par->nx;
        if(iind >= par->nx) iind -= par->nx;

        if(jind < 0) jind += par->ny;
        if(jind >= par->ny) jind -= par->ny;

        if(kind < 0) kind += par->nz;
        if(kind >= par->nz) kind -= par->nz;

        b1 += weight[i][j][k] * fld->b1[iind][jind][kind];
        b2 += weight[i][j][k] * fld->b2[iind][jind][kind];
        b3 += weight[i][j][k] * fld->b3[iind][jind][kind];

        u1 += weight[i][j][k] * mom->u1[iind][jind][kind];
        u2 += weight[i][j][k] * mom->u2[iind][jind][kind];
        u3 += weight[i][j][k] * mom->u3[iind][jind][kind];

        totwei += weight[i][j][k];
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;
  
  /* ensure interpolated E dot interpolated B = interpolated E dot B */
  totwei1 = 1.0/totwei;
  *bfld1 = b1*totwei1; *bfld2 = b2*totwei1; *bfld3 = b3*totwei1;
  *ufld1 = u1*totwei1; *ufld2 = u2*totwei1; *ufld3 = u3*totwei1;
  
  return 0;
}

/*----------------------------------------------------------------------------*/
/*! \fn void getwei_linear(GridS *pG, Real x1, Real x2, Real x3,
 *           Real3Vect cell1, Real weight[3][3][3], int *is, int *js, int *ks)
 *  \brief Get weight using linear interpolation
 *
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */

void getwei_linear(const struct Parameters *par, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c;    /* grid coordinate for the position (x1,x2,x3) */
  Real wei1[2], wei2[2], wei3[2];/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  i = celli(par, x1, cell1.x1, &i1, &a);  /* x1 index */
  i1 = i+i1-1;  *is = i1;    /* starting x1 index */
  wei1[1] = a - i1 - 0.5;      /* one direction weight */
  wei1[0] = 1.0 - wei1[1];      /* 0: left; 1: right */

  /* x2 direction */
  j = cellj(par, x2, cell1.x2, &j1, &b);  /* x2 index */
  j1 = j+j1-1;  *js = j1;    /* starting x2 index */
  wei2[1] = b - j1 - 0.5;      /* one direction weight */
  wei2[0] = 1.0 - wei2[1];      /* 0: left; 1: right */

  /* x3 direction */
  k = cellk(par, x3, cell1.x3, &k1, &c);  /* x3 index */
  k1 = k+k1-1;  *ks = k1;    /* starting x3 index */
  wei3[1] = c - k1 - 0.5;      /* one direction weight */
  wei3[0] = 1.0 - wei3[1];      /* 0: left; 1: right */
  

  /* calculate 3D weight */
  for (k=0; k<2; k++)
    for (j=0; j<2; j++)
      for (i=0; i<2; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void getwei_TSC(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
 *                      Real weight[3][3][3], int *is, int *js, int *ks)
 *
 *  \brief Get weight using Triangular Shaped Cloud (TSC) interpolation 
 * Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
 * Output: weight: weight function; is,js,ks: starting cell indices in the grid.
 * Note: this interpolation works in any 1-3 dimensions.
 */
void getwei_TSC(const struct Parameters *par, Real x1, Real x2, Real x3, Real3Vect cell1,
                           Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k;
  Real a, b, c, d;  /* grid coordinate for the position (x1,x2,x3) */
  Real wei1[3], wei2[3], wei3[3];/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  celli(par, x1, cell1.x1, &i, &a);    /* x1 index */
  *is = i - 1;        /* starting x1 index, wei[0] */
  d = a - i;
  wei1[0] = 0.5*SQR(1.0-d);      /* 0: left; 2: right */
  wei1[1] = 0.75-SQR(d-0.5);      /* one direction weight */
  wei1[2] = 0.5*SQR(d);

  /* x2 direction */
  cellj(par, x2, cell1.x2, &j, &b);    /* x2 index */
  *js = j - 1;        /* starting x2 index */
  d = b - j;
  wei2[0] = 0.5*SQR(1.0-d);      /* 0: left; 2: right */
  wei2[1] = 0.75-SQR(d-0.5);      /* one direction weight */
  wei2[2] = 0.5*SQR(d);

  /* x3 direction */
  cellk(par, x3, cell1.x3, &k, &c);    /* x3 index */
  *ks = k - 1;        /* starting x3 index */
  d = c - k;
  wei3[0] = 0.5*SQR(1.0-d);      /* 0: left; 2: right */
  wei3[1] = 0.75-SQR(d-0.5);      /* one direction weight */
  wei3[2] = 0.5*SQR(d);

  /* calculate 3D weight */
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}

/*----------------------------------------------------------------------------*/
/* Input: pGrid: grid; x: global x coordinate;
 *        dx1_1: 1/dx1 (to improve performance)
 * Output: i: i-index containing x; a: grid index coordinate of x;
 * Return: 0: x is on the left of the ith cell;
 *         1: x is on the right of the ith cell;
 */

/*! \fn int celli(const GridS* pG, const Real x, const Real dx1_1, 
 *      int *i, Real *a)
 *  \brief given x, returns containing cell first index.  */
int celli(const struct Parameters* par, const Real x, const Real dx1_1, int *i, Real *a)
{
  *a = (x - par->x0) * dx1_1; // + (Real)pG->is;
  *i = (int)(*a);
  if (((*a)-(*i)) < 0.5) return 0;  /* in the left half of the cell*/
  else return 1;      /* in the right half of the cell*/
}

/*! \fn cellj(const GridS* pG, const Real y, const Real dx2_1, 
 *        int *j, Real *b)
 *  \brief given y, returns containing cell first index.  */
int cellj(const struct Parameters* par, const Real y, const Real dx2_1, int *j, Real *b)
{
  *b = (y - par->y0) * dx2_1;// + (Real)pG->js;
  *j = (int)(*b);
  if (((*b)-(*j)) < 0.5) return 0;  /* in the left half of the cell*/
  else return 1;      /* in the right half of the cell*/
}

int cellk(const struct Parameters* par, const Real z, const Real dx3_1, int *k, Real *c)
{
  *c = (z - par->z0) * dx3_1;// + (Real)pG->ks;
  *k = (int)(*c);
  if (((*c)-(*k)) < 0.5) return 0;  /* in the left half of the cell*/
  else return 1;      /* in the right half of the cell*/
}
