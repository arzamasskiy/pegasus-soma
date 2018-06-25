#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "read_vtk.h"
#include "math.h"
#include "utils.h"
#include "time.h"

#include <mpi.h>


#define COMP 1
#define LOG_BIN

#define NBINS_T 9
#define NBINS_L 25

struct Triplet{
  short indx;
  short indy;
  short indz;
};


int rank = 0, numtasks;

void shuffle_triplets(struct Triplet *trips, long size, long int idum);

static double ran2(long int *idum);

/*------------------------------------------------------------------------------
 *  This file calculates the structure function of the magnetic field.
 *  It does so by creating a list of triplets (all points on the magnetic field 
 *  grid) and shuffles it randomly. By using this as the top loop, it gradually 
 *  builds the ensemble average, and can be trivially parallelized in a logical
 *  way.
 *  
 *  While a total calculation of the full ensemble average is unfeasable 
 *  (quadratic in Ncell), the shuffle should hopefully provide enough statistics.
 *
 */
int main(int argc, char* argv[]){
  struct Fields fld;
  struct Moments mom;
  struct Parameters p;
  struct Triplet *trips;
  
  char fname[50];

  double histo[NBINS_L*NBINS_T]; 
  double histo_par[NBINS_L*NBINS_T]; 
  double histo_perp[NBINS_L*NBINS_T]; 

  double rhisto[NBINS_L*NBINS_T]; 
  double rhisto_par[NBINS_L*NBINS_T]; 
  double rhisto_perp[NBINS_L*NBINS_T]; 
  
  double histou[NBINS_L*NBINS_T]; 
  double histou_par[NBINS_L*NBINS_T]; 
  double histou_perp[NBINS_L*NBINS_T]; 

  double rhistou[NBINS_L*NBINS_T]; 
  double rhistou_par[NBINS_L*NBINS_T]; 
  double rhistou_perp[NBINS_L*NBINS_T]; 

  double comp_t[NBINS_T];
  double comp_l[NBINS_L];
  long num, nlim, ncell, counts_tot;
  long int idum;

  int i,j,k;
  int dind_x, dind_y, dind_z;
  int nx,ny,nz,indx,indy,indz;
  int li, ti;

  double dx,dy,dz, lx, ly, lz, lmax,lmin,opp; 
  float ***l; //precompute lengths (saves a sqrt() call)
  double lenx, leny, lenz;

  double l_curr, l_par;

  double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12;

  double Bavg_x, Bavg_y, Bavg_z, Bavg_m, inv_Bavg_m;
  double theta,dl,dtheta,idtheta,idl;

  double B1x,B1y,B1z,B2x,B2y,B2z;
  double dBx,dBy,dBz,dBsqr;
  double dB_par, dBsqr_par, dBsqr_perp;

  double U1x,U1y,U1z,U2x,U2y,U2z;
  double dUx,dUy,dUz,dUsqr;
  double dU_par, dUsqr_par, dUsqr_perp;

  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == 0){
    if(argc != 3 && argc != 4){
      printf("Usage: %s [ensembles] [fld_vtk] [mom_vtk optional]\n",argv[0]);
      exit(0);
    }

    if(argc == 4){
      read_mom(argv[3],0,&p,&mom);
      free_contiguous_3dArray(mom.n);
      mom.n = NULL;
    }

    nlim=atoi(argv[1]);
    read_fld(argv[2],0,&p,&fld);
    idum = (long int) time(NULL);
  } 

  MPI_Bcast(&p.nx,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.ny,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.nz,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.dx,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.dy,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&p.dz,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nlim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&idum,1,MPI_LONG,0,MPI_COMM_WORLD); //make sure use the same random number

  nx = p.nx; ny = p.ny; nz = p.nz;
  dx = p.dx; dy = p.dy; dz = p.dz;

  if(rank != 0){
    fld.b1   = new_contiguous_3dArray(nz,ny,nx);
    fld.b2   = new_contiguous_3dArray(nz,ny,nx);
    fld.b3   = new_contiguous_3dArray(nz,ny,nx);
    if(argc == 4){
      mom.u1   = new_contiguous_3dArray(nz,ny,nx);
      mom.u2   = new_contiguous_3dArray(nz,ny,nx);
      mom.u3   = new_contiguous_3dArray(nz,ny,nx);
    }
  }
  

  MPI_Bcast(&(fld.b1[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(fld.b2[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&(fld.b3[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  if(argc==4){
    MPI_Bcast(&(mom.u1[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&(mom.u2[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&(mom.u3[0][0][0]),nx*ny*nz,MPI_FLOAT,0,MPI_COMM_WORLD);
  }

  ncell = nx*ny*nz;
  if(ncell % numtasks != 0){
    printf("Please choose a divisor of %ld for the number of processes.\n",ncell);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(1);
  }

 
  lenx = dx*nx;
  leny = dy*ny;
  lenz = dz*nz;

  lmax = (0.5+0.0001)*sqrt(lenx*lenx + leny*leny + lenz*lenz);
  //lmax = (0.5 + 0.0001)*lenx;
  lmin = (1.0-0.0001)*dx;

#ifdef LOG_BIN
  double llmin = log10(lmin);
  dl = (log10(lmax) - llmin)/NBINS_L;
  idl = 1.0/dl;
#else
  double llmin = 0;
  dl = (lmax-lmin)/NBINS_L;   
  idl = 1.0/dl;
#endif

  idtheta = 2.0*NBINS_T/M_PI; // 0 to pi/2
  dtheta = 1.0/idtheta;
  

  l = new_contiguous_3dArray(nz,ny,nx);

  trips = malloc(sizeof(struct Triplet)*nx*ny*nz);

  // contruct triplets and precompute l
  for(k = 0; k < nz; k++){
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        struct Triplet t;
        t.indx = i;
        t.indy = j;
        t.indz = k;
        trips[k*(ny*nx)+j*nx+i]=t;
        
        lx = i > nx/2 ? dx*(i-nx) : dx*i;
        ly = j > ny/2 ? dy*(j-ny) : dy*j;
        lz = k > nz/2 ? dz*(k-nz) : dz*k;
        
        l[k][j][i] = sqrt(lx*lx + ly*ly + lz*lz);
      }
    }
  }

  //random sort here
  shuffle_triplets(trips,ncell,idum);

  if(COMP){
    for(i = 0; i < NBINS_T; i++) comp_t[i] = 1.0/(cos(i*dtheta) - cos((i+1)*dtheta));//average over solid angle
    for(i = 0; i < NBINS_L; i++) comp_l[i] = 0.0;
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
    //for(k = 0; k < 1; k++){
    //  for(j = 0; j < 1; j++){
        for(i = 0; i < nx; i++){
          if( i == 0 && j == 0 && k == 0) continue;
          l_curr = l[k][j][i];
#ifdef LOG_BIN
          li = (int)((log10(l_curr)-llmin)*idl);
          if(li == NBINS_L) li--;
          comp_l[li]++;
#else
          li = (int)((l_curr-lmin)*idl);
          if(li == NBINS_L) li--;
          comp_l[li]++;
#endif
        }
      }
    }
    for(i = 0; i < NBINS_L; i++) comp_l[i] = comp_l[i] == 0.0 ? 1.0 : 1.0/comp_l[i];
  } else {
    for(i = 0; i < NBINS_T; i++) comp_t[i] = 1.0;
    for(i = 0; i < NBINS_L; i++) comp_l[i] = 1.0;
  }



  for(i =0; i < NBINS_L; i++){
    for(j =0; j < NBINS_T; j++){
      histo[i*NBINS_T +j]   = 0.0;
      histo_par[i*NBINS_T +j]  = 0.0;
      histo_perp[i*NBINS_T +j] = 0.0;
      histou[i*NBINS_T +j]   = 0.0;
      histou_par[i*NBINS_T +j]  = 0.0;
      histou_perp[i*NBINS_T +j] = 0.0;
      rhistou[i*NBINS_T +j]   = 0.0;
      rhistou_par[i*NBINS_T +j]  = 0.0;
      rhistou_perp[i*NBINS_T +j] = 0.0;
    }
  }  
  counts_tot = 0;

  //main loop
  long offset = (ncell/numtasks)*rank;
  for(num = 0; num <  ncell/numtasks && num < nlim; num++){
    indx = trips[num + offset].indx;
    indy = trips[num + offset].indy;
    indz = trips[num + offset].indz;

    B1x = fld.b1[indz][indy][indx];
    B1y = fld.b2[indz][indy][indx];
    B1z = fld.b3[indz][indy][indx];

    if(argc == 4){
      U1x = mom.u1[indz][indy][indx];
      U1y = mom.u2[indz][indy][indx];
      U1z = mom.u3[indz][indy][indx];
    }
        
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
    //for(k = indz; k < indz+1; k++){
    //  for(j = indy; j < indy+1; j++){
        for(i = 0; i < nx; i++){
          if(i == indx && j == indy && k == indz) continue;

          B2x = fld.b1[k][j][i];
          B2y = fld.b2[k][j][i];
          B2z = fld.b3[k][j][i];

          Bavg_x = 0.5*(B1x + B2x);
          Bavg_y = 0.5*(B1y + B2y);
          Bavg_z = 0.5*(B1z + B2z);
          Bavg_m = sqrt(Bavg_x*Bavg_x + Bavg_y*Bavg_y + Bavg_z*Bavg_z);
          inv_Bavg_m = 1.0 / Bavg_m;
      
          dind_x = i - indx;    
          dind_y = j - indy;    
          dind_z = k - indz;    

          if(dind_x < 0) dind_x += nx;
          if(dind_y < 0) dind_y += ny;
          if(dind_z < 0) dind_z += nz;

          lx = dind_x > nx/2 ? dx*(dind_x-nx) : dx*dind_x;
          ly = dind_y > ny/2 ? dy*(dind_y-ny) : dy*dind_y;
          lz = dind_z > nz/2 ? dz*(dind_z-nz) : dz*dind_z;

          l_curr = (double) l[dind_z][dind_y][dind_x];
          //l_curr = sqrt(lx*lx + ly*ly + lz*lz);
          l_par = inv_Bavg_m*(lx*Bavg_x + ly*Bavg_y + lz*Bavg_z);
  
          if(l_par < 0) l_par = -l_par;
          opp = l_par/l_curr;
          if(opp > 1.0){//round-off?
            if(opp > 1.1){//much too large!
              printf("Arccos error with %f %f %f %f %f.\n", l_curr,sqrt(lx*lx+ly*ly+lz*lz),inv_Bavg_m, l_par, opp-1.0);
              printf("error with %f %f %f %d %d %d.\n", lx, ly, lz, dind_x, dind_y, dind_z);
              exit(1);
            }
            opp = 1.0;
          } 
          theta = acos(opp);
#ifdef LOG_BIN
          li = (int)((log10(l_curr)-llmin)*idl);
#else
          li = (int)((l_curr-lmin)*idl);
#endif
          ti = (int)(theta*idtheta);

          if(ti >= NBINS_T || ti < 0){
            if(ti == NBINS_T )
              ti--;
            else
              printf("T binning error %d %d %e %e %e %e %e\n",li,ti,l_par, l_curr,idl,theta,idtheta);
          }

          if(li >= NBINS_L || li < 0){
            if(li == NBINS_L )
              li--;
            else
              printf("L binning error %d %d %e %e %e %e %e\n",li,ti,l_par, l_curr,idl,theta,idtheta);
          }

          dBx = B2x - B1x;
          dBy = B2y - B1y;
          dBz = B2z - B1z;
    
          dBsqr = dBx*dBx + dBy*dBy + dBz*dBz;
          dB_par = inv_Bavg_m*(dBx*Bavg_x + dBy*Bavg_y + dBz*Bavg_z);
          dBsqr_par = dB_par*dB_par;
          dBsqr_perp = dBsqr - dBsqr_par;    

          histo[li*NBINS_T + ti]      += comp_l[li]*comp_t[ti]*dBsqr;
          histo_par[li*NBINS_T + ti]  += comp_l[li]*comp_t[ti]*dBsqr_par;
          histo_perp[li*NBINS_T + ti] += comp_l[li]*comp_t[ti]*dBsqr_perp * 0.5;

          if(argc == 4){
            U2x = mom.u1[k][j][i];
            U2y = mom.u2[k][j][i];
            U2z = mom.u3[k][j][i];

            dUx = U2x - U1x;
            dUy = U2y - U1y;
            dUz = U2z - U1z;
    
            dUsqr = dUx*dUx + dUy*dUy + dUz*dUz;
            dU_par = inv_Bavg_m*(dUx*Bavg_x + dUy*Bavg_y + dUz*Bavg_z);
            dUsqr_par = dU_par*dU_par;
            dUsqr_perp = dUsqr - dUsqr_par;    

            histou[li*NBINS_T + ti]      += comp_l[li]*comp_t[ti]*dUsqr;
            histou_par[li*NBINS_T + ti]  += comp_l[li]*comp_t[ti]*dUsqr_par;
            histou_perp[li*NBINS_T + ti] += comp_l[li]*comp_t[ti]*dUsqr_perp * 0.5;
          }

          counts_tot++;
        }
      }
    } 
    if(num % 100 == 0 && rank == 0){
      printf("Step %ld\n",num);
    }
    if((num-200)% 1000 == 0 || num == nlim - 1){

      MPI_Reduce(histo,rhisto,NBINS_T*NBINS_L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(histo_par,rhisto_par,NBINS_T*NBINS_L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(histo_perp,rhisto_perp,NBINS_T*NBINS_L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      if(argc == 4){
        MPI_Reduce(histou,rhistou,NBINS_T*NBINS_L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(histou_par,rhistou_par,NBINS_T*NBINS_L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(histou_perp,rhistou_perp,NBINS_T*NBINS_L,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      }

      if(rank==0){
        printf("Dumping histogram at step %ld with %ld points\n",num, counts_tot*numtasks);
    
        float inum = 1.0/((double)(num*numtasks));
        sprintf(fname,"histo_%.0f.dat",p.time);
        FILE* fp = fopen(fname,"w");

        fprintf(fp,"[1] l  [2] dB  [3] dB_par [4] dB_perp  [5] dU  [6] dU_par  [7] dU_perp\n");
        for(i =0; i < NBINS_L; i++){
          for(j =0; j < NBINS_T; j++){
            h1 = rhisto[i*NBINS_T + j]*inum;
            h2 = rhisto_par[i*NBINS_T + j]*inum;
            h3 = rhisto_perp[i*NBINS_T + j]*inum;
            h4 = rhistou[i*NBINS_T + j]*inum;
            h5 = rhistou_par[i*NBINS_T + j]*inum;
            h6 = rhistou_perp[i*NBINS_T + j]*inum;
#ifdef LOG_BIN
            fprintf(fp,"%le %le %le %le %le %le %le %le\n",i*dl+llmin,j*dtheta,h1,h2,h3,h4,h5,h6);
#else
            fprintf(fp,"%le %le %le %le %le %le %le %le\n",i*dl+llmin,j*dtheta,h1,h2,h3,h4,h5,h6);
#endif

          }
          fprintf(fp,"\n");
        } 
        fclose(fp);

        sprintf(fname,"dump_%.0f.dat",p.time);
        fp = fopen(fname,"w");
        fprintf(fp,"[1] l  [2] k [3] dB(k_par) [4] dB(k_perp) [5] dB_par(k_par) [6] dB_par(k_perp) [7] dB_perp(k_par) [8] dB_perp(k_perp) [9] dU(k_par) [10] dU(k_perp) [11] dU_par(k_par) [12] dU_par(k_perp) [13] dU_perp(k_par) [14] dU_perp(k_perp)\n");  
        for(i =1; i < NBINS_L; i++){
          h1 = inum*rhisto[i*NBINS_T];
          h2 = inum*rhisto[i*NBINS_T + NBINS_T - 1];
          h3 = inum*rhisto_par[i*NBINS_T];
          h4 = inum*rhisto_par[i*NBINS_T + NBINS_T - 1];
          h5 = inum*rhisto_perp[i*NBINS_T];
          h6 = inum*rhisto_perp[i*NBINS_T + NBINS_T - 1];
          h7 = inum*rhistou[i*NBINS_T];
          h8 = inum*rhistou[i*NBINS_T + NBINS_T - 1];
          h9 = inum*rhistou_par[i*NBINS_T];
          h10= inum*rhistou_par[i*NBINS_T + NBINS_T - 1];
          h11= inum*rhistou_perp[i*NBINS_T];
          h12= inum*rhistou_perp[i*NBINS_T + NBINS_T - 1];
#ifdef LOG_BIN
          float l = pow(10,i*dl+llmin);
          fprintf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",l,2*M_PI/l,
                  h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12); 
#else
          fprintf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",i*dl+llmin,2*M_PI/(i*dl),
                  h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12); 
#endif

        }
        fclose(fp);
      }
    }
  } 

  freeFieldStruct(&fld);
  free_contiguous_3dArray(l);
  free(trips);
  
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
