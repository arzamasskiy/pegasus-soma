#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "read_vtk.h"
#include "math.h"
#include "utils.h"
#include "time.h"

#include <mpi.h>


#define COMP 0
#define LOG_BIN 1

#define NBINS_PAR  50
#define NBINS_PERP 50

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
 *  This version bins in l_par and l_perp, rather than l and theta
 *
 */
int main(int argc, char* argv[]){
  struct Fields fld;
  struct Moments mom;
  struct Parameters p;
  struct Triplet *trips;
  
  char fname[50];


  long  counts[NBINS_PAR*NBINS_PERP]; 
  double histo[NBINS_PAR*NBINS_PERP]; 
  double histo_par[NBINS_PAR*NBINS_PERP]; 
  double histo_perp[NBINS_PAR*NBINS_PERP]; 

  long   rcounts[NBINS_PAR*NBINS_PERP]; 
  double rhisto[NBINS_PAR*NBINS_PERP]; 
  double rhisto_par[NBINS_PAR*NBINS_PERP]; 
  double rhisto_perp[NBINS_PAR*NBINS_PERP]; 
  
  double histou[NBINS_PAR*NBINS_PERP]; 
  double histou_par[NBINS_PAR*NBINS_PERP]; 
  double histou_perp[NBINS_PAR*NBINS_PERP]; 

  double rhistou[NBINS_PAR*NBINS_PERP]; 
  double rhistou_par[NBINS_PAR*NBINS_PERP]; 
  double rhistou_perp[NBINS_PAR*NBINS_PERP]; 

  double comp_par[NBINS_PAR];
  double comp_perp[NBINS_PERP];

  double lpar_widths[NBINS_PAR];
  double lperp_widths[NBINS_PERP];

  long num,tnum,kmod,kmod2, nlim, ncell, counts_tot;
  long int idum;

  int i,j,k,i2,j2,k2;
  int dind_x, dind_y, dind_z;
  int nx,ny,nz,indx,indy,indz;
  double parf, perpf;
  int pari, perpi;

  double dx,dy,dz, lx, ly, lz, lmax,lmin,opp; 
  double lpx,lpy,lpz; //l_perp vector
  float ***l; //precompute lengths (saves a sqrt() call)
  double lenx, leny, lenz;
  int remain,extra,off;

  double l_curr, l_par, l_perp;

  double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,cc,ccp;

  double Bavg_x, Bavg_y, Bavg_z, Bavg_m, inv_Bavg_m;
  double theta,dl_par,dl_perp,idl_par,idl_perp;

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

  printf("%d %d %d %e %e %e %e %e %e\n",nx,ny,nz,dx,dy,dz,nx*dx,ny*dy,nz*dz);

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

  kmod = (long)ceil(pow(nlim*numtasks,1.0/3.0));
  kmod2 =kmod*kmod;

 
  lenx = dx*nx;
  leny = dy*ny;
  lenz = dz*nz;

  lmax = (0.5+0.0001)*sqrt(lenx*lenx + leny*leny + lenz*lenz);
#if LOG_BIN
  lmin = dx;
#else
  lmin = 0;
#endif

#if LOG_BIN
  double llmin = log10(lmin);
  dl_perp = log10(lmax/lmin)/NBINS_PERP;
  dl_par  = log10(lmax/lmin)/NBINS_PAR;
  idl_par = 1.0/dl_par;
  idl_perp = 1.0/dl_perp;
  for(i = 0; i < NBINS_PAR ; i++)  lpar_widths[i] = pow(10,dl_par*(i+1)  + llmin) - pow(10,dl_par*i +llmin);
  for(i = 0; i < NBINS_PERP; i++) lperp_widths[i] = pow(10,dl_perp*(i+1) + llmin) - pow(10,dl_perp*i+llmin);
#else
  double llmin = 0;
  dl_par = (lmax-lmin)/NBINS_PAR;   
  dl_perp= (lmax-lmin)/NBINS_PERP;   
  idl_par = 1.0/dl_par;
  idl_perp= 1.0/dl_perp;
  for(i = 0; i < NBINS_PAR; i++)   lpar_widths[i] = dl_par;
  for(i = 0; i < NBINS_PERP; i++) lperp_widths[i] = dl_perp;
#endif

  

  l = new_contiguous_3dArray(nz,ny,nx);

  trips = malloc(sizeof(struct Triplet)*nx*ny*nz);

  double brms=0;
  // contruct triplets and precompute l
  for(k = 0; k < nz; k++){
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        struct Triplet t;
        t.indx = i;
        t.indy = j;
        t.indz = k;
        trips[k*(ny*nx)+j*nx+i]=t;
        
        //lx = dx*(i-nx);
        //ly = dy*(j-ny);
        //lz = dz*(k-nz);
        lx = i > nx/2 ? dx*(i-nx) : dx*i;
        ly = j > ny/2 ? dy*(j-ny) : dy*j;
        lz = k > nz/2 ? dz*(k-nz) : dz*k;
        
        brms +=  fld.b1[k][j][i]*fld.b1[k][j][i]
               + fld.b2[k][j][i]*fld.b2[k][j][i]
               + fld.b3[k][j][i]*fld.b3[k][j][i];
        l[k][j][i] = sqrt(lx*lx + ly*ly + lz*lz);
      }
    }
  }
  brms = sqrt(brms/((double)(nx*ny*nz)));

  //random sort here
  shuffle_triplets(trips,ncell,idum);

  for(i = 0; i < NBINS_PERP; i++) comp_perp[i] = 1.0;
  for(i = 0; i < NBINS_PAR; i++) comp_par[i] = 1.0;



  for(i =0; i < NBINS_PAR; i++){
    for(j =0; j < NBINS_PERP; j++){
      counts[i*NBINS_PERP +j]   = 0;
      histo[i*NBINS_PERP +j]   = 0.0;
      histo_par[i*NBINS_PERP +j]  = 0.0;
      histo_perp[i*NBINS_PERP +j] = 0.0;
      histou[i*NBINS_PERP +j]   = 0.0;
      histou_par[i*NBINS_PERP +j]  = 0.0;
      histou_perp[i*NBINS_PERP +j] = 0.0;
      rhistou[i*NBINS_PERP +j]   = 0.0;
      rhistou_par[i*NBINS_PERP +j]  = 0.0;
      rhistou_perp[i*NBINS_PERP +j] = 0.0;
    }
  }  
  counts_tot = 0;

  //main loop
  remain = ncell % numtasks;
  extra = rank < remain ? 1 : 0;
  off = rank < extra ? rank : extra;
  long offset = (ncell/numtasks)*rank + off;
  for(num = 0; num <  (ncell/numtasks+extra) && num < nlim; num++){
    indx = trips[num + offset].indx;
    indy = trips[num + offset].indy;
    indz = trips[num + offset].indz;

/*
 *  LOOP OVER POINTS AND FIND CORRELATIONS
 */  

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

          //lx = dx*dind_x;
          //ly = dy*dind_y;
          //lz = dz*dind_z;

          if(dind_x < 0) dind_x += nx;
          if(dind_y < 0) dind_y += ny;
          if(dind_z < 0) dind_z += nz;
          lx = dind_x > nx/2 ? dx*(dind_x-nx) : dx*dind_x;
          ly = dind_y > ny/2 ? dy*(dind_y-ny) : dy*dind_y;
          lz = dind_z > nz/2 ? dz*(dind_z-nz) : dz*dind_z;

          //l_curr = (double) l[dind_z][dind_y][dind_x];
          l_curr = sqrt(lx*lx + ly*ly + lz*lz);
          l_par = inv_Bavg_m*(lx*Bavg_x + ly*Bavg_y + lz*Bavg_z);
          lpx = (lx - l_par*Bavg_x*inv_Bavg_m);
          lpy = (ly - l_par*Bavg_y*inv_Bavg_m);
          lpz = (lz - l_par*Bavg_z*inv_Bavg_m);
          l_perp = sqrt(lpx*lpx + lpy*lpy + lpz*lpz);

  
          if(l_par < 0) l_par = -l_par;


#if LOG_BIN
          parf  = (log10(l_par)-llmin)*idl_par;
          perpf = (log10(l_perp)-llmin)*idl_perp;
          if(parf < 0 || perpf < 0) continue;
          pari  = (int)parf;
          perpi = (int)perpf;
#else
          pari = (int)((l_par-lmin)*idl_par);
          perpi = (int)((l_perp-lmin)*idl_perp);
#endif


          if(pari >= NBINS_PAR || pari < 0){
#if LOG_BIN
            continue;
#else
            if(pari == NBINS_PAR )
              pari--;
            else{
              printf("Lpar binning error %d %d %d %d %d %e %e %e %e %e\n",i,j,k,pari,perpi,l_par, l_perp,idl_par,idl_perp,l_curr);
              return -1;
            }
#endif
          }
          if(perpi >= NBINS_PERP || perpi < 0){
#if LOG_BIN
            continue;
#else
            if(perpi == NBINS_PERP )
              perpi--;
            else{
              printf("Lperp binning error %d %d %d %d %d %e %e %e %e %e\n",i,j,k,pari,perpi,l_par, l_perp,idl_par,idl_perp,l_curr);
              return -1;
            }
#endif
          }

          dBx = B2x - B1x;
          dBy = B2y - B1y;
          dBz = B2z - B1z;
    
          dBsqr = dBx*dBx + dBy*dBy + dBz*dBz;
          dB_par = inv_Bavg_m*(dBx*Bavg_x + dBy*Bavg_y + dBz*Bavg_z);
          dBsqr_par = dB_par*dB_par;
          dBsqr_perp = dBsqr - dBsqr_par;    

          counts[pari*NBINS_PERP + perpi]++;
          histo[pari*NBINS_PERP + perpi]      += comp_par[pari]*comp_perp[perpi]*dBsqr;
          histo_par[pari*NBINS_PERP + perpi]  += comp_par[pari]*comp_perp[perpi]*dBsqr_par;
          histo_perp[pari*NBINS_PERP + perpi] += comp_par[pari]*comp_perp[perpi]*dBsqr_perp * 0.5;

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

            histou[pari*NBINS_PERP + perpi]      += comp_par[pari]*comp_perp[perpi]*dUsqr;
            histou_par[pari*NBINS_PERP + perpi]  += comp_par[pari]*comp_perp[perpi]*dUsqr_par;
            histou_perp[pari*NBINS_PERP + perpi] += comp_par[pari]*comp_perp[perpi]*dUsqr_perp * 0.5;
          }

          counts_tot++;
        }
      }
    } 

/*
 * LOOP OVER LENGTHS (for better statistics ;) ) 
 */

    tnum= num + nlim*rank;
    indx=tnum/kmod2;
    indy=(tnum%kmod2)/kmod;
    indz=tnum%kmod;

    lx = indx > nx/2 ? dx*(indx-nx): dx*indx;    
    ly = indy > ny/2 ? dy*(indy-ny): dy*indy;    
    lz = indz > nz/2 ? dz*(indz-nz): dz*indz;    
    l_curr = sqrt(lx*lx + ly*ly + lz*lz);
        
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          if(indx == 0 && indy == 0 && indz == 0) continue;

          
          i2 = i + indx;
          j2 = j + indy;
          k2 = k + indz;

          if(i2 >= nx) i2 -= nx;
          if(j2 >= ny) j2 -= ny;
          if(k2 >= nz) k2 -= nz;

          B1x = fld.b1[k][j][i];
          B1y = fld.b2[k][j][i];
          B1z = fld.b3[k][j][i];
          B2x = fld.b1[k2][j2][i2];
          B2y = fld.b2[k2][j2][i2];
          B2z = fld.b3[k2][j2][i2];

          if(argc == 4){
            U1x = mom.u1[k][j][i];
            U1y = mom.u2[k][j][i];
            U1z = mom.u3[k][j][i];
            U2x = mom.u1[k2][j2][i2];
            U2y = mom.u2[k2][j2][i2];
            U2z = mom.u3[k2][j2][i2];
          }


          Bavg_x = 0.5*(B1x + B2x);
          Bavg_y = 0.5*(B1y + B2y);
          Bavg_z = 0.5*(B1z + B2z);
          Bavg_m = sqrt(Bavg_x*Bavg_x + Bavg_y*Bavg_y + Bavg_z*Bavg_z);
          inv_Bavg_m = 1.0 / Bavg_m;
      

          //l_curr = (double) l[dind_z][dind_y][dind_x];
          l_par = inv_Bavg_m*(lx*Bavg_x + ly*Bavg_y + lz*Bavg_z);
          lpx = (lx - l_par*Bavg_x*inv_Bavg_m);
          lpy = (ly - l_par*Bavg_y*inv_Bavg_m);
          lpz = (lz - l_par*Bavg_z*inv_Bavg_m);
          l_perp = sqrt(lpx*lpx + lpy*lpy + lpz*lpz);

  
          if(l_par < 0) l_par = -l_par;


#if LOG_BIN
          parf  = (log10(l_par)-llmin)*idl_par;
          perpf = (log10(l_perp)-llmin)*idl_perp;
          if(parf < 0 || perpf < 0) continue;
          pari  = (int)parf;
          perpi = (int)perpf;
#else
          pari = (int)((l_par-lmin)*idl_par);
          perpi = (int)((l_perp-lmin)*idl_perp);
#endif


          if(pari >= NBINS_PAR || pari < 0){
#if LOG_BIN
            continue;
#else
            if(pari == NBINS_PAR )
              pari--;
            else{
              printf("Lpar binning error %d %d %d %d %d %e %e %e %e %e\n",i,j,k,pari,perpi,l_par, l_perp,idl_par,idl_perp,l_curr);
              return -1;
            }
#endif
          }
          if(perpi >= NBINS_PERP || perpi < 0){
#if LOG_BIN
            continue;
#else
            if(perpi == NBINS_PERP )
              perpi--;
            else{
              printf("Lperp binning error %d %d %d %d %d %e %e %e %e %e\n",i,j,k,pari,perpi,l_par, l_perp,idl_par,idl_perp,l_curr);
              return -1;
            }
#endif
          }

          dBx = B2x - B1x;
          dBy = B2y - B1y;
          dBz = B2z - B1z;
    
          dBsqr = dBx*dBx + dBy*dBy + dBz*dBz;
          dB_par = inv_Bavg_m*(dBx*Bavg_x + dBy*Bavg_y + dBz*Bavg_z);
          dBsqr_par = dB_par*dB_par;
          dBsqr_perp = dBsqr - dBsqr_par;    

          counts[pari*NBINS_PERP + perpi]++;
          histo[pari*NBINS_PERP + perpi]      += comp_par[pari]*comp_perp[perpi]*dBsqr;
          histo_par[pari*NBINS_PERP + perpi]  += comp_par[pari]*comp_perp[perpi]*dBsqr_par;
          histo_perp[pari*NBINS_PERP + perpi] += comp_par[pari]*comp_perp[perpi]*dBsqr_perp * 0.5;

          if(argc == 4){

            dUx = U2x - U1x;
            dUy = U2y - U1y;
            dUz = U2z - U1z;
    
            dUsqr = dUx*dUx + dUy*dUy + dUz*dUz;
            dU_par = inv_Bavg_m*(dUx*Bavg_x + dUy*Bavg_y + dUz*Bavg_z);
            dUsqr_par = dU_par*dU_par;
            dUsqr_perp = dUsqr - dUsqr_par;    

            histou[pari*NBINS_PERP + perpi]      += comp_par[pari]*comp_perp[perpi]*dUsqr;
            histou_par[pari*NBINS_PERP + perpi]  += comp_par[pari]*comp_perp[perpi]*dUsqr_par;
            histou_perp[pari*NBINS_PERP + perpi] += comp_par[pari]*comp_perp[perpi]*dUsqr_perp * 0.5;
          }

          counts_tot++;
        }
      }
    }

    printf("\rStep %ld",num);
    fflush(stdout);
    if((num-200)% 1000 == 0 || num == nlim - 1){

      MPI_Reduce(counts,rcounts,NBINS_PERP*NBINS_PAR,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(histo,rhisto,NBINS_PERP*NBINS_PAR,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(histo_par,rhisto_par,NBINS_PERP*NBINS_PAR,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(histo_perp,rhisto_perp,NBINS_PERP*NBINS_PAR,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      if(argc == 4){
        MPI_Reduce(histou,rhistou,NBINS_PERP*NBINS_PAR,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(histou_par,rhistou_par,NBINS_PERP*NBINS_PAR,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(histou_perp,rhistou_perp,NBINS_PERP*NBINS_PAR,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      }

      if(rank==0){
        printf("Dumping histogram at step %ld with %ld points\n",num, counts_tot*numtasks);
    
        float inum = 1.0/((double)(num*numtasks));
        sprintf(fname,"histo_%.0f.dat",p.time);
        FILE* fp = fopen(fname,"w");

        fprintf(fp,"#Histogram values are squared. Brms=%e Count=%e\n",brms,inum);
        fprintf(fp,"#[1] l_par  [2] l_perp [3] dB  [4] dB_par [5] dB_perp  [6] dU  [7] dU_par  [8] dU_perp\n");
        for(i =0; i < NBINS_PAR; i++){
          for(j =0; j < NBINS_PERP; j++){
            cc = 1.0/rcounts[i*NBINS_PERP + j];
            h1 = rhisto[i*NBINS_PERP + j]*cc;
            h2 = rhisto_par[i*NBINS_PERP + j]*cc;
            h3 = rhisto_perp[i*NBINS_PERP + j]*cc;
            h4 = rhistou[i*NBINS_PERP + j]*cc;
            h5 = rhistou_par[i*NBINS_PERP + j]*cc;
            h6 = rhistou_perp[i*NBINS_PERP + j]*cc;
            fprintf(fp,"%le %le %le %le %le %le %le %le %ld %le\n",i*dl_par+llmin,j*dl_perp+llmin,
                      h1,h2,h3,h4,h5,h6,
                      rcounts[i*NBINS_PERP + j],lpar_widths[i]*lperp_widths[j]);

          }
          fprintf(fp,"\n");
        } 
        fclose(fp);

        
        sprintf(fname,"hslice_%.0f.dat",p.time);
        fp = fopen(fname,"w");
        fprintf(fp,"[1] l_par  [2] l_perp [3] dB(k_par) [4] dB(k_perp) [5] dB_par(k_par) [6] dB_par(k_perp) [7] dB_perp(k_par) [8] dB_perp(k_perp) [9] dU(k_par) [10] dU(k_perp) [11] dU_par(k_par) [12] dU_par(k_perp) [13] dU_perp(k_par) [14] dU_perp(k_perp)\n");  
        for(i =1; i < NBINS_PAR; i++){
          cc = 1.0/rcounts[i*NBINS_PERP];
          ccp = 1.0/rcounts[i];
          h1 = cc*rhisto[i*NBINS_PERP];
          h2 = ccp*rhisto[i];
          h3 = cc*rhisto_par[i*NBINS_PERP];
          h4 = ccp*rhisto_par[i];
          h5 = cc*rhisto_perp[i*NBINS_PERP];
          h6 = ccp*rhisto_perp[i];
          h7 = cc*rhistou[i*NBINS_PERP];
          h8 = ccp*rhistou[i];
          h9 = cc*rhistou_par[i*NBINS_PERP];
          h10= ccp*rhistou_par[i];
          h11= cc*rhistou_perp[i*NBINS_PERP];
          h12= ccp*rhistou_perp[i];
#if LOG_BIN
          fprintf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",pow(10,i*dl_par+llmin),pow(10,i*dl_perp+llmin),
                  h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12); 
#else
          fprintf(fp,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",i*dl_par+llmin,i*dl_perp+llmin,
                  h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12); 
#endif

        }
        fclose(fp);
        
      }
/*
      printf("Dumping histogram at step %ld with %ld points\n",num, counts_tot*numtasks);
    
      float inum = 1.0/((double)(num));
      sprintf(fname,"histo_id%d_%.0f.dat",rank,p.time);
      FILE* fp = fopen(fname,"w");

      fprintf(fp,"#Histogram values are squared. Brms=%e\n",brms);
      fprintf(fp,"#[1] l_par  [2] l_perp [3] dB  [4] dB_par [5] dB_perp  [6] dU  [7] dU_par  [8] dU_perp\n");
        for(i =0; i < NBINS_PAR; i++){
          for(j =0; j < NBINS_PERP; j++){
            h1 = histo[i*NBINS_PERP + j]*inum;
            h2 = histo_par[i*NBINS_PERP + j]*inum;
            h3 = histo_perp[i*NBINS_PERP + j]*inum;
            h4 = histou[i*NBINS_PERP + j]*inum;
            h5 = histou_par[i*NBINS_PERP + j]*inum;
            h6 = histou_perp[i*NBINS_PERP + j]*inum;
            fprintf(fp,"%le %le %le %le %le %le %le %le %ld %le\n",i*dl_par+llmin,j*dl_perp+llmin,
                      h1/(lpar_widths[i]*lperp_widths[j]),
                      h2/(lpar_widths[i]*lperp_widths[j]),
                      h3/(lpar_widths[i]*lperp_widths[j]),
                      h4/(lpar_widths[i]*lperp_widths[j]),
                      h5/(lpar_widths[i]*lperp_widths[j]),
                      h6/(lpar_widths[i]*lperp_widths[j]),
                      counts[i*NBINS_PERP + j],lpar_widths[i]*lperp_widths[j]);

        }
        fprintf(fp,"\n");
      } 
      fclose(fp);
      */
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
