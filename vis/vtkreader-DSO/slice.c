#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "read_vtk.h"
#define SQR(x) ((x)*(x))

int main(int argc, char* argv[]){
  struct Fields fld;
  struct Moments mom;
  struct Parameters p,p1;
  char fname[200];

  int cull = 10;
  int dslice=1000000;

  int i,j,k;
  double mag, bxhat,byhat,bzhat,umag,upar,uperpsq;
  double P, P_par, P_perp, Delt;


  if(argc < 4){
    fprintf(stderr, "Usage: %s fld mom dslice (cull)\n",argv[0]);
    exit(1);
  }

  dslice = atoi(argv[3]);
  if(argc == 5) cull = atoi(argv[4]);

  read_fld(argv[1],0,&p,&fld);
  read_mom(argv[2],1,&p1,&mom);
  assert(p.time == p1.time);

  for(k = 0; k < p.nz; k+=dslice){

    sprintf(fname,"Bslice_t%.0f_z%d.dat",p.time,k);
    FILE* fx   = fopen(fname,"w");
    fprintf(fx,"#[1] x [2] y [3] b1 [4] b2 [5] b3 [6] bmag [7] delt [8] u [9] upar [10] uperp\n");

    sprintf(fname,"Bcull_t%.0f_z%d.dat",p.time,k);
    FILE* fy   = fopen(fname,"w");
    fprintf(fy,"#[1] x [2] y [3] b1 [4] b2 [5] b3 [6] bmag [7] delt [8] u [9] upar [10] uperp\n");

    for(j = 0; j < p.ny; j++){
      for(i = 0; i < p.nx; i++){
        mag = sqrt(SQR(fld.b1[k][j][i])+SQR(fld.b2[k][j][i])+SQR(fld.b3[k][j][i]));
        bxhat = fld.b1[k][j][i]/mag;
        byhat = fld.b2[k][j][i]/mag;
        bzhat = fld.b3[k][j][i]/mag;

        umag = sqrt(SQR(mom.u1[k][j][i])+SQR(mom.u2[k][j][i])+SQR(mom.u3[k][j][i]));
        upar = bxhat*mom.u1[k][j][i] + byhat*mom.u2[k][j][i] + bzhat*mom.u3[k][j][i];
        uperpsq = SQR(mom.u1[k][j][i] - upar*bxhat) 
                + SQR(mom.u2[k][j][i] - upar*byhat) 
                + SQR(mom.u3[k][j][i] - upar*bzhat); 

        P = (mom.p11[k][j][i] + mom.p22[k][j][i] + mom.p33[k][j][i])/3.0;
        P_par = bxhat*(mom.p11[k][j][i]*bxhat + mom.p12[k][j][i]*byhat + mom.p13[k][j][i]*bzhat)
              + byhat*(mom.p12[k][j][i]*bxhat + mom.p22[k][j][i]*byhat + mom.p23[k][j][i]*bzhat)
              + bzhat*(mom.p13[k][j][i]*bxhat + mom.p23[k][j][i]*byhat + mom.p33[k][j][i]*bzhat);
        P_perp = 1.5*P - 0.5*P_par;
        Delt = P_perp/P_par - 1.0;


        fprintf(fx,"%e %e %e %e %e %e %e %e %e %e\n", i*p.dx,j*p.dy,fld.b1[k][j][i], fld.b2[k][j][i],
                                                         fld.b3[k][j][i], mag, Delt,umag,upar,sqrt(uperpsq));

        if((i % cull == 0) && (j % cull == 0)){
          fprintf(fy,"%e %e %e %e %e %e %e %e %e %e\n", i*p.dx,j*p.dy,fld.b1[k][j][i],
                                            fld.b2[k][j][i],fld.b3[k][j][i], mag, Delt,umag,upar,sqrt(uperpsq));
        }
      }
      fprintf(fx,"\n");
      fprintf(fy,"\n");
    }

    fclose(fx);
    fclose(fy);
  }

  freeMomentStruct(&mom);
  freeFieldStruct(&fld);

  return 0;
}
