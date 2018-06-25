#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "read_vtk.h"
#define SQR(x) ((x)*(x))

int main(int argc, char* argv[]){
  struct Fields fld;
  struct Parameters p;
  char fname[200];

  int dslice=1000000;

  int i,j,k;
  double mag;


  if(argc < 3){
    fprintf(stderr, "Usage: %s fld dslice\n",argv[0]);
    exit(1);
  }

  dslice = atoi(argv[2]);

  read_fld(argv[1],0,&p,&fld);

  for(k = 0; k < p.nz; k+=dslice){

    sprintf(fname,"Bslice_t%.0f_z%d.dat",p.time,k);
    FILE* fx   = fopen(fname,"w");
    fprintf(fx,"[1] x [2] y [3] b1 [4] b2 [5] b3 [6] bmag\n");


    for(j = 0; j < p.ny; j++){
      for(i = 0; i < p.nx; i++){
        mag = sqrt(SQR(fld.b1[k][j][i])+SQR(fld.b2[k][j][i])+SQR(fld.b3[k][j][i]));
        fprintf(fx,"%e %e %e %e %e %e \n", j*p.dy,i*p.dx,fld.b1[k][j][i], fld.b2[k][j][i],
                                                         fld.b3[k][j][i], mag);
      }
      fprintf(fx,"\n");
    }

    fclose(fx);
  }

  freeFieldStruct(&fld);

  return 0;
}
