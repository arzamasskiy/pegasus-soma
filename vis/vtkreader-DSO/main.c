#include <stdlib.h>
#include <stdio.h>
#include "read_vtk.h"

void main(int argc, char* argv[]){
  struct Fields fld;
  struct Moments mom;
  struct Parameters p;
  read_fld("combined.0200.fld.vtk",0,&p,&fld);
  read_mom("combined.0100.mom.vtk",0,&p,&mom);

  int i,j,k;

  freeFieldStruct(&fld);
  freeMomentStruct(&mom);
}
