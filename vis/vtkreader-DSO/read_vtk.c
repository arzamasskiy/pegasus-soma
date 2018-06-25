#include "read_vtk.h"
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "assert.h"

void parse_header(FILE *fp, struct Parameters *pOut);
void parse_fld(FILE *fp, const short readEfield, struct Parameters *pOut, struct Fields *fOut);
void parse_mom(FILE *fp, const short readPressure, struct Parameters* pOut, struct Moments *mOut);

float fix_endian(float in_number){
  static int doneTest = 0;
  static int shouldSwap = 0;

  if(!doneTest){
    int tmp1 = 1;
    unsigned char *tmp2 = (unsigned char *) &tmp1;
    if (*tmp2 != 0)
      shouldSwap = 1;
    doneTest = 1;
  }

  if(shouldSwap){
    unsigned char *bytes = (unsigned char*) &in_number;
    unsigned char tmp;
    tmp = bytes[0];    
    bytes[0] = bytes[3];
    bytes[3] = tmp;
    tmp = bytes[1];    
    bytes[1] = bytes[2];
    bytes[2] = tmp;
  }
  return(in_number);
}

void read_all(const char *fname, struct Parameters *pOut, 
              struct Fields *fOut, struct Moments *mOut){

  FILE *fp = fopen(fname,"r");
  assert(fp != NULL);

  parse_header(fp,pOut);
  parse_fld(fp,1,pOut,fOut);
  parse_mom(fp,1,pOut,mOut);
}

void read_fld(const char *fname, const short readEfield, struct Parameters *pOut, struct Fields *fOut){
  FILE *fp = fopen(fname,"r");
  assert(fp != NULL);

  parse_header(fp,pOut);
  parse_fld(fp,readEfield,pOut,fOut);
}

void read_mom(const char *fname, const short readPressure, struct Parameters* pOut, struct Moments *mOut){
  FILE *fp = fopen(fname,"r");
  assert(fp != NULL);

  parse_header(fp,pOut);
  parse_mom(fp,readPressure,pOut,mOut);

}


void parse_header(FILE *fp, struct Parameters *pOut){
  int nx,ny,nz;
  char buffer[200];
  fgets(buffer, 200, fp);  // # vtk Datafile
  fgets(buffer, 200, fp);  // # **** vars at time=...

  // Get time
  sscanf(buffer,"%*s %*s %*s %*s %f",&(pOut->time));

  fgets(buffer, 200, fp);  // # BINARY
  fgets(buffer, 200, fp);  // # DATASET STRUCTURED POINTS
  fgets(buffer, 200, fp);  // # DIMENSIONS %d %d %d \n
    
  // Get dimensions
  sscanf(buffer,"%*s %d %d %d",&nx,&ny,&nz);
  nx--; ny--; nz--;

  pOut->nx = nx;
  pOut->ny = ny;
  pOut->nz = nz;
  

  fgets(buffer, 200, fp);  // # ORIGIN
  // Get Origin
  sscanf(buffer,"%*s %f %f %f",&(pOut->x0),&(pOut->y0),&(pOut->z0));

  fgets(buffer, 200, fp);  // # SPACING
  // Get Spacing
  sscanf(buffer,"%*s %f %f %f",&(pOut->dx),&(pOut->dy),&(pOut->dz));

  fgets(buffer, 200, fp);  // # CELL_DATA

}

void parse_fld(FILE *fp, const short readEfield, struct Parameters *pOut, struct Fields *fOut){
  int nx,ny,nz;
  int i, j, k;
  float data;

  char buffer[200];
  
  nx = pOut->nx;
  ny = pOut->ny;
  nz = pOut->nz;

  fOut->b1   = new_contiguous_3dArray(nz,ny,nx);
  fOut->b2   = new_contiguous_3dArray(nz,ny,nx);
  fOut->b3   = new_contiguous_3dArray(nz,ny,nx);


  fgets(buffer, 200, fp);  // # VECTOR cell_centered_B
  while(buffer[0] != 'V') fgets(buffer, 200, fp); 

  // Read in file, accounting for endianess.
  for(k = 0; k < nz; k++){
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        fread(&data,sizeof(float),1,fp); 
        fOut->b1[k][j][i] = fix_endian(data);
        fread(&data,sizeof(float),1,fp); 
        fOut->b2[k][j][i] = fix_endian(data);
        fread(&data,sizeof(float),1,fp); 
        fOut->b3[k][j][i] = fix_endian(data);
      }
    }
  }

  if(readEfield){
    fOut->e1   = new_contiguous_3dArray(nz,ny,nx);
    fOut->e2   = new_contiguous_3dArray(nz,ny,nx);
    fOut->e3   = new_contiguous_3dArray(nz,ny,nx);

    fgets(buffer, 200, fp);  // # VECTOR cell_centered_B

    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          fOut->e1[k][j][i] = fix_endian(data);
          fread(&data,sizeof(float),1,fp); 
          fOut->e2[k][j][i] = fix_endian(data);
          fread(&data,sizeof(float),1,fp); 
          fOut->e3[k][j][i] = fix_endian(data);
        }
      }
    }
  } else {
    fOut->e1   = NULL;
    fOut->e2   = NULL;
    fOut->e3   = NULL;
  }
  fclose(fp);
}


void parse_mom(FILE *fp, const short readPressure, struct Parameters* pOut, struct Moments *mOut){
  int nx,ny,nz;
  int i, j, k;
  float data;

  char buffer[200];

  nx = pOut->nx;
  ny = pOut->ny;
  nz = pOut->nz;

  mOut->n   = new_contiguous_3dArray(nz,ny,nx);

  mOut->u1   = new_contiguous_3dArray(nz,ny,nx);
  mOut->u2   = new_contiguous_3dArray(nz,ny,nx);
  mOut->u3   = new_contiguous_3dArray(nz,ny,nx);

  float ninv; // to go from momentum to velocity

  // Read in file, accounting for endianess.
  fgets(buffer, 200, fp);  // # SCALARS particle_density float
  fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
  for(k = 0; k < nz; k++){
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        fread(&data,sizeof(float),1,fp); 
        mOut->n[k][j][i] = fix_endian(data);
      }
    }
  }
  fgets(buffer, 200, fp);  // # VECTORS particle_momentum float
  for(k = 0; k < nz; k++){
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        ninv = 1.0/mOut->n[k][j][i]; 
        fread(&data,sizeof(float),1,fp); 
        mOut->u1[k][j][i] = ninv*fix_endian(data);
        fread(&data,sizeof(float),1,fp); 
        mOut->u2[k][j][i] = ninv*fix_endian(data);
        fread(&data,sizeof(float),1,fp); 
        mOut->u3[k][j][i] = ninv*fix_endian(data);
      }
    }
  }

  if(readPressure){
    mOut->p11   = new_contiguous_3dArray(nz,ny,nx);
    mOut->p12   = new_contiguous_3dArray(nz,ny,nx);
    mOut->p13   = new_contiguous_3dArray(nz,ny,nx);
    mOut->p22   = new_contiguous_3dArray(nz,ny,nx);
    mOut->p23   = new_contiguous_3dArray(nz,ny,nx);
    mOut->p33   = new_contiguous_3dArray(nz,ny,nx);

    fgets(buffer, 200, fp);  // # SCALARS particle_pressure_11 float
    fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          mOut->p11[k][j][i] = fix_endian(data);
        }
      }
    }
    fgets(buffer, 200, fp);  // # SCALARS particle_pressure_12 float
    fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          mOut->p12[k][j][i] = fix_endian(data);
        }
      }
    }
    fgets(buffer, 200, fp);  // # SCALARS particle_pressure_13 float
    fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          mOut->p13[k][j][i] = fix_endian(data);
        }
      }
    }
    fgets(buffer, 200, fp);  // # SCALARS particle_pressure_22 float
    fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          mOut->p22[k][j][i] = fix_endian(data);
        }
      }
    }
    fgets(buffer, 200, fp);  // # SCALARS particle_pressure_23 float
    fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          mOut->p23[k][j][i] = fix_endian(data);
        }
      }
    }
    fgets(buffer, 200, fp);  // # SCALARS particle_pressure_33 float
    fgets(buffer, 200, fp);  // # LOOKUP_TABLE default
    for(k = 0; k < nz; k++){
      for(j = 0; j < ny; j++){
        for(i = 0; i < nx; i++){
          fread(&data,sizeof(float),1,fp); 
          mOut->p33[k][j][i] = fix_endian(data);
        }
      }
    }
  } else {
    mOut->p11   = NULL;
    mOut->p12   = NULL;
    mOut->p13   = NULL;
    mOut->p22   = NULL;
    mOut->p23   = NULL;
    mOut->p33   = NULL;
  }
  fclose(fp);
}




//-----------------------------------------------------------//
// Helper functions
//-----------------------------------------------------------//
void freeFieldStruct(struct Fields *fld){
  free_contiguous_3dArray(fld->b1);
  free_contiguous_3dArray(fld->b2);
  free_contiguous_3dArray(fld->b3);
  
  free_contiguous_3dArray(fld->e1);
  free_contiguous_3dArray(fld->e2);
  free_contiguous_3dArray(fld->e3);
}

void freeMomentStruct(struct Moments *mom){
  free_contiguous_3dArray(mom->n);
  free_contiguous_3dArray(mom->u1);
  free_contiguous_3dArray(mom->u2);
  free_contiguous_3dArray(mom->u3);

  free_contiguous_3dArray(mom->p11);
  free_contiguous_3dArray(mom->p12);
  free_contiguous_3dArray(mom->p13);
  free_contiguous_3dArray(mom->p22);
  free_contiguous_3dArray(mom->p23);
  free_contiguous_3dArray(mom->p33);
}
