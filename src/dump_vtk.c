#include "copyright.h"
/*============================================================================*/
/*! \file dump_vtk.c 
 *  \brief Function to write a dump in VTK "legacy" format.
 *
 * PURPOSE: Function to write a dump in VTK "legacy" format.  With SMR,
 *   dumps are made for all levels and domains, unless nlevel and ndomain are
 *   specified in <output> block.  Works for BOTH conserved and primitives.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - dump_vtk() - writes VTK dump (all variables).			      */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "pegasus.h"
#include "prototypes.h"
#include "particles/particle.h"

/*----------------------------------------------------------------------------*/
/*! \fn void dump_vtk(MeshS *pM, OutputS *pOut)
 *  \brief Writes VTK dump (all variables).				      */

void dump_vtk(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd;
  int big_end = peg_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3;

  /* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pGrid = pM->Domain[nl][nd].Grid;

        il = pGrid->is, iu = pGrid->ie;
        jl = pGrid->js, ju = pGrid->je;
        kl = pGrid->ks, ku = pGrid->ke;

#ifdef WRITE_GHOST_CELLS
        iu = pGrid->ie + nghost;
        il = pGrid->is - nghost;

        if(pGrid->Nx[1] > 1) {
          ju = pGrid->je + nghost;
          jl = pGrid->js - nghost;
        }

        if(pGrid->Nx[2] > 1) {
          ku = pGrid->ke + nghost;
          kl = pGrid->ks - nghost;
        }
#endif /* WRITE_GHOST_CELLS */

        ndata0 = iu-il+1;
        ndata1 = ju-jl+1;
        ndata2 = ku-kl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
        }
        if((fname = peg_fname(plev,pM->outfilename,plev,pdom,num_digit,
			      pOut->num,pOut->id,"vtk")) == NULL){
          peg_error("[dump_vtk]: Error constructing filename\n");
        }

        if((pfile = fopen(fname,"w")) == NULL){
          peg_error("[dump_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(3*ndata0*sizeof(float))) == NULL){
//        if((data = (float *)malloc(9*ndata0*sizeof(float))) == NULL){
          peg_error("[dump_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */
        
        if (strcmp(pOut->out,"all") == 0){
          fprintf(pfile,"ALL vars at time= %e, level= %i, domain= %i\n",
                  pGrid->time,nl,nd);
        } else if(strcmp(pOut->out,"field") == 0) {
          fprintf(pfile,"FIELD vars at time= %e, level= %i, domain= %i\n",
                  pGrid->time,nl,nd);
        } else if(strcmp(pOut->out,"moment") == 0) {
          fprintf(pfile,"MOMENT vars at time= %e, level= %i, domain= %i\n",
                  pGrid->time,nl,nd);
        }

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */

        fc_pos(pGrid, il, jl, kl, &x1, &x2, &x3);;

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
        if (pGrid->Nx[1] == 1) {
          fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,1,1);
        } else {
          if (pGrid->Nx[2] == 1) {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,1);
          } else {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,ku-kl+2);
          }
        }
        fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pGrid->dx1,pGrid->dx2,pGrid->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ju-jl+1)*(ku-kl+1));

	if ((strcmp(pOut->out,"all") == 0) || (strcmp(pOut->out,"field") == 0)) {

/* Write cell-centered B */

	fprintf(pfile,"\nVECTORS cell_centered_B float\n");
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              data[3*(i-il)  ] = (float)pGrid->U[k][j][i].B1c;
              data[3*(i-il)+1] = (float)pGrid->U[k][j][i].B2c;
              data[3*(i-il)+2] = (float)pGrid->U[k][j][i].B3c;
            }
            if(!big_end) peg_bswap(data,sizeof(float),3*(iu-il+1));
            fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
          }
        }

/* Write cell-centered E */
          
        fprintf(pfile,"\nVECTORS cell_centered_E float\n");
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              data[3*(i-il)  ] = (float)pGrid->U[k][j][i].E1c;
              data[3*(i-il)+1] = (float)pGrid->U[k][j][i].E2c;
              data[3*(i-il)+2] = (float)pGrid->U[k][j][i].E3c;
            }
            if(!big_end) peg_bswap(data,sizeof(float),3*(iu-il+1));
            fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
          }
        }

        }
     
              
/* Write binned particle grid: density, momentum, pressure tensor */

        if ((strcmp(pOut->out,"all") == 0) || (strcmp(pOut->out,"moment") == 0)) {
        
        if (pOut->out_pargrid) {
          fprintf(pfile,"\nSCALARS particle_density float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_d;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nVECTORS particle_momentum float\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[3*(i-il)] = pGrid->Coup[k][j][i].grid_M1;
                data[3*(i-il)+1] = pGrid->Coup[k][j][i].grid_M2;
                data[3*(i-il)+2] = pGrid->Coup[k][j][i].grid_M3;
              }
              if(!big_end) peg_bswap(data,sizeof(float),3*(iu-il+1));
              fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
            }
          }
          fprintf(pfile,"\nSCALARS particle_pressure_11 float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_p11;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nSCALARS particle_pressure_12 float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_p12;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nSCALARS particle_pressure_13 float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_p13;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nSCALARS particle_pressure_22 float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_p22;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nSCALARS particle_pressure_23 float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_p23;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nSCALARS particle_pressure_33 float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_p33;
              }
              if(!big_end) peg_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }

          /*
          fprintf(pfile,"\nTENSORS particle_pressure float\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[9*(i-il)]   = pGrid->Coup[k][j][i].grid_p11;
                data[9*(i-il)+1] = pGrid->Coup[k][j][i].grid_p12;
                data[9*(i-il)+2] = pGrid->Coup[k][j][i].grid_p13;
                data[9*(i-il)+3] = pGrid->Coup[k][j][i].grid_p12;
                data[9*(i-il)+4] = pGrid->Coup[k][j][i].grid_p22;
                data[9*(i-il)+5] = pGrid->Coup[k][j][i].grid_p23;
                data[9*(i-il)+6] = pGrid->Coup[k][j][i].grid_p13;
                data[9*(i-il)+7] = pGrid->Coup[k][j][i].grid_p23;
                data[9*(i-il)+8] = pGrid->Coup[k][j][i].grid_p33;
              }
              if(!big_end) peg_bswap(data,sizeof(float),9*(iu-il+1));
              fwrite(data,sizeof(float),(size_t)(9*ndata0),pfile);
            }
          }
          */
          
        }
        
        }

/* close file and free memory */

        fclose(pfile);
        free(data);
      }}
    }
  }

  return;
}
