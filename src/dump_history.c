#include "copyright.h"
/*============================================================================*/
/*! \file dump_history.c
 *  \brief Functions to write dumps of scalar "history" variables in a
 *  formatted table.
 *
 * PURPOSE: Functions to write dumps of scalar "history" variables in a
 *   formatted table.  "History" variables are mostly VOLUME AVERAGES, e.g.
 *   - S = \Sum_{ijk} q[k][j][i]*dx1[i]*dx2[j]*dx3[k] /
 *             \Sum_{ijk} dx1[i]*dx2[j]*dx3[k],
 *
 *   although other quantities like time and dt are also output.  To compute
 *   TOTAL, VOLUME INTEGRATED quantities, just multiply by the Domain volume,
 *   which is output in the first line of the history file.  History data is
 *   written periodically, giving the time-evolution of that quantitity.
 *   Default (hardwired) values are:
 *   - scal[0] = time
 *   - scal[1] = dt
 *   - scal[2] = mass
 *   - scal[3] = total energy
 *   - scal[4] = d*v1
 *   - scal[5] = d*v2
 *   - scal[6] = d*v3
 *   - scal[7] = 0.5*d*v1**2
 *   - scal[8] = 0.5*d*v2**2
 *   - scal[9] = 0.5*d*v3**2
 *   - scal[10] = 0.5*b1**2
 *   - scal[11] = 0.5*b2**2
 *   - scal[12] = 0.5*b3**2
 *   - scal[13] = 0.5*e1**2
 *   - scal[14] = 0.5*e2**2
 *   - scal[15] = 0.5*e3**2
 *
 * More variables can be hardwired by increasing NSCAL=number of variables, and
 * adding calculation of desired quantities below.
 *
 * Alternatively, up to MAX_USR_H_COUNT new history variables can be added using
 * dump_history_enroll() in the problem generator.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - dump_history()        - Writes variables as formatted table
 * - dump_history_enroll() - Adds new user-defined history variable           */
/*============================================================================*/

#include <stdio.h>
#include "defs.h"
#include "pegasus.h"
#include "globals.h"
#include "prototypes.h"

/* Maximum Number of default history dump columns. */
#define NSCAL 16

/* Maximum number of history dump columns that the user routine can add. */
#define MAX_USR_H_COUNT 30

/* Array of strings / labels for the user added history column. */
static char *usr_label[MAX_USR_H_COUNT];

/* Array of history dump function pointers for user added history columns. */
static ConsFun_t phst_fun[MAX_USR_H_COUNT];

static int usr_hst_cnt = 0; /* User History Counter <= MAX_USR_H_COUNT */

/*----------------------------------------------------------------------------*/
/*! \fn void dump_history(MeshS *pM, OutputS *pOut)
 *  \brief Function to write dumps of scalar "history" variables in a
 *  formatted table. */

void dump_history(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  GPCouple *pq;
  DomainS *pD;
  int i,j,k,is,ie,js,je,ks,ke,nl,nd;
  double dVol, scal[NSCAL + MAX_USR_H_COUNT], d1;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL,*pdir=NULL,fmt[80];
  char levstr[8],domstr[8],dirstr[20];
  int n, total_hst_cnt, mhst, myID_Comm_Domain=1;
#ifdef MPI_PARALLEL
  double my_scal[NSCAL + MAX_USR_H_COUNT]; /* My Volume averages */
  int ierr;
#endif
  
  total_hst_cnt = 16 + usr_hst_cnt;

/* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%14.6e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

/* store time and dt in first two elements of output vector */

  scal[0] = pM->time;
  scal[1] = pM->dt;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){

    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        is = pG->is, ie = pG->ie;
        js = pG->js, je = pG->je;
        ks = pG->ks, ke = pG->ke;

        for (i=2; i<total_hst_cnt; i++) {
          scal[i] = 0.0;
        }
 
/* Compute history variables */

        for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
              
              pq = (&(pG->Coup[k][j][i]));
              
              dVol = 1.0; 
              if (pG->dx1 > 0.0) dVol *= pG->dx1;
              if (pG->dx2 > 0.0) dVol *= pG->dx2;
              if (pG->dx3 > 0.0) dVol *= pG->dx3;

              mhst = 2;
              scal[mhst] += dVol*pq->grid_d;
              d1 = 1.0/pq->grid_d;

              mhst++;
              scal[mhst] += dVol*0.5*(SQR(pG->U[k][j][i].B1c)+
                                      SQR(pG->U[k][j][i].B2c)+
                                      SQR(pG->U[k][j][i].B3c)+
                                      SQR(pG->U[k][j][i].E1c)+
                                      SQR(pG->U[k][j][i].E2c)+
                                      SQR(pG->U[k][j][i].E3c)+
                                      SQR(pq->grid_M1)*d1+
                                      SQR(pq->grid_M2)*d1+
                                      SQR(pq->grid_M3)*d1);

              mhst++;
              scal[mhst] += dVol*pq->grid_M1;
              mhst++;
              scal[mhst] += dVol*pq->grid_M2;
              mhst++;
              scal[mhst] += dVol*pq->grid_M3;
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pq->grid_M1)*d1;
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pq->grid_M2)*d1;
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pq->grid_M3)*d1;

              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B1c);
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B2c);
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B3c);
              
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].E1c);
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].E2c);
              mhst++;
              scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].E3c);

/* Calculate the user defined history variables */
              for(n=0; n<usr_hst_cnt; n++){
                mhst++;
                scal[mhst] += dVol*(*phst_fun[n])(pG, i, j, k);
              }
            }
          }
        }

/* Compute the sum over all Grids in Domain */

#ifdef MPI_PARALLEL 
        for(i=2; i<total_hst_cnt; i++){
          my_scal[i] = scal[i];
        }
        ierr = MPI_Reduce(&(my_scal[2]), &(scal[2]), (total_hst_cnt - 2),
          MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#endif

/* Only the parent (rank=0) process computes the average and writes output.
 * For single-processor jobs, myID_Comm_world is always zero. */

#ifdef MPI_PARALLEL
        ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
#endif
        if((myID_Comm_Domain==0) || (myID_Comm_world==0)){  /* I'm the parent */

/* Compute volume averages */

          dVol = pD->MaxX[0] - pD->MinX[0];
          if (pD->Nx[1] > 1) dVol *= (pD->MaxX[1] - pD->MinX[1]);
          if (pD->Nx[2] > 1) dVol *= (pD->MaxX[2] - pD->MinX[2]);
          for(i=2; i<total_hst_cnt; i++){
            scal[i] /= dVol;
          }

/* Create filename and open file.  History files are always written in lev#
 * directories of root process (rank=0 in MPI_COMM_WORLD) */
#ifdef MPI_PARALLEL
          if (nl>0) {
            plev = &levstr[0];
            sprintf(plev,"lev%d",nl);
            pdir = &dirstr[0];
            sprintf(pdir,"../id0/lev%d",nl);
          }
#else
          if (nl>0) {
            plev = &levstr[0];
            sprintf(plev,"lev%d",nl);
            pdir = &dirstr[0];
            sprintf(pdir,"lev%d",nl);
          }
#endif

          if (nd>0) {
            pdom = &domstr[0];
            sprintf(pdom,"dom%d",nd);
          }

          fname = peg_fname(pdir,pM->outfilename,plev,pdom,0,0,NULL,"hst");
          if(fname == NULL){
            peg_perr(-1,"[dump_history]: Unable to create history filename\n");
          }
          pfile = fopen(fname,"a");
          if(pfile == NULL){
            peg_perr(-1,"[dump_history]: Unable to open the history file\n");
          }

/* Write out column headers, but only for first dump */

          mhst = 0;
          if(pOut->num == 0){
            fprintf(pfile,
         "# Pegasus history dump for level=%i domain=%i volume=%e\n",nl,nd,dVol);
            mhst++;
            fprintf(pfile,"#   [%i]=time   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=dt      ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=mass    ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=total E ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3 Mom. ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-KE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1-ME   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-ME   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-ME   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x1-EE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x2-EE   ",mhst);
            mhst++;
            fprintf(pfile,"   [%i]=x3-EE   ",mhst);

            for(n=0; n<usr_hst_cnt; n++){
              mhst++;
              fprintf(pfile,"  [%i]=%s",mhst,usr_label[n]);
            }
            fprintf(pfile,"\n#\n");
          }

/* Write out data, and close file */

          for (i=0; i<total_hst_cnt; i++) {
            fprintf(pfile,fmt,scal[i]);
          }
          fprintf(pfile,"\n");
          fclose(pfile);
  
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_history_enroll(const ConsFun_t pfun, const char *label)
 *  \brief Adds new user-defined history variables                            */

void dump_history_enroll(const ConsFun_t pfun, const char *label){

  if(usr_hst_cnt >= MAX_USR_H_COUNT)
    peg_error("[dump_history_enroll]: MAX_USR_H_COUNT = %d exceeded\n",
              MAX_USR_H_COUNT);

/* Copy the label string */
  if((usr_label[usr_hst_cnt] = peg_strdup(label)) == NULL)
    peg_error("[dump_history_enroll]: Error on sim_strdup(\"%s\")\n",label);

/* Store the function pointer */
  phst_fun[usr_hst_cnt] = pfun;

  usr_hst_cnt++;

  return;

}

#undef NSCAL
#undef MAX_USR_H_COUNT
