#include "copyright.h"
/*============================================================================*/
/*! \file restart.c
 *  \brief Functions for writing and reading restart files.
 *
 * PURPOSE: Functions for writing and reading restart files.  Restart files
 *   begin with the input parameter file (peginput.XX) used to start the run in
 *   text format, followed by binary format data for selected quantities (time,
 *   dt, etc) for each of the Grid structures being updated by this processor.
 *   Lastly, any problem-specific data in binary format is included, which must
 *   be written by a user function in the problem generator.  Since all of the
 *   mesh data (NLevels, number of Domains, size of Grids, etc.) can be
 *   recalculated from the peginput file, only the dependent variables (time,
 *   dt, nstep, U, B, etc.) are dumped in restart files.
 *
 *   The peginput file included at the start of the restart file is parsed by
 *   main() in Step 2, using the standard par_* functions, with the parameters
 *   used by init_mesh() and init_grid() to re-initialize the grid.  Parameter
 *   values read from the peginput file contained in the resfile can be
 *   superceded by input from the command line, or another input file.
 *
 * MPI parallel jobs must be restarted on the same number of processors as they
 * were run originally.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - restart_grids() - reads nstep,time,dt,ConsS and B from restart file 
 * - dump_restart()  - writes a restart file
 *									      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "pegasus.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

/*----------------------------------------------------------------------------*/
/*! \fn void restart_grids(char *res_file, MeshS *pM)
 *  \brief Reads nstep, time, dt, and arrays of ConsS and interface B
 *   for each of the Grid structures in the restart file.  
 *
 *    By the time this
 *   function is called (in Step 6 of main()), the Mesh hierarchy has already
 *   been re-initialized by init_mesh() and init_grid() in Step 4 of main()
 *   using parameters in the peginput file at the start of this restart file,
 *   the command line, or from a new input file.
 */

void restart_grids(char *res_file, MeshS *pM)
{
  GridS *pG;
  FILE *fp;
  char line[MAXLEN];
  int i,j,k,is,ie,js,je,ks,ke,nl,nd;
  int ib=0,jb=0,kb=0;
  long p;


/* Open the restart file */

  if((fp = fopen(res_file,"r")) == NULL)
    peg_error("[restart_grids]: Error opening the restart file\nIf this is a MPI job, make sure each file from each processor is in the same directory.\n");

/* Skip over the parameter file at the start of the restart file */

  do{
    fgets(line,MAXLEN,fp);
  }while(strncmp(line,"<par_end>",9) != 0);

/* read nstep */

  fgets(line,MAXLEN,fp);
  if(strncmp(line,"N_STEP",6) != 0)
    peg_error("[restart_grids]: Expected N_STEP, found %s",line);
  fread(&(pM->nstep),sizeof(int),1,fp);

/* read time */

  fgets(line,MAXLEN,fp);   /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"TIME",4) != 0)
    peg_error("[restart_grids]: Expected TIME, found %s",line);
  fread(&(pM->time),sizeof(Real),1,fp);

/* read dt */

  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"TIME_STEP",9) != 0)
    peg_error("[restart_grids]: Expected TIME_STEP, found %s",line);
  fread(&(pM->dt),sizeof(Real),1,fp);
#ifdef STS
  fread(&(pM->STS_dt),sizeof(Real),1,fp);
#endif

/* Now loop over all Domains containing a Grid on this processor */

  for (nl=0; nl<=(pM->NLevels)-1; nl++){
  for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pG=pM->Domain[nl][nd].Grid;
      is = pG->is;
      ie = pG->ie;
      js = pG->js;
      je = pG->je;
      ks = pG->ks;
      ke = pG->ke;

/* propagate time and dt to all Grids */

      pG->time = pM->time;
      pG->dt   = pM->dt;

/* if there is more than one cell in each dimension, need to read one more
 * face-centered field component than the number of cells.  [ijk]b is
 * the number of extra cells to be read  */

      if (ie > is) ib=1;
      if (je > js) jb=1;
      if (ke > ks) kb=1;

/* Read the face-centered x1 B-field */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"1-FIELD",7) != 0)
        peg_error("[restart_grids]: Expected 1-FIELD, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie+ib; i++) {
            fread(&(pG->B1i[k][j][i]),sizeof(Real),1,fp);
          }
        }
      }

/* Read the face-centered x2 B-field */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"2-FIELD",7) != 0)
        peg_error("[restart_grids]: Expected 2-FIELD, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je+jb; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->B2i[k][j][i]),sizeof(Real),1,fp);
          }
        }
      }

/* Read the face-centered x3 B-field */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"3-FIELD",7) != 0)
        peg_error("[restart_grids]: Expected 3-FIELD, found %s",line);
      for (k=ks; k<=ke+kb; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->B3i[k][j][i]),sizeof(Real),1,fp);
          }
        }
      }

/* initialize the cell center magnetic fields as either the average of the face
 * centered field if there is more than one cell in that dimension, or just
 * the face centered field if not  */

      for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = pG->B1i[k][j][i];
        pG->U[k][j][i].B2c = pG->B2i[k][j][i];
        pG->U[k][j][i].B3c = pG->B3i[k][j][i];
        if(ib==1) pG->U[k][j][i].B1c=0.5*(pG->B1i[k][j][i] +pG->B1i[k][j][i+1]);
        if(jb==1) pG->U[k][j][i].B2c=0.5*(pG->B2i[k][j][i] +pG->B2i[k][j+1][i]);
        if(kb==1) pG->U[k][j][i].B3c=0.5*(pG->B3i[k][j][i] +pG->B3i[k+1][j][i]);
      }}}
      
/* Read particle properties and the complete particle list */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE LIST",13) != 0)
        peg_error("[restart_grids]: Expected PARTICLE LIST, found %s",line);
      fread(&(pG->nparticle),sizeof(long),1,fp);

      if (pG->nparticle > pG->arrsize-2)
        particle_realloc(pG, pG->nparticle+2);

      fread(&(npartypes),sizeof(int),1,fp);
      for (i=0; i<npartypes; i++) {          /* particle property list */
        fread(&(grproperty[i].m),sizeof(Real),1,fp);
        fread(&(grproperty[i].qomc),sizeof(Real),1,fp);
      }

      for (i=0; i<npartypes; i++)
        fread(&(grproperty[i].integrator),sizeof(short),1,fp);
      
      int nfpassFAIL=0;
      fread(&(nfpassFAIL),sizeof(int),1,fp); /* number of filter passes */

/* Read the x1-positions */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE X1",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE X1, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].x1),sizeof(Real),1,fp);
      }

/* Read the x2-positions */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE X2",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE X2, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].x2),sizeof(Real),1,fp);
      }

/* Read the x3-positions */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE X3",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE X3, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].x3),sizeof(Real),1,fp);
      }

/* Read the v1 velocity */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE V1",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE V1, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].v1),sizeof(Real),1,fp);
      }

/* Read the v2 velocity */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE V2",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE V2, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].v2),sizeof(Real),1,fp);
      }

/* Read the v3 velocity */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE V3",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE V3, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].v3),sizeof(Real),1,fp);
      }
      
#ifdef DELTA_F
      /* Read f(t=0) */
      
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE F0",11) != 0)
        peg_error("[restart_grids]: Expected PARTICLE F0, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].f_0),sizeof(Real),1,fp);
      }
      
#endif /* DELTA_F */
      
/* Read particle properties */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE PROPERTY",17) != 0)
        peg_error("[restart_grids]: Expected PARTICLE PROPERTY, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].property),sizeof(int),1,fp);
        pG->particle[p].pos = 1;	/* grid particle */
      }

/* Read particle my_id */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE MY_ID",14) != 0)
        peg_error("[restart_grids]: Expected PARTICLE MY_ID, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].my_id),sizeof(long),1,fp);
      }

#ifdef MPI_PARALLEL
/* Read particle init_id */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE INIT_ID",16) != 0)
        peg_error("[restart_grids]: Expected PARTICLE INIT_ID, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].init_id),sizeof(int),1,fp);
      }
#endif

/* count the number of particles with different types */

      for (i=0; i<npartypes; i++)
        grproperty[i].num = 0;
      for (p=0; p<pG->nparticle; p++)
        grproperty[pG->particle[p].property].num += 1;

    }
  }} /* End loop over all Domains --------------------------------------------*/

/* Call a user function to read his/her problem-specific data! */

  fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"USER_DATA",9) != 0)
    peg_error("[restart_grids]: Expected USER_DATA, found %s",line);
  problem_read_restart(pM, fp);

  fclose(fp);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_restart(MeshS *pM, OutputS *pout)
 *  \brief Writes a restart file, including problem-specific data from
 *   a user defined function  */

void dump_restart(MeshS *pM, OutputS *pout)
{
  GridS *pG;
  FILE *fp;
  char *fname;
  int i,j,k,is,ie,js,je,ks,ke,nl,nd;
  int ib=0,jb=0,kb=0;
  int nprop, *ibuf = NULL, ibufsize, lbufsize, sbufsize, nibuf, nlbuf, nsbuf;
  long np, p, *lbuf = NULL;
  short *sbuf = NULL;
  nibuf = 0;	nsbuf = 0;	nlbuf = 0;
  int bufsize, nbuf = 0;
  Real *buf = NULL;

#ifdef DRIVING
  struct RNG_State curState;
#endif

/* Allocate memory for buffer */
  bufsize = 262144 / sizeof(Real);  /* 256 KB worth of Reals */
  if ((buf = (Real*)calloc_1d_array(bufsize, sizeof(Real))) == NULL) {
    peg_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
  ibufsize = 262144 / sizeof(int);  /* 256 KB worth of ints */
  if ((ibuf = (int*)calloc_1d_array(ibufsize, sizeof(int))) == NULL) {
    peg_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
  lbufsize = 262144 / sizeof(long);  /* 256 KB worth of longs */
  if ((lbuf = (long*)calloc_1d_array(lbufsize, sizeof(long))) == NULL) {
    peg_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
  sbufsize = 262144 / sizeof(short);  /* 256 KB worth of longs */
  if ((sbuf = (short*)calloc_1d_array(sbufsize, sizeof(short))) == NULL) {
    peg_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }

/* Create filename and Open the output file */

  if((fname = peg_fname(NULL,pM->outfilename,NULL,NULL,num_digit,
      pout->num,NULL,"rst")) == NULL){
    peg_error("[dump_restart]: Error constructing filename\n");
  }

  if((fp = fopen(fname,"wb")) == NULL){
    peg_error("[dump_restart]: Unable to open restart file\n");
    return;
  }

/* Add the current time & nstep to the parameter file */

  par_setd("time","time","%e",pM->time,"Current Simulation Time");
  par_seti("time","nstep","%d",pM->nstep,"Current Simulation Time Step");

/* Write the current state of the parameter file */

  par_dump(2,fp);

/* Write out the current simulation step number */

  fprintf(fp,"N_STEP\n");
  if(fwrite(&(pM->nstep),sizeof(int),1,fp) != 1)
    peg_error("[dump_restart]: fwrite() error\n");

/* Write out the current simulation time */

  fprintf(fp,"\nTIME\n");
  if(fwrite(&(pM->time),sizeof(Real),1,fp) != 1)
    peg_error("[dump_restart]: fwrite() error\n");

/* Write out the current simulation time step */

  fprintf(fp,"\nTIME_STEP\n");
  if(fwrite(&(pM->dt),sizeof(Real),1,fp) != 1)
    peg_error("[dump_restart]: fwrite() error\n");
#ifdef STS
  if (fwrite(&(pM->STS_dt),sizeof(Real),1,fp) != 1)
    peg_error("[dump_restart]: fwrite() error\n");
#endif

/* Now loop over all Domains containing a Grid on this processor */

  for (nl=0; nl<=(pM->NLevels)-1; nl++){
  for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pG=pM->Domain[nl][nd].Grid;
      is = pG->is;
      ie = pG->ie;
      js = pG->js;
      je = pG->je;
      ks = pG->ks;
      ke = pG->ke;

/* see comments in restart_grid_block() for use of [ijk]b */

  if (ie > is) ib = 1;
  if (je > js) jb = 1;
  if (ke > ks) kb = 1;

/* Write the x1-field */

      fprintf(fp,"\n1-FIELD\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie+ib; i++) {
            buf[nbuf++] = pG->B1i[k][j][i];
            if ((nbuf+1) > bufsize) {
              fwrite(buf,sizeof(Real),nbuf,fp);
              nbuf = 0;
            }
          }
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }

/* Write the x2-field */

      fprintf(fp,"\n2-FIELD\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je+jb; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->B2i[k][j][i];
            if ((nbuf+1) > bufsize) {
              fwrite(buf,sizeof(Real),nbuf,fp);
              nbuf = 0;
            }
          }
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }

/* Write the x3-field */

      fprintf(fp,"\n3-FIELD\n");
      for (k=ks; k<=ke+kb; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->B3i[k][j][i];
            if ((nbuf+1) > bufsize) {
              fwrite(buf,sizeof(Real),nbuf,fp);
              nbuf = 0;
            }
          }
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
      

#ifdef DRIVING
/* Write the x1-force */
      
      fprintf(fp,"\n1-FORCE\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->force[k][j][i].x1;
            if ((nbuf+1) > bufsize) {
              fwrite(buf,sizeof(Real),nbuf,fp);
              nbuf = 0;
            }
          }
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
      
/* Write the x2-force */
      
      fprintf(fp,"\n2-FORCE\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->force[k][j][i].x2;
            if ((nbuf+1) > bufsize) {
              fwrite(buf,sizeof(Real),nbuf,fp);
              nbuf = 0;
            }
          }
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
      
/* Write the x3-force */
      
      fprintf(fp,"\n3-FORCE\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->force[k][j][i].x3;
            if ((nbuf+1) > bufsize) {
              fwrite(buf,sizeof(Real),nbuf,fp);
              nbuf = 0;
            }
          }
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }

/* Write the current random state */
      
	  curState=getRNGState();
      fprintf(fp,"\nRSTATE\n");
      fwrite(&(curState),sizeof(struct RNG_State),1,fp);

#endif /* DRIVING */


/* Write out the number of particles */

      fprintf(fp,"\nPARTICLE LIST\n");
      np = 0;
      for (p=0; p<pG->nparticle; p++)
        if (pG->particle[p].pos == 1) np += 1;
      fwrite(&(np),sizeof(long),1,fp);
    
/* Write out the particle properties */
    
      nprop = 2;

      fwrite(&(npartypes),sizeof(int),1,fp); /* number of particle types */
      for (i=0; i<npartypes; i++) {          /* particle property list */
        buf[nbuf++] = grproperty[i].m;
        buf[nbuf++] = grproperty[i].qomc;
        if ((nbuf+nprop) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
    
      for (i=0; i<npartypes; i++) {         /* particle integrator type */
        sbuf[nsbuf++] = grproperty[i].integrator;
        if ((nsbuf+1) > sbufsize) {
          fwrite(sbuf,sizeof(short),nsbuf,fp);
          nsbuf = 0;
        }
      }
      if (nsbuf > 0) {
        fwrite(sbuf,sizeof(short),nsbuf,fp);
        nsbuf = 0;
      }
      
      fwrite(&(nfpass),sizeof(int),1,fp); /* number of filter passes */
    
/* Write x1 */
    
      fprintf(fp,"\nPARTICLE X1\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        buf[nbuf++] = pG->particle[p].x1;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
    
/* Write x2 */
    
      fprintf(fp,"\nPARTICLE X2\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        buf[nbuf++] = pG->particle[p].x2;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
    
/* Write x3 */
    
      fprintf(fp,"\nPARTICLE X3\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        buf[nbuf++] = pG->particle[p].x3;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
    
/* Write v1 */
    
      fprintf(fp,"\nPARTICLE V1\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        buf[nbuf++] = pG->particle[p].v1;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
    
/* Write v2 */
    
      fprintf(fp,"\nPARTICLE V2\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        buf[nbuf++] = pG->particle[p].v2;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
    
/* Write v3 */
    
      fprintf(fp,"\nPARTICLE V3\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        buf[nbuf++] = pG->particle[p].v3;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
      
#ifdef DELTA_F
      
/* Write f_0 */
      
      fprintf(fp,"\nPARTICLE F0\n");
      for (p=0;p<pG->nparticle;p++)
        if (pG->particle[p].pos == 1){
          buf[nbuf++] = pG->particle[p].f_0;
          if ((nbuf+1) > bufsize) {
            fwrite(buf,sizeof(Real),nbuf,fp);
            nbuf = 0;
          }
        }
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }
      
#endif /*DELTA_F */
      
/* Write properties */
    
      fprintf(fp,"\nPARTICLE PROPERTY\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        ibuf[nibuf++] = pG->particle[p].property;
        if ((nibuf+1) > ibufsize) {
          fwrite(ibuf,sizeof(int),nibuf,fp);
          nibuf = 0;
        }
      }
      if (nibuf > 0) {
        fwrite(ibuf,sizeof(int),nibuf,fp);
        nibuf = 0;
      }
    
/* Write my_id */
    
      fprintf(fp,"\nPARTICLE MY_ID\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        lbuf[nlbuf++] = pG->particle[p].my_id;
        if ((nlbuf+1) > lbufsize) {
          fwrite(lbuf,sizeof(long),nlbuf,fp);
          nlbuf = 0;
        }
      }
      if (nlbuf > 0) {
        fwrite(lbuf,sizeof(long),nlbuf,fp);
        nlbuf = 0;
      }
    
#ifdef MPI_PARALLEL
/* Write init_id */
    
      fprintf(fp,"\nPARTICLE INIT_ID\n");
      for (p=0;p<pG->nparticle;p++)
      if (pG->particle[p].pos == 1){
        ibuf[nibuf++] = pG->particle[p].init_id;
        if ((nibuf+1) > ibufsize) {
          fwrite(ibuf,sizeof(int),nibuf,fp);
          nibuf = 0;
        }
      }
      if (nibuf > 0) {
        fwrite(ibuf,sizeof(int),nibuf,fp);
        nibuf = 0;
      }
#endif

    }
  }}  /*---------- End loop over all Domains ---------------------------------*/
    
/* call a user function to write his/her problem-specific data! */
    
  fprintf(fp,"\nUSER_DATA\n");
  problem_write_restart(pM, fp);

  fclose(fp);

  free_1d_array(buf);

  return;
}
