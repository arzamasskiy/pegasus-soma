/*==============================================================================
 * FILE: spec_to_gnuplot.c
 *
 * PURPOSE: This is go through a .spec file dumped out by pegasus and output
 *          a new file that's easier to read by gnuplot (my plotter of choice)
 *
 *       Code will output initial header, but will not output any subsequent
 *       comments.
 *
 * COMPILE USING: gcc -o spec_to_gnuplot spec_to_gnuplot.c
 *
 * USAGE: ./spec_to_gnuplot <infile>
 *                
 *
 * WRITTEN BY: Denis St-Onge, March 2018
 *============================================================================*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

# define LENGTH 10000000
# define MAX_LINES 5000000

static void usage(const char *arg);

/* ========================================================================== */

int main(int argc, char* argv[])
{
  /* file variables */
  FILE *fidin, *fidout;
  char *in_name;
  char *tok;
  char line[LENGTH];
  char fname[1000];
  /* data variables */
  float times[MAX_LINES];
  double perp_m, perp_M, par_m,par_M;
  double dperp,dpar;
  int Nperp, Npar;
  int n,i,j;

  /* Read Arguments */

  if(argc != 2) usage(argv[0]);

  in_name = argv[1];

  fidin = fopen(in_name,"r");
  if (fidin == NULL){
    fprintf(stderr,"Failed to open input file %s.\n",in_name);
    exit(1);
  }

/* Grab header and determine size of histogram */
 
  fgets(line,LENGTH,fidin);
  if(line[0] == '\0') exit(0);
  

  sscanf(line,"%*s %*s %*s %*s \
               %d %*s %*s      \
               %lf %*s         \
               %lf %*s         \
               %d %*s %*s      \
               %lf %*s         \
               %lf",&Nperp,&perp_m,&perp_M, &Npar,&par_m,&par_M);

  dpar = (par_M-par_m)/((double)Npar);
  dperp = (perp_M-perp_m)/((double)Nperp);

  //more header
  fgets(line,LENGTH,fidin);
  fgets(line,LENGTH,fidin);

  //data starts here
  fgets(line,LENGTH,fidin);
  n=0;
  while(line != NULL && ! feof(fidin)){
    i=0; j=0;
    sprintf(fname,"test_%d",n);
    fidout = fopen(fname,"w");

    tok=strtok(line," \t"); //time
    fprintf(fidout,"# time = %s\n",tok);
    
    tok=strtok(NULL," \t");    
    while(tok != NULL ){
      fprintf(fidout,"%e %e %s\n",i*dperp + perp_m,j*dpar + par_m,tok);
      tok=strtok(NULL," ");  
      i++;
      if(i % Nperp == 0){
        i=0;
        j++;
        fprintf(fidout,"\n");
      }
    }
    fclose(fidout);
    n++;
    fgets(line,LENGTH,fidin);
  }
  return 0;
}


/* ========================================================================== */


static void usage(const char *arg)
{
  fprintf(stderr,"Required Usage: %s infile\n", arg);
  exit(0);
}
