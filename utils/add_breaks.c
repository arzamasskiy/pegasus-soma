/*==============================================================================
 * FILE: add_breaks.c
 *
 * PURPOSE: For period boundary conditions, adds breaks into the data by 
 *          inserting a line of "inf"'s between steps where particle travels
 *          across the boundary.
 *
 *          Assumes a square box. Should be generalized to rectangular boxes.
 *
 * COMPILE USING: gcc -o add_breaks add_breaks.c -lm
 *
 * USAGE: ./add_breaks -s <box size> -o <basename-out> -i <basename-in> -s <post-name>
 *                
 *
 * WRITTEN BY: Denis St-Onge, June 2016
 *============================================================================*/

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

# define LENGTH 1000
# define N_LINES 50000
# define HEADER 3
# define COLUMNS 9
# define THRESHOLD 0.8

static void usage(const char *arg);

/* ========================================================================== */

int main(int argc, char* argv[])
{
  /* file variables */
  FILE *fidin,*fidout;
  char *in_name, *out_name;
  char line[LENGTH],extra[LENGTH];
  /* data variables */
  float t,box_size, thres;
  float x0,y0,z0,x1,y1,z1;
  float times[N_LINES], last_good_time;
  int i,expand = 0;
  long data_pos, nx=0,ny=0,nz=0,nlines=0;

  /* Read Arguments */

  for (i=1; i<argc; i++) {
/* If argv[i] is a 2 character string of the form "-?" then: */
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                                /* -d <basename-in>   */
        in_name = argv[++i];
        break;
      case 's':                                /* -s <post-name> */
        sscanf(argv[++i],"%f",&box_size);
        break;
      case 'o':                                /* -o <basename-out>   */
        out_name = argv[++i];
        break;
      case 'e':                                /* Expand box  */
        expand=1;
        break;
      case 'h':                                /* -h */
        usage(argv[0]);
        break;
      default:
        usage(argv[0]);
        break;
      }
    }
  }

  fidin = fopen(in_name,"r");
  if (fidin == NULL){
    fprintf(stderr,"Fail to open input file %s!\n",in_name);
    exit(1);
  }

  fidout = fopen(out_name,"w");
  if (fidout == NULL){
    fprintf(stderr,"Fail to open output file %s!\n",out_name);
    exit(1);
  }

/* Grad header and write into output */ 
  for(i = 0; i < HEADER; i++){
  	fgets(line,LENGTH,fidin);
	if(line == NULL){
  	    fprintf(stderr,"Particle tracking not of the right format.\n");
	    exit(1);
 	}
	fputs(line,fidout);
  }
  data_pos=ftell(fidin);

/* if the run has been restarted, clean the extraneous times*/
  fgets(line,LENGTH,fidin);
  while(line != NULL && !feof(fidin)){
    sscanf(line,"%f", &times[nlines++]);
    fgets(line,LENGTH,fidin);
  }

  fseek(fidin,data_pos,SEEK_SET);

  last_good_time=times[nlines-1];
  for(i = nlines-1; i >= 0; i--){
      if(times[i-1] < last_good_time)
          last_good_time=times[i-1];
      else
  	times[i-1] = -1; 
  }  
 
  i=0;
  x0=0.25*box_size;
  y0=0.25*box_size;
  z0=0.25*box_size;
  thres = box_size * THRESHOLD;

  fgets(line,LENGTH,fidin);
  while(line != NULL && !feof(fidin)){
	if(times[i] != -1){
    	    sscanf(line,"%f %f %f %f %2000[^ \n]s", &t, &x1,&y1,&z1,&extra);
	    if(!expand){
	        if(abs(x1-x0) > thres || abs(y1-y0) > thres || abs(z1-z0) > thres)
                fputs("inf inf inf inf inf inf inf inf inf \n",fidout);
	            fputs(line,fidout);
	    } else {
	        if(abs(x1-x0)>thres) (x1-x0) > 0 ? nx-- : nx++; 
	        if(abs(y1-y0)>thres) (y1-y0) > 0 ? ny-- : ny++; 
	        if(abs(z1-z0)>thres) (z1-z0) > 0 ? nz-- : nz++; 
	        fprintf(fidout,"%f %f %f %f %s\n",t,x1 + nx*box_size, y1 + ny*box_size, z1 +nz*box_size,extra);
	    }
	    x0=x1;
	    y0=y1;
	    z0=z1;
	}
	i++;
	fgets(line,LENGTH,fidin);
  }

  return 0;
}


/* ========================================================================== */


static void usage(const char *arg)
{
  fprintf(stderr,"\nRequired Usage: %s -i infile -o outfile -s box_size\n", arg);
  fprintf(stderr,"  -e         extend the domain instead of breaking lines?\n");

  exit(0);
}
