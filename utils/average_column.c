#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "assert.h"

#define LENGTH 2000

int main(int argc, char* argv[]){
  char buf[LENGTH], *token, *fname;
  char ch;
  int nCols = 0, count = 0, lines = 0;
  int i;
  double t,t0, *ave, *stddev, val, dc;
  FILE* in;


  if(argc != 3){
    printf("Usage: %s t0\n",argv[0]);
    return 1; 
  }

  fname = argv[1];
  t0 = atof(argv[2]);

  in = fopen(fname,"r");
  assert(in != NULL);
  
  fgets(buf,LENGTH,in);
  //if(buf[0] == '\0') exit(0);

  while(buf[0] != '\0' && !feof(in)){
    if(buf[0] != '#'){
      lines++;
      if(lines == 1){
        token = strtok(buf,"   ");
        while( token != NULL){
          nCols++;
          token = strtok(NULL,"   ");
        }
      }
    }
    fgets(buf,LENGTH,in);
  }

  rewind(in);
  fgets(buf,LENGTH,in);

  ave    = calloc(nCols-1,sizeof(double));
  stddev = calloc(nCols-1,sizeof(double));

  while(buf[0] != '\0' && !feof(in)){
    if(buf[0] != '#'){
      i=0;
      token = strtok(buf,"   ");
      t=atof(token);
      if(t >= t0){
        count++;
        token = strtok(NULL,"   ");
        while( token != NULL){
          val = atof(token);
          ave[i] += val;
          stddev[i] += val*val;
          token = strtok(NULL,"   ");
          i++;
        }
      }
    }
    fgets(buf,LENGTH,in);
  }


  printf("Columns: %d  Lines: %d\n",nCols,count);

  dc = (double)count;
  for(i = 0; i < nCols-1; i++){
    ave[i];
    stddev[i] = sqrt((dc*stddev[i] - ave[i]*ave[i])/(dc*(dc - 1.0))); 
    printf("%e    ",ave[i]);
  }
  printf("\n");
  for(i = 0; i < nCols-1; i++){
    ave[i]/=((double)count);
    printf("%e    ",ave[i]);
  }
  printf("\n");

  for(i = 0; i < nCols-1; i++){
    printf("%e    ",stddev[i]);
  }
  printf("\n");

  fclose(in);
  free(ave);
  free(stddev);

  return 0;
}
