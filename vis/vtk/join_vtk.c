/*==============================================================================
 * FILE: join_vtk.c
 *
 * PURPOSE: Joins together multiple vtk files generated by an MPI job into one
 *   file for visualization and analysis.  The code is not optimal and can be
 *   a memory hog, as it loads all data into memory before writing out the 
 *   joined file.
 *
 * Modified July 2011: Reads NGrid_y times NGrid_x files at one time, to avoid
 * exceeding file count (1024) limits in big simulations.
 *
 * COMPILE USING: gcc -Wall -W -o join_vtk join_vtk.c -lm
 *
 * USAGE: ./join_vtk -o <outfile.vtk> infile1.vtk infile2.vtk ...
 *
 * WRITTEN BY: Tom Gardiner, November 2004
 *============================================================================*/

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

/* This code is not bulletproof.  In a variety of places I assume a
   simplified version of the VTK file format.  It should work for our
   VTK data dumps though.  --  T. A. Gardiner  --  Nov. 17, 2004 */

/* NOTE: I copied various routines such as strdup() and the 3d_array()
   routines into this file so as to make it self contained.  Including
   prototypes.h can lead to funny dependencies such as MPI,
   etc. depending on the athena configuration.
   --  T. A. Gardiner -- Nov. 17, 2004 */

/* Added fseek() to write_joined_vtk() to prevent leading data bytes
   that happen to correspond to "white space" from being consumed by
   fscanf().  -- Nicole Lemaster -- Feb. 16, 2006 */

/* Added tensor field capability. -- Matthew Kunz -- Dec. 17, 2012 */


static void join_error(const char *fmt, ...);
static void init_domain_1d(void);
static void sort_domain_1d(void);
static void write_joined_vtk(const char *out_name);
static char *my_strdup(const char *in);
static void free_3d_array(void ***array);
static void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);


/* This stores the domain information of each vtk file */
typedef struct Domain_s{
  char *fname;
  char *comment;
  FILE *fp;
  long pos;
  int Nx, Ny, Nz;    /* Grid dimensions */
  double ox, oy, oz; /* Origin of this particular domain */
  double dx, dy, dz; /* grid cell size */
}VTK_Domain;


static int file_count; /* Number of input vtk files */
static VTK_Domain *domain_1d;

/* NGrid_{x,y,z} -> Number of grids in x-, y- and z-direction */
static int NGrid_x, NGrid_y, NGrid_z;
static VTK_Domain ***domain_3d=NULL;


/* ========================================================================== */


int main(int argc, char* argv[]){

  int i, j, k, n;
  char *out_name=NULL;

  if(argc < 5)
    join_error("Usage: %s -o <out_name.vtk> file1.vtk file2.vtk ...\n",argv[0]);

  file_count = argc-3; /* Number of input vtk files */

  /* Parse the command line for the output filename */
  for(i=1; i<argc-1; i++){
    if(strcmp(argv[i],"-o") == 0){
      i++; /* increment to the filename */
      if((out_name = my_strdup(argv[i])) == NULL)
	join_error("out_name = my_strdup(\"%s\") failed\n",argv[i]);
      break;
    }
  }

  /* An output filename is required */
  if(out_name == NULL)
    join_error("Usage: %s -o <out_name.vtk> file1.vtk file2.vtk ...\n",argv[0]);

  printf("Output filename is \"%s\"\n",out_name);
  printf("Found %d files on the command line\n",file_count);

  /* ====================================================================== */

  domain_1d = (VTK_Domain*)calloc(file_count,sizeof(VTK_Domain));
  if(domain_1d == NULL)
    join_error("calloc returned a NULL pointer for domain_1d\n");

  /* Populate the 1d domain array with the filenames */
  for(n=0, i=1; i<argc; i++){
    if(strcmp(argv[i],"-o") == 0){
      i++; /* increment to the filename */
    }
    else{
      if((domain_1d[n].fname = my_strdup(argv[i])) == NULL){
        join_error("domain_1d[%d].fname = my_strdup(\"%s\") failed\n",
		   n,argv[i]);
      }
      /* printf("domain_1d[%d].fname = \"%s\"\n",n,domain_1d[n].fname); */
      n++; /* increment the domain counter */
    }
  }

  init_domain_1d(); /* Read in the header information */

  sort_domain_1d(); /* Sort the VTK_Domain elements in [k][j][i] order */

  printf("NGrid_z %d NGrid_y %d NGrid_x %d\n",NGrid_z,NGrid_y, NGrid_x);

  /* Allocate the domain_3d[][][] array */
    domain_3d = (VTK_Domain***)
      calloc_3d_array(NGrid_z, NGrid_y, NGrid_x, sizeof(VTK_Domain));
    if(domain_3d == NULL)
      join_error("calloc_3d_array() returned a NULL pointer\n");

  /* Copy the contents of the domain_1d[] array to the domain_3d[][][] array */
  n=0;
  for(k=0; k<NGrid_z; k++){
    for(j=0; j<NGrid_y; j++){
      for(i=0; i<NGrid_x; i++){
        domain_3d[k][j][i] = domain_1d[n++];
      }
    }
  }

  /* For debugging purposes, write out the origin for all of the grid
     domains in what should now be an ascending order for a [k][j][i] array. */
  /* This mimics the output in sort_domain_1d() and a diff should show
     that they are exactly the same. */
  /* for(k=0; k<NGrid_z; k++){
     for(j=0; j<NGrid_y; j++){
     for(i=0; i<NGrid_x; i++){
     printf("[%d]: ox=%e  oy=%e  oz=%e\n",i + NGrid_x*(j + k*NGrid_y),
     domain_3d[k][j][i].ox, domain_3d[k][j][i].oy, 
     domain_3d[k][j][i].oz);
     }
     }
     } */

  /* TO MAKE THIS CODE MORE BULLETPROOF ADD A CALL TO SOME FUNCTION TO
     CHECK THAT THIS DOMAIN ARRAY SATISFIES A SET OF CONSISTENCY
     CONDITIONS LIKE Nx IS CONSTANT ALONG Y- AND Z-DIRECTIONS, ETC. */

  /* ====================================================================== */

  write_joined_vtk(out_name);

  /* ====================================================================== */

  /* Now we need to do some cleanup */
  free_3d_array((void ***)domain_3d);
  domain_3d = NULL;

  /* Now close the input files and free the input file names */
  /*for(i=0; i<file_count; i++){
    fclose(domain_1d[i].fp); 
    free(domain_1d[i].fname);
    domain_1d[i].fname = NULL;
    free(domain_1d[i].comment);
    domain_1d[i].comment = NULL;
  }

  free(domain_1d);
  domain_1d = NULL;

  free(out_name);
  out_name = NULL; */

  return(0) ;
}


/* ========================================================================== */


/* Write an error message and terminate the simulation with an error status. */
static void join_error(const char *fmt, ...){
  va_list ap;

  va_start(ap, fmt);         /* ap starts after the fmt parameter */
  vfprintf(stderr, fmt, ap); /* print the error message to stderr */
  va_end(ap);                /* end stdargs (clean up the va_list ap) */

  fflush(stderr);            /* flush it NOW */
  exit(1);                   /* clean up and exit */
}


/* ========================================================================== */


/* Input character pointer is to the start of a NUL terminated string */
static void strip_trail_white(char *pc){
  char *cp = pc;

  /* iterate down to the NUL terminator */
  while(*cp != '\0') cp++;

  while(cp > pc){
    cp--;
    if(isspace(*cp)) *cp = '\0';
    else break;
  }

  return;
}


static void init_domain_1d(void){
  FILE *fp; /* A temporary copy to make the code cleaner */
  int i, ndat, cell_dat;
  char line[256];

  for(i=0; i<file_count; i++){
    if((fp = domain_1d[i].fp = fopen(domain_1d[i].fname,"r")) == NULL)
      join_error("Error opening file \"%s\"\n",domain_1d[i].fname);

    printf("\ndomain_1d[%d].fname = \"%s\"\n",i,domain_1d[i].fname);

    /* get header */
    fgets(line,256,fp);
    strip_trail_white(line);
    if(strcmp(line,"# vtk DataFile Version 3.0") != 0 /* mymhd  */ &&
       strcmp(line,"# vtk DataFile Version 2.0") != 0 /* athena */ )
      join_error("First Line is \"%s\"\n",line);

    /* get comment field */
    fgets(line,256,fp);
    strip_trail_white(line);
    printf("Comment Field: \"%s\"\n",line);
    /* store the comment field */
    if((domain_1d[i].comment = my_strdup(line)) == NULL){
      join_error("domain_1d[%d].comment = my_strdup(\"%s\") failed\n",
		 i,line);
    }

    /* get BINARY or ASCII */
    fgets(line,256,fp);
    strip_trail_white(line);
    if(strcmp(line,"BINARY") != 0)
      join_error("Unsupported file type: %s",line);

    /* get DATASET STRUCTURED_POINTS */
    fgets(line,256,fp);
    strip_trail_white(line);
    if(strcmp(line,"DATASET STRUCTURED_POINTS") != 0)
      join_error("Unsupported file type: %s",line);

    /* I'm assuming from this point on that the header is in good shape */

    /* Dimensions */
    fscanf(fp,"DIMENSIONS %d %d %d\n",
	   &(domain_1d[i].Nx), &(domain_1d[i].Ny), &(domain_1d[i].Nz));
    printf("DIMENSIONS %d %d %d\n",
	   domain_1d[i].Nx, domain_1d[i].Ny, domain_1d[i].Nz);

    /* We want to store the number of grid cells, not the number of grid
       cell corners */
    if(domain_1d[i].Nx > 1) domain_1d[i].Nx--;
    if(domain_1d[i].Ny > 1) domain_1d[i].Ny--;
    if(domain_1d[i].Nz > 1) domain_1d[i].Nz--;

    /* Origin */
    fscanf(fp,"ORIGIN %le %le %le\n",
	   &(domain_1d[i].ox), &(domain_1d[i].oy), &(domain_1d[i].oz));
    printf("ORIGIN %e %e %e\n",
	   domain_1d[i].ox, domain_1d[i].oy, domain_1d[i].oz);

    /* Spacing, dx, dy, dz */
    fscanf(fp,"SPACING %le %le %le\n",
	   &(domain_1d[i].dx), &(domain_1d[i].dy), &(domain_1d[i].dz));
    printf("SPACING %e %e %e\n",
	   domain_1d[i].dx, domain_1d[i].dy, domain_1d[i].dz);

    /* Cell Data = Nx*Ny*Nz */
    fscanf(fp,"CELL_DATA %d\n",&cell_dat);
    printf("CELL_DATA %d\n",cell_dat);
    ndat = (domain_1d[i].Nx)*(domain_1d[i].Ny)*(domain_1d[i].Nz);
    if(cell_dat != ndat)
      join_error("Nx*Ny*Nz = %d\n",ndat);

    domain_1d[i].pos = ftell(fp);
    fclose(fp);
  }

  return;
}


/* ========================================================================== */


static int compare_ox(const void *p1, const void *p2){

  if(((VTK_Domain*)p1)->ox < ((VTK_Domain*)p2)->ox) return -1;
  else if(((VTK_Domain*)p1)->ox > ((VTK_Domain*)p2)->ox) return 1;
  else return 0;
}


static int compare_oy(const void *p1, const void *p2){

  if(((VTK_Domain*)p1)->oy < ((VTK_Domain*)p2)->oy) return -1;
  else if(((VTK_Domain*)p1)->oy > ((VTK_Domain*)p2)->oy) return 1;
  else return 0;
}


static int compare_oz(const void *p1, const void *p2){

  if(((VTK_Domain*)p1)->oz < ((VTK_Domain*)p2)->oz) return -1;
  else if(((VTK_Domain*)p1)->oz > ((VTK_Domain*)p2)->oz) return 1;
  else return 0;
}


static void sort_domain_1d(void){
  int i, j, xcount, ycount, zcount, xy_cnt, yz_cnt;
  double oy, oz;
  div_t cnt_div;

  /* Sort the domains in ascending order of the z-origin in each domain */
  qsort(domain_1d, file_count, sizeof(VTK_Domain), compare_oz);

  /* Now count the number of Grid domains in the z-direction */
  oz = domain_1d[0].oz;
  zcount=1;
  printf("zcount = 1 @ oz = %e\n",oz);
  for(i=1; i<file_count; i++){
    if(domain_1d[i].oz > oz){
      oz = domain_1d[i].oz;
      zcount++;
      printf("zcount = %d @ oz = %e\n",zcount,oz);
    }
  }

  cnt_div = div(file_count, zcount);

  if(cnt_div.rem != 0)
    join_error("file_count%%zcount = %d\n",cnt_div.rem);

  xy_cnt = cnt_div.quot;

  /* Sort each group of domains with the same z-origin in order of
     ascending y-origin. */
  for(i=0; i<zcount; i++){
    qsort(&(domain_1d[i*xy_cnt]), xy_cnt, sizeof(VTK_Domain), compare_oy);

    /* Count the number of grid domains in the y-direction */
    if(i == 0){
      oy = domain_1d[0].oy;
      ycount=1;
      printf("ycount = 1 @ oy = %e\n",oy);
      for(j=1; j<xy_cnt; j++){
	if(domain_1d[j].oy > oy){
	  oy = domain_1d[j].oy;
	  ycount++;
	  printf("ycount = %d @ oy = %e\n",ycount,oy);
	}
      }
    }
    /* For i != 0 we will make sure this is uniformly consistent later */
  }

  cnt_div = div(xy_cnt, ycount);

  if(cnt_div.rem != 0)
    join_error("xy_cnt%%ycount = %d\n",cnt_div.rem);

  xcount = cnt_div.quot;
  printf("xcount = %d\n",xcount);
  yz_cnt = ycount*zcount;

  /* Sort each group of domains with the same y-origin and z-origin in order of
     ascending x-origin. */
  for(i=0; i<yz_cnt; i++){
    qsort(&(domain_1d[i*xcount]), xcount, sizeof(VTK_Domain), compare_ox);
  }

  /* For debugging purposes, write out the origin for all of the grid
     domains in what should now be an ascending order for a [k][j][i] array. */
  /* for(i=0; i<file_count; i++){
     printf("[%d]: ox=%e  oy=%e  oz=%e\n",i,
     domain_1d[i].ox, domain_1d[i].oy, domain_1d[i].oz);
     } */

  /* Initialize NGrid_{x,y,z} */
  NGrid_x = xcount;
  NGrid_y = ycount;
  NGrid_z = zcount;

  return;
}


/* ========================================================================== */


static void read_write_scalar(FILE *fp_out){

  FILE *fp;
  int i, j, k, ig, jg, kg, nread;
  float fdat;

  for(kg=0; kg<NGrid_z; kg++){
    for(jg=0; jg<NGrid_y; jg++){
      for(ig=0; ig<NGrid_x; ig++){
        if((domain_3d[kg][jg][ig].fp = fopen(domain_3d[kg][jg][ig].fname,"r")) == NULL)
           join_error("Error opening file \"%s\"\n",domain_3d[kg][jg][ig].fname);
        fseek(domain_3d[kg][jg][ig].fp,domain_3d[kg][jg][ig].pos, SEEK_SET);
      }
    }
    for(k=0; k<domain_3d[kg][0][0].Nz; k++){
      for(jg=0; jg<NGrid_y; jg++){
	for(j=0; j<domain_3d[0][jg][0].Ny; j++){
	  for(ig=0; ig<NGrid_x; ig++){

            fp = domain_3d[kg][jg][ig].fp;

	    for(i=0; i<domain_3d[0][0][ig].Nx; i++){

	      if((nread = fread(&fdat, sizeof(float), 1, fp)) != 1)
		join_error("read_write_scalar error\n");

	      fwrite(&fdat, sizeof(float), 1, fp_out);
	    }
	  }
	}
      }
    }
    for(jg=0; jg<NGrid_y; jg++){
      for(ig=0; ig<NGrid_x; ig++){
        domain_3d[kg][jg][ig].pos=ftell(domain_3d[kg][jg][ig].fp);
        fclose(domain_3d[kg][jg][ig].fp);
      }
    }
  }

  return;
}


/* ========================================================================== */


static void read_write_vector(FILE *fp_out){

  FILE *fp;
  int i, j, k, ig, jg, kg, nread;
  float fvec[3];


  for(kg=0; kg<NGrid_z; kg++){
    for(jg=0; jg<NGrid_y; jg++){
      for(ig=0; ig<NGrid_x; ig++){
        if((domain_3d[kg][jg][ig].fp = fopen(domain_3d[kg][jg][ig].fname,"r")) == NULL)
           join_error("Error opening file \"%s\"\n",domain_3d[kg][jg][ig].fname);
        fseek(domain_3d[kg][jg][ig].fp,domain_3d[kg][jg][ig].pos, SEEK_SET);
      }
    }
    for(k=0; k<domain_3d[kg][0][0].Nz; k++){
      for(jg=0; jg<NGrid_y; jg++){
	for(j=0; j<domain_3d[0][jg][0].Ny; j++){
	  for(ig=0; ig<NGrid_x; ig++){

	    fp = domain_3d[kg][jg][ig].fp;

	    for(i=0; i<domain_3d[0][0][ig].Nx; i++){

	      if((nread = fread(fvec, sizeof(float), 3, fp)) != 3)
		join_error("read_write_vector error\n");

	      fwrite(fvec, sizeof(float), 3, fp_out);
	    }
	  }
	}
      }
    }
    for(jg=0; jg<NGrid_y; jg++){
      for(ig=0; ig<NGrid_x; ig++){
        domain_3d[kg][jg][ig].pos=ftell(domain_3d[kg][jg][ig].fp);
        fclose(domain_3d[kg][jg][ig].fp);
      }
    }
  }

  return;
}

/* ========================================================================== */


static void read_write_tensor(FILE *fp_out){
  
  FILE *fp;
  int i, j, k, ig, jg, kg, nread;
  float ften[9];
  
  
  for(kg=0; kg<NGrid_z; kg++){
    for(jg=0; jg<NGrid_y; jg++){
      for(ig=0; ig<NGrid_x; ig++){
        if((domain_3d[kg][jg][ig].fp = fopen(domain_3d[kg][jg][ig].fname,"r")) == NULL)
          join_error("Error opening file \"%s\"\n",domain_3d[kg][jg][ig].fname);
        fseek(domain_3d[kg][jg][ig].fp,domain_3d[kg][jg][ig].pos, SEEK_SET);
      }
    }
    for(k=0; k<domain_3d[kg][0][0].Nz; k++){
      for(jg=0; jg<NGrid_y; jg++){
        for(j=0; j<domain_3d[0][jg][0].Ny; j++){
          for(ig=0; ig<NGrid_x; ig++){
            
            fp = domain_3d[kg][jg][ig].fp;
            
            for(i=0; i<domain_3d[0][0][ig].Nx; i++){
              
              if((nread = fread(ften, sizeof(float), 9, fp)) != 9)
                join_error("read_write_tensor error\n");
              
              fwrite(ften, sizeof(float), 9, fp_out);
            }
          }
        }
      }
    }
    for(jg=0; jg<NGrid_y; jg++){
      for(ig=0; ig<NGrid_x; ig++){
        domain_3d[kg][jg][ig].pos=ftell(domain_3d[kg][jg][ig].fp);
        fclose(domain_3d[kg][jg][ig].fp);
      }
    }
  }
  
  return;
}

/* ========================================================================== */


static void write_joined_vtk(const char *out_name){
  FILE *fp_out;
  int nxt, nyt, nzt; /* Total number of grid cells in each dir. */
  int nxp, nyp, nzp;
  int i, j, k;
  double ox, oy, oz, dx, dy, dz;
  char type[128], variable[128], format[128];
  char t_type[128], t_variable[128], t_format[128]; /* Temporary versions */
  int retval;

  /* Count the total number of grid cells in each direction */
  nxt = nyt = nzt = 0;
  for(i=0; i<NGrid_x; i++)
    nxt += domain_3d[0][0][i].Nx;

  for(j=0; j<NGrid_y; j++)
    nyt += domain_3d[0][j][0].Ny;

  for(k=0; k<NGrid_z; k++)
    nzt += domain_3d[k][0][0].Nz;

  /* Initialize ox, oy, oz */
  ox = domain_3d[0][0][0].ox;
  oy = domain_3d[0][0][0].oy;
  oz = domain_3d[0][0][0].oz;

  /* Initialize dx, dy, dz */
  dx = domain_3d[0][0][0].dx;
  dy = domain_3d[0][0][0].dy;
  dz = domain_3d[0][0][0].dz;

  /* Count the number of grid cell corners */
  if(nxt >= 1 && dx > 0.0) nxp = nxt+1;
  else nxp = nxt; /* dx = 0.0 */

  if(nyt >= 1 && dy > 0.0) nyp = nyt+1;
  else nyp = nyt; /* dy = 0.0 */

  if(nzt >= 1 && dz > 0.0) nzp = nzt+1;
  else nzp = nzt; /* dz = 0.0 */

  /* Open the output file */
  if((fp_out = fopen(out_name,"w")) == NULL)
    join_error("Error opening the output file \"%s\"\n",out_name);

  /* Write out some header information */
  fprintf(fp_out,"# vtk DataFile Version 3.0\n");
  /* Save the comment field from the [0][0][0] vtk domain file */
  /* fprintf(fp_out,"Joined VTK files\n"); */
  fprintf(fp_out,"%s\n",domain_3d[0][0][0].comment);
  fprintf(fp_out,"BINARY\n");
  fprintf(fp_out,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp_out,"DIMENSIONS %d %d %d\n", nxp, nyp, nzp);
  fprintf(fp_out,"ORIGIN %e %e %e\n", ox, oy, oz);
  fprintf(fp_out,"SPACING %e %e %e\n", dx, dy, dz);
  fprintf(fp_out,"CELL_DATA %d\n",nxt*nyt*nzt);

  while(1){
    for(k=0; k<NGrid_z; k++){
      for(j=0; j<NGrid_y; j++){
	for(i=0; i<NGrid_x; i++){

          if((domain_3d[k][j][i].fp = fopen(domain_3d[k][j][i].fname,"r")) == NULL)
            join_error("Error opening file \"%s\"\n",domain_3d[k][j][i].fname);
          fseek(domain_3d[k][j][i].fp,domain_3d[k][j][i].pos, SEEK_SET);

	  if(i == 0 && j == 0 && k == 0){
	    retval = fscanf(domain_3d[k][j][i].fp,"%s %s %s\n",
			    type, variable, format);

	    if(retval == EOF){ /* Assuming no errors, we are done... */
	      fclose(fp_out);
	      return;
	    }

	    if(strcmp(format, "float") != 0){
	      fclose(fp_out);
	      join_error("Expected \"float\" format, found \"%s\"\n",format);
	    }
	  }
	  else{
	    retval = fscanf(domain_3d[k][j][i].fp,"%s %s %s\n",
			    t_type, t_variable, t_format);

	    if(retval == EOF){ /* This shouldn't occur... */
	      fclose(fp_out);
	      fprintf(stderr,"[join_vtk_dump]: EOF returned for file \"%s\"\n",
		      domain_3d[k][j][i].fname);
	      return;
	    }

	    if(strcmp(type, t_type) != 0 ||
	       strcmp(variable, t_variable) != 0 ||
	       strcmp(format, t_format) != 0 ){
	      fclose(fp_out);
	      join_error("mismatch in input file positions\n");
	    }
	  }

	  if(strcmp(type, "SCALARS") == 0){
	    /* Read in the LOOKUP_TABLE (only default supported for now) */
	    fscanf(domain_3d[k][j][i].fp,"%s %s\n", t_type, t_format);
	    if(strcmp(t_type, "LOOKUP_TABLE") != 0 ||
	       strcmp(t_format, "default") != 0 ){
	      fclose(fp_out);
	      fprintf(stderr,"Expected \"LOOKUP_TABLE default\"\n");
	      join_error("Found \"%s %s\"\n",t_type,t_format);
	    }
	  }

	  /* Prevent leading data bytes that correspond to "white space"
	     from being consumed by the fscanf's above -- MNL 2/6/06 */
          assert(fseek(domain_3d[k][j][i].fp, -1, SEEK_CUR) == 0);
	  while (isspace(fgetc(domain_3d[k][j][i].fp)))
            assert(fseek(domain_3d[k][j][i].fp, -2, SEEK_CUR) == 0);
	  assert(fgetc(domain_3d[k][j][i].fp) == '\n');

          domain_3d[k][j][i].pos=ftell(domain_3d[k][j][i].fp);
          fclose(domain_3d[k][j][i].fp);
	}
      }
    }

    printf("Reading: \"%s %s %s\"\n",type,variable,format);

    /* Now, every file should agree that we either have SCALARS or
       VECTORS or TENSORS data */
    fprintf(fp_out,"\n%s %s %s\n",type,variable,format);
    if(strcmp(type, "SCALARS") == 0){
      fprintf(fp_out,"LOOKUP_TABLE default\n");
      read_write_scalar(fp_out);
    }
    else if(strcmp(type, "VECTORS") == 0){
      read_write_vector(fp_out);
    }
    else if (strcmp(type, "TENSORS") == 0){
      read_write_tensor(fp_out);
    }
    else
      join_error("Input type = \"%s\"\n",type);
  }

  return;
}


/* ========================================================================== */


static char *my_strdup(const char *in){
  char *out = (char *)malloc((1+strlen(in))*sizeof(char));
  if(out == NULL) {
    fprintf(stderr,"my_strdup: failed to allocate %d\n",(1+strlen(in)));
    return NULL; /* malloc failed */
  }
  return strcpy(out,in);
}


/* ========================================================================== */


/* Free the 3-D array: array[nt][nr][nc]
   Usage: free_2d_array((void ***)array);
*/
static void free_3d_array(void ***array){

  free(array[0][0]);
  free(array[0]);
  free(array);
}


/* Construct a 3-D array: array[nt][nr][nc] 
   Usage: array = (double ***)calloc_3d_array(nt,nr,nc,sizeof(double));
*/
static void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size){

  void ***array;
  size_t i,j;

  if((array = (void ***)calloc(nt,sizeof(void**))) == NULL){
    fprintf(stderr,"failed calloc_3d 1(%d)",nt);
    return NULL;
  }

  if((array[0] = (void **)calloc(nt*nr,sizeof(void*))) == NULL){
    fprintf(stderr,"failed calloc_3d 2(%d)",nt*nr);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<nt; i++){
    array[i] = (void **)((unsigned char *)array[0] + i*nr*sizeof(void*));
  }

  if((array[0][0] = (void *)calloc(nt*nr*nc,size)) == NULL){
    fprintf(stderr,"failed calloc_3d(%d,%d,%d,%d)",nt,nr,nr,size);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(j=1; j<nr; j++){
    array[0][j] = (void **)((unsigned char *)array[0][j-1] + nc*size);
  }

  for(i=1; i<nt; i++){
    array[i][0] = (void **)((unsigned char *)array[i-1][0] + nr*nc*size);
    for(j=1; j<nr; j++){
      array[i][j] = (void **)((unsigned char *)array[i][j-1] + nc*size);
    }
  }

  return array;
}


/* ========================================================================== */
