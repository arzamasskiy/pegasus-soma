#include "copyright.h"
#define MAIN_C
/*============================================================================*/
/*! \file main.c
 *  \brief Pegasus main program file.
 *
 * //////////////////////////// PEGASUS Main Program \\\\\\\\\\\\\\\\\\\\\\\\\\\
 *
 *  Pegasus - developed by MW Kunz, JM Stone, X Bai.
 *  SOMA2018 version modified by L Arzamasskiy 
 *
 *  Infrastructure based upon Athena v4
 *
 *  History:
 * - v1.0 [Nov 2012] - 3D
 *
 * See the GNU General Public License for usage restrictions. 
 *									        
 * PRIVATE FUNCTION PROTOTYPES:
 * - change_rundir() - creates and outputs data to new directory
 * - usage()         - outputs help message and terminates execution	      */
/*============================================================================*/
static char *pegasus_version = "version 1.1 - 29-5-2014";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include "defs.h"
#include "pegasus.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   change_rundir - creates and outputs data to new directory
 *   usage         - outputs help message and terminates execution
 *============================================================================*/

static void change_rundir(const char *name);
static void usage(const char *prog);

/* Maximum number of mkdir() and chdir() file operations that will be executed
 * at once in the change_rundir() function when running in parallel, passed to
 * baton_start() and baton_end().
 */ 
#define MAX_FILE_OP 256

/*----------------------------------------------------------------------------*/
/*! \fn int main(int argc, char *argv[]) 
 *  \brief Pegasus main program  
 *
 * Steps in main:
 * - 1 - check for command line options and respond
 * - 2 - read input file and parse command line for changes
 * - 3 - set up diagnostic log files
 * - 4 - initialize Mesh, Domains, and Grids
 * - 5 - set initial conditions
 * - 6 - set boundary condition function pointers, and use to set BCs
 * - 7 - set function pointers for desired algorithms and physics
 * - 8 - write initial conditions to output file(s)
 * - 9 - main loop
 * - 10 - finish by writing final output(s), diagnostics, and free memory     */

int main(int argc, char *argv[])
{
  MeshS Mesh;             /* the entire mesh hierarchy, see pegasus.h */
  VDFun_t Integrate;      /* function pointer to integrator, set at runtime */
  int nl,nd;
  char *definput = "peginput";  /* default input filename */
  char *peginput = definput;
  int ires=0;             /* restart flag, set to 1 if -r argument on cmdline */
  char *res_file = NULL;  /* restart filename */
  char *rundir = NULL;    /* directory to which output is directed */
  int nstep_start=0;      /* number of cycles already completed on restart */
  char *name = NULL;
  char level_dir[80];     /* names of directories for levels > 0 */
  FILE *fp;               /* file pointer for data outputs */
  int nflag=0;            /* set to 1 if -n argument is given on command line */
  int i,nlim;             /* cycle index and limit */
  Real tlim;              /* time limit (in code units) */
  int dtncheck;		  /*!< number of cycles between checking if dt < dt_CFL */
  
  int out_level, err_level, lazy; /* diagnostic output & error log levels */
  int iflush, nflush;             /* flush buffers every iflush cycles */

  int iquit=0;  /* quit signal sent to peg_sig_act, our system signal handler */

/* local variables used for timing and performance measures */

  time_t start, stop;
  int have_time = time(&start);  /* Is current calendar time (UTC) available? */
  int zones;
  double cpu_time, step_time, zcs;
  long clk_tck = sysconf(_SC_CLK_TCK);
  struct tms tbuf;
  clock_t time0,time1, have_times;
  struct timeval tvs, tve;
  struct timeval tv_curr, tv_prev;
  struct timeval tv_init, tvd; 
  Real dt_done;

#ifdef MPI_PARALLEL
  char *pc, *suffix, new_name[MAXLEN];
  int len, h, m, s, err, use_wtlim=0;
  double wtend;
  if(MPI_SUCCESS != MPI_Init(&argc, &argv))
    peg_error("[main]: Error on calling MPI_Init\n");
#endif /* MPI_PARALLEL */

/*----------------------------------------------------------------------------*/
/* Steps in main:
 *  1 - check for command line options and respond
 *  2 - read input file and parse command line for changes
 *  3 - set up diagnostic log files
 *  4 - initialize Mesh, Domains, and Grids
 *  5 - set initial conditions
 *  6 - set boundary condition function pointers, and use to set BCs
 *  7 - set function pointers for desired algorithms and physics
 *  8 - write initial conditions to output file(s)
 *  9 - main loop
 *  10 - finish by writing final output(s), diagnostics, and free memory
 */

/*--- Step 1. ----------------------------------------------------------------*/
/* Check for command line options and respond.  See comments in usage()
 * for description of options.  */

  for (i=1; i<argc; i++) {
/* If argv[i] is a 2 character string of the form "-?" then: */
    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                      /* -i <file>   */
	peginput = argv[++i];
	break;
      case 'r':                      /* -r <file>   */
	ires = 1;
	res_file = argv[++i];
/* If input file is not set on command line, use the restart file */
	if(peginput == definput) peginput = res_file;
	break;
      case 'd':                      /* -d <directory>   */
	rundir = argv[++i];
	break;
      case 'n':                      /* -n */
	nflag = 1;
	break;
      case 'h':                      /* -h */
	usage(argv[0]);
	break;
      case 'c':                      /* -c */
	show_config();
	exit(0);
	break;
#ifdef MPI_PARALLEL
      case 't':                      /* -t hh:mm:ss */
	use_wtlim = 1; /* Logical to use a wall time limit */
	sscanf(argv[++i],"%d:%d:%d",&h,&m,&s);
	wtend = MPI_Wtime() + s + 60*(m + 60*h);
	printf("Wall time limit: %d hrs, %d min, %d sec\n",h,m,s);
	break;
#else
      default:
	usage(argv[0]);
	break;
#endif /* MPI_PARALLEL */
      }
    }
  }

/*--- Step 2. ----------------------------------------------------------------*/
/* Read input file and parse command line.  For MPI_PARALLEL jobs, parent reads
 * input file and distributes information to children  */

#ifdef MPI_PARALLEL
/* Get proc id (rank) in MPI_COMM_WORLD, store as global variable */

  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &myID_Comm_world))
    peg_error("[main]: Error on calling MPI_Comm_rank\n");

/* Only rank=0 processor reads input parameter file, parses command line,
 * broadcasts the contents of the (updated) parameter file to the children. */

  if(myID_Comm_world == 0){
    par_open(peginput);   /* for restarts, default is peginput=resfile */ 
    par_cmdline(argc,argv);
  }
  par_dist_mpi(myID_Comm_world,MPI_COMM_WORLD);

/* Modify the problem_id name in the <job> block to include information about
 * processor ids, so that all output filenames constructed from this name will
 * include myID_Comm_world as an identifier.  Only child processes modify
 * name, rank=0 process does not have myID_Comm_world in the filename */

  if(myID_Comm_world != 0){
    name = par_gets("job","problem_id");
    sprintf(new_name,"%s-id%d",name,myID_Comm_world);
    free(name);
    par_sets("job","problem_id",new_name,NULL);
  }

  show_config_par(); /* Add the configure block to the parameter database */

/* Share the restart flag with the children */

  if(MPI_SUCCESS != MPI_Bcast(&ires, 1, MPI_INT, 0, MPI_COMM_WORLD))
    peg_error("[main]: Error on calling MPI_Bcast\n");

/* rank=0 needs to send the restart file name to the children.  This requires 
 * sending the length of the restart filename string, the string, and then
 * having each child add my_id to the name so it opens the appropriate file */

/* Parent finds length of restart filename */

  if(ires){ 
    if(myID_Comm_world == 0)
      len = 1 + (int)strlen(res_file);

/* Share this length with the children */

    if(MPI_SUCCESS != MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD))
      peg_error("[main]: Error on calling MPI_Bcast\n");

    if(len + 10 > MAXLEN)
      peg_error("[main]: Restart filename length = %d is too large\n",len);

/* Share the restart filename with the children */

    if(myID_Comm_world == 0) strcpy(new_name, res_file);
    if(MPI_SUCCESS != MPI_Bcast(new_name, len, MPI_CHAR, 0, MPI_COMM_WORLD))
      peg_error("[main]: Error on calling MPI_Bcast\n");

/* Assume the restart file name is of the form
 *  [/some/dir/]basename.0000.rst and search for the periods in the name. */

    pc = &(new_name[len - 5]);
    if(*pc != '.'){
      peg_error("[main]: Bad Restart filename: %s\n",new_name);
    }

    do{ /* Position the char pointer at the first period */
      pc--;
      if(pc == new_name)
	peg_error("[main]: Bad Restart filename: %s\n",new_name);
    }while(*pc != '.');

/* Only children add myID_Comm_world to the filename */

    if(myID_Comm_world == 0) {
      strcpy(new_name, res_file);
    } else {       
      suffix = peg_strdup(pc);
      sprintf(pc,"-id%d%s",myID_Comm_world,suffix);
      free(suffix);
      res_file = new_name;
    }
  }

/* Quit MPI_PARALLEL job if code was run with -n option. */

  if(nflag){          
    par_dump(0,stdout);   
    par_close();
    MPI_Finalize();
    return 0;
  }

#else
/* For serial (single processor) job, there is only one process to open and
 * read input file  */

  myID_Comm_world = 0;
  par_open(peginput);   /* opens AND reads */
  par_cmdline(argc,argv);
  show_config_par();   /* Add the configure block to the parameter database */

/* Quit non-MPI_PARALLEL job if code was run with -n option. */

  if(nflag){
    par_dump(0,stdout);
    par_close();
    return 0;
  }
#endif /* MPI_PARALLEL */

/*--- Step 3. ----------------------------------------------------------------*/
/* set up the simulation log files */

/* Open <problem_id>.out and <problem_id>.err files if file_open=1 in the 
 * <log> block of the input file.  Otherwise, diagnositic output will go to
 * stdout and stderr streams. */

  if(par_geti_def("log","file_open",0)){
    iflush = par_geti_def("log","iflush",0);
    name = par_gets("job","problem_id");
    lazy = par_geti_def("log","lazy",1);
    /* On restart we use mode "a", otherwise we use mode "w". */
    peg_log_open(name, lazy, (ires ? "a" : "w"));
    free(name);  name = NULL;
  }
  else{
    iflush = par_geti_def("log","iflush",1);
  }
  iflush = iflush > 0 ? iflush : 0; /* make iflush non-negative */

/* Set the peg_log output and error logging levels */
  out_level = par_geti_def("log","out_level",0);
  err_level = par_geti_def("log","err_level",0);
#ifdef MPI_PARALLEL
    if(myID_Comm_world > 0){   /* Children may use different log levels */
    out_level = par_geti_def("log","child_out_level",-1);
    err_level = par_geti_def("log","child_err_level",-1);
  }
#endif /* MPI_PARALLEL */
  peg_log_set_level(out_level, err_level);

  if(have_time > 0) /* current calendar time (UTC) is available */
    peg_pout(0,"Simulation started on %s\n",ctime(&start));

/*--- Step 4. ----------------------------------------------------------------*/
/* Initialize nested mesh hierarchy. */

  init_mesh(&Mesh);
  init_grid(&Mesh);
  init_particle(&Mesh);

/*--- Step 5. ----------------------------------------------------------------*/
/* Set initial conditions, either by reading from restart or calling problem
 * generator.  But first start by setting variables in <time> block (these
 * control execution time), and reading EOS parameters from <problem> block.  */

  CourNo = par_getd("time","cour_no");
#ifdef VARIABLE_DT
  Safety = par_getd_def("time","safety",1.0);
#endif
  nlim   = par_geti_def("time","nlim",-1);
  tlim   = par_getd("time","tlim");
  
  dtncheck = par_geti_def("time","dtncheck",1);

  beta     = par_getd("problem","beta");
  beta_prp = par_getd_def("problem","beta_prp",beta);
  beta_prl = par_getd_def("problem","beta_prl",beta);
  ZTeTi    = par_getd_def("problem","ZTeTi",1.0);
#ifdef ADIABATIC
  Gamma   = par_getd("problem","gamma");
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;
#endif

  if(ires) {
    restart_grids(res_file, &Mesh);  /*  Restart */
    nstep_start = Mesh.nstep;
  } else {                           /* New problem */
    Mesh.dt = par_getd("time","timestep");
    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL) {
          problem(&(Mesh.Domain[nl][nd]));
          Mesh.Domain[nl][nd].Grid->dt = Mesh.dt;
        }
      }
    }
  }

/* Initialize the first nstep value to flush the output and error logs. */
  nflush = nstep_start + iflush;

/*--- Step 6. ----------------------------------------------------------------*/
/* set boundary value function pointers using BC flags in <grid> blocks, then
 * set boundary conditions for initial conditions. This includes particle
 * crossings, particle injection, deposition onto the grid, computation of 
 * the driving force (if desired), and MHD boundary conditions. */

  bvals_mhd_init(&Mesh);

  bvals_particle_init(&Mesh);
  exchange_gpcouple_init(&Mesh);

  for (nl=0; nl<(Mesh.NLevels); nl++){ 
    for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
      if (Mesh.Domain[nl][nd].Grid != NULL){
        bvals_particle(&(Mesh.Domain[nl][nd]));

        bvals_mhd(&(Mesh.Domain[nl][nd]));
      }
    }
  }


/*--- Step 7. ----------------------------------------------------------------*/
/* Set function pointers for integrator;
 * Allocate temporary arrays */

  init_output(&Mesh); 
  Integrate = integrate_init(&Mesh);

/*--- Step 8. ----------------------------------------------------------------*/
/* Setup complete, output initial conditions */

  if(out_level >= 0){
    fp = pegout_fp();
    par_dump(0,fp);      /* Dump a copy of the parsed information to pegout */
  }
  change_rundir(rundir); /* Change to run directory */
  peg_sig_init();        /* Install a signal handler */
  for (nl=1; nl<(Mesh.NLevels); nl++){
    sprintf(level_dir,"lev%d",nl);
    mkdir(level_dir, 0775); /* Create directories for levels > 0 */
  }

  gettimeofday(&tvs,NULL);
  tv_prev=tvs;
  if((have_times = times(&tbuf)) > 0)
    time0 = tbuf.tms_utime + tbuf.tms_stime;
  else
    time0 = clock();

/* Force output of everything (by passing last argument of data_output = 1) */

  if (ires==0) data_output(&Mesh, 0);

  peg_pout(0,"\nSetup complete, entering main loop...\n\n");
  peg_pout(0,"cycle=%i time=%e next dt=%e\n\n",Mesh.nstep, Mesh.time, Mesh.dt);

/*--- Step 9. ----------------------------------------------------------------*/
/* START OF MAIN INTEGRATION LOOP ==============================================
 * Steps are: (a) Check for data ouput
 *            (b) Integrate all Grids over Mesh hierarchy
 *            (c) Userwork
 *            (d) Update time, set new timestep
 *            (e) Set boundary values
 *            (f) check for stopping criteria
 */

  while (Mesh.time < tlim && (nlim < 0 || Mesh.nstep < nlim)) {

/* Check if timestep is small enough */
if (Mesh.nstep % dtncheck == 0){
#ifdef VARIABLE_DT
      new_dt(&Mesh,dtncheck);
#else
      peg_dt(&Mesh);
#endif
}

/*--- Step 9a. ---------------------------------------------------------------*/
/* Only write output's with t_out>t (last argument of data_output = 0) */
    data_output(&Mesh, 0);
/*--- Step 9c. ---------------------------------------------------------------*/
/* User work (defined in problem()) */
    Userwork_in_loop(&Mesh);

/*--- Step 9b. ---------------------------------------------------------------*/
/* Loop over all Domains and call Integrator */
    gettimeofday(&tvd,NULL);
    tv_init=tvd;
    for (nl=0; nl<(Mesh.NLevels); nl++){
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){
        deposit(&(Mesh.Domain[nl][nd]),1);
      }
    }
    gettimeofday(&tvd,NULL);
    step_time = (double)(tvd.tv_sec - tv_init.tv_sec) +
      1.0e-6*(double)(tvd.tv_usec - tv_init.tv_usec);
    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL){
          (*Integrate)(&(Mesh.Domain[nl][nd]));

        }
      }
    }

/*--- Step 9d. ---------------------------------------------------------------*/
/* Apply boundary conditions to exiting particles, inject particles if desired,
 * deposit particles onto grid and set boundary conditions on moments, compute
 * new driving force if desired, then set boundary conditions for magnetic 
 * field (and driving force, if desired)
 *
 * Note: Boundary values must be set after time is updated for t-dependent BCs.
 */
    gettimeofday(&tvd,NULL);
    tv_init=tvd;
    for (nl=0; nl<(Mesh.NLevels); nl++){ 
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
        if (Mesh.Domain[nl][nd].Grid != NULL){
          bvals_particle(&(Mesh.Domain[nl][nd]));

          bvals_mhd(&(Mesh.Domain[nl][nd]));
        }
      }
    }
    gettimeofday(&tvd,NULL);
    step_time = (double)(tvd.tv_sec - tv_init.tv_sec) +
      1.0e-6*(double)(tvd.tv_usec - tv_init.tv_usec);


/*--- Step 9e. ---------------------------------------------------------------*/
/* Update Mesh time, and time in all Grid's. */

    Mesh.nstep++;
    Mesh.time += Mesh.dt;
    for (nl=0; nl<(Mesh.NLevels); nl++){
      for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){
        if (Mesh.Domain[nl][nd].Grid != NULL){
          Mesh.Domain[nl][nd].Grid->time = Mesh.time;
          Mesh.Domain[nl][nd].Grid->nstep = Mesh.nstep;
        }
      }
    }

/*--- Step 9f. ---------------------------------------------------------------*/
/* Force quit if wall time limit reached.  Check signals from system */

#ifdef MPI_PARALLEL
    if(use_wtlim && (MPI_Wtime() > wtend))
      iquit = 103; /* an arbitrary, unused signal number */
#endif /* MPI_PARALLEL */
    if(peg_sig_act(&iquit) != 0) break;

/* Print diagnostic message, flush message buffers, and continue... */

    gettimeofday(&tv_curr,NULL);

    step_time = (double)(tv_curr.tv_sec - tv_prev.tv_sec) +
    1.0e-6*(double)(tv_curr.tv_usec - tv_prev.tv_usec);

    peg_pout(0,"cycle=%i time=%e next dt=%e step time=%e s\n",
	     Mesh.nstep,Mesh.time,Mesh.dt, step_time);

    tv_prev = tv_curr;

    if(nflush == Mesh.nstep){
      peg_flush_out();
      peg_flush_err();
      nflush += iflush;
    }
  } /* END OF MAIN INTEGRATION LOOP ==========================================*/

/*--- Step 10. ---------------------------------------------------------------*/
/* Finish up by computing zc/sec, dumping data, and deallocate memory */

/* Print diagnostic message as to why run terminated */

  if (Mesh.nstep == nlim)
    peg_pout(0,"\nterminating on cycle limit\n");
#ifdef MPI_PARALLEL
  else if(use_wtlim && iquit == 103)
    peg_pout(0,"\nterminating on wall-time limit\n");
#endif /* MPI_PARALLEL */
  else
    peg_pout(0,"\nterminating on time limit\n");

/* Get time used */

  gettimeofday(&tve,NULL);
  if(have_times > 0) {
    times(&tbuf);
    time1 = tbuf.tms_utime + tbuf.tms_stime;
    cpu_time = (time1 > time0 ? (double)(time1 - time0) : 1.0)/
      (double)clk_tck;
  } else {
    time1 = clock();
    cpu_time = (time1 > time0 ? (double)(time1 - time0) : 1.0)/
      (double)CLOCKS_PER_SEC;
  }

/* Calculate and print the zone-cycles/cpu-second on this processor */

  zones = 0;
  for (nl=0; nl<(Mesh.NLevels); nl++){ 
  for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
    if (Mesh.Domain[nl][nd].Grid != NULL) {
      zones += (Mesh.Domain[nl][nd].Grid->Nx[0])*
               (Mesh.Domain[nl][nd].Grid->Nx[1])*
               (Mesh.Domain[nl][nd].Grid->Nx[2]);
    }
  }}
  zcs = (double)zones*(double)((Mesh.nstep) - nstep_start)/cpu_time;

  peg_pout(0,"  tlim= %e   nlim= %i\n",tlim,nlim);
  peg_pout(0,"  time= %e  cycle= %i\n",Mesh.time,Mesh.nstep);
  peg_pout(0,"\nzone-cycles/cpu-second = %e\n",zcs);

/* Calculate and print the zone-cycles/wall-second on this processor */

  cpu_time = (double)(tve.tv_sec - tvs.tv_sec) +
    1.0e-6*(double)(tve.tv_usec - tvs.tv_usec);
  zcs = (double)zones*(double)((Mesh.nstep) - nstep_start)/cpu_time;
  peg_pout(0,"\nelapsed wall time = %e sec.\n",cpu_time);
  peg_pout(0,"\nzone-cycles/wall-second = %e\n",zcs);

/* Calculate and print total zone-cycles/wall-second on all processors */
#ifdef MPI_PARALLEL
  zones = 0;
  for (nl=0; nl<(Mesh.NLevels); nl++){ 
  for (nd=0; nd<(Mesh.DomainsPerLevel[nl]); nd++){  
    zones += (Mesh.Domain[nl][nd].Nx[0])*
               (Mesh.Domain[nl][nd].Nx[1])*
               (Mesh.Domain[nl][nd].Nx[2]);
  }}
  zcs = (double)zones*(double)(Mesh.nstep - nstep_start)/cpu_time;
  peg_pout(0,"\ntotal zone-cycles/wall-second = %e\n",zcs);
#endif /* MPI_PARALLEL */

/* complete any final User work */

  Userwork_after_loop(&Mesh);

/* Final output everything (last argument of data_output = 1) */

  data_output(&Mesh, 1);
  
/* Free all memory */

  integrate_destruct();
  data_output_destruct();
  particle_destruct(&Mesh);
  bvals_particle_destruct(&Mesh);
  exchange_gpcouple_destruct(&Mesh);

  par_close();

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif /* MPI_PARALLEL */

  if(time(&stop)>0) /* current calendar time (UTC) is available */
    peg_pout(0,"\nSimulation terminated on %s",ctime(&stop));

  peg_log_close(); /* close the simulation log files */

  return EXIT_SUCCESS;
}

/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*! \fn void change_rundir(const char *name) 
 *  \brief Change run directory;  create it if it does not exist yet
 */

void change_rundir(const char *name)
{
#ifdef MPI_PARALLEL

  int err=0;
  int rerr, gerr, my_id, status;
  char mydir[80];

  status = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  if(status != MPI_SUCCESS)
    peg_error("[change_rundir]: MPI_Comm_rank error = %d\n",status);

  if(name != NULL && *name != '\0'){

    if(my_id == 0)
      mkdir(name, 0775); /* May return an error, e.g. the directory exists */

    MPI_Barrier(MPI_COMM_WORLD); /* Wait for rank 0 to mkdir() */

    baton_start(MAX_FILE_OP, ch_rundir0_tag);

    if(chdir(name)){
      peg_perr(-1,"[change_rundir]: Cannot change directory to \"%s\"\n",name);
      err = 1;
    }

    baton_stop(MAX_FILE_OP, ch_rundir0_tag);

    /* Did anyone fail to make and change to the run directory? */
    rerr = MPI_Allreduce(&err, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(rerr) peg_perr(-1,"[change_rundir]: MPI_Allreduce error = %d\n",rerr);

    if(rerr || gerr){
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(EXIT_FAILURE);
    }
  }

/* Next, change to the local run directory */

  sprintf(mydir, "id%d", my_id);

  baton_start(MAX_FILE_OP, ch_rundir1_tag);

  mkdir(mydir, 0775); /* May return an error, e.g. the directory exists */
  if(chdir(mydir)){
    peg_perr(-1,"[change_rundir]: Cannot change directory to \"%s\"\n",mydir);
    err = 1;
  }

  baton_stop(MAX_FILE_OP, ch_rundir1_tag);

  /* Did anyone fail to make and change to the local run directory? */
  rerr = MPI_Allreduce(&err, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(rerr) peg_perr(-1,"[change_rundir]: MPI_Allreduce error = %d\n",rerr);

  if(rerr || gerr){
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(EXIT_FAILURE);
  }

#else /* Serial job */

  if(name == NULL || *name == '\0') return;

  mkdir(name, 0775); /* May return an error, e.g. the directory exists */
  if(chdir(name))
    peg_error("[change_rundir]: Cannot change directory to \"%s\"\n",name);

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void usage(const char *prog)
 *  \brief Outputs help
 *
 *    pegasus_version is hardwired at beginning of this file
 *    CONFIGURE_DATE is macro set when configure script runs
 */

static void usage(const char *prog)
{
  peg_perr(-1,"Pegasus %s\n",pegasus_version);
  peg_perr(-1,"  Last configure: %s\n",CONFIGURE_DATE);
  peg_perr(-1,"\nUsage: %s [options] [block/par=value ...]\n",prog);
  peg_perr(-1,"\nOptions:\n");
  peg_perr(-1,"  -i <file>       Alternate input file [peginput]\n");
  peg_perr(-1,"  -r <file>       Restart a simulation with this file\n");
  peg_perr(-1,"  -d <directory>  Alternate run dir [current dir]\n");
  peg_perr(-1,"  -h              This Help, and configuration settings\n");
  peg_perr(-1,"  -n              Parse input, but don't run program\n");
  peg_perr(-1,"  -c              Show Configuration details and quit\n");
  peg_perr(-1,"  -t hh:mm:ss     With MPI, wall time limit for final output\n");
  show_config();
  exit(0);
}
