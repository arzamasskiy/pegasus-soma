#-------------------------------------------------------------------------------
# configure.ac: Pegasus configuration script template file.  Run 'autoconf'
# to generate a new 'configure' script from this file.
#
# When configure is run, it uses the command line options and default settings
# to create custom versions of the following files:
#    Makefile, Makeoption, src/Makefile, src/defs.h, src/*/Makefile,
#
# If new configure options are added to this file, they should also be
# added in show_config() in src/show_config.c, and to the diagnostic message
# output at the end of this file.
#
# This template adds the following options to the 'configure' script:
#
# PHYSICS "packages":
#   --with-problem='name'                 (use prob/name.c as problem generator)
#   --with-eos=[isothermal, adiabatic]                       (equation of state)
#   --with-driving='name'	      (use src/particles/name.c driving routine)
#
# PHYSICS "features":
#   --enable-resistivity                                  (explicit resistivity)
#
# ALGORITHM "packages":
#   --with-cflags=[opt,debug,profile]                       (set compiler flags)
#
# ALGORITHM "features":
#   --enable-fargo                                      (enable FARGO algorithm)
#   --enable-ghost                      (write out ghost cells in outputs/dumps)
#   --enable-mpi                                          (parallelize with MPI)
#   --enable-shearing box                    (include shearing box source terms)
#   --enable-single                                 (double or single precision)
#   --enable-deltaf                                 (employ delta-f formulation)
#   --enable-fft	        (compile and link with FFTW block decomposition)
#   --enable-variable-dt	 (enable variable timestep)
#
#-------------------------------------------------------------------------------
# generic things

AC_INIT(src/pegasus.h)
AC_CONFIG_HEADER(src/config.h)
AC_PROG_CC
AC_PROG_INSTALL

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: name of problem generator from src/prob/name.c 
#   --with-problem=name (default is linear_wave1d)

AC_SUBST(PROBLEM)
AC_ARG_WITH(problem,
	[--with-problem=PROB  Problem generator from src/prob/name.c],
	with_problem=$withval, with_problem="cpaw1d")
PROBLEM=$with_problem
if test -e src/prob/$PROBLEM.c; then
   rm -f src/problem.c
   ln -s prob/$PROBLEM.c src/problem.c
else
  ls src/prob
  AC_MSG_ERROR([Invalid problem filename, valid names listed above])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: equation of state
#   --with-eos=[adiabatic,isothermal] (default is isothermal)

AC_SUBST(EOS)
AC_ARG_WITH(eos,
	[--with-eos=EOS  Equation of state (adiabatic,isothermal)],
	with_eos=$withval, with_eos=isothermal)
if test   "$with_eos" = "adiabatic";  then
  EOS="ADIABATIC"
elif test "$with_eos" = "isothermal"; then
  EOS="ISOTHERMAL"
else
  AC_MSG_ERROR([expected --with-eos=adiabatic or isothermal])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: turn on driving
#   --with-driving=name (default is none)

AC_SUBST(DRIVING_MODE)
AC_ARG_WITH(driving,
	[--with-driving=DRIV  Driving Function defined in  src/particles/driving/name.c],
	with_driving=$withval, with_driving="none")
DRIVING_ROUTINE=$with_driving
if test "$with_driving" = "none"; then
  DRIVING_MODE="NO_DRIVING"
  DRIVING_MODE_USER="OFF"
else
  DRIVING_MODE="DRIVING"
  DRIVING_MODE_USER="ON"
fi
if test -e src/particles/driving/$DRIVING_ROUTINE.c; then
    rm -f src/particles/driving.c
    ln -s driving/$DRIVING_ROUTINE.c src/particles/driving.c
else
  ls src/particles/driving
  AC_MSG_ERROR([Invalid driving filename, valid names listed above])
fi

#-------------------------------------------------------------------------------
# PHYSICS FEATURE: explicit resistivity
#  --enable-resistivity
 
AC_SUBST(RESISTIVITY_MODE)
AC_ARG_ENABLE(resistivity,
	[--enable-resistivity  enable resistivity (default is no)],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  RESISTIVITY_MODE="RESISTIVITY"
  RESISTIVITY_MODE_USER="ON"
else
  RESISTIVITY_MODE="NO_RESISTIVITY"
  RESISTIVITY_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM PACKAGE: set compiler options.
#   --with-cflags=[opt,debug,profile] (default is opt)

AC_SUBST(DEBUG_MODE)
AC_SUBST(COMPILER_OPTS)

AC_ARG_WITH(cflags,
	[--with-cflags=DEBUG_MODE  Set debugging mode (opt,debug,profile)],
	with_debug=$withval, with_debug="opt")
if test "$with_debug" = "debug"; then
  COMPILER_OPTS="-g -Wall -W -ansi -pedantic"
  DEBUG_MODE="DEBUG"
elif test "$with_debug" = "profile"; then
  COMPILER_OPTS="-O3 -pg"
  DEBUG_MODE="PROFILE"
elif test "$with_debug" = "opt"; then
  COMPILER_OPTS="-O3"
  DEBUG_MODE="OPTIMIZE"
else
  AC_MSG_ERROR([expected --with-cflags=opt,debug, or profile])
fi


#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: precision of floating point arithmetic
#   --enable-single (default is double)

AC_SUBST(PRECISION)
AC_ARG_ENABLE(single,
	[--enable-single  single-precision (default is double)], 
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  PRECISION="SINGLE_PREC"
else
  PRECISION="DOUBLE_PREC"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: write ghost cells in outputs/dumps
#   --enable-ghost

AC_SUBST(WRITE_GHOST_MODE)
AC_ARG_ENABLE(ghost,
	[--enable-ghost  write ghost zones],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  WRITE_GHOST_MODE="WRITE_GHOST_CELLS" 
  WRITE_GHOST_MODE_USER="ON" 
else
  WRITE_GHOST_MODE="NO_WRITE_GHOST_CELLS"
  WRITE_GHOST_MODE_USER="OFF" 
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: parallelize with MPI, --enable-mpi (default is no MPI)

AC_SUBST(MPI_MODE)
AC_ARG_ENABLE(mpi,
	[--enable-mpi  enable MPI parellelization],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  MPI_MODE="MPI_PARALLEL"
  MPI_MODE_USER="ON"
else
  MPI_MODE="NO_MPI_PARALLEL"
  MPI_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on shearing box evolution
#   --enable-shearing-box
# Note shearing box only works with CTU+CT and HYBRID integrators.  Prints error message if
# shearing box is enabled with VL+CT integrator.

AC_SUBST(SHEARING_BOX_MODE)
AC_ARG_ENABLE(shearing-box,
	[--enable-shearing-box  turn on shearing-box],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
SHEARING_BOX_MODE="SHEARING_BOX"
SHEARING_BOX_MODE_USER="ON"
else
  SHEARING_BOX_MODE="NO_SHEARING_BOX"
  SHEARING_BOX_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on FARGO algorithm
#   --enable-fargo

AC_SUBST(FARGO_MODE)
AC_ARG_ENABLE(fargo,
	[--enable-fargo  turn on fargo],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  FARGO_MODE="FARGO"
  FARGO_MODE_USER="ON"
else
  FARGO_MODE="NO_FARGO"
  FARGO_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------                         
# ALGORITHM FEATURE: turn on FFT libraries                                                               
#   --enable-fft                                                                                         

AC_SUBST(FFT_MODE)
AC_ARG_ENABLE(fft,
        [--enable-fft  compile and link FFT interface code (requires FFTW)],
        ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  FFT_MODE="FFT_ENABLED"
  FFT_MODE_USER="ON"
else
  FFT_MODE="NO_FFT"
  FFT_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on delta-f algorithm
#   --enable-deltaf

AC_SUBST(DIST_FUNC_MODE)
AC_ARG_ENABLE(deltaf,
	[--enable-deltaf  turn on deltaf],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  DIST_FUNC_MODE="DELTA_F"
  DIST_FUNC_MODE_USER="DELTA-F"
else
  DIST_FUNC_MODE="FULL_F"
  DIST_FUNC_MODE_USER="FULL-F"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on variable timestep
#   --enable-variable-dt

AC_SUBST(DT_MODE)
AC_ARG_ENABLE(variable-dt,
	[--enable-variable-dt  variable timestep],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  DT_MODE="VARIABLE_DT"
  DT_MODE_USER="VARIABLE-DT"
else
  DT_MODE="FIXED_DT"
  DT_MODE_USER="FIXED-DT"
fi

#-------------------------------------------------------------------------------
# check for compatibility of various options

#-------------------------------------------------------------------------------
# check for compatibility of various options

#-------------------------------------------------------------------------------
# check for various library functions

AC_PROG_GCC_TRADITIONAL
AC_CHECK_FUNCS(strdup)
AC_C_BIGENDIAN

#-------------------------------------------------------------------------------
# date and time of configure, becomes a macro in defs.h

A_CONFIGURE_DATE="`date`"
AC_SUBST(A_CONFIGURE_DATE)

#-------------------------------------------------------------------------------
# write final diagnostic output

AC_SUBST(WARNING1)
AC_SUBST(WARNING2)
WARNING1="WARNING! This file has been automatically generated by configure."
WARNING2="Any changes to it will be overwritten the next time configure is run."

AC_CONFIG_FILES(Makefile Makeoptions src/Makefile src/defs.h src/fftsrc/Makefile src/integrators/Makefile src/microphysics/Makefile src/particles/Makefile)
AC_OUTPUT

echo ""
echo "Your Pegasus distribution has now been configured:"
echo ""
echo "--- Physics Options ---"
echo ""
echo "Problem:                 $PROBLEM"
echo "Electron Eqn of State:   $EOS"
echo "Resistivity:             $RESISTIVITY_MODE_USER"
echo "Driving:                 $DRIVING_ROUTINE"
echo ""
echo "--- Algorithm Options ----"
echo ""
echo "Precision:               $PRECISION"
echo "Compiler options:        $COMPILER_OPTS"
echo "Ghost cell output:       $WRITE_GHOST_MODE_USER"
echo "Parallel modes: MPI      $MPI_MODE_USER"
echo "Shearing-box:            $SHEARING_BOX_MODE_USER"
echo "FARGO:                   $FARGO_MODE_USER"
echo "Distribution Function:   $DIST_FUNC_MODE_USER"
echo "FFT:                     $FFT_MODE_USER"
echo "Variable dt:             $DT_MODE_USER"
echo ""
