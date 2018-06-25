# ===========================================================================
#            http://www.nongnu.org/autoconf-archive/ax_openmp.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_OPENMP([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use OpenMP a
#   standard API and set of compiler directives for parallel programming
#   (see http://www-unix.mcs/)
#
#   On success, it sets the OPENMP_CFLAGS/OPENMP_CXXFLAGS/OPENMP_F77FLAGS
#   output variable to the flag (e.g. -omp) used both to compile *and* link
#   OpenMP programs in the current language.
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well.
#
#   If you want to compile everything with OpenMP, you should set:
#
#       CFLAGS="$CFLAGS $OPENMP_CFLAGS"
#       #OR#  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
#       #OR#  FFLAGS="$FFLAGS $OPENMP_FFLAGS"
#
#   (depending on the selected language).
#
#   The user can override the default choice by setting the corresponding
#   environment variable (e.g. OPENMP_CFLAGS).
#
#   ACTION-IF-FOUND is a list of shell commands to run if an OpenMP flag is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_OPENMP.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_OPENMP], [
AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp=unknown
# Flags to try:  -fopenmp (gcc), -openmp (icc), -mp (SGI & PGI),
#                -xopenmp (Sun), -omp (Tru64), -qsmp=omp (AIX), none
ax_openmp_flags="-fopenmp -openmp -mp -xopenmp -omp -qsmp=omp none"
if test "x$OPENMP_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_flags="$OPENMP_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flags"
fi
for ax_openmp_flag in $ax_openmp_flags; do
  case $ax_openmp_flag in
    none) []_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[] ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flag" ;;
  esac
  AC_TRY_LINK_FUNC(omp_set_num_threads,
	[ax_cv_[]_AC_LANG_ABBREV[]_openmp=$ax_openmp_flag; break])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xnone"; then
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp
  fi
  m4_default([$1], [AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
])dnl AX_OPENMP

# ===========================================================================
#          http://www.nongnu.org/autoconf-archive/ax_cc_maxopt.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CC_MAXOPT
#
# DESCRIPTION
#
#   Try to turn on "good" C optimization flags for various compilers and
#   architectures, for some definition of "good". (In our case, good for
#   FFTW and hopefully for other scientific codes. Modify as needed.)
#
#   The user can override the flags by setting the CFLAGS environment
#   variable. The user can also specify --enable-portable-binary in order to
#   disable any optimization flags that might result in a binary that only
#   runs on the host architecture.
#
#   Note also that the flags assume that ANSI C aliasing rules are followed
#   by the code (e.g. for gcc's -fstrict-aliasing), and that floating-point
#   computations can be re-ordered as needed.
#
#	G. Lesur Modif:
#   Includes icc11 option -xHOST test for newer Intel processors
#
#   Requires macros: AX_CHECK_COMPILER_FLAGS, AX_COMPILER_VENDOR,
#   AX_GCC_ARCHFLAG, AX_GCC_X86_CPUID.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#   Copyright (c) 2009 Geoffroy Lesur
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_CC_MAXOPT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_ENABLE(portable-binary, [AC_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])],
	acx_maxopt_portable=$withval, acx_maxopt_portable=no)

# Try to determine "good" native compiler flags if none specified via CFLAGS
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS=""
  case $ax_cv_c_compiler_vendor in
    dec) CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
	 if test "x$acx_maxopt_portable" = xno; then
           CFLAGS="$CFLAGS -arch host"
         fi;;

    sun) CFLAGS="-native -fast -xO5 -dalign"
	 if test "x$acx_maxopt_portable" = xyes; then
	   CFLAGS="$CFLAGS -xarch=generic"
         fi;;

    hp)  CFLAGS="+Oall +Optrs_ansi +DSnative"
	 if test "x$acx_maxopt_portable" = xyes; then
	   CFLAGS="$CFLAGS +DAportable"
	 fi;;

    ibm) if test "x$acx_maxopt_portable" = xno; then
           xlc_opt="-qarch=auto -qtune=auto"
	 else
           xlc_opt="-qtune=auto"
	 fi
         AX_CHECK_COMPILER_FLAGS($xlc_opt,
         	CFLAGS="-O3 -qansialias -w $xlc_opt",
               [CFLAGS="-O3 -qansialias -w"
                echo "******************************************************"
                echo "*  You seem to have the IBM  C compiler.  It is      *"
                echo "*  recommended for best performance that you use:    *"
                echo "*                                                    *"
                echo "*    CFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
                echo "*                      ^^^        ^^^                *"
                echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
                echo "*  CPU you have.  (Set the CFLAGS environment var.   *"
                echo "*  and re-run configure.)  For more info, man cc.    *"
                echo "******************************************************"])
         ;;

    intel) CFLAGS="-O3 -ansi_alias"
	if test "x$acx_maxopt_portable" = xno; then
	  icc_archflag=unknown
	  icc_flags=""
	  case $host_cpu in
	    i686*|x86_64*)
              # icc accepts gcc assembly syntax, so these should work:
	      AX_GCC_X86_CPUID(0)
              AX_GCC_X86_CPUID(1)
	      case $ax_cv_gcc_x86_cpuid_0 in # see AX_GCC_ARCHFLAG
                *:756e6547:*:*) # Intel
                  case $ax_cv_gcc_x86_cpuid_1 in
                    *6a?:*[[234]]:*:*|*6[[789b]]?:*:*:*) icc_flags="-xK";;
                    *f3[[347]]:*:*:*|*f4[1347]:*:*:*) icc_flags="-xP -xN -xW -xK";;
                    *f??:*:*:*) icc_flags="-xN -xW -xK";;
                  esac ;;
              esac ;;
          esac
          if test "x$icc_flags" != x; then
            for flag in $icc_flags; do
              AX_CHECK_COMPILER_FLAGS($flag, [icc_archflag=$flag; break])
            done
          fi
          AC_MSG_CHECKING([for icc architecture flag])
	  AC_MSG_RESULT($icc_archflag)
          if test "x$icc_archflag" != xunknown; then
            CFLAGS="$CFLAGS $icc_archflag"
		  else
		    # Test -xHOST for icc11 if we have not found anything else...
			AX_CHECK_COMPILER_FLAGS(-xHOST, CFLAGS="$CFLAGS -xHOST")
          fi
        fi
	;;

    gnu)
     # default optimization flags for gcc on all systems
     CFLAGS="-O3 -fomit-frame-pointer"

     # -malign-double for x86 systems
     AX_CHECK_COMPILER_FLAGS(-malign-double, CFLAGS="$CFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILER_FLAGS(-fstrict-aliasing,
	CFLAGS="$CFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILER_FLAGS(-ffast-math, CFLAGS="$CFLAGS -ffast-math")

     AX_GCC_ARCHFLAG($acx_maxopt_portable)
     ;;
  esac

  if test -z "$CFLAGS"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
	echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
	echo "********************************************************"
	echo ""
        CFLAGS="-O3"
  fi

  AX_CHECK_COMPILER_FLAGS($CFLAGS, [], [
	echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "********************************************************"
        echo ""
        CFLAGS=""
  ])

fi
])

# ===========================================================================
#     http://www.nongnu.org/autoconf-archive/ax_check_compiler_flags.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_COMPILER_FLAGS(FLAGS, [ACTION-SUCCESS], [ACTION-FAILURE])
#
# DESCRIPTION
#
#   Check whether the given compiler FLAGS work with the current language's
#   compiler, or whether they give an error. (Warnings, however, are
#   ignored.)
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.
#
# LICENSE
#
#   Copyright (c) 2009 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2009 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_CHECK_COMPILER_FLAGS],
[AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX
AC_MSG_CHECKING([whether _AC_LANG compiler accepts $1])
dnl Some hackery here since AC_CACHE_VAL can't handle a non-literal varname:
AS_LITERAL_IF([$1],
  [AC_CACHE_VAL(AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_[$1]), [
      ax_save_FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      _AC_LANG_PREFIX[]FLAGS="$1"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
        AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_[$1])=yes,
        AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_[$1])=no)
      _AC_LANG_PREFIX[]FLAGS=$ax_save_FLAGS])],
  [ax_save_FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
   _AC_LANG_PREFIX[]FLAGS="$1"
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
     eval AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_[$1])=yes,
     eval AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_[$1])=no)
   _AC_LANG_PREFIX[]FLAGS=$ax_save_FLAGS])
eval ax_check_compiler_flags=$AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_[$1])
AC_MSG_RESULT($ax_check_compiler_flags)
if test "x$ax_check_compiler_flags" = xyes; then
	m4_default([$2], :)
else
	m4_default([$3], :)
fi
])dnl AX_CHECK_COMPILER_FLAGS

# ===========================================================================
#       http://www.nongnu.org/autoconf-archive/ax_compiler_vendor.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_COMPILER_VENDOR
#
# DESCRIPTION
#
#   Determine the vendor of the C/C++ compiler, e.g., gnu, intel, ibm, sun,
#   hp, borland, comeau, dec, cray, kai, lcc, metrowerks, sgi, microsoft,
#   watcom, etc. The vendor is returned in the cache variable
#   $ax_cv_c_compiler_vendor for C and $ax_cv_cxx_compiler_vendor for C++.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_COMPILER_VENDOR],
[
AC_CACHE_CHECK([for _AC_LANG compiler vendor], ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor,
 [ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor=unknown
  # note: don't check for gcc first since some other compilers define __GNUC__
  for ventest in intel:__ICC,__ECC,__INTEL_COMPILER ibm:__xlc__,__xlC__,__IBMC__,__IBMCPP__ pathscale:__PATHCC__,__PATHSCALE__ gnu:__GNUC__ sun:__SUNPRO_C,__SUNPRO_CC hp:__HP_cc,__HP_aCC dec:__DECC,__DECCXX,__DECC_VER,__DECCXX_VER borland:__BORLANDC__,__TURBOC__ comeau:__COMO__ cray:_CRAYC kai:__KCC lcc:__LCC__ metrowerks:__MWERKS__ sgi:__sgi,sgi microsoft:_MSC_VER watcom:__WATCOMC__ portland:__PGI; do
    vencpp="defined("`echo $ventest | cut -d: -f2 | sed 's/,/) || defined(/g'`")"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
#if !($vencpp)
      thisisanerror;
#endif
])], [ax_cv_]_AC_LANG_ABBREV[_compiler_vendor=`echo $ventest | cut -d: -f1`; break])
  done
 ])
])

# ===========================================================================
#         http://www.nongnu.org/autoconf-archive/ax_gcc_archflag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_GCC_ARCHFLAG([PORTABLE?], [ACTION-SUCCESS], [ACTION-FAILURE])
#
# DESCRIPTION
#
#   This macro tries to guess the "native" arch corresponding to the target
#   architecture for use with gcc's -march=arch or -mtune=arch flags. If
#   found, the cache variable $ax_cv_gcc_archflag is set to this flag and
#   ACTION-SUCCESS is executed; otherwise $ax_cv_gcc_archflag is is set to
#   "unknown" and ACTION-FAILURE is executed. The default ACTION-SUCCESS is
#   to add $ax_cv_gcc_archflag to the end of $CFLAGS.
#
#   PORTABLE? should be either [yes] (default) or [no]. In the former case,
#   the flag is set to -mtune (or equivalent) so that the architecture is
#   only used for tuning, but the instruction set used is still portable. In
#   the latter case, the flag is set to -march (or equivalent) so that
#   architecture-specific instructions are enabled.
#
#   The user can specify --with-gcc-arch=<arch> in order to override the
#   macro's choice of architecture, or --without-gcc-arch to disable this.
#
#   When cross-compiling, or if $CC is not gcc, then ACTION-FAILURE is
#   called unless the user specified --with-gcc-arch manually.
#
#   Requires macros: AX_CHECK_COMPILER_FLAGS, AX_GCC_X86_CPUID
#
#   (The main emphasis here is on recent CPUs, on the principle that doing
#   high-performance computing on old hardware is uncommon.)
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_GCC_ARCHFLAG],
[AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_WITH(gcc-arch, [AC_HELP_STRING([--with-gcc-arch=<arch>], [use architecture <arch> for gcc -march/-mtune, instead of guessing])],
	ax_gcc_arch=$withval, ax_gcc_arch=yes)

AC_MSG_CHECKING([for gcc architecture flag])
AC_MSG_RESULT([])
AC_CACHE_VAL(ax_cv_gcc_archflag,
[
ax_cv_gcc_archflag="unknown"

if test "$GCC" = yes; then

if test "x$ax_gcc_arch" = xyes; then
ax_gcc_arch=""
if test "$cross_compiling" = no; then
case $host_cpu in
  i[[3456]]86*|x86_64*) # use cpuid codes, in part from x86info-1.7 by D. Jones
     AX_GCC_X86_CPUID(0)
     AX_GCC_X86_CPUID(1)
     case $ax_cv_gcc_x86_cpuid_0 in
       *:756e6547:*:*) # Intel
          case $ax_cv_gcc_x86_cpuid_1 in
	    *5[[48]]?:*:*:*) ax_gcc_arch="pentium-mmx pentium" ;;
	    *5??:*:*:*) ax_gcc_arch=pentium ;;
	    *6[[3456]]?:*:*:*) ax_gcc_arch="pentium2 pentiumpro" ;;
	    *6a?:*[[01]]:*:*) ax_gcc_arch="pentium2 pentiumpro" ;;
	    *6a?:*[[234]]:*:*) ax_gcc_arch="pentium3 pentiumpro" ;;
	    *6[[9d]]?:*:*:*) ax_gcc_arch="pentium-m pentium3 pentiumpro" ;;
	    *6[[78b]]?:*:*:*) ax_gcc_arch="pentium3 pentiumpro" ;;
	    *6??:*:*:*) ax_gcc_arch=pentiumpro ;;
            *f3[[347]]:*:*:*|*f4[1347]:*:*:*)
		case $host_cpu in
                  x86_64*) ax_gcc_arch="nocona pentium4 pentiumpro" ;;
                  *) ax_gcc_arch="prescott pentium4 pentiumpro" ;;
                esac ;;
            *f??:*:*:*) ax_gcc_arch="pentium4 pentiumpro";;
          esac ;;
       *:68747541:*:*) # AMD
          case $ax_cv_gcc_x86_cpuid_1 in
	    *5[[67]]?:*:*:*) ax_gcc_arch=k6 ;;
	    *5[[8d]]?:*:*:*) ax_gcc_arch="k6-2 k6" ;;
	    *5[[9]]?:*:*:*) ax_gcc_arch="k6-3 k6" ;;
	    *60?:*:*:*) ax_gcc_arch=k7 ;;
	    *6[[12]]?:*:*:*) ax_gcc_arch="athlon k7" ;;
	    *6[[34]]?:*:*:*) ax_gcc_arch="athlon-tbird k7" ;;
	    *67?:*:*:*) ax_gcc_arch="athlon-4 athlon k7" ;;
	    *6[[68a]]?:*:*:*)
	       AX_GCC_X86_CPUID(0x80000006) # L2 cache size
	       case $ax_cv_gcc_x86_cpuid_0x80000006 in
                 *:*:*[[1-9a-f]]??????:*) # (L2 = ecx >> 16) >= 256
			ax_gcc_arch="athlon-xp athlon-4 athlon k7" ;;
                 *) ax_gcc_arch="athlon-4 athlon k7" ;;
	       esac ;;
	    *f[[4cef8b]]?:*:*:*) ax_gcc_arch="athlon64 k8" ;;
	    *f5?:*:*:*) ax_gcc_arch="opteron k8" ;;
	    *f7?:*:*:*) ax_gcc_arch="athlon-fx opteron k8" ;;
	    *f??:*:*:*) ax_gcc_arch="k8" ;;
          esac ;;
	*:746e6543:*:*) # IDT
	   case $ax_cv_gcc_x86_cpuid_1 in
	     *54?:*:*:*) ax_gcc_arch=winchip-c6 ;;
	     *58?:*:*:*) ax_gcc_arch=winchip2 ;;
	     *6[[78]]?:*:*:*) ax_gcc_arch=c3 ;;
	     *69?:*:*:*) ax_gcc_arch="c3-2 c3" ;;
	   esac ;;
     esac
     if test x"$ax_gcc_arch" = x; then # fallback
	case $host_cpu in
	  i586*) ax_gcc_arch=pentium ;;
	  i686*) ax_gcc_arch=pentiumpro ;;
        esac
     fi
     ;;

  sparc*)
     AC_PATH_PROG([PRTDIAG], [prtdiag], [prtdiag], [$PATH:/usr/platform/`uname -i`/sbin/:/usr/platform/`uname -m`/sbin/])
     cputype=`(((grep cpu /proc/cpuinfo | cut -d: -f2) ; ($PRTDIAG -v |grep -i sparc) ; grep -i cpu /var/run/dmesg.boot ) | head -n 1) 2> /dev/null`
     cputype=`echo "$cputype" | tr -d ' -' |tr $as_cr_LETTERS $as_cr_letters`
     case $cputype in
         *ultrasparciv*) ax_gcc_arch="ultrasparc4 ultrasparc3 ultrasparc v9" ;;
         *ultrasparciii*) ax_gcc_arch="ultrasparc3 ultrasparc v9" ;;
         *ultrasparc*) ax_gcc_arch="ultrasparc v9" ;;
         *supersparc*|*tms390z5[[05]]*) ax_gcc_arch="supersparc v8" ;;
         *hypersparc*|*rt62[[056]]*) ax_gcc_arch="hypersparc v8" ;;
         *cypress*) ax_gcc_arch=cypress ;;
     esac ;;

  alphaev5) ax_gcc_arch=ev5 ;;
  alphaev56) ax_gcc_arch=ev56 ;;
  alphapca56) ax_gcc_arch="pca56 ev56" ;;
  alphapca57) ax_gcc_arch="pca57 pca56 ev56" ;;
  alphaev6) ax_gcc_arch=ev6 ;;
  alphaev67) ax_gcc_arch=ev67 ;;
  alphaev68) ax_gcc_arch="ev68 ev67" ;;
  alphaev69) ax_gcc_arch="ev69 ev68 ev67" ;;
  alphaev7) ax_gcc_arch="ev7 ev69 ev68 ev67" ;;
  alphaev79) ax_gcc_arch="ev79 ev7 ev69 ev68 ev67" ;;

  powerpc*)
     cputype=`((grep cpu /proc/cpuinfo | head -n 1 | cut -d: -f2 | cut -d, -f1 | sed 's/ //g') ; /usr/bin/machine ; /bin/machine; grep CPU /var/run/dmesg.boot | head -n 1 | cut -d" " -f2) 2> /dev/null`
     cputype=`echo $cputype | sed -e 's/ppc//g;s/ *//g'`
     case $cputype in
       *750*) ax_gcc_arch="750 G3" ;;
       *740[[0-9]]*) ax_gcc_arch="$cputype 7400 G4" ;;
       *74[[4-5]][[0-9]]*) ax_gcc_arch="$cputype 7450 G4" ;;
       *74[[0-9]][[0-9]]*) ax_gcc_arch="$cputype G4" ;;
       *970*) ax_gcc_arch="970 G5 power4";;
       *POWER4*|*power4*|*gq*) ax_gcc_arch="power4 970";;
       *POWER5*|*power5*|*gr*|*gs*) ax_gcc_arch="power5 power4 970";;
       603ev|8240) ax_gcc_arch="$cputype 603e 603";;
       *) ax_gcc_arch=$cputype ;;
     esac
     ax_gcc_arch="$ax_gcc_arch powerpc"
     ;;
esac
fi # not cross-compiling
fi # guess arch

if test "x$ax_gcc_arch" != x -a "x$ax_gcc_arch" != xno; then
for arch in $ax_gcc_arch; do
  if test "x[]m4_default([$1],yes)" = xyes; then # if we require portable code
    flags="-mtune=$arch"
    # -mcpu=$arch and m$arch generate nonportable code on every arch except
    # x86.  And some other arches (e.g. Alpha) don't accept -mtune.  Grrr.
    case $host_cpu in i*86|x86_64*) flags="$flags -mcpu=$arch -m$arch";; esac
  else
    flags="-march=$arch -mcpu=$arch -m$arch"
  fi
  for flag in $flags; do
    AX_CHECK_COMPILER_FLAGS($flag, [ax_cv_gcc_archflag=$flag; break])
  done
  test "x$ax_cv_gcc_archflag" = xunknown || break
done
fi

fi # $GCC=yes
])
AC_MSG_CHECKING([for gcc architecture flag])
AC_MSG_RESULT($ax_cv_gcc_archflag)
if test "x$ax_cv_gcc_archflag" = xunknown; then
  m4_default([$3],:)
else
  m4_default([$2], [CFLAGS="$CFLAGS $ax_cv_gcc_archflag"])
fi
])

# ===========================================================================
#        http://www.nongnu.org/autoconf-archive/ax_gcc_x86_cpuid.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_GCC_X86_CPUID(OP)
#
# DESCRIPTION
#
#   On Pentium and later x86 processors, with gcc or a compiler that has a
#   compatible syntax for inline assembly instructions, run a small program
#   that executes the cpuid instruction with input OP. This can be used to
#   detect the CPU type.
#
#   On output, the values of the eax, ebx, ecx, and edx registers are stored
#   as hexadecimal strings as "eax:ebx:ecx:edx" in the cache variable
#   ax_cv_gcc_x86_cpuid_OP.
#
#   If the cpuid instruction fails (because you are running a
#   cross-compiler, or because you are not using gcc, or because you are on
#   a processor that doesn't have this instruction), ax_cv_gcc_x86_cpuid_OP
#   is set to the string "unknown".
#
#   This macro mainly exists to be used in AX_GCC_ARCHFLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_GCC_X86_CPUID],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for x86 cpuid $1 output, ax_cv_gcc_x86_cpuid_$1,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>], [
     int op = $1, eax, ebx, ecx, edx;
     FILE *f;
      __asm__("cpuid"
        : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
        : "a" (op));
     f = fopen("conftest_cpuid", "w"); if (!f) return 1;
     fprintf(f, "%x:%x:%x:%x\n", eax, ebx, ecx, edx);
     fclose(f);
     return 0;
])],
     [ax_cv_gcc_x86_cpuid_$1=`cat conftest_cpuid`; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown])])
AC_LANG_POP([C])
])

# ===========================================================================
#             http://www.nongnu.org/autoconf-archive/acx_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   ACX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
#   On success, it sets the MPICC, MPICXX, MPIF77, or MPIFC output variable
#   to the name of the MPI compiler, depending upon the current language.
#   (This may just be $CC/$CXX/$F77/$FC, but is more often something like
#   mpicc/mpiCC/mpif77/mpif90.) It also sets MPILIBS to any libraries that
#   are needed for linking MPI (e.g. -lmpi or -lfmpi, if a special
#   MPICC/MPICXX/MPIF77/MPIFC was not found).
#
#   If you want to compile everything with MPI, you should set:
#
#       CC="MPICC" #OR# CXX="MPICXX" #OR# F77="MPIF77" #OR# FC="MPIFC"
#       LIBS="$MPILIBS $LIBS"
#
#   NOTE: The above assumes that you will use $CC (or whatever) for linking
#   as well as for compiling. (This is the default for automake and most
#   Makefiles.)
#
#   The user can force a particular library/compiler by setting the
#   MPICC/MPICXX/MPIF77/MPIFC and/or MPILIBS environment variables.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an MPI library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run if it is not
#   found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_MPI.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	acx_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	acx_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)
],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf_r mpxlf mpf77 cmpifc, $F77)
	acx_mpi_save_F77="$F77"
	F77="$MPIF77"
	AC_SUBST(MPIF77)
],
[Fortran], [
	AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
	acx_mpi_save_FC="$FC"
	FC="$MPIFC"
	AC_SUBST(MPIFC)
])

if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
		[Fortran], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
AC_LANG_CASE([Fortran 77], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpich, MPI_Init, [MPILIBS="-lfmpich"])
	fi
],
[Fortran], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
	fi
])
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran 77], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
	[C++], [CXX="$acx_mpi_save_CXX"],
	[Fortran 77], [F77="$acx_mpi_save_F77"],
	[Fortran], [FC="$acx_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI

