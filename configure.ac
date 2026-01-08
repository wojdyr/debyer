
AC_PREREQ([2.71])
# keep version the same in debyer/debyer.ggo and debyer/version.h
AC_INIT([debyer],[0.4],[])
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_SRCDIR([debyer/debyer.c])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([1.11 foreign subdir-objects silent-rules -Wall])

AM_MAINTAINER_MODE(enable)

### MPI or not
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--enable-mpi],
                                  [compile for parallel runs using MPI]))

AC_ARG_ENABLE(single, AS_HELP_STRING([--enable-single],
                                  [use single precision for atom coordinates]))

# optional zlib and bzlib libraries
AC_ARG_WITH(zlib, AS_HELP_STRING([--without-zlib],
                                 [disable zlib support (reading .gz files)]))

AC_ARG_WITH(bzlib, AS_HELP_STRING([--without-bzlib],
                               [disable bzlib support (reading .bz2 files)]))

# We are checking for compilers in the same way as GROMACS does. They say:
# We cannot set that [MPICC] in the PROG_CC macro [...], 
# since autoconf checks that the created file can be executed. This would
# fail on platforms where MPI executables can only be run through a batchqueue.
#
# This creates a problem when only mpi compilers are installed. In this case
# use CC and CXX options for configure.

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR

AS_IF(test x$enable_mpi = xyes, [
  AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)
  AC_CHECK_PROGS(MPICXX, mpic++ mpiCC mpicxx mpCC hcp mpxlC mpxlC_r cmpic++, $CXX)
  # now change the normal cc to the MPI one - see the comment above.
  CC=$MPICC
  CXX=$MPICXX
  AC_MSG_CHECKING([whether the MPI cc command works]) # be paranoid
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[int argc; char **argv; MPI_Init(&argc,&argv);]])],[AC_MSG_RESULT(yes)],[AC_MSG_ERROR(Cannot compile and link MPI code with $CC)])
  AC_DEFINE(USE_MPI, 1, [Make a parallel version of Debyer using MPI])
],[
  AC_OPENMP
])

if test x$enable_single = xyes; then
  AC_DEFINE(USE_SINGLE, 1, [Use single precision for atom coordinates])
fi

AC_SEARCH_LIBS(cos, m, , AC_MSG_ERROR([no math library]))

if test x$with_zlib != xno; then
  AC_CHECK_LIB(z, gzopen, , AC_MSG_ERROR([
   zlib library was not found.
   Either use flag --without-zlib or install the library.]))
  AC_CHECK_HEADER(zlib.h, [dummy=yes], AC_MSG_ERROR([
   zlib.h header was not found. Either use flag --without-zlib or install
   the zlib library including header files.]))
fi

if test x$with_bzlib != xno; then
  AC_CHECK_LIB(bz2, BZ2_bzReadOpen, , AC_MSG_ERROR([
   bzlib library was not found.
   Either use flag --without-bzlib or install the library.]))
  AC_CHECK_HEADER(bzlib.h, [dummy=yes], AC_MSG_ERROR([
   bzlib.h header was not found. Either use flag --without-bzlib or install
   bzlib library including header files.]))
fi

AC_CONFIG_FILES([Makefile debyer/Makefile])
AC_OUTPUT
