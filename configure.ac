
AC_PREREQ(2.59)
AC_INIT([debyer], [0.3], [])
AC_CONFIG_SRCDIR([debyer/debyer.c])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([1.9 dist-bzip2 foreign])

AM_MAINTAINER_MODE

### MPI or not
AC_ARG_ENABLE(mpi,
 [  --enable-mpi                  compile for parallel runs using MPI],,enable_mpi=no)

AC_ARG_ENABLE(double,
 [  --enable-double               use double precision for atom coordinates],,enable_double=no)

# optional zlib and bzlib libraries
AC_ARG_WITH(zlib, 
 [  --without-zlib                disable zlib support (reading .gz files)])
test -z "$with_zlib" && with_zlib="yes"  

AC_ARG_WITH(bzlib, 
 [  --without-bzlib               disable bzlib support (reading .bz2 files)])
test -z "$with_bzlib" && with_bzlib="yes"  

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

if test "$enable_mpi" = "yes"; then

  AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)
  AC_CHECK_PROGS(MPICXX, mpic++ mpiCC mpicxx mpCC hcp mpxlC mpxlC_r cmpic++, $CXX)
  # now change the normal cc to the MPI one - see the comment above.
  CC=$MPICC
  CXX=$MPICXX
  AC_MSG_CHECKING([whether the MPI cc command works]) # be paranoid
  AC_TRY_LINK([#include <mpi.h>],
	      [int argc; char **argv; MPI_Init(&argc,&argv);],
	      AC_MSG_RESULT([yes]),
	      AC_MSG_ERROR([Cannot compile and link MPI code with $CC]))
  AC_DEFINE(USE_MPI, 1, [Make a parallel version of Debyer using MPI])
fi

if test "$enable_double" = "yes"; then
  AC_DEFINE(USE_DOUBLE, 1, [Use double precision for atom coordinates])
fi

# Checks for libraries.
if test "$with_zlib" = "yes"; then
  AC_CHECK_LIB([z], [gzopen],
               [AC_DEFINE([HAVE_ZLIB], [1], [Define if zlib is available])],
 	       [AC_MSG_ERROR(
  [zlib library was not found. 
   Either use flag --without-zlib or install the library.
 ])])
  LIBS="$LIBS -lz"
  AC_CHECK_HEADER([zlib.h], [], [AC_MSG_ERROR([
   zlib.h header was not found. Either use flag --without-zlib or install 
   zlib library with development files (some Linux distros have development
   files in a separate you may find them in a separate -dev or -devel package).
 ])])
fi

if test "$with_bzlib" = "yes"; then
  AC_CHECK_LIB([bz2], [BZ2_bzReadOpen], 
	       [AC_DEFINE([HAVE_BZLIB], [1], [Define if bzlib is available])],
	       [AC_MSG_ERROR(
  [bzlib library was not found. 
   Either use flag --without-bzlib or install the library.
 ])])
  LIBS="$LIBS -lbz2"
  AC_CHECK_HEADER([bzlib.h], [], [AC_MSG_ERROR([
   bzlib.h header was not found. Either use flag --without-bzlib or install 
   bzlib library with development files (some Linux distros have development
   files in a separate you may find them in a separate -dev or -devel package).
 ])])
fi

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_STRTOD
#AC_FUNC_VPRINTF
AC_CHECK_FUNCS([strchr strtol])
AC_CHECK_LIB([m], [pow], [], [AC_MSG_ERROR([pow function not found.])])
AC_CHECK_LIB([m], [floor], [], [AC_MSG_ERROR([floor function not found.])])
AC_CHECK_LIB([m], [sqrt], [], [AC_MSG_ERROR([sqrt function not found.])])

AC_OUTPUT([
Makefile
debyer/Makefile
])

