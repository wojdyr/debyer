Debyer and companion programs analyze and manipulate atomistic models.
In particular, debyer can calculate powder diffraction pattern of virtual
sample using the Debye scattering formula.
This project used to be hosted at google code.

[![Travis Status](https://travis-ci.org/wojdyr/debyer.svg?branch=master)](https://travis-ci.org/wojdyr/debyer)
[![Coverity Status](https://scan.coverity.com/projects/4820/badge.svg)](https://scan.coverity.com/projects/4820)

Documentation
=============

is kept in doc/ and can be read at http://debyer.readthedocs.org/

Installation
============

When installing from git, you must have `autoconf`, `automake` and
`gengetopt`, and run first:

    autoreconf -i

The rest of the build procedure is:

    ./configure CFLAGS="-O3 -ffast-math" [MORE OPTIONS]
    make
    sudo make install

Zlib and bzlib libraries are optional prerequisites.
If you don't have them installed, use `--without-zlib`
or `--without-bzlib` options, respectively -- the programs will work fine,
just won't be able to read `.gz` and `.bz2` files.

The `make` command builds `debyer` and a few programs with
names staring with `dbr_`.

`make install` copies the programs to *PREFIX*/bin
(*PREFIX* is `/usr/local` unless set explicitly with
`configure --prefix=...`)

A few configure options have been introduced to speed up Debyer.
Defaults are: OpenMP, double precision and `-O2 -g` compiler flags.
If your compiler has option `-ffast-math`, adding it to `CFLAGS`
is safe and makes the program notably faster.

* `--disable-openmp` builds serial version
* `--enable-mpi` builds parallel version using MPI (not OpenMP)
  and LDFLAGS turns on OpenMP-based parallelization
* `--enable-single` uses single precision (double precision is default)
* `CFLAGS` -- Debyer includes both C++ and C code. The computationally
  intensive part happens to be written in C, so `CFLAGS` matters more than
  `CXXFLAGS`. `CFLAGS="-O3 -ffast-math -march=native"` should give
  good results.

If Debyer is compiled with OpenMP, environment variable `OMP_NUM_THREADS`
can be used to control the number of threads.


Alternatively,
--------------

you may compile Debyer without autotools and make,
directly calling your compiler:

    cd debyer

    gengetopt -i debyer.ggo
    c++ -O3 -ffast-math -fopenmp -o debyer cmdline.c debyer.c atomtables.c \
        main.cc lineio.cc fileio.cc -lm

and the same for other utilities:

    c++ -O3 -o dbr_conv conv.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

    gengetopt -i extend_cmd.ggo
    c++ -O3 -DVERSION='"0.4"' -o dbr_extend dbr_extend.cc extend_cmd.c \
        cells.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

