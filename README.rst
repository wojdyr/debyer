Debyer and companion programs analyze and manipulate atomistic models.
In particular, debyer can calculate powder diffraction pattern of virtual
sample using the Debye scattering formula.
This project used to be hosted at google code.

|travis-status|_

.. _travis-status: https://travis-ci.org/wojdyr/debyer/
.. |travis-status| image:: https://api.travis-ci.org/wojdyr/debyer.png

Documentation
=============

is kept in doc/ and can be read at http://debyer.readthedocs.org/

Installation
============

When installing from git, you must have ``autoconf``, ``automake`` and
``gengetopt``, and run first::

    autoreconf -i

The rest of the build procedure is::

    ./configure CFLAGS="-O3 -ffast-math" [--prefix=...] [MORE OPTIONS]
    make
    sudo make install

Zlib and bzlib libraries are optional prerequisites.
If you don't have them installed, use ``--without-zlib``
or ``--without-bzlib`` options, respectively -- the programs will work fine,
just won't be able to read ``.gz`` and ``.bz2`` files.

The ``make`` command builds ``debyer`` and a few programs with
names staring with ``dbr_``.

``make install`` copies the programs to *PREFIX*/bin
(*PREFIX* is ``/usr/local`` unless set explicitly with ``--prefix``).

A few configure options have been introduced to speed up Debyer:

* ``--enable-mpi`` builds parallel (MPI) version
* ``-fopenmp`` flag (or equivalent for your compiler) added to CFLAGS
  and LDFLAGS turns on OpenMP-based parallelization
* ``--enable-single`` uses single precision (double precision is default)
* ``CFLAGS`` - the debyer program can be safely compiled with the GCC/Clang
  ``-ffast-math`` option. You may also use ``-O3 -march=native``.


Alternatively,
--------------

you may compile Debyer without autotools and make,
directly calling your compiler::

    cd debyer

    gengetopt -i debyer.ggo
    c++ -O3 -ffast-math -fopenmp -DVERSION='"0.3"' -o debyer cmdline.c \
        main.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

and the same for other utilities::

    c++ -O3 -o dbr_conv conv.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

    gengetopt -i extend_cmd.ggo
    c++ -O3 -DVERSION='"0.3"' -o dbr_extend dbr_extend.cc extend_cmd.c \
        cells.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

