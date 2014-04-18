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

zlib and bzlib libraries are optional prerequisites.
Building is typical::

    ./configure [OPTIONS]
    make
    sudo make install

If you don't have zlib or bzlib libraries installed, use ``--without-zlib``
or ``--without-bzlib`` options, respectively -- the programs will work fine,
just won't be able to read ``.gz`` and ``.bz2`` files.

The ``make`` command builds debyer and a few programs with the prefix ``dbr_``.

The ``make install`` command installs the programs, except for a couple of
undocumented and unfinished ones.

A few configure options have been introduced to speed up the debyer program:

* ``--enable-mpi`` builds parallel (MPI) version of debyer
* ``--enable-single`` uses single precision (double precision is default)
* ``CFLAGS`` - the debyer program can be safely compiled with GCC ``-ffast-math`` option.


Alternative compilation
-----------------------

If you really want you may build the programs you need manually instead
of using autotools and make. For example::

    cd debyer

    gengetopt -i debyer.ggo
    c++ -O3 -ffast-math -DVERSION='"0.3"' -o debyer cmdline.c cmdline.h \
        main.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

    gengetopt -i extend_cmd.ggo
    c++ -O3 -DVERSION='"0.3"' -o dbr_extend dbr_extend.cc extend_cmd.c \
        cells.cc debyer.c lineio.cc fileio.cc atomtables.c -lm

    c++ -O3 -o dbr_conv conv.cc debyer.c lineio.cc fileio.cc atomtables.c -lm
