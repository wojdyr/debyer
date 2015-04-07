
You can find here a few tools to analyze and manipulate atomistic models.
The main tool, a program called :ref:`Debyer <debyer>`, is well-documented
and used in several labs. Other programs were written for very specific
tasks and may not be ever re-used -- but just in case,
they are briefly described in the :ref:`companions` section.

The `source code (C/C++) and installation instructions are on GitHub`__.

__ https://github.com/wojdyr/debyer

.. _debyer:

Debyer
######

Debyer takes as an **input**
a file with positions of all the atoms in the virtual sample
(up to tens of millions of atoms, perhaps even more) and can **output**

* x-ray and neutron powder diffraction pattern,
  total scattering structure function,
* pair distribution/correlation function (PDF/PCF)
  and related functions (RDF, reduced PDF).

A few use cases for debyer can be found `in the papers citing its website`__.
Usually, the analyzed model

* is a result of molecular dynamics simulation,
* or is constructed using simple geometrical rules,
  with no interatomic potentials and no forces.

__ http://scholar.google.com/scholar?q=unipress.waw.pl%2Fdebyer+OR+code.google.com%2Fp%2Fdebyer+OR+github.com%2Fwojdyr%2Fdebyer

Debyer is **not** designed to calculate a diffraction pattern of a perfect
crystal. It does not make sense to employ the Debye's formula for a unit cell:
you can easily get indexed pattern using other programs, such as PowderCell_
(sadly it's not maintained since 1990s;
no idea what program should be recommended nowadays).

.. _PowderCell: http://www.iucr.org/resources/other-directories/software/powder-cell

.. note::

    If the math below doesn't render properly, right-click on it and change
    ``Math Renderer`` in ``Math Settings``.
    If the math is not displayed at all and there is ``https://`` in the
    address bar, change it to ``http://`` (no ``s``).

Direct-space patterns
=====================

Unlike in experiments, the pair distribution function
(pair correlation function)
is obtained directly from the atomic coordinates,
not through the diffraction pattern.
In other words, we do not emulate errors of experimental PDFs.

The definitions of PDF/RDF/PCF vary. Here, we stick to names
and symbols from the *Underneath the Bragg Peaks* book
(`ISBN 0-08-042698-0 <http://books.google.pl/books?id=ek2ymu7_NfgC>`_):

* *g(r)* -- atomic pair distribution function (PDF),
  a.k.a. pair correlation function (PCF), it converges to 1 as *r* increases.
* *G(r)* -- reduced PDF (rPDF), which oscillates around 0,
* *R(r)* -- radial distribution function (RDF), which goes up like parabola,
* *ρ(r)* -- atomic pair density function, proportional to *g(r)*, rarely
  used but listed here for completeness.

.. math::

 R(r)    &= \frac{1}{N}\sum_{\nu}\sum_{\mu} \frac{b_{\nu}b_{\mu}}
              {\left\langle b\right\rangle ^2} \delta(r-r_{\nu\mu}) \\
 \rho(r) &= \frac{1}{4\pi r^{2}} R(r) \\
 g(r)    &= \rho(r) / \rho_0 = \frac{R(r)}{4\pi\rho_0 r^2} \\
 G(r)    &= 4\pi r\rho_0 \left[ g(r)-1 \right]

The average density `\rho_0` is either defined by the user or determined
automatically.

Average density is well-defined for a continuous system in the PBC (periodic
boundary conditions), but not for a single grain in the vacuum.
In the latter case *G(r)* will be distored.
Although corrections for particular grain shapes has been proposed
in the literature, they are implemented in this program.

Weights *b* can be set as neutron scattering lengths, x-ray *f*\ (0)
or equal 1 (see the section :ref:`scat`).

Debye scattering formula
========================

Reciprocal space patterns (X-ray and neutron powder diffraction patterns)
are calculated using the Debye scattering equation:

.. math:: I(Q)=\sum_i \sum_j f_i f_j \frac{\sin(Qr_{ij})}{Qr_{ij}}
 
where

* `Q` is the scattering vector, called also momentum transfer vector
  (`Q = \left|\boldsymbol{Q}\right| = 4 \pi \sin \theta / \lambda` , where
  `\theta` is diffraction half-angle and `\lambda` is the wavelength),
* `r_{ij}=\left|\boldsymbol{r}_i - \boldsymbol{r}_j\right|` is the distance
  between atoms *i* and *j*,
* and `f_i` is the atomic scattering factor of *i*-th atom;
  in general it depends on `Q` and we should write it as `f(Q)` ,
  but we don't, to keep the notation simple.

.. note::

    Different symbols are used in the literature. In particular,
    the `q` symbol is sometimes used instead of `Q` (with the same meaning).
    `q` may also be defined as `q = 2\sin\theta / \lambda = Q / 2\pi` .
    Or `s` may be used for the latter.

The best way to understand what this equation
calculates, what assumptions and  approximations are used, is to derive it.

Derivation
----------

Derivation of the Debye formula can be found in many places. The standard
citation is a German-language paper, `P. Debye, Ann. Physik., 46 (1915) 809
<http://onlinelibrary.wiley.com/doi/10.1002/andp.19153510606/abstract>`_.
Re-derivation can be found for example in `Farrow & Billinge, Acta Cryst A65
(2009) 232 <http://dx.doi.org/10.1107/S0108767309009714>`_
(`pre-print <http://arxiv.org/pdf/0811.1140.pdf>`_).

Let's re-derive it once more. We start with the amplitude of the
scattered wave. If we ignore non-elastic scattering and assume that each photon
is scattered only once, the amplitude is a simple sum over all atoms:

.. math::

    \Psi(\boldsymbol{Q}) = \sum_i \psi_i
    =\sum_i f_i \exp(-i \boldsymbol{Q} \cdot \boldsymbol{r}_i)

(*i* has two unrelated meanings here, index and the imaginary unit,
hopefully it's not too confusing).

The intensity of the scattered wave is given as

.. math::

    I(\boldsymbol{Q}) = \left| \Psi (\boldsymbol{Q}) \right|^2
    = \Psi(\boldsymbol{Q}) \Psi^{*}(\boldsymbol{Q})

where the complex conjugate

.. math::

    \Psi^{*}(\boldsymbol{Q}) = \sum_i \psi^{*}_i
    = \sum_i f_i^{*} \exp(i \boldsymbol{Q} \cdot \boldsymbol{r}_i)

Therefore, with the assumption that atomic form factors are real and `f^{*}=f` ,

.. math::

    I(\boldsymbol{Q})
    = \sum_i \sum_j \psi_i \psi^{*}_j
    = \sum_i \sum_j f_i f_j \exp \left[ -i \boldsymbol{Q} \cdot
    \left( \boldsymbol{r}_i-\boldsymbol{r}_j \right) \right]


The Debye scattering equation gives spherically avaraged intensity.
The averaging is similar to calculating the surface area of sphere.
It is more elegant to use double integral,
but knowing the circumference formula
we can take a shortcut and use one integral:

.. math::

    A_{sphere} = \int_0^\pi 2\pi r\sin\theta \, r \mathrm{d} \theta
    = 2\pi r^2 [ - \cos \theta ]_0^\pi
    = 4 \pi r^2

Now, let `\gamma` be the angle between `\boldsymbol{Q}` and
`\boldsymbol{r}_{ij} \, (\equiv \boldsymbol{r}_i - \boldsymbol{r}_j)` .

.. math::

    I(\boldsymbol{Q})
    = \sum_i \sum_j f_i f_j \exp \left( -i Q r_{ij} \cos \gamma \right)


Avaraging the exponential part gives

.. math::
    \left\langle \exp( -iQr_{ij} \cos \gamma ) \right\rangle
    & = \frac{1}{4\pi r_{ij}^2} \int_0^\pi \exp\left( -iQr_{ij} \cos\gamma
    \right) \, 2\pi r_{ij}^2 \sin\gamma \, \mathrm{d}\gamma
    \\
    & = \frac{1}{2} \int_0^\pi \exp\left( -iQr_{ij} \cos\gamma \right)
          \sin\gamma \, \mathrm{d}\gamma
    \\
    & = \frac{1}{2} \left[ \frac{\exp\left( -iQr_{ij} \cos\gamma \right)}
          {iQr_{ij}} \right]_0^\pi
    \\
    & = \frac{\exp(iQr_{ij})-\exp(-iQr_{ij})}{2iQr_{ij}}
    \\
    & = \frac{\sin(Qr_{ij})}{Qr_{ij}}

which proves the Debye formula.


Modifications
-------------

Usually this formula is normalized by `\frac{1}{N}` or `\frac{1}{N<f>^2}` .
TBC

include Debye-Waller factor?
the correlated broadening factor for the atom pair (as mentioned in
Farrow 2009)?



Cut-off
-------

In this section, to simplify notation, we will consider monoatomic system.
It can be easily generalized to multiple species.

The Debye formula allows us to directly calculate the intensity
from an isolated particle.
But in "infinite" bulk material (simulated in periodic boundary conditions)
we must somehow limit the number of considered atomic pairs.

The simplest idea could be to pick a cut-off distance and limit the Debye
formula to atomic pairs not further apart than this distance.
But the termination effect would create a large sinusoid.
So it is a common practice to add corrections for this effect.

Further calculations will be easier if we write the Debye formula
as an integral,

.. math:: I(Q) = f^{2} \int_0^\infty n(r) \frac{\sin(Qr)}{Qr} \, \mathrm{d}r ,
    :label: debye-integral

where

.. math:: n(r) \equiv \sum_{i,j} \delta(r-r_{ij}) .

Compensation
^^^^^^^^^^^^

Let us compensate for the missing atomic pairs,
approximating the structure beyond the cut-off distance `r_c` with a continuum:

.. math:: I(Q) \approx I^{r<r_c}(Q) + I_{cont}^{r>r_c}(Q)

The density of the continuum `\rho` is set to the average density of
the structure, and

.. math:: n_{cont}(r) = N \, 4\pi r^{2} \, \rho,

where *N* is the total number of atoms.
We need to calculate the second addend.

.. math::
    I_{cont}^{r>r_c}(Q)
    & = f^{2} \int_{r_c}^\infty n_{cont}(r) \frac{\sin(Qr)}{Qr} \, \mathrm{d}r
    \\
    & = \frac{4\pi N \rho f^2}{Q} \int_{r_c}^\infty r \sin(Qr) \mathrm{d}r
    \\
    & = \frac{4\pi N \rho f^2}{Q}
          \left[ \frac{\sin(Qr) - Qr \cos(Qr)}{Q^2} \right]_{r_c}^\infty

Oops, it doesn't converge.

It will converge if we subtract from *I(Q)* the intensity diffracted from
a continuus system with density `\rho`.

Since continuum does not add to diffraction (at non-zero *Q*),
it should not harm to subtract `I_{cont}` from the right side of the
approximation above.
(I haven't seen it explained like this in the literature, but this
operation is quite common).

So now the correction is `I_{cont}^{r>r_c} - I_{cont} = - I_{cont}^{r<r_c}` ,

.. math::

    I_{cont}^{r<r_c}(Q) = 
    \frac{4\pi N \rho f^2}{Q} \left[ \frac{\sin(Qr) - Qr \cos(Qr)}{Q^2} \right]_0^{r_c}
    = \frac{4\pi N \rho f^2}{Q^3} \left[ \sin(Qr_c) - Qr_c \cos(Qr_c) \right]


Finally,

.. math::

    I(Q)/N = f^2\left[\frac{1}{N}
    \left( \sum_i \sum_{j,r_{ij}<r_c} \frac{\sin(Qr_{ij})}{Qr_{ij}} \right)
    + \frac{4\pi\rho}{Q^{3}}\left(Qr_{c}\cos(Qr_{c})-\sin(Qr_{c})\right)\right].
 
That's the simplest correction.
It can be also applied without the analytical form above.
If the summation is using eq. :eq:`debye-integral`
(with histogram approximation, i.e. *n(r)* is counted in finite intervals),
the alternative is to subtract *n*:sub:`cont`\ *(r)* from *n(r)* in each
interval.

TODO: introduce structure factor S(Q). Would using S(Q) instead of I(Q)
make things simpler?

Damping
^^^^^^^

The correction above works well enough for polycrystalline systems,
but may not work for a single crystal.
It should work fine if the pair correlation function is flat at the cut-off
distance. If it is not flat, it is necessary to smooth it
using damping function.

In a few papers
the `sinc function <http://en.wikipedia.org/wiki/Sinc_function>`_
is used for this purpose:

* E. Lorch in
  `J. Phys. C in 1969 <http://dx.doi.org/10.1088/0022-3719/2/2/305>`_.
  Actually, he was going the opposite way -- from *I(Q)* to *g(r)*,
  but the truncation effect is similar,

* G. Gutiérrez *et al.* in
  `PRB in 2002 <http://link.aps.org/doi/10.1103/PhysRevB.65.104202>`_
  (`copy <http://www.nucleo-milenio.cl/interior/publications/PRB04202.pdf>`__)
  -- the "window function" in eq. (2) there,

* Z. Lin & L. Zhigilei in
  `PRB in 2006 <http://link.aps.org/doi/10.1103/PhysRevB.73.184113>`_
  (`copy <http://www.dtic.mil/dtic/tr/fulltext/u2/a465173.pdf>`__)
  -- the "damping function", eq. (8) there.

The continuous density approximation from the previous section,
in histogram approximation, replaced *n(r)* in eq. :eq:`debye-integral`
with

.. math:: [n(r) - n_{cont}(r)]

With the sinc damping, *n(r)* is replaced by:

.. math:: [n(r) - n_{cont}(r)] \frac{\sin(\pi r / r_c)}{\pi r / r_c}


Computational approach
----------------------

.. _scat:

Scattering factors
^^^^^^^^^^^^^^^^^^

Atomic scattering factor are assigned automatically, by searching atom names in
`built-in tables <https://raw.github.com/wojdyr/fityk/master/wxgui/atomtables.c>`_.
These tables are based on

* International Tables for Crystallography, Volume C 1992,
  Table 6.1.1.4 (pp. 500-502),
  *Coefficients for analytical approximation to the scattering factors*
  (for x-rays)
* `Neutron scattering lengths and cross-sections`__ (for neutrons)

__ http://www.ncnr.nist.gov/resources/n-lengths/list.html

Histogram approximation
^^^^^^^^^^^^^^^^^^^^^^^

Distance-histogram approximation avoids calculation of expensive sine function
for each atomic pair. It splits computations into two steps.

* calculate a histogram of distances for each pair of atomic species
  (the most computationally intensive step),
* the Debye formula is applied treating distances in single histogram bin
  together.

For monoatomic system it can be written as:

.. math::

  I(Q) &= f^2 \sum_i^N \sum_j^N \frac{\sin(Qr_{ij})}{Qr_{ij}} \\
       &= f^2 \left( N + 2\sum_i^N \sum_{j>i}^N
                     \frac{\sin(Qr_{ij})}{Qr_{ij}} \right) \\
       &\approx f^2 \left( N + 2\sum_k^{N_{bins}} n_k
                     \frac{\sin(Qr_k)}{Qr_k} \right)

where *n*:sub:`k` and *r*:sub:`k` are the number of pairs and the distance
corresponding to the *k*-th bin.

The errors from this approximation (when using default histogram bin size)
are negligible.

(TODO: describe not used alternatives, such as spherical harmonics approx.)

Cell lists
^^^^^^^^^^

When cut-off *r*:sub:`c` is small `cell lists`__ are used to speed up
computations (currently it's implemented only for systems in PBC).

__ http://en.wikipedia.org/wiki/Cell_lists


Parallelization
^^^^^^^^^^^^^^^

If you are working with large configurations (millions of atoms),
you may build the parallel version of the program. It uses the MPI library.
Also note that compiler options related to floating point arithmetics
may notably improve performance.

Only calculation of atomic distances is parallelized and optimized.


Usage
=====

``debyer --help`` shows this summary:

.. highlight:: none

::

 Usage: debyer [OPTIONS]... [FILE]...
 
   -h, --help                    Print help and exit
       --full-help               Print help, including hidden options, and exit
   -V, --version                 Print version and exit
   -q, --quiet                   silent mode
   -v, --verbose                 increase verbosity level (can be used twice)
       --bench=ATOM-COUNT        benchmark - reports time of ID calculation for
                                   ATOM-COUNT atoms.
 
 Stage 1: calculation of ID (Interatomic Distances):
   -r, --cutoff=FLOAT            cut-off distance for ID calculation
       --quanta=FLOAT            ID discretization quanta  (default=`0.001')
   -a, --pbc-a=FLOAT             PBC box length in x direction
   -b, --pbc-b=FLOAT             PBC box length in y direction
   -c, --pbc-c=FLOAT             PBC box length in z direction
       --sample=INT              calculate ID by random sampling n atoms
   -d, --save-id[=FILENAME]      save ID to file
 
 Stage 2: calculation of the final result from ID:
 
  Group: mode
   what to calculate - pick one
   -x, --xray                    x-ray powder pattern
   -n, --neutron                 neutron powder pattern
   -S, --sf                      scattering factor (a.k.a total scattering
                                   structure function)
   -R, --RDF                     radial distribution function (RDF)
   -g, --PDF                     pair distribution function (PDF)
   -G, --rPDF                    reduced PDF
 
  Output range (for example -f20 -t100 -s0.1):
   -f, --from=FLOAT              start of calculated pattern
   -t, --to=FLOAT                end of calculated pattern
   -s, --step=FLOAT              step of calculated pattern
 
  Options for real space patterns:
   -w, --weight=STRING           weighting by x-ray f(0) or neutron b  (possible
                                   values="x", "n", "1" default=`1')
   -p, --partials                include partials as extra columns
 
  Options for reciprocal space patterns:
   -l, --lambda=FLOAT            wavelength (omit for a pattern in Q)
       --sinc                    use sinc as damping function (requires cut-off)
 
  Options valid for both real and reciprocal space patterns:
       --ro=FLOAT                numeric density, required for (r)PDF and
                                   diffractogram with cutoff
 
   -o, --output=FILENAME         output filename. If not given, will be
                                   auto-generated.
 
**IDs** (interatomic distances) can be saved with ``-d``.
It is useful only when calculating multiple patterns of the same sample.
The saved ID can be then used as an input file.

**Supported file formats**:
AtomEye extended CFG, DL_POLY CONFIG/REVCON,
LAMMPS input file, XMOL XYZ,
and plain text with *x y z symbol* or *symbol x y z* in each line.

**PBC**:
if the input file contains the box size (e.g., in AtomEye CFG, LAMMPS and
DL_POLY file formats) and the options ``-a``, ``-b``, ``-c`` are not given,
the size from the file is used.
To discard PBC use options ``-a0``, ``-b0``, ``-c0``.
PBC in only one or two dimensions are not supported.

**Sampling** large systems: when the ``--sample=N`` option
is specified, *N* atoms are randomly chosen and only atomic distance
between these *N* atoms and all other atoms are calculated.
The only point of this option is to make computations faster (but less exact).

**Units of length** should be used consistently. The program is unit agnostic.
If the coordinates in the input file are in Angstroms,
the ``--lambda`` option should be also in Å,
the values of `Q` will be in Å\ :sup:`-1`, and so on.

**Other units**: diffraction angle in the output is in degrees of `2\theta`.
(anything else?)

Examples
========

Calculate diffraction pattern calculation::

 debyer -x -f1 -t20 -s0.01 -l0.1 -o t.dat zns.xyz

TBC


.. _companions:

Companions
##########

dbr_extend
==========

Various modifications of atomistic system in orthorhombic PBC.
This program was written with bicrystal geometry in mind.
Options ``--help`` and ``--show-examples`` display basic documentation.

::

 $ dbr_extend --help

 dbr_extend 0.3

 Usage: dbr_extend [OPTIONS]... [FILE]...

   -h, --help                    Print help and exit
   -V, --version                 Print version and exit
       --show-examples           show examples, can be more useful than --help
   -v, --verbose                 increase verbosity level (can be used 3x)
   -s, --min-cell=FLOAT          (internal) size (lower limit) of cells used for
                                   searching atoms [A]  (default=`2.0')
   -r, --reduced                 use reduced coordinates (between 0 and 1)

 Most of the actions in this program requires a defined `slab'.
 The slab is constructed by a bounding plane (which must be normal
 to one of the x, y, z axes) and width.
   -x, --x=FLOAT                 defines a plane x=FLOAT
   -y, --y=FLOAT                 defines a plane y=FLOAT
   -z, --z=FLOAT                 defines a plane z=FLOAT
   -b, --bound=ENUM              for system finite in given direction: set a
                                   plane to the lower (x,y,z) or upper (X,Y,Z)
                                   bound  (possible values="x", "X", "y",
                                   "Y", "z", "Z")
   -w, --width=FLOAT             width of the slab that will be processed (or
                                   |t|)

 If width is not given, it will be equal to the shortest translation vector
 found by the program. The following options affect searching of this vector.
   -e, --epsilon=FLOAT           epsilon used to compare coordinates
                                   (default=`0.1')
   -m, --min-delta=FLOAT         min length of translation vector
   -M, --max-delta=FLOAT         max length of translation vector
   -p, --periodic                require translational symmetry of all the
                                   system

 Actions to be performed on atomic configuration file.
   -d, --delete                  delete atoms in the slab
   -C, --cut                     cut out the slab and change PBC accordingly
   -u, --add-vacuum              add vacuum 'slab'; extends PBC
   -c, --add-copy[=N]            multiplicate the slab N times (extend PBC and
                                   add atoms)  (default=`1')
   -a, --add=WIDTH               extend PBC by WIDTH and fill the new space
                                   using translation symmetry found in the
                                   defined slab
   -N, --multiply=NxNxN          multiply configuration in x y and z
   -S, --shift=x,y,z             shift the system (all atoms) under PBC
   -R, --resize=x,y,z            resize the system, changing interatomic
                                   distances; the argument is either new size (0
                                   means no change) or another filename
   -F, --find-trans              find translation symmetries
   -U, --make-cubic[=a1,a2,a3,b1,b2,b3,c1,c2,c3]
                                 remove small distortions from perfect cubic
                                   lattice
       --merge                   merge atoms that are closer to each other than
                                   epsilon
       --t1                      transformation defined in the code as
                                   tranform1()

   -D, --density                 calculate numeric density of the slab

 Output file. Default is a dry run.
   -i, --in-place                replace input file with output
   -o, --output=FILENAME         output file

 For usage examples, invoke the program with --show-examples option.


 > dbr_extend --show-examples

                               USAGE EXAMPLES

 dbr_extend -z41.5 -e0.2 -vvv file.cfg
   Try to find periodicity of the structure in the z direction starting from
   z0=41.5. It tries to find z1 such that each atom with z0 < z < z1 has
   a periodic image with x'=x, y'=y, z'=z+delta, delta=z1-z0.
   Comparisons of coordinates are done with epsilon=0.2.

 dbr_extend -c10 -z41.5 -e0.2 -i file.cfg
   Extend PBC box by in the z direction by 10 times the value reported
   from the command above, copy atoms to the newly created space,
   write the configuration back to the same file.

 dbr_extend -bz -w3 -d -i file.cfg
 dbr_extend -bZ -w-3 -d -i file.cfg
   (Initially, file.cfg contained a slab with surfaces normal to z axis.)
   Delete surfaces (3A deep) of the slab.

 dbr_extend -v -bz -w2 -a3 -o tmp2.cfg tmp1.cfg
 dbr_extend -v -bZ -w-2 -a3 -o tmp3.cfg tmp2.cfg
   (Initially, tmp1.cfg contained a slab with surfaces normal to z axis.)
   Extend the slab in the z direction, 3A from each surface.

 dbr_extend -S0,0,0.5 -r -i file.cfg
   Shift object under PBC, by half of the PBC box, in the z direction,
   write the configuration back to file.cfg.

 dbr_extend -N1x2x1 -o out.cfg file.cfg
   Duplicate the system in the y direction (create a supercell).

 dbr_extend -v -z-7 -w12. --density file.cfg
   Calculate (in a smart way) numeric density of the slab defined by planes
   z=-7 and z=5 (its in PBC, so it's continuus region that includes z=0).
   Designed to calculate density of a GB in bicrystal geometry.

 dbr_extend --resize=ref.cfg -o output.cfg input.cfg
   Resize the PBC box, make it the same as the size of the file ref.cfg.
   Atomic positions are scaled with the box.


dbr_conv
========

dbr_conv converts between atomistic model file formats. It supports:

* AtomEye extended CFG,
* DL_POLY CONFIG/REVCON,
* LAMMPS data format,
* XMOL XYZ,
* PDB (write only)
* plain format (xyza): each line contains either x y z symbol or symbol x y z.

Compressed files (.gz, .bz2) can be read, but not written.

DL_POLY supports polarizable shell models and the CONFIG/REVCON format contains
positions of both atoms and shells. Debyer relies on the convention that shells
have names ending with one of the following strings: _sh, _shl, _shel, _shell,
-sh, -shl, -shel, -shell. The shells are ignored when reading files.

LAMMPS data file format does not contain contain atom types, only numbers
that are associated to types in an input script to LAMMPS.
Here we use a convention (both when writing and reading the file)
that atom types are given in a comment after the "atom types" line, e.g.
``2 atom types # C Si`` means that type 1 is C and type 2 is Si.

Another program for converting atomistic file formats is
mdfile.py from `gosam`__ (Python).

__ https://github.com/wojdyr/gosam/

::

 $ ./dbr_conv -h
 Usage: dbr_conv [OPTIONS]... INPUT_FILE OUTPUT_FILE
        dbr_conv [OPTIONS...] -t TO_FORMAT -m INPUT_FILE1 [INPUT_FILE2]...
 
   -h          Print help and exit
   -V          Print version and exit
   -q          Silent mode
   -v          Increase verbosity level (can be used twice)
   -s          Sort atoms by atomic symbols.
   -f FORMAT   Convert file from format.
   -t FORMAT   Convert file to format.
   -m          Convert multiple files. Output files have only file extension
               changed.
 
 Supported format names: atomeye, dlpoly, lammps, pdb, xyz, xyza.
 Compressed files (.gz, .bz2) can be read, but not written.

others
======

* dbr_bonds -- show some statistics about bonds in the system.
  The program takes two arguments: input file and maximum bond length, e.g:
  ``dbr_bonds file.cfg 2.1``.

* dbr_sic -- utility to calculate coordination numbers of atoms, so-called
  ring distribution and other features of zinc-blende structure.
  Named 'sic' because it was used to study SiC structure.

* dbr_symmetry -- obsolete

* dbr_volume -- undocumented

* dbr_diffus -- undocumented

gosam
=====

It is actually a separate set of programs, available at
https://github.com/wojdyr/gosam/ .

Gosam is a bunch of Python scripts that can:

* create monocrystal in PBC box (atomistic model),
* create bicrystals with coincidence site lattice (CSL) boundaries in PBC box,
* create crystalline grains of given shape, with vacancies,
  thermal vibrations, etc.
* read/write several file formats (AtomEye cfg, VASP POSCAR, LAMMPS, DL_POLY,
  XMOL XYZ).
