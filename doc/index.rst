
Draft of new Debyer's documentation
###################################

I just started writting down some missing pieces of docs here.
It's not organized yet. Some points are just notes to myself.
It's meant for the next version of Debyer, if this ever happens.

**For now, visit**: http://code.google.com/p/debyer

.. toctree::
   :maxdepth: 2

.. note::

    If the math doesn't render properly, right-click on it and change
    ``Math Renderer`` in ``Math Settings``.

Introduction
============

Debyer takes as an input a file with positions of all the atoms in the
virtual sample (up to tens of millions of atoms, perhaps even more) and can
output x-ray and neutron powder diffraction pattern, total scattering structure
function, pair distribution function (PDF) and related functions (RDF, reduced
PDF).

A few use cases for debyer can be found `in papers citing its website`__.
Note that Debyer is not designed to calculate a diffraction pattern
of a perfect crystal.

__ http://scholar.google.com/scholar?q=unipress.waw.pl%2Fdebyer+OR+code.google.com%2Fp%2Fdebyer

Debye scattering formula
========================

X-ray and neutron powder diffraction patterns are calculated using
the Debye scattering equation:

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
the correlated broadening factor for the atom pair (as mentioned in Farrow 2009)



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
It would not work, though.
The termination effect would create a large sinusoid.
So it needs to be more complicated.

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
 
That's the simplest correction (and the only one implemented in debyer for now).

This correction can be also applied without the analytical form above.
If the summation is using eq. :eq:`debye-integral`
(with histogram approximation, i.e. *n(r)* is counted in finite intervals),
just subtract *n*:sub:`cont`\ *(r)* from *n(r)* in each interval.

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

Together with the continuous density approximation from the previous
section, the sinc damping function can be added by replacing *n(r)*
in eq. :eq:`debye-integral` with

.. math:: [n(r) - n_{cont}(r)] \frac{\sin(\pi r / r_c)}{\pi r / r_c}

TBC

Computational approach
----------------------

TODO: histogram approximation, spherical harmonics approx., ...

