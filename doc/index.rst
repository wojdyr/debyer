
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
   (`Q = \left|\boldsymbol{Q}\right| = 4 \pi \sin \theta / \lambda` where
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

Let's re-derive it once more. We start with formula for the amplitude of the
scattered wave. If we ignore non-elastic scattering and assume that each photon
is scattered only once, the formula is a simple sum over all atoms:

.. math::

    \Psi(\boldsymbol{Q}) = \sum_i \psi_i
    =\sum_i f_i \exp(-i \boldsymbol{Q} \cdot \boldsymbol{r}_i)

(*i* has two meanings here, hopefully it's not confusing).

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
`\boldsymbol{r}_{ij}`,
`(\boldsymbol{r}_{ij} = \boldsymbol{r}_i - \boldsymbol{r}_j)`.

.. math::

    I(\boldsymbol{Q})
    = \sum_i \sum_j f_i f_j \exp \left( -i Q r_{ij} \cos \gamma \right)


Avaraging the exponential part gives

.. math::

    \begin{eqnarray}
    \left\langle \exp( -iQr_{ij} \cos \gamma ) \right\rangle
    & = &
    \frac{1}{4\pi r_{ij}^2} \int_0^\pi \exp\left( -iQr_{ij} \cos\gamma \right)
    \, 2\pi r_{ij}^2 \sin\gamma \, \mathrm{d}\gamma \\
    & = &
    \frac{1}{2} \int_0^\pi \exp\left( -iQr_{ij}\cos\gamma \right)
    \sin\gamma \, \mathrm{d}\gamma \\
    & = &
    \frac{1}{2} \left[ \frac{\exp\left( -iQr_{ij} \cos\gamma \right)}
    {iQr_{ij}} \right]_0^\pi \\
    & = &
    \frac{\exp(iQr_{ij})-\exp(-iQr_{ij})}{2iQr_{ij}} \\
    & = &
    \frac{\sin(Qr_{ij})}{Qr_{ij}}
    \end{eqnarray}

which proves the Debye formula.


Modifications
-------------

Usually this formula is normalized by `\frac{1}{N}` or `\frac{1}{N<f>^2}` .
TBC

include Debye-Waller factor?
the correlated broadening factor for the atom pair (as mentioned in Farrow 2009)



Cut-off
-------

For isolated particles we can directly calculate the intensity.
In case of bulk material (simulated in periodic boundary conditions)
we must use a cut-off. Simply discarding all atomic pairs longer
than a choosen cut-off distance would give a pattern that looks like sinusoid.

To have simpler notation, let's consider monoatomic system.
(*TODO: extend it to multiple species later*).

It is instructive to write the Debye formula as a function of
`n(r) \equiv \sum_{i,j} \delta(r-r_{ij})` ,

.. math:: I(Q) = f^{2} \int_0^\infty n(r) \frac{\sin(Qr)}{Qr} \, \mathrm{d}r.

To use a simple cut-off, we need to estimate and compensate for the atomic
pairs beyond the cut-off. Let us approximate the structure beyond the
cut-off distance with a continuum of the density equal to the average density
of the structure,

.. math:: n(r) = N \, 4\pi r^{2} \, \rho,

where *N* is the total number of atoms and `\rho` is the numeric density.
The contribution of the discared pairs would be

.. math::
    \begin{eqnarray}
    I_{cont}^{r>r_c}(Q)
    & = &
    f^{2} \int_{r_c}^\infty n(r) \frac{\sin(Qr)}{Qr} \, \mathrm{d}r \\
    & = &
    \frac{4\pi N \rho f^2}{Q} \int_{r_c}^\infty r \sin(Qr) \mathrm{d}r \\
    & = &
    \frac{4\pi N \rho f^2}{Q} \left[ \frac{\sin(Qr) - Qr cos(Qr)}{Q^2} \right]_{r_c}^\infty 
    \end{eqnarray}

Oops, it doesn't converge.

It will converge if we subtract from `I(Q)` the intensity diffracted from
a continuus system with density `\rho`. Since continuum does not
add to diffraction (at non-zero `Q`), this operation doesn't change `I(Q)`.
(I haven't seen it explained like this in the literature, but subtracting
the average density is quite common).

So now the correction is `I_{cont}^{r>r_c} - I_{cont} = - I_{cont}^{r<r_c}` ,

.. math::

    I_{cont}^{r<r_c}(Q) = 
    \frac{4\pi N \rho f^2}{Q} \left[ \frac{\sin(Qr) - Qr \cos(Qr)}{Q^2} \right]_0^{r_c}
    = \frac{4\pi N \rho f^2}{Q^3} \left[ \sin(Qr_c) - Qr \cos(Qr_c) \right]


Finally,

.. math::

    I(Q)/N = f^2\left[\frac{1}{N}
    \left( \sum_i \sum_{j,r_{ij}<r_c} \frac{\sin(Qr_{ij})}{Qr_{ij}} \right)
    + \frac{4\pi\rho}{Q^{3}}\left(Qr_{c}\cos(Qr_{c})-\sin(Qr_{c})\right)\right].
 
That's the simplest correction (and the only one implemented in debyer for now).
It works well enough for polycrystalline systems, but not for a single crystal.
Probably a damping function would work better for single crystals, for example
`Lin & Zhigilei in PRB <http://dx.doi.org/10.1103/PhysRevB.73.184113>`_
(eq. (8) there and the surrounding text)

TODO: introduce structure factor S(Q). Would using S(Q) instead of I(Q)
make things simpler?

Computational approaches
------------------------

TODO: histogram approximation, spherical harmonics approx., ...

