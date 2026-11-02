"""\
=====================================================================
Start configurations --- :mod:`polyphys.simulate.build`
=====================================================================

The :mod:`~polyphys.simulate.build` module builds LAMMPS start
configurations of a bead-spring polymer in a predictable, unphysical shape,
a helix or a closed plectoneme, which a LAMMPS run then equilibrates. It
places beads along a curve a fixed bond length apart, adds the chain
topology, checks that LAMMPS can start from the result, and hands the
configuration to :func:`~polyphys.simulate.data.write_data` as a
:class:`~polyphys.simulate.data.DataFile`. The design follows moltemplate's
confined-DNA example (Jewett et al., 2021): a curve, beads placed along it,
the topology, then the file. Nothing is copied from moltemplate.

.. versionadded:: 0.6.0

Contents
========

:func:`helix` and :func:`plectoneme_ring` place beads on the two built-in
shapes; :func:`resample_curve` places them along any curve given as points.
:func:`chain_bonds` and :func:`chain_angles` give the topology of a linear
or ring chain, and :func:`chain_data` wraps positions and topology in a
:class:`~polyphys.simulate.data.DataFile`. :func:`check_conformation` tests
a data file against the conditions under which LAMMPS stops or warns, and
raises :class:`ConformationError` or returns a :class:`ConformationReport`.

Notation
========

:math:`N`
    The number of beads.
:math:`b`
    The bond length, the distance between bonded bead centres.
:math:`\\sigma`, :math:`\\epsilon`
    The bead diameter and the energy unit of the Lennard-Jones (LJ) units of
    the run; all lengths here are in :math:`\\sigma`.
:math:`R`, :math:`p`
    The radius and the pitch (the axial rise per turn) of a helix;
    :math:`c = p / 2\\pi` is the rise per radian.
:math:`\\Delta t`
    The angle between neighbouring beads on a helix.
:math:`\\theta_A`
    The helical angle, between a strand's tangent and the axis:
    :math:`\\tan\\theta_A = 2\\pi R / p`.
:math:`L_A`
    The arclength of one strand of a plectoneme.
:math:`K`, :math:`R_0`
    The spring constant and the maximum extension of a FENE bond.
:math:`d_{\\min}`
    The smallest distance :func:`check_conformation` allows between two
    beads that are not bonded.
:math:`D`, :math:`a`
    The diameter of a cylindrical confinement and the bead diameter;
    :math:`D - a` is the effective diameter, the diameter the bead centres
    can reach.

Physical background
===================

The chains are Kremer-Grest bead-spring chains. Every pair of beads repels
through the Weeks-Chandler-Andersen (WCA) potential, the LJ potential cut
at its minimum :math:`2^{1/6}\\sigma` and shifted up by :math:`\\epsilon`
(Weeks, Chandler, and Andersen, 1971, Eq. 4). Bonded beads add the finitely
extensible nonlinear elastic (FENE) spring, whose force
:math:`K r / (1 - (r/R_0)^2)` diverges at :math:`R_0` (Warner, 1972, Eq.
4). Grest and Kremer (1986, Eqs. 8 and 9) set :math:`K = 30\\,
\\epsilon/\\sigma^2` and :math:`R_0 = 1.5\\,\\sigma` so that chains cannot
cross; with these values at :math:`k_B T = \\epsilon` the bond never
exceeds :math:`1.2\\,\\sigma` (Kremer and Grest, 1990, Appendix). The LAMMPS
``bond_style fene`` holds both terms of the bond,

.. math::

    E = -\\tfrac{1}{2} K R_0^2 \\ln\\left[1 - \\left(\\frac{r}{R_0}\\right)^2
    \\right] + 4\\epsilon\\left[\\left(\\frac{\\sigma}{r}\\right)^{12} -
    \\left(\\frac{\\sigma}{r}\\right)^{6}\\right] + \\epsilon,

the second part cut at :math:`2^{1/6}\\sigma`, and ``special_bonds fene``
switches the pair interaction off between bonded beads only, so beads two
apart along the chain repel like any other pair.

*The start bond length.* Kremer and Grest (1990, Sec. II.B) built their
start configurations as random walks with a bond length of :math:`0.97
\\sigma`, kept beads two apart more than :math:`1.02\\,\\sigma` from each
other, and measured a mean bond length of :math:`0.97\\,\\sigma` in the
equilibrated melt. The force of the FENE and WCA terms together is zero at
:math:`0.961\\,\\sigma`, so a bond of :math:`0.97\\,\\sigma` stores
:math:`0.04\\,\\epsilon`; the same :math:`0.97` is the default of Pizza.py's
``chain.py``, of LAMMPS's ``tools/def.chain``, and of
``BondFENE::equilibrium_distance``. The builders here use it as their
default, `bond_length` ``= 0.97``.

*What stops a LAMMPS run.* LAMMPS warns ``FENE bond too long`` and clamps
the bond force when :math:`1 - (r/R_0)^2 < 0.1`, that is for :math:`r >
\\sqrt{0.9}\\,R_0 = 1.423\\,\\sigma` at :math:`R_0 = 1.5\\,\\sigma`, and
stops with ``Bad FENE bond`` at :math:`r \\ge 2R_0` (``bond_fene.cpp``). A
bond longer than the ghost-atom cutoff, the pair cutoff plus the neighbour
skin (:math:`1.42\\,\\sigma` for a WCA cutoff and a skin of :math:`0.3`),
stops a run with ``Bond atoms missing`` or a wrong periodic image as soon
as it crosses a periodic face or a processor boundary. ``read_data`` loses
an atom outside the box on a non-periodic axis, and ``fix wall/region``
stops when a bead sits on or outside its surface. Two beads at one position
give a non-finite force under every integrator. A bond of :math:`1.425\\,
\\sigma`, that is :math:`0.95\\,R_0`, is above both the warning length and
the default ghost cutoff; in test runs of this model (:math:`K = 30`,
:math:`R_0 = 1.5`, skin 0.3, :math:`k_B T = \\epsilon`) starts whose bonds
were all between :math:`0.97` and :math:`1.10\\,\\sigma` ran under plain
``nve``, under Langevin dynamics, and under ``nve/limit`` alike, while a
start with every bond at :math:`1.40\\,\\sigma` stopped a plain ``nve`` run
within 70 steps. :func:`check_conformation` tests each of these conditions.

Shapes
======

*Helix.* Bead :math:`k` sits at :math:`(R\\cos k\\Delta t, R\\sin k\\Delta
t, c k\\Delta t)`, centred on :math:`z = 0`. The distance between
neighbouring beads is the chord

.. math::

    b^2 = 4R^2 \\sin^2(\\Delta t / 2) + c^2 \\Delta t^2,

whose right side increases on :math:`(0, \\pi]`, so :math:`\\Delta t` is
its one root there whenever :math:`b^2 \\le 4R^2 + c^2\\pi^2`, that is
whenever a bond spans at most half a turn. Stepping by arclength instead
would give bonds shorter than the arc, by an amount that depends on the
curvature. Adjacent turns are :math:`p\\sin\\theta_A` apart, so a helix
needs :math:`p \\ge d_{\\min} / \\sin\\theta_A` to pass the pair check.

*Plectoneme ring.* The closed curve is the plasmid of Neukirch and
Starostin (2008, Appendix A, Eqs. A1 to A4): a right-handed helix of radius
:math:`R` and helical angle :math:`\\theta_A` (strand A), an end loop
(B), the same helix rotated by :math:`\\pi` about the axis and run back
down (strand C), and a second end loop (D). The loops pass through the axis
and keep the curve inside the cylinder of radius :math:`R`; the tangent is
continuous everywhere and the curvature jumps at the four joints. It is the
interwound superhelix of Boles, White, and Cozzarelli (1990, Fig. 3), with
their semicircular end caps replaced by the smooth loops. Beads are placed
along the curve by chord stepping, each at the first point of the curve a
distance :math:`b` from the one before, and the strand arclength
:math:`L_A` is solved so that the :math:`N` chords close the ring exactly,
the bond between the last bead and the first included. :math:`R`,
:math:`p`, :math:`b`, and :math:`N` are therefore exactly as given, and only
the height of the ring is a result. The two strands are half a pitch apart
along the axis, so a ring needs :math:`p \\ge 2 d_{\\min} / \\sin\\theta_A`
(about :math:`2.0` to :math:`2.1\\,\\sigma` for :math:`R \\ge 2\\sigma` and
:math:`d_{\\min} = \\sigma`) to pass the pair check.

*Any curve.* :func:`resample_curve` applies the same chord stepping to a
polyline. A closed polyline has no height to solve for, so it is scaled
uniformly about its centroid until :math:`N` chords close it; the shape is
kept and the size changes by about :math:`(N b - \\ell) / \\ell` for a
polyline of length :math:`\\ell`. Closing needs a curve whose features are
large compared with :math:`b`; a curve with smaller corners or teeth is
refused rather than closed with one odd bond.

Equilibration
=============

A start configuration only has to be one LAMMPS runs from; the run removes
its memory. Kremer and Grest (1990, Sec. II.B) relaxed their starts with a
soft pair potential, :math:`A[1 + \\cos(\\pi r / 2^{1/6}\\sigma)]`, raising
:math:`A` from 1 to about 60, and moltemplate's confined-DNA example uses
``fix nve/limit`` with a Langevin thermostat and recommends soft potentials
for a long confined chain; LAMMPS provides both (``pair_style soft`` with
``fix adapt``, and ``fix nve/limit``). Neither is part of this module.

Complexity
==========

:func:`helix` takes :math:`O(N)` time and memory. :func:`plectoneme_ring`
and a closed :func:`resample_curve` take :math:`O(N)` time and memory per
march along the curve, and about ten marches to close the ring.
:func:`check_conformation` takes :math:`O(N^2)` time, in row chunks, so
its memory is :math:`O(N)` beyond the positions.

References
==========

- Kremer, K. and Grest, G. S., "Dynamics of entangled linear polymer melts:
  A molecular-dynamics simulation", J. Chem. Phys. 92, 5057 (1990),
  `doi:10.1063/1.458541 <https://doi.org/10.1063/1.458541>`_: Sec. II.B,
  Eqs. (2.3) and (2.4), p. 5060, and the Appendix, Eqs. (A1) and (A2).
- Grest, G. S. and Kremer, K., "Molecular dynamics simulation for polymers
  in the presence of a heat bath", Phys. Rev. A 33, 3628 (1986),
  `doi:10.1103/PhysRevA.33.3628
  <https://doi.org/10.1103/PhysRevA.33.3628>`_: Eqs. (8) and (9).
- Weeks, J. D., Chandler, D., and Andersen, H. C., "Role of repulsive
  forces in determining the equilibrium structure of simple liquids",
  J. Chem. Phys. 54, 5237 (1971), `doi:10.1063/1.1674820
  <https://doi.org/10.1063/1.1674820>`_: Eq. (4).
- Warner, H. R., "Kinetic theory and rheology of dilute suspensions of
  finitely extendible dumbbells", Ind. Eng. Chem. Fundam. 11, 379 (1972),
  `doi:10.1021/i160043a017 <https://doi.org/10.1021/i160043a017>`_: Eq. (4).
- Neukirch, S. and Starostin, E. L., "Writhe formulas and antipodal points
  in plectonemic DNA configurations", Phys. Rev. E 78, 041912 (2008),
  `doi:10.1103/PhysRevE.78.041912
  <https://doi.org/10.1103/PhysRevE.78.041912>`_: Appendix A.
- Boles, T. C., White, J. H., and Cozzarelli, N. R., "Structure of
  plectonemically supercoiled DNA", J. Mol. Biol. 213, 931 (1990),
  `doi:10.1016/s0022-2836(05)80272-4
  <https://doi.org/10.1016/s0022-2836(05)80272-4>`_: Fig. 3 and Sec. 2(f).
- Jewett, A. I. et al., "Moltemplate: A tool for coarse-grained modeling of
  complex biological matter and soft condensed matter physics", J. Mol.
  Biol. 433, 166841 (2021), `doi:10.1016/j.jmb.2021.166841
  <https://doi.org/10.1016/j.jmb.2021.166841>`_, and its example
  `confined_viral_DNA
  <https://github.com/jewettaij/moltemplate/tree/031ecd4efc7abdad8d49238404de2d3ef388b37b/examples/coarse_grained/DNA_models/dsDNA_only/2strands/3bp_2particles/confined_viral_DNA>`_.
- The LAMMPS documentation of `bond_style fene
  <https://docs.lammps.org/bond_fene.html>`_, `special_bonds
  <https://docs.lammps.org/special_bonds.html>`_, `read_data
  <https://docs.lammps.org/read_data.html>`_, `comm_modify
  <https://docs.lammps.org/comm_modify.html>`_, `fix wall/region
  <https://docs.lammps.org/fix_wall_region.html>`_, and `fix nve/limit
  <https://docs.lammps.org/fix_nve_limit.html>`_, and the LAMMPS source at
  release ``stable_30Sep2026_update1``: `bond_fene.cpp, lines 89 to 94
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/MOLECULE/bond_fene.cpp#L89-L94>`_
  (the warning and the stop) and `bond_fene.cpp, lines 192 to 195
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/MOLECULE/bond_fene.cpp#L192-L195>`_
  (the equilibrium distance :math:`0.97\\,\\sigma`).
- Plimpton, S., `Pizza.py chain.py
  <https://github.com/thatgeeman/pizza/blob/master/src/chain.py>`_ (2005)
  and LAMMPS `tools/def.chain
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/tools/def.chain>`_:
  bond 0.97 and 1-3 distance 1.02.

Dependencies
============
- Python standard library: `dataclasses`, `logging`, `math`
- Third-party: `numpy`, `scipy` (:func:`scipy.optimize.brentq`)
- PolyPhys: :mod:`polyphys.simulate.core`, :mod:`polyphys.simulate.data`
"""

import logging
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq

from polyphys.simulate.core import Box, FloatArray, IntArray
from polyphys.simulate.data import DataFile

__all__ = [
    "FENE_WARNING_FRACTION",
    "ConformationError",
    "ConformationReport",
    "chain_angles",
    "chain_bonds",
    "chain_data",
    "check_conformation",
    "helix",
    "plectoneme_ring",
    "resample_curve",
]

logger = logging.getLogger(__name__)

#: LAMMPS warns ``FENE bond too long`` when :math:`1 - (r/R_0)^2 < 0.1`,
#: that is for a bond longer than :math:`\sqrt{0.9}\,R_0`
#: (``bond_fene.cpp``, lines 89 to 94).
FENE_WARNING_FRACTION = math.sqrt(0.9)

_Curve = Callable[[FloatArray], FloatArray]
_SAMPLES_PER_BOND = 64
_COINCIDENT = 1e-9
_TINY = 1e-9


class ConformationError(ValueError):
    """A start configuration would make LAMMPS stop or warn.

    Raised by :func:`check_conformation`; the message names every failed
    condition with its worst value and the atoms involved.

    .. versionadded:: 0.6.0

    Examples
    --------
    >>> box = Box(lo=np.full(3, -5.0), hi=np.full(3, 5.0))
    >>> pair = chain_data([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], box)
    >>> try:
    ...     check_conformation(pair, r0=1.5, ghost_cutoff=1.42,
    ...                        periodic=(True, True, True))
    ... except ConformationError as err:
    ...     print(err)
    1 bond at or above sqrt(0.9) R0 = 1.423: longest 1.5 (atoms 1 2); \
1 bond at or above the ghost cutoff 1.42: longest 1.5 (atoms 1 2)
    """


@dataclass(frozen=True)
class ConformationReport:
    """The extremes :func:`check_conformation` measured on a passing file.

    .. versionadded:: 0.6.0

    Attributes
    ----------
    n_atoms : int
        The number of atoms.
    min_bond, max_bond : float
        The shortest and the longest bond, in the length unit of the file;
        NaN when there are no bonds.
    min_pair : float
        The shortest distance between two atoms that are not bonded, with
        the minimum image on the periodic axes; NaN when there is no such
        pair.
    max_radial : float or None
        The largest distance of an atom from the cylinder axis, or None
        when no cylinder was given.
    """

    n_atoms: int
    min_bond: float
    max_bond: float
    min_pair: float
    max_radial: float | None


def helix(
    n_beads: int, *, radius: float, pitch: float, bond_length: float = 0.97
) -> FloatArray:
    r"""Place beads a fixed distance apart along a helix.

    Bead :math:`k` sits at :math:`(R\cos k\Delta t, R\sin k\Delta t, c k
    \Delta t)` on a right-handed helix of radius :math:`R` and pitch
    :math:`p` about the :math:`z` axis, with :math:`c = p/2\pi`, centred on
    :math:`z = 0`. The angle step :math:`\Delta t` is the root of the chord
    equation :math:`b^2 = 4R^2\sin^2(\Delta t/2) + c^2\Delta t^2`, so every
    bond has the length :math:`b` exactly.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    n_beads : int
        The number of beads :math:`N`, 1 or more.
    radius : float
        The helix radius :math:`R`, positive, in the length unit of the run.
    pitch : float
        The pitch :math:`p`, the rise per turn, positive. Adjacent turns are
        :math:`p\sin\theta_A` apart, with :math:`\tan\theta_A = 2\pi R/p`.
    bond_length : float, optional
        The bond length :math:`b`, positive; by default :math:`0.97`, the
        start bond length of Kremer and Grest (1990) in units of the bead
        diameter.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_beads, 3)``, float64: the bead positions, the first at
        angle zero.

    Raises
    ------
    ValueError
        If `n_beads` is below 1, a length is not positive, or the bond does
        not fit in half a turn, :math:`b^2 > 4R^2 + c^2\pi^2`.

    Notes
    -----
    The chord equation's right side increases on :math:`(0, \pi]`, so
    :func:`scipy.optimize.brentq` finds its one root there. Time and memory
    are :math:`O(N)`.

    Examples
    --------
    A bond of :math:`\sqrt 2` on a helix of radius 1 and pitch 6 steps by
    :math:`\pi/3` and rises by 1 per bead:

    >>> beads = helix(4, radius=1.0, pitch=6.0, bond_length=2 ** 0.5)
    >>> np.round(beads, 6) + 0.0
    array([[ 1.      ,  0.      , -1.5     ],
           [ 0.5     ,  0.866025, -0.5     ],
           [-0.5     ,  0.866025,  0.5     ],
           [-1.      ,  0.      ,  1.5     ]])
    """
    _check_count("n_beads", n_beads, 1)
    for name, value in (("radius", radius), ("pitch", pitch),
                        ("bond_length", bond_length)):
        _check_positive(name, value)
    rise = pitch / (2.0 * math.pi)
    step = _chord_angle(bond_length, radius, rise)
    angles = np.arange(n_beads, dtype=np.float64) * step
    z = rise * angles
    z -= 0.5 * (z[0] + z[-1])
    positions: FloatArray = np.column_stack(
        [radius * np.cos(angles), radius * np.sin(angles), z])
    return positions


def plectoneme_ring(
    n_beads: int, *, radius: float, pitch: float, bond_length: float = 0.97
) -> FloatArray:
    r"""Place beads a fixed distance apart along a closed plectoneme.

    The curve is the plasmid of Neukirch and Starostin (2008, Appendix A):
    two coaxial right-handed helices of radius :math:`R` and pitch
    :math:`p`, half a turn apart, one run up and one run down, joined at
    each end by a smooth loop through the axis. The beads are placed by
    chord stepping, each at the first point of the curve a distance
    :math:`b` from the one before, and the strand arclength :math:`L_A` is
    solved so that :math:`N` chords close the curve: every bond, the one
    between the last bead and the first included, has the length :math:`b`,
    and :math:`R` and :math:`p` are exactly as given. The ring is centred on
    the origin with its axis along :math:`z`; no bead lies farther than
    :math:`R` from the axis.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    n_beads : int
        The number of beads :math:`N`. The two end loops alone take about
        :math:`7.7 R / b` beads each; fewer raises.
    radius : float
        The helix radius :math:`R`, positive, in the length unit of the run.
    pitch : float
        The pitch :math:`p` of each strand, positive. The two strands are
        :math:`(p/2)\sin\theta_A` apart and adjacent turns of one strand
        :math:`p\sin\theta_A`, with :math:`\tan\theta_A = 2\pi R/p`; a
        pitch of about :math:`2\sigma` or more keeps the strands a bead
        diameter apart.
    bond_length : float, optional
        The bond length :math:`b`, positive; by default :math:`0.97`, as in
        :func:`helix`.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_beads, 3)``, float64: the bead positions, from the foot
        of the first strand up and around.

    Raises
    ------
    ValueError
        If `n_beads` is too small for the loops, a length is not positive,
        or the bond does not fit in half a turn of the helix.

    Notes
    -----
    The closing gap, the distance from the last bead back to the first
    minus :math:`b`, grows with :math:`L_A` and is continuous, so
    :func:`scipy.optimize.brentq` finds its root in a bracket around the
    estimate :math:`(N s - \ell_{loops}) / 2`, with :math:`s` the arc a
    chord of :math:`b` spans on the helix and :math:`\ell_{loops}` the
    length of the two loops; the root lies within a fraction of a bond of
    the estimate. Each evaluation is one march along the curve, :math:`O(N)`
    in time and memory, and the solve takes about ten marches.

    Examples
    --------
    >>> ring = plectoneme_ring(60, radius=2.0, pitch=2.5)
    >>> ring.shape
    (60, 3)
    >>> bonds = np.linalg.norm(np.roll(ring, -1, axis=0) - ring, axis=1)
    >>> bool(np.allclose(bonds, 0.97, rtol=1e-12, atol=0.0))
    True
    >>> round(float(np.hypot(ring[:, 0], ring[:, 1]).max()), 6)
    2.0
    """
    _check_count("n_beads", n_beads, 3)
    for name, value in (("radius", radius), ("pitch", pitch),
                        ("bond_length", bond_length)):
        _check_positive(name, value)
    theta = math.atan2(2.0 * math.pi * radius, pitch)
    bond = bond_length / radius                     # in units of R
    cot = 1.0 / math.tan(theta)
    arc_per_chord = _chord_angle(bond, 1.0, cot) / math.sin(theta)
    sample_step = min(bond, 1.0) / _SAMPLES_PER_BOND
    loops = _plectoneme_loop_length(theta, sample_step)
    estimate = (n_beads * arc_per_chord - loops) / 2.0

    def closing_gap(strand: float) -> float:
        beads, gap = _march_curve(
            lambda u: _plectoneme_curve(u, strand, theta),
            2.0 * strand + 2.0 * math.pi, n_beads, bond, sample_step)
        return -bond if gap is None else gap - bond

    too_few = ValueError(
        f"'n_beads' must be at least about "
        f"{math.ceil((loops + 2 * arc_per_chord) / arc_per_chord)} for a "
        f"plectoneme ring of radius {radius} and pitch {pitch} with bonds "
        f"of {bond_length}, got {n_beads}")
    if estimate <= 0.0:
        raise too_few
    width = arc_per_chord
    lo, hi = max(estimate - width, _TINY), estimate + width
    while closing_gap(lo) >= 0.0:
        if lo <= _TINY:
            raise too_few
        width *= 2.0
        lo = max(lo - width, _TINY)
    while closing_gap(hi) <= 0.0:
        width *= 2.0
        hi += width
    strand = brentq(closing_gap, lo, hi, xtol=1e-14, maxiter=200)
    beads, gap = _march_curve(
        lambda u: _plectoneme_curve(u, strand, theta),
        2.0 * strand + 2.0 * math.pi, n_beads, bond, sample_step)
    _check_closed(gap, bond)
    beads *= radius
    beads[:, 2] -= 0.5 * (beads[:, 2].min() + beads[:, 2].max())
    logger.debug("plectoneme_ring: strand arclength %.6f, height %.6f",
                 strand * radius, strand * radius * math.cos(theta))
    return beads


def resample_curve(
    points: npt.ArrayLike, n_beads: int, bond_length: float = 0.97, *,
    closed: bool = False,
) -> FloatArray:
    r"""Place beads a fixed distance apart along a curve given as points.

    The curve is the polyline through `points`. The first bead sits on the
    first point, and each next bead at the first point of the polyline a
    distance :math:`b` from the bead before (the first crossing of the
    sphere of radius :math:`b` about it), so every bond has the length
    :math:`b` exactly. On an open curve the beads cover the first
    :math:`N - 1` chords and the rest of the curve is unused. A closed
    curve is scaled uniformly about the centroid of its points until
    :math:`N` chords close it, the bond from the last bead to the first
    included; the shape is kept and the size changes by about :math:`(Nb -
    \ell)/\ell` for a polyline of length :math:`\ell`.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    points : array_like
        Shape ``(n_points, 3)``: the vertices of the polyline, 2 or more,
        in the length unit of the run. Repeated consecutive points are
        dropped. A closed curve may or may not repeat its first point last.
    n_beads : int
        The number of beads :math:`N`: 1 or more on an open curve, 3 or more
        on a closed one.
    bond_length : float, optional
        The bond length :math:`b`, positive; by default :math:`0.97`, as in
        :func:`helix`.
    closed : bool, optional
        Whether the curve is closed and the ring must close with the bond
        length too; False by default.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_beads, 3)``, float64: the bead positions, in the order of
        the curve.

    Raises
    ------
    ValueError
        If `points` is not ``(n_points, 3)`` with 2 or more distinct
        consecutive points, a value is not finite, `n_beads` is too small,
        `bond_length` is not positive, an open curve is too short for
        `n_beads` beads, or a closed curve cannot be closed with equal
        bonds (see Notes).

    Notes
    -----
    Each chord is a line-sphere intersection, exact to round-off. A march
    takes :math:`O(N + n_{points})` time and memory; closing a curve takes
    about ten marches, the scale found by :func:`scipy.optimize.brentq`.
    The polyline, not the smooth curve it may stand for, carries the beads:
    sample a smooth curve finely compared with :math:`b`. Closing needs the
    beads to move continuously as the curve is scaled, which holds when the
    curve's features, its radius of curvature and the spacing of its
    corners, are large compared with :math:`b`; on a curve with corners or
    teeth smaller than a bond the closing gap jumps, no scale closes the
    curve exactly, and ``ValueError`` is raised rather than a ring with one
    odd bond.

    Examples
    --------
    A corner: the second bead is where the sphere of radius :math:`\sqrt 2`
    about the origin first meets the polyline, past the corner.

    >>> corner = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 5.0, 0.0]]
    >>> np.round(resample_curve(corner, 3, 2 ** 0.5), 6)
    array([[0.      , 0.      , 0.      ],
           [1.      , 1.      , 0.      ],
           [1.      , 2.414214, 0.      ]])

    A closed square of side 4 holds 16 bonds of length 1 as it is:

    >>> square = [[0, 0, 0], [4, 0, 0], [4, 4, 0], [0, 4, 0]]
    >>> beads = resample_curve(square, 16, 1.0, closed=True)
    >>> np.round(beads[3:6], 6)
    array([[3., 0., 0.],
           [4., 0., 0.],
           [4., 1., 0.]])
    """
    vertices = np.array(points, dtype=np.float64)
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError(
            f"'points' must have shape (n_points, 3), got {vertices.shape}")
    if not np.all(np.isfinite(vertices)):
        raise ValueError("'points' must be finite")
    _check_count("n_beads", n_beads, 3 if closed else 1)
    _check_positive("bond_length", bond_length)
    keep = np.ones(vertices.shape[0], dtype=bool)
    keep[1:] = np.any(np.diff(vertices, axis=0) != 0.0, axis=1)
    vertices = vertices[keep]
    if closed and vertices.shape[0] > 1 and np.array_equal(
        vertices[0], vertices[-1]
    ):
        vertices = vertices[:-1]
    if vertices.shape[0] < 2:
        raise ValueError("'points' must hold 2 or more distinct points")
    if not closed:
        beads, _ = _march_polyline(vertices, n_beads, bond_length)
        if beads.shape[0] < n_beads:
            raise ValueError(
                f"the curve holds {beads.shape[0]} beads {bond_length} "
                f"apart, not {n_beads}")
        return beads
    centroid = vertices.mean(axis=0)
    polyline = np.vstack([vertices, vertices[:1]])
    length = float(np.linalg.norm(np.diff(polyline, axis=0), axis=1).sum())

    def closing_gap(scale: float) -> float:
        beads, gap = _march_polyline(
            centroid + scale * (polyline - centroid), n_beads, bond_length)
        return -bond_length if gap is None else gap - bond_length

    estimate = n_beads * bond_length / length
    width = 2.0 * bond_length / length
    lo, hi = max(estimate - width, _TINY), estimate + width
    while closing_gap(lo) >= 0.0:
        width *= 2.0
        lo = max(lo - width, _TINY)
    while closing_gap(hi) <= 0.0:
        width *= 2.0
        hi += width
    scale = brentq(closing_gap, lo, hi, xtol=1e-15, maxiter=200)
    beads, gap = _march_polyline(
        centroid + scale * (polyline - centroid), n_beads, bond_length)
    _check_closed(gap, bond_length)
    logger.debug("resample_curve: scaled the closed curve by %.9f", scale)
    return beads


def chain_bonds(n_beads: int, *, ring: bool = False,
                bond_type: int = 1) -> IntArray:
    """Return the bonds of a linear or ring chain of beads numbered from 1.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    n_beads : int
        The number of beads, 2 or more for a linear chain and 3 or more for
        a ring.
    ring : bool, optional
        Whether to add the bond from the last bead to the first; False by
        default.
    bond_type : int, optional
        The bond type of every bond, 1 or more; 1 by default.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_bonds, 3)``, int64: ``type, atom1, atom2`` per bond, the
        layout of :attr:`~polyphys.simulate.data.DataFile.bonds`, with
        ``n_bonds`` equal to `n_beads` minus one, or to `n_beads` for a
        ring. Time and memory are :math:`O(N)`.

    Raises
    ------
    ValueError
        If `n_beads` or `bond_type` is too small.

    Examples
    --------
    >>> chain_bonds(4, ring=True)
    array([[1, 1, 2],
           [1, 2, 3],
           [1, 3, 4],
           [1, 4, 1]])
    """
    _check_count("n_beads", n_beads, 3 if ring else 2)
    _check_count("bond_type", bond_type, 1)
    first = np.arange(1, n_beads + 1 if ring else n_beads, dtype=np.int64)
    second = first % n_beads + 1
    bonds: IntArray = np.column_stack(
        [np.full_like(first, bond_type), first, second])
    return bonds


def chain_angles(n_beads: int, *, ring: bool = False,
                 angle_type: int = 1) -> IntArray:
    """Return the angles of a linear or ring chain of beads numbered from 1.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    n_beads : int
        The number of beads, 3 or more.
    ring : bool, optional
        Whether to add the two angles that span the closing bond; False by
        default.
    angle_type : int, optional
        The angle type of every angle, 1 or more; 1 by default.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_angles, 4)``, int64: ``type, atom1, atom2, atom3`` per
        angle with ``atom2`` the vertex, the layout of
        :attr:`~polyphys.simulate.data.DataFile.angles`, with ``n_angles``
        equal to `n_beads` minus two, or to `n_beads` for a ring. Time and
        memory are :math:`O(N)`.

    Raises
    ------
    ValueError
        If `n_beads` or `angle_type` is too small.

    Examples
    --------
    >>> chain_angles(4, ring=True)
    array([[1, 1, 2, 3],
           [1, 2, 3, 4],
           [1, 3, 4, 1],
           [1, 4, 1, 2]])
    """
    _check_count("n_beads", n_beads, 3)
    _check_count("angle_type", angle_type, 1)
    first = np.arange(1, n_beads + 1 if ring else n_beads - 1,
                      dtype=np.int64)
    angles: IntArray = np.column_stack(
        [np.full_like(first, angle_type), first, first % n_beads + 1,
         (first + 1) % n_beads + 1])
    return angles


def chain_data(
    positions: npt.ArrayLike, box: Box, *, ring: bool = False,
    angles: bool = False, types: npt.ArrayLike | None = None,
) -> DataFile:
    """Return a data file of one chain from bead positions.

    The atoms get the IDs 1 to :math:`N` in the order of `positions`,
    molecule 1, and the bonds of :func:`chain_bonds` and, if asked, the
    angles of :func:`chain_angles`, all of type 1. The masses and the
    velocities are left to the input script.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    positions : array_like
        Shape ``(n_beads, 3)``: the bead positions along the chain, as
        :func:`helix`, :func:`plectoneme_ring`, or :func:`resample_curve`
        return them.
    box : Box
        The simulation box; it is not checked against the positions here,
        :func:`check_conformation` does that.
    ring : bool, optional
        Whether the chain is a ring; False by default.
    angles : bool, optional
        Whether to add the angles and use the atom style ``angle`` instead
        of ``bond``; False by default.
    types : array_like, optional
        Shape ``(n_beads,)``: the atom types, from 1; all 1 by default.

    Returns
    -------
    DataFile
        With ``n_atom_types`` the largest type, ``n_bond_types`` 1, and
        ``n_angle_types`` 1 if `angles` else 0. Time and memory are
        :math:`O(N)`.

    Raises
    ------
    ValueError
        If `positions` is not ``(n_beads, 3)`` with enough beads for the
        topology, or `types` has the wrong shape or a type below 1.
    TypeError
        If `box` is not a :class:`~polyphys.simulate.core.Box`.

    Examples
    --------
    >>> box = Box(lo=np.full(3, -5.0), hi=np.full(3, 5.0))
    >>> data = chain_data(helix(6, radius=1.0, pitch=2.0), box, angles=True)
    >>> data.atom_style, data.n_atoms, data.n_bonds, data.n_angles
    ('angle', 6, 5, 4)
    >>> data.masses is None, data.mol_ids.tolist()
    (True, [1, 1, 1, 1, 1, 1])
    """
    coordinates = np.array(positions, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            f"'positions' must have shape (n_beads, 3), got "
            f"{coordinates.shape}")
    n_beads = coordinates.shape[0]
    if types is None:
        atom_types = np.ones(n_beads, dtype=np.int64)
    else:
        atom_types = np.array(types)
        if atom_types.shape != (n_beads,):
            raise ValueError(
                f"'types' must have shape ({n_beads},), got "
                f"{atom_types.shape}")
    return DataFile(
        atom_style="angle" if angles else "bond",
        box=box,
        ids=np.arange(1, n_beads + 1),
        types=atom_types,
        positions=coordinates,
        n_atom_types=int(atom_types.max()) if n_beads else 1,
        mol_ids=np.ones(n_beads, dtype=np.int64),
        bonds=chain_bonds(n_beads, ring=ring),
        angles=chain_angles(n_beads, ring=ring) if angles else np.empty(
            (0, 4), dtype=np.int64),
        n_bond_types=1,
        n_angle_types=1 if angles else 0,
    )


def check_conformation(
    data: DataFile, *, r0: float, ghost_cutoff: float,
    periodic: Sequence[bool], min_distance: float = 1.0,
    dcyl: float | None = None, dmon: float | None = None,
) -> ConformationReport:
    r"""Check that LAMMPS can start a Kremer-Grest run from a data file.

    Each condition is one that makes LAMMPS warn or stop (see the module
    documentation). The file fails if any holds:

    1. A bond is at or above the FENE warning length :math:`\sqrt{0.9}\,
       R_0`, where LAMMPS warns ``FENE bond too long`` and clamps the force
       (and stops at :math:`2R_0`).
    2. A bond is at or above `ghost_cutoff`, the pair cutoff plus the
       neighbour skin, beyond which a bond across a periodic face or a
       processor boundary stops the run.
    3. Two atoms, bonded or not, are closer than :math:`10^{-9}`, which
       gives a non-finite force.
    4. Two atoms that are not bonded are closer than `min_distance`;
       beads two apart along the chain count, as ``special_bonds fene``
       treats them.
    5. An atom lies outside the box, or on its upper face, along a
       non-periodic axis, where ``read_data`` loses it. Along a periodic
       axis ``read_data`` remaps it, and the distances above use the
       minimum image.
    6. With a cylinder, an atom lies farther than :math:`(D - a)/2` from
       the :math:`z` axis through the box centre, the effective radius the
       bead centres can reach; a bead there is one bead diameter from a
       wall of diameter :math:`D` built from beads of diameter :math:`a`.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    data : DataFile
        The start configuration; its positions as written, in the length
        unit of the run.
    r0 : float
        The FENE maximum extension :math:`R_0` of the run, positive; 1.5 in
        the Kremer-Grest model.
    ghost_cutoff : float
        The ghost-atom cutoff of the run, positive: the largest pair cutoff
        plus the neighbour skin, unless ``comm_modify cutoff`` sets more;
        :math:`2^{1/6} + 0.3 = 1.42` for a WCA run with the LJ default skin.
    periodic : sequence of bool
        Which of the x, y, and z axes are periodic in the run; a data file
        does not record it.
    min_distance : float, optional
        :math:`d_{\min}`, the smallest allowed distance between two atoms
        that are not bonded, 0 or more; by default 1.0, one bead diameter,
        where the WCA energy of the pair is :math:`\epsilon`.
    dcyl, dmon : float, optional
        The cylinder diameter :math:`D` and the bead diameter :math:`a`,
        given together or not at all, with :math:`D > a`.

    Returns
    -------
    ConformationReport
        The extremes measured, when every condition passes.

    Raises
    ------
    ConformationError
        If a condition fails; the message lists every failed condition
        with its count, its worst value, and the atoms involved.
    ValueError
        If an argument is out of range, `periodic` does not have three
        entries, only one of `dcyl` and `dmon` is given, or the box is
        triclinic and an axis is periodic.

    Notes
    -----
    Bonds are checked for length only; the builders never make a short
    bond, and a bond below about :math:`0.9\,\sigma` stores several
    :math:`\epsilon` in the WCA term of the FENE bond that this check does
    not see. Time is :math:`O(N^2)`, in row chunks, so memory is
    :math:`O(N)` beyond the positions.

    Examples
    --------
    >>> box = Box(lo=np.array([-4.0, -4.0, -8.0]),
    ...           hi=np.array([4.0, 4.0, 8.0]))
    >>> data = chain_data(helix(12, radius=2.0, pitch=2.5), box)
    >>> report = check_conformation(data, r0=1.5, ghost_cutoff=1.42,
    ...                             periodic=(False, False, True),
    ...                             dcyl=6.0, dmon=1.0)
    >>> round(report.max_bond, 6), round(report.max_radial, 6)
    (0.97, 2.0)
    """
    for name, value in (("r0", r0), ("ghost_cutoff", ghost_cutoff)):
        _check_positive(name, value)
    if min_distance < 0.0:
        raise ValueError(f"'min_distance' must be 0 or more, got "
                         f"{min_distance}")
    axes = tuple(bool(axis) for axis in periodic)
    if len(axes) != 3:
        raise ValueError(f"'periodic' must have three entries, got "
                         f"{len(axes)}")
    if (dcyl is None) != (dmon is None):
        raise ValueError("'dcyl' and 'dmon' must be given together")
    if dcyl is not None and dmon is not None and not 0.0 < dmon < dcyl:
        raise ValueError(f"'dcyl' must exceed 'dmon' and both must be "
                         f"positive, got {dcyl} and {dmon}")
    box = data.box
    if box.triclinic and any(axes):
        raise ValueError("a periodic axis needs an orthogonal box")
    positions = data.positions
    n_atoms = data.n_atoms
    lengths = box.lengths
    periodic_axes = np.array(axes)
    failures: list[str] = []

    def minimum_image(delta: FloatArray) -> FloatArray:
        if periodic_axes.any():
            span = lengths[periodic_axes]
            delta[..., periodic_axes] -= span * np.round(
                delta[..., periodic_axes] / span)
        return delta

    order = np.argsort(data.ids, kind="stable")
    rows = order[np.searchsorted(data.ids[order], data.bonds[:, 1:])]
    bond_vectors = minimum_image(
        positions[rows[:, 1]] - positions[rows[:, 0]])
    bond_lengths = np.linalg.norm(bond_vectors, axis=1)
    warning_length = FENE_WARNING_FRACTION * r0
    for limit, label in ((warning_length,
                          f"sqrt(0.9) R0 = {warning_length:.4g}"),
                         (ghost_cutoff,
                          f"the ghost cutoff {ghost_cutoff:.4g}")):
        over = bond_lengths >= limit
        if over.any():
            worst = int(np.argmax(bond_lengths))
            failures.append(
                f"{_count(int(over.sum()), 'bond', 'bonds')} at or above "
                f"{label}: longest {bond_lengths[worst]:.4g} (atoms "
                f"{data.bonds[worst, 1]} {data.bonds[worst, 2]})")
    bonded_keys = np.unique(rows.min(axis=1) * n_atoms + rows.max(axis=1))
    min_pair, pair_atoms = math.inf, (0, 0)
    min_any, any_atoms = math.inf, (0, 0)
    n_close = n_coincident = 0
    chunk = max(1, 2_000_000 // max(n_atoms, 1))
    for start in range(0, n_atoms, chunk):
        stop = min(start + chunk, n_atoms)
        delta = minimum_image(
            positions[None, :, :] - positions[start:stop, None, :])
        squared = np.einsum("ijk,ijk->ij", delta, delta)
        i_rows = np.arange(start, stop)[:, None]
        j_cols = np.arange(n_atoms)[None, :]
        squared = np.where(j_cols > i_rows, squared, math.inf)
        n_coincident += int(np.count_nonzero(squared < _COINCIDENT ** 2))
        flat = int(np.argmin(squared))
        if squared.flat[flat] < min_any ** 2:
            min_any = math.sqrt(squared.flat[flat])
            any_atoms = (start + flat // n_atoms, flat % n_atoms)
        bonded = np.isin(i_rows * n_atoms + j_cols, bonded_keys)
        squared = np.where(bonded, math.inf, squared)
        n_close += int(np.count_nonzero(squared < min_distance ** 2))
        flat = int(np.argmin(squared))
        if squared.flat[flat] < min_pair ** 2:
            min_pair = math.sqrt(squared.flat[flat])
            pair_atoms = (start + flat // n_atoms, flat % n_atoms)
    if n_coincident:
        failures.append(
            f"{_count(n_coincident, 'pair of atoms', 'pairs of atoms')} "
            f"closer than "
            f"{_COINCIDENT:g}: {min_any:.3g} (atoms "
            f"{data.ids[any_atoms[0]]} {data.ids[any_atoms[1]]})")
    if n_close:
        failures.append(
            f"{_count(n_close, 'non-bonded pair', 'non-bonded pairs')} "
            f"closer than "
            f"{min_distance:.4g}: closest {min_pair:.4g} (atoms "
            f"{data.ids[pair_atoms[0]]} {data.ids[pair_atoms[1]]})")
    if box.triclinic:
        fractions = (positions - box.lo) @ np.linalg.inv(box.edges())
    else:
        fractions = (positions - box.lo) / lengths
    outside = ((fractions < 0.0) | (fractions >= 1.0))[:, ~periodic_axes]
    out_rows = np.flatnonzero(outside.any(axis=1))
    if out_rows.size:
        failures.append(
            f"{_count(out_rows.size, 'atom', 'atoms')} outside the box on a "
            f"non-periodic axis (atoms "
            f"{' '.join(map(str, data.ids[out_rows[:3]].tolist()))}"
            f"{' ...' if out_rows.size > 3 else ''})")
    max_radial: float | None = None
    if dcyl is not None and dmon is not None:
        centre = 0.5 * (box.lo + box.hi)
        radial = np.hypot(positions[:, 0] - centre[0],
                          positions[:, 1] - centre[1])
        max_radial = float(radial.max()) if n_atoms else math.nan
        limit = 0.5 * (dcyl - dmon)
        beyond = radial > limit
        if beyond.any():
            worst = int(np.argmax(radial))
            failures.append(
                f"{_count(int(beyond.sum()), 'atom', 'atoms')} farther than "
                f"(D - a)/2 = {limit:.4g} from the cylinder axis: farthest "
                f"{radial[worst]:.4g} (atom {data.ids[worst]})")
    if failures:
        raise ConformationError("; ".join(failures))
    return ConformationReport(
        n_atoms=n_atoms,
        min_bond=float(bond_lengths.min()) if bond_lengths.size else math.nan,
        max_bond=float(bond_lengths.max()) if bond_lengths.size else math.nan,
        min_pair=min_pair if math.isfinite(min_pair) else math.nan,
        max_radial=max_radial,
    )


def _count(n: int, singular: str, plural: str) -> str:
    """Return `n` with the fitting one of `singular` and `plural`."""
    return f"{n} {singular if n == 1 else plural}"


def _check_count(name: str, value: int, minimum: int) -> None:
    """Raise unless `value` is an integer of at least `minimum`."""
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"'{name}' must be an integer, got {value!r}")
    if value < minimum:
        raise ValueError(f"'{name}' must be {minimum} or more, got {value}")


def _check_positive(name: str, value: float) -> None:
    """Raise unless `value` is a finite positive number."""
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"'{name}' must be positive, got {value}")


def _check_closed(gap: float | None, bond: float) -> None:
    """Raise unless the closing chord `gap` equals `bond` to round-off; the
    root finder lands on a jump of the closing gap otherwise, which happens
    when a curve has features smaller than the bond."""
    if gap is None or abs(gap - bond) > 1e-9 * bond:
        raise ValueError(
            "the curve cannot be closed with equal bonds: its corners or "
            "features are smaller than the bond length")


def _chord_angle(bond: float, radius: float, rise: float) -> float:
    """Return the angle step at which the chord of the helix of `radius`
    and rise per radian `rise` equals `bond`; the root of the chord
    equation on (0, pi]."""
    def chord_excess(angle: float) -> float:
        return (4.0 * radius * radius * math.sin(0.5 * angle) ** 2
                + rise * rise * angle * angle - bond * bond)
    if chord_excess(math.pi) < 0.0:
        raise ValueError(
            f"a bond of {bond} does not fit in half a turn of a helix of "
            f"radius {radius} and pitch {2 * math.pi * rise}")
    angle: float = brentq(chord_excess, 0.0, math.pi, xtol=1e-14,
                          maxiter=200)
    return angle


def _plectoneme_curve(u: FloatArray, strand: float,
                      theta: float) -> FloatArray:
    """Return points of the unit-radius Neukirch-Starostin plasmid at the
    parameters `u` in [0, 2 strand + 2 pi]: strand A for u in [0, strand),
    loop B for the next pi, strand C for the next `strand`, loop D for the
    last pi (Neukirch and Starostin, 2008, Eqs. A1 to A4)."""
    u = np.asarray(u, dtype=np.float64)
    sin_t, cos_t = math.sin(theta), math.cos(theta)
    turn = 1.0 + 2.0 * theta / math.pi
    twist = strand * sin_t                   # the angle strand A turns
    points = np.empty(u.shape + (3,), dtype=np.float64)
    in_a = u < strand
    s = u[in_a]
    points[in_a] = np.stack(
        [np.sin(s * sin_t), -np.cos(s * sin_t), s * cos_t], axis=-1)
    in_b = (u >= strand) & (u < strand + math.pi)
    s = u[in_b] - strand
    theta_b = theta - s * turn
    points[in_b] = np.stack([
        math.sin(twist) * np.cos(s)
        + math.cos(twist) * sin_t * np.sin(s) * np.cos(s),
        -math.cos(twist) * np.cos(s)
        + math.sin(twist) * sin_t * np.sin(s) * np.cos(s),
        strand * cos_t + (sin_t - np.sin(theta_b)) / turn], axis=-1)
    in_c = (u >= strand + math.pi) & (u < 2.0 * strand + math.pi)
    s = u[in_c] - strand - math.pi
    psi = twist - s * sin_t
    points[in_c] = np.stack(
        [-np.sin(psi), np.cos(psi), (strand - s) * cos_t], axis=-1)
    in_d = u >= 2.0 * strand + math.pi
    s = np.minimum(u[in_d] - 2.0 * strand - math.pi, math.pi)
    theta_d = -math.pi - theta + s * turn
    points[in_d] = np.stack([
        sin_t * np.sin(s) * np.cos(s), np.cos(s),
        (np.sin(theta_d) - sin_t) / turn], axis=-1)
    return points


def _plectoneme_loop_length(theta: float, step: float) -> float:
    """Return the arclength of the two end loops of the unit-radius
    plasmid, the whole curve when the strands have zero length."""
    u = np.arange(0.0, 2.0 * math.pi + step, step)
    points = _plectoneme_curve(u, 0.0, theta)
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def _march_curve(
    curve: _Curve, u_end: float, n_beads: int, bond: float, step: float
) -> tuple[FloatArray, float | None]:
    """Place up to `n_beads` beads `bond` apart along `curve`, sampled at
    parameter steps of `step` from 0 to `u_end`, each bead at the first
    crossing of the sphere about the bead before; return the beads and the
    distance from the last bead back to the first, or None when the curve
    ran out. O(n_beads + u_end / step) time and memory."""
    u = np.arange(0.0, u_end, step)
    if u[-1] < u_end:
        u = np.append(u, u_end)
    samples = curve(u)
    beads = np.empty((n_beads, 3), dtype=np.float64)
    beads[0] = samples[0]
    window = int(3.0 * bond / step) + 16
    index = 0
    for k in range(1, n_beads):
        last = beads[k - 1]
        beyond = _first_beyond(samples, index, window, last, bond)
        if beyond is None:
            return beads[:k], None
        beyond = max(beyond, index + 1)

        def excess(parameter: float) -> float:
            point = curve(np.array([parameter]))[0]
            return float(np.dot(point - last, point - last)) - bond * bond

        lo, hi = u[beyond - 1], u[beyond]
        if excess(lo) >= 0.0:
            crossing = lo
        else:
            crossing = brentq(excess, lo, hi, xtol=1e-14, maxiter=200)
        beads[k] = curve(np.array([crossing]))[0]
        index = beyond - 1
    return beads, float(np.linalg.norm(beads[-1] - beads[0]))


def _first_beyond(
    samples: FloatArray, start: int, window: int, centre: FloatArray,
    radius: float,
) -> int | None:
    """Return the index of the first sample from `start` on that is farther
    than `radius` from `centre`, scanning `window` samples at a time."""
    while True:
        block = samples[start:start + window]
        delta = block - centre
        outside = np.flatnonzero(
            np.einsum("ij,ij->i", delta, delta) > radius * radius)
        if outside.size:
            return start + int(outside[0])
        if start + window >= samples.shape[0]:
            return None
        start += window - 1


def _march_polyline(
    vertices: FloatArray, n_beads: int, bond: float
) -> tuple[FloatArray, float | None]:
    """Place up to `n_beads` beads `bond` apart along the polyline through
    `vertices`, each at the first crossing of the sphere about the bead
    before, by exact line-sphere intersection; return the beads and the
    distance from the last bead back to the first vertex, or None when the
    polyline ran out. O(n_beads + n_vertices) time and memory."""
    segments = np.diff(vertices, axis=0)
    beads = np.empty((n_beads, 3), dtype=np.float64)
    beads[0] = vertices[0]
    segment = 0
    for k in range(1, n_beads):
        last = beads[k - 1]
        end = _first_beyond_vertex(vertices, segment + 1, last, bond)
        if end is None:
            return beads[:k], None
        segment = end - 1
        offset = vertices[segment] - last
        along = segments[segment]
        alpha = float(np.dot(along, along))
        beta = float(np.dot(offset, along))
        gamma = float(np.dot(offset, offset)) - bond * bond
        root = math.sqrt(max(beta * beta - alpha * gamma, 0.0))
        if beta <= 0.0:
            fraction = (root - beta) / alpha
        else:
            fraction = -gamma / (beta + root)
        beads[k] = vertices[segment] + fraction * along
    return beads, float(np.linalg.norm(beads[-1] - vertices[0]))


def _first_beyond_vertex(
    vertices: FloatArray, start: int, centre: FloatArray, radius: float
) -> int | None:
    """Return the index of the first vertex from `start` on that is at
    least `radius` from `centre`, scanning 256 vertices at a time."""
    while start < vertices.shape[0]:
        block = vertices[start:start + 256]
        delta = block - centre
        outside = np.flatnonzero(
            np.einsum("ij,ij->i", delta, delta) >= radius * radius)
        if outside.size:
            return start + int(outside[0])
        start += 256
    return None
