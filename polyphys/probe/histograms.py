"""
==========================================================
:mod:`polyphys.probe.histograms`
==========================================================

The :mod:`polyphys.probe.histograms` module bins the positions of every
species of a study in every frame of a trajectory: one-dimensional
histograms along the directions of the geometry, kept per frame, and
two-dimensional maps in the three Cartesian planes, kept as the mean map
of each block of consecutive frames. :class:`SpatialHistograms` is an
MDAnalysis :class:`~MDAnalysis.analysis.base.AnalysisBase`, so it runs over
a slice of frames in serial or, split across workers, in parallel;
:func:`histogram_edges` builds the bin edges of one whole from its probe
spec and parser, with the research code's bin sizes and the studies'
geometry.

The research code binned the same species in one loop over the frames of
each study's probe function, after shifting every position by half the
parser's box, and summed the counts over the frames. Here the one-dimensional
counts are kept per frame, so that the time dependence and the block
statistics of a profile remain available, and the two-dimensional maps,
which would not fit in memory per frame, are averaged over blocks of
``block_frames`` frames. The sum over frames of either is what the research
code saved.

Contents
========

- :class:`SpatialHistograms`, the per-frame and per-block binning.
- :func:`histogram_edges`, the bin edges of a whole.

Dependencies
============
- Python standard library: `math`, `pathlib`, `typing`.
- `NumPy <https://numpy.org/>`_.
- `MDAnalysis <https://docs.mdanalysis.org/>`_ 2.8 or later, through the
  ``mda`` extra: :class:`~MDAnalysis.analysis.base.AnalysisBase` and
  :class:`~MDAnalysis.analysis.results.ResultsGroup`.
- :mod:`polyphys.analyze.measurer`: the histogram helpers
  :func:`~polyphys.analyze.measurer.fixedsize_bins`,
  :func:`~polyphys.analyze.measurer.radial_histogram`,
  :func:`~polyphys.analyze.measurer.radial_cyl_histogram`,
  :func:`~polyphys.analyze.measurer.axial_histogram`,
  :func:`~polyphys.analyze.measurer.azimuth_cyl_histogram`, and
  :func:`~polyphys.analyze.measurer.planar_cartesian_histogram`.
- :mod:`polyphys.manage.projects`:
  :class:`~polyphys.manage.projects.ProbeSpec` and
  :data:`~polyphys.manage.projects.RecenterT`.
- :mod:`polyphys.manage.parser`:
  :class:`~polyphys.manage.parser.ParserBase`.
- :mod:`polyphys.probe.core`: ``_compute_part``, which computes each
  part of a split run on its own copy of the analysis.

Notation
========
:math:`F`, :math:`b`, :math:`K = \\lceil F / b \\rceil`
    The frames of a run, the block length, and the number of blocks; block
    :math:`k` holds the frames :math:`kb, \\ldots, \\min((k + 1) b, F) - 1`
    in the order the run visits them.
:math:`\\mathbf{r}_i = (x_i, y_i, z_i)`
    The position of particle :math:`i` after recentering, in the dump's
    length unit, the monomer diameter in the studies.
:math:`r`, :math:`\\theta`, :math:`z`
    In a cylinder whose axis is :math:`z`: :math:`\\sqrt{x^2 + y^2}`,
    :math:`\\operatorname{atan2}(y, x)` in :math:`[-\\pi, \\pi]`, and the
    axial coordinate. In a cube, :math:`r = |\\mathbf{r}|` from the box
    center.
:math:`D`, :math:`Z`, :math:`L`, :math:`a_w`
    The parser's ``dcyl``, ``lcyl``, and ``lcube``, and the wall-bead size
    :math:`a_w = 1` of the cylindrical studies, so that the LAMMPS box is
    :math:`(D + a_w) \\times (D + a_w) \\times Z` or :math:`L^3`.

Examples
========
Two particles in a periodic cube of side 20, three frames, with the dump
reader's frame already translated to the origin, so that ``'half'`` moves
the box center to the origin:

>>> import MDAnalysis as mda
>>> import numpy as np
>>> frames = np.array([
...     [[10.0, 10.0, 10.0], [12.5, 10.0, 10.0]],
...     [[10.0, 10.0, 10.0], [10.0, 10.0, 11.0]],
...     [[10.0, 10.0, 10.0], [10.0, 13.5, 10.0]]])
>>> universe = mda.Universe.empty(2, trajectory=True)
>>> universe.load_new(frames, format='MEMORY', dimensions=np.tile(
...     [20.0, 20.0, 20.0, 90.0, 90.0, 90.0], (3, 1)))
<Universe with 2 atoms>
>>> edges = {'r': np.arange(0.0, 17.5, 1.0)}
>>> for key in 'xyz':
...     edges[key] = np.linspace(-10.0, 10.0, 5)
>>> analysis = SpatialHistograms(
...     {'Mon': universe.atoms}, edges, geometry='cubic', block_frames=2)
>>> analysis = analysis.run()
>>> analysis.results['rHistTMon'][:, :4]
array([[1, 0, 1, 0],
       [1, 1, 0, 0],
       [1, 0, 0, 1]])
>>> analysis.results['blockFrames']
array([2, 1])

.. versionadded:: 0.7.0
"""
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Self

import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe

from ..analyze import measurer
from ..manage.parser import ParserBase
from ..manage.projects import ProbeSpec, RecenterT
from ..manage.types import (
    AxisT,
    BinT,
    DirectionT,
    EntityT,
    GeometryT,
    Pathish,
)
from .core import _compute_part

# The size of the beads whose centers form the cylindrical wall, which the
# parsers document: the LAMMPS box spans the wall-bead centers, so its
# transverse length is dcyl + _WALL_BEAD.
_WALL_BEAD = 1.0

# The directions of the one-dimensional histograms in each geometry, in the
# research code's order, and the edges each geometry's two-dimensional maps
# need. 'z' edges serve both the axial histogram and the maps of a cylinder.
_DIRECTIONS: Mapping[GeometryT, tuple[DirectionT, ...]] = {
    'cubic': ('r',),
    'cylindrical': ('r', 'z', 'theta'),
}
_CARTESIAN: tuple[DirectionT, ...] = ('x', 'y', 'z')
_PLANES: tuple[str, ...] = ('xy', 'xz', 'yz')
_AXIS: AxisT = 2  # the cylinder axis, z, in every study

# The directions each geometry makes periodic, along which a position is
# wrapped into the box before binning: every one in a periodic cube, the
# axis only in a cylinder (LAMMPS boundary ``ff ff pp``).
_PERIODIC: Mapping[GeometryT, tuple[int, ...]] = {
    'cubic': (0, 1, 2),
    'cylindrical': (_AXIS,),
}

# One frame's one-dimensional histogram along a direction.
_Histogram = Callable[[np.ndarray, np.ndarray], np.ndarray]


def _radial_cyl(positions: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return measurer.radial_cyl_histogram(positions, edges, _AXIS)


def _axial(positions: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return measurer.axial_histogram(positions, edges, _AXIS)


def _azimuth(positions: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return measurer.azimuth_cyl_histogram(positions, edges, _AXIS)


_HISTOGRAMS: Mapping[GeometryT, Mapping[DirectionT, _Histogram]] = {
    'cubic': {'r': measurer.radial_histogram},
    'cylindrical': {'r': _radial_cyl, 'z': _axial, 'theta': _azimuth},
}

_BACKENDS = ('serial', 'multiprocessing', 'dask')


def _required_edges(geometry: GeometryT) -> tuple[DirectionT, ...]:
    """The edge keys a geometry's histograms and maps need, in order."""
    return tuple(dict.fromkeys((*_DIRECTIONS[geometry], *_CARTESIAN)))


def _check_recenter(recenter: Any) -> RecenterT:
    """Validate a recentering and return it with a triple as floats."""
    if isinstance(recenter, str) and recenter == 'none':
        return 'none'
    if isinstance(recenter, str) and recenter == 'half':
        return 'half'
    if (isinstance(recenter, (tuple, list)) and len(recenter) == 3
            and all(isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(value) for value in recenter)):
        return (float(recenter[0]), float(recenter[1]), float(recenter[2]))
    raise ValueError(
        f"'recenter' must be 'none', 'half', or a triple of finite numbers, "
        f"not {recenter!r}.")


def _planar(
    positions: np.ndarray, plane: str, edges: Mapping[str, np.ndarray]
) -> np.ndarray:
    """One frame's map in `plane`, rows along its first axis."""
    if plane == 'xy':
        return measurer.planar_cartesian_histogram(
            positions, [edges['x'], edges['y']], 2)
    if plane == 'yz':
        return measurer.planar_cartesian_histogram(
            positions, [edges['y'], edges['z']], 0)
    # dim=1 bins (z, x); transpose to the research code's (x, z).
    return measurer.planar_cartesian_histogram(
        positions, [edges['z'], edges['x']], 1).T


def histogram_edges(
    spec: ProbeSpec, parser: ParserBase
) -> dict[str, np.ndarray]:
    """
    The bin edges of one whole's histograms and maps.

    The bin sizes are the spec's, resolved for this whole by
    :meth:`~polyphys.manage.projects.ProbeSpec.resolve`; the ranges follow
    the geometry, so that every position inside the confinement falls
    between the edges:

    - Cylinder: ``'r'`` over :math:`[0, (D + a_w) / 2]`, the wall surface,
      with 'nonnegative' bins; ``'z'`` over :math:`[-Z / 2, Z / 2]` and
      ``'x'``, ``'y'`` over :math:`[-(D + a_w) / 2, (D + a_w) / 2]` with
      'ordinary' bins; ``'theta'`` over :math:`[-\\pi, \\pi]` with
      'periodic' bins.
    - Cube: ``'r'`` over :math:`[0, \\sqrt{3} L / 2]`, the corner, with
      'nonnegative' bins; ``'x'``, ``'y'``, ``'z'`` over
      :math:`[-L / 2, L / 2]` with 'ordinary' bins.

    The research code ranged ``'r'``, ``'x'``, and ``'y'`` of a cylinder
    over :math:`D / 2` and ``'r'`` of a cube over :math:`L / 2`, and
    silently dropped the positions beyond them; the maintainer extended the
    ranges on 2026-11-04, since :class:`SpatialHistograms` raises on such a
    position instead.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    spec : :class:`~polyphys.manage.projects.ProbeSpec`
        The study's probe spec, whose ``bin_sizes`` holds every key that
        ``geometry`` needs: ``'r'``, ``'z'``, ``'theta'``, ``'x'``, ``'y'``
        for a cylinder; ``'r'``, ``'x'``, ``'y'``, ``'z'`` for a cube.
    parser : :class:`~polyphys.manage.parser.ParserBase`
        The whole's parser, with ``dcyl`` and ``lcyl`` or ``lcube``.

    Returns
    -------
    dict of (str -> numpy.ndarray)
        The edges keyed by direction, each float64 and increasing, with
        ``n_bins + 1`` entries, in the order the geometry needs them.

    Raises
    ------
    ValueError
        If ``spec.bin_sizes`` lacks a key the geometry needs, or if a
        resolved bin size is not a finite positive number.

    Notes
    -----
    :math:`O(B)` time and space for :math:`B` edges in all.

    Examples
    --------
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> from polyphys.manage.projects import ProbeSpec, ProjectSpec
    >>> project = ProjectSpec(
    ...     parser=SumRuleCyl, group='bug', divisor=0.025,
    ...     physical_attrs=('space', 'nmon', 'dcyl', 'dcrowd'),
    ...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
    ...     tseries_properties=('gyrTMon',))
    >>> spec = ProbeSpec(
    ...     project, {'Mon': 'resid 1', 'Crd': 'resid 0'}, 'Mon', ('gyrT',),
    ...     unwrap=False, histogram_species=('Mon', 'Crd'), block_frames=100,
    ...     bin_sizes={'r': 0.5, 'z': 1.0, 'theta': np.pi / 2, 'x': 0.5,
    ...                'y': 0.5})
    >>> whole = SumRuleCyl(
    ...     'N20epsilon5r2.5lz5sig1nc0dt0.005bdump1adump1ens1', 'whole',
    ...     'all')
    >>> edges = histogram_edges(spec, whole)
    >>> list(edges)
    ['r', 'z', 'theta', 'x', 'y']
    >>> edges['r']
    array([0. , 0.5, 1. , 1.5, 2. , 2.5])
    >>> edges['z'][[0, -1]], edges['theta'].size
    (array([-5.,  5.]), 5)
    """
    geometry = spec.geometry
    missing = [key for key in _required_edges(geometry)
               if key not in spec.bin_sizes]
    if missing:
        raise ValueError(
            f"{spec.name}: 'bin_sizes' lacks {missing}, which the "
            f"{geometry} histograms need.")
    sizes = {key: ProbeSpec.resolve(spec.bin_sizes[key], parser)
             for key in _required_edges(geometry)}
    ranges: dict[DirectionT, tuple[float, float, BinT]]
    if geometry == 'cylindrical':
        half_box = 0.5 * (float(getattr(parser, 'dcyl')) + _WALL_BEAD)
        half_length = 0.5 * float(getattr(parser, 'lcyl'))
        ranges = {
            'r': (0.0, half_box, 'nonnegative'),
            'z': (-half_length, half_length, 'ordinary'),
            'theta': (-np.pi, np.pi, 'periodic'),
            'x': (-half_box, half_box, 'ordinary'),
            'y': (-half_box, half_box, 'ordinary'),
        }
    else:
        half_side = 0.5 * float(getattr(parser, 'lcube'))
        ranges = {
            'r': (0.0, math.sqrt(3.0) * half_side, 'nonnegative'),
            'x': (-half_side, half_side, 'ordinary'),
            'y': (-half_side, half_side, 'ordinary'),
            'z': (-half_side, half_side, 'ordinary'),
        }
    return {
        key: np.asarray(
            measurer.fixedsize_bins(sizes[key], lmin, lmax,
                                    bin_type)['bin_edges'],
            dtype=np.float64)
        for key, (lmin, lmax, bin_type) in ranges.items()}


class SpatialHistograms(AnalysisBase):
    """
    The positions of each species binned in every frame: one-dimensional
    histograms per frame and two-dimensional maps per block of frames.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    groups : mapping of (:py:data:`EntityT` -> AtomGroup)
        Each species label mapped to its atoms, all in one universe, in the
        order ``results`` and :meth:`save` list the species. A group may be
        empty, for a species a whole lacks; its histograms are then zero.
    edges : mapping of (str -> numpy.ndarray)
        The bin edges, increasing, keyed as :func:`histogram_edges` returns
        them: the directions of `geometry` and ``'x'``, ``'y'``, ``'z'``.
    geometry : :py:data:`GeometryT`
        ``'cylindrical'``, with the histograms along ``'r'``, ``'z'``, and
        ``'theta'`` about the ``z`` axis, or ``'cubic'``, with the
        histogram along ``'r'`` from the box center.
    block_frames : int
        The block length :math:`b`, positive: consecutive frames of the
        run whose maps are averaged into one stored map.
    recenter : :py:data:`~polyphys.manage.projects.RecenterT`, default 'half'
        How a frame's positions are moved before binning: ``'half'``
        subtracts half the frame's box lengths, ``'none'`` nothing, a
        triple its values. The positions in the trajectory are untouched.
    verbose : bool, default False
        Whether to show a progress bar during :meth:`run`.

    Attributes
    ----------
    groups : mapping of (:py:data:`EntityT` -> AtomGroup)
        The species, a copy of the mapping given.
    edges : mapping of (str -> numpy.ndarray)
        The edges, a copy of the mapping given, each float64.
    geometry : :py:data:`GeometryT`
        The geometry.
    block_frames : int
        The block length.
    recenter : :py:data:`~polyphys.manage.projects.RecenterT`
        The recentering, a triple as floats.
    results : :class:`~MDAnalysis.analysis.results.Results`
        After :meth:`run`, for each species ``S`` and direction ``d`` of
        the geometry, ``dHistTS``, int64 ``(n_frames, n_bins)``, the counts
        of frame ``i`` in row ``i``, and ``dHistS``, int64 ``(n_bins,)``,
        their sum over the frames; for each plane ``p`` in ``xy``, ``xz``,
        ``yz``, ``pHistBlockS``, float64 ``(n_blocks, n_a, n_b)``, the mean
        map of each block, rows along the plane's first axis, and
        ``pHistS``, float64 ``(n_a, n_b)``, the counts summed over the
        frames; and ``blockFrames``, int64 ``(n_blocks,)``, the frames in
        each block.
    n_frames : int
        After :meth:`run`, the number of frames analyzed, from
        :class:`~MDAnalysis.analysis.base.AnalysisBase`.
    frames : numpy.ndarray
        After :meth:`run`, the indices of the frames analyzed, shape
        ``(n_frames,)``.
    times : numpy.ndarray
        After :meth:`run`, the times of the frames analyzed, shape
        ``(n_frames,)``, in the trajectory's time unit.

    Raises
    ------
    ValueError
        If `groups` is empty, has an empty label, or spans two universes;
        if `edges` lacks a key the geometry needs or has edges that are
        not a one-dimensional increasing array of at least two values; if
        `geometry` is unknown; if `block_frames` is not a positive integer;
        if `recenter` is neither ``'none'``, ``'half'``, nor a triple of
        finite numbers; or, during :meth:`run`, if a frame carries no box,
        or if a wrapped and recentered position falls outside the edges of
        a direction or a Cartesian axis.

    Notes
    -----
    Per frame and species, the positions are copied, wrapped, recentered,
    and binned. Wrapping maps each position into the box of the reader's
    frame, :math:`[0, L)`, along the directions the geometry makes
    periodic, every one in a cube and ``z`` in a cylinder: LAMMPS enforces
    the periodic boundaries only when it rebuilds the neighbor lists, so a
    dump can place an atom slightly outside the box [1]_, as the HnsCub
    dumps do by up to :math:`0.04` monomer diameters. The box of the
    reader's frame is MDAnalysis's primary unit cell, :math:`[0, L)`, where
    its LAMMPS dump reader puts the box by subtracting the lower bounds; a
    universe built in memory must place its positions there too. A position
    past a wall is not wrapped. The positions are binned with
    :func:`~polyphys.analyze.measurer.radial_cyl_histogram`,
    :func:`~polyphys.analyze.measurer.axial_histogram`, and
    :func:`~polyphys.analyze.measurer.azimuth_cyl_histogram` in a cylinder,
    :func:`~polyphys.analyze.measurer.radial_histogram` in a cube, and
    :func:`~polyphys.analyze.measurer.planar_cartesian_histogram` in each
    plane, with the ``xz`` map transposed to rows along ``x`` as the
    research code stored it. A position outside the edges after wrapping
    raises, so every frame's counts sum to the species' particle number:
    the research code's ``np.histogram`` dropped such a position silently,
    and did not wrap.

    Blocks are numbered over the frames of the whole run in the order
    :meth:`run` visits them, whatever ``start``, ``stop``, ``step``, or
    ``frames`` select, so that a run split across workers gives the same
    blocks as a serial one: :meth:`_setup_frames` records the run's frame
    count before the split and :meth:`_compute` each worker's offset into
    it. Each worker accumulates the block sums and frame counts of its own
    frames, which :meth:`_get_aggregator` adds
    (``ResultsGroup.ndarray_sum``), and its per-frame rows, which it
    concatenates (``ResultsGroup.ndarray_vstack``); :meth:`_conclude` then
    divides each block's sum by its frame count and sums the counts over
    the frames.

    Cost per frame: :math:`O(S N \\log B)` time for :math:`S` species of
    :math:`N` particles and :math:`B` bins along a direction, from the
    sorted-edge search of ``np.histogram`` and ``np.histogram2d``, and
    :math:`O(N)` extra space for the recentered positions. Over a run, each
    one-dimensional histogram holds :math:`8 F B` bytes and each map
    :math:`8 K B_a B_b` bytes, for :math:`F` frames, :math:`K = \\lceil F /
    b \\rceil` blocks, and :math:`B_a \\times B_b` cells: with
    :math:`b = 100` and :math:`200 \\times 200` cells, about 1 MB per
    block per species for the three planes. Wrapping adds :math:`O(N)`
    time per frame and species.

    References
    ----------
    .. [1] LAMMPS Developers, "dump command", LAMMPS documentation, note
       on periodic boundaries: "the coordinates of an atom written to a dump
       file may be slightly outside the simulation box".
       https://docs.lammps.org/dump.html

    Examples
    --------
    One particle in a cylinder of radius 2 and length 4, two frames in the
    reader's frame, the box from 0 to 4, recentered by a given shift to
    (1.5, 0, 1) and (0, 0.5, -1); the frame's counts sum to its one
    particle along every direction:

    >>> import MDAnalysis as mda
    >>> import numpy as np
    >>> frames = np.array([[[3.5, 2.0, 3.0]], [[2.0, 2.5, 1.0]]])
    >>> universe = mda.Universe.empty(1, trajectory=True)
    >>> universe.load_new(frames, format='MEMORY', dimensions=np.tile(
    ...     [4.0, 4.0, 4.0, 90.0, 90.0, 90.0], (2, 1)))
    <Universe with 1 atoms>
    >>> edges = {'r': np.array([0.0, 1.0, 2.0]),
    ...          'z': np.array([-2.0, 0.0, 2.0]),
    ...          'theta': np.linspace(-np.pi, np.pi, 5),
    ...          'x': np.array([-2.0, 0.0, 2.0]),
    ...          'y': np.array([-2.0, 0.0, 2.0])}
    >>> analysis = SpatialHistograms(
    ...     {'Crd': universe.atoms}, edges, geometry='cylindrical',
    ...     block_frames=2, recenter=(2.0, 2.0, 2.0)).run()
    >>> analysis.results['rHistTCrd'], analysis.results['zHistTCrd']
    (array([[0, 1],
           [1, 0]]), array([[0, 1],
           [1, 0]]))
    >>> analysis.results['xyHistBlockCrd']
    array([[[0., 0.],
            [0., 1.]]])
    >>> float(analysis.results['xyHistCrd'].sum())
    2.0
    """

    _analysis_algorithm_is_parallelizable = True

    def __init__(
        self,
        groups: Mapping[EntityT, AtomGroup],
        edges: Mapping[str, np.ndarray],
        *,
        geometry: GeometryT,
        block_frames: int,
        recenter: RecenterT = 'half',
        verbose: bool = False,
    ) -> None:
        groups = dict(groups)
        if not groups or any(not (isinstance(label, str) and label)
                             for label in groups):
            raise ValueError(
                f"'groups' must map each non-empty species label to its "
                f"atoms, not {list(groups)!r}.")
        universes = {id(group.universe) for group in groups.values()}
        if len(universes) != 1:
            raise ValueError("'groups' must all belong to one universe.")
        if geometry not in _DIRECTIONS:
            raise ValueError(
                f"'geometry' must be one of {list(_DIRECTIONS)}, not "
                f"{geometry!r}.")
        required = _required_edges(geometry)
        missing = [key for key in required if key not in edges]
        if missing:
            raise ValueError(
                f"'edges' lacks {missing}, which the {geometry} histograms "
                f"need.")
        checked: dict[str, np.ndarray] = {}
        for key in required:
            array = np.asarray(edges[key], dtype=np.float64)
            if array.ndim != 1 or array.size < 2 or np.any(
                    np.diff(array) <= 0):
                raise ValueError(
                    f"The edges along {key!r} must be a one-dimensional "
                    f"increasing array of at least two values.")
            checked[key] = array
        if not (isinstance(block_frames, int)
                and not isinstance(block_frames, bool) and block_frames > 0):
            raise ValueError(
                f"'block_frames' must be a positive integer, not "
                f"{block_frames!r}.")
        universe = next(iter(groups.values())).universe
        super().__init__(universe.trajectory, verbose=verbose)
        # Plain dicts, not read-only proxies: the multiprocessing backend
        # pickles the analysis, and a mappingproxy cannot be pickled.
        self.groups: Mapping[EntityT, AtomGroup] = groups
        self.edges: Mapping[str, np.ndarray] = checked
        self.geometry: GeometryT = geometry
        self.block_frames = block_frames
        self.recenter = _check_recenter(recenter)
        self._directions = _DIRECTIONS[geometry]
        self._histograms = _HISTOGRAMS[geometry]
        self._run_n_frames = 0
        self._run_offset = 0

    @classmethod
    def from_spec(
        cls,
        universe: Universe,
        spec: ProbeSpec,
        parser: ParserBase,
        *,
        verbose: bool = False,
    ) -> Self:
        """
        Build the analysis of a whole from its study's probe spec.

        The species are ``spec.histogram_species``, each selected with its
        string in ``spec.selections``; the edges are
        :func:`histogram_edges`'s; the geometry, block length, and
        recentering are the spec's.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        universe : :class:`~MDAnalysis.core.universe.Universe`
            The trajectory, as :func:`~polyphys.probe.core.make_universe`
            opens it.
        spec : :class:`~polyphys.manage.projects.ProbeSpec`
            The study's probe spec, with ``histogram_species``,
            ``bin_sizes``, and ``block_frames`` set.
        parser : :class:`~polyphys.manage.parser.ParserBase`
            The whole's parser, for the edges.
        verbose : bool, default False
            Whether to show a progress bar during :meth:`run`.

        Returns
        -------
        SpatialHistograms
            The analysis, not yet run.

        Raises
        ------
        ValueError
            If the spec names no histogram species, or from
            :func:`histogram_edges`.

        Examples
        --------
        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> from polyphys.manage.parser import SumRuleCyl
        >>> from polyphys.manage.projects import ProbeSpec, ProjectSpec
        >>> project = ProjectSpec(
        ...     parser=SumRuleCyl, group='bug', divisor=0.025,
        ...     physical_attrs=('space', 'nmon', 'dcyl', 'dcrowd'),
        ...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
        ...     tseries_properties=('gyrTMon',))
        >>> spec = ProbeSpec(
        ...     project, {'Mon': 'resid 1', 'Crd': 'resid 0'}, 'Mon',
        ...     ('gyrT',), unwrap=False, histogram_species=('Mon', 'Crd'),
        ...     block_frames=100, bin_sizes={'r': 0.5, 'z': 1.0,
        ...     'theta': np.pi / 2, 'x': 0.5, 'y': 0.5})
        >>> whole = SumRuleCyl(
        ...     'N2epsilon5r2.5lz5sig1nc1dt0.005bdump1adump1ens1', 'whole',
        ...     'all')
        >>> universe = mda.Universe.empty(
        ...     3, n_residues=2, atom_resindex=[0, 1, 1], trajectory=True)
        >>> universe.add_TopologyAttr('resid', [0, 1])
        >>> universe.load_new(np.array([[[1.0, 2.0, 3.0], [2.0, 2.0, 5.0],
        ...     [3.0, 2.0, 5.0]]]), format='MEMORY', dimensions=np.array(
        ...     [[5.0, 5.0, 10.0, 90.0, 90.0, 90.0]]))
        <Universe with 3 atoms>
        >>> analysis = SpatialHistograms.from_spec(universe, spec, whole)
        >>> [group.n_atoms for group in analysis.groups.values()]
        [2, 1]
        >>> analysis.run().results['rHistCrd']
        array([0, 0, 0, 1, 0])
        """
        if not spec.histogram_species or spec.block_frames is None:
            raise ValueError(
                f"{spec.name}: the probe spec names no histogram species.")
        groups = {species: universe.select_atoms(spec.selections[species])
                  for species in spec.histogram_species}
        return cls(groups, histogram_edges(spec, parser),
                   geometry=spec.geometry, block_frames=spec.block_frames,
                   recenter=spec.recenter, verbose=verbose)

    @classmethod
    def get_supported_backends(cls) -> tuple[str, ...]:
        """
        The backends :meth:`run` accepts: ``'serial'``,
        ``'multiprocessing'``, and ``'dask'``.

        Examples
        --------
        >>> SpatialHistograms.get_supported_backends()
        ('serial', 'multiprocessing', 'dask')
        """
        return _BACKENDS

    def _setup_frames(
        self, trajectory: Any, start: Any = None, stop: Any = None,
        step: Any = None, frames: Any = None,
    ) -> None:
        """
        Record the frame count of the whole run before it is split across
        workers, which fixes the number of blocks.
        """
        super()._setup_frames(trajectory, start=start, stop=stop, step=step,
                              frames=frames)
        self._run_n_frames = int(self.n_frames)
        self._run_offset = 0

    def _compute(
        self, indexed_frames: np.ndarray, verbose: Any = None, *,
        progressbar_kwargs: Any = None,
    ) -> Any:
        """
        Record where this part's frames start in the run, then compute the
        part on a copy of this analysis, so that serial parts do not
        overwrite each other.
        """
        self._run_offset = (int(indexed_frames[0, 0])
                            if len(indexed_frames) else 0)
        return _compute_part(self, indexed_frames, verbose,
                             progressbar_kwargs=progressbar_kwargs)

    def _get_aggregator(self) -> ResultsGroup:
        """
        Concatenate each worker's per-frame rows and add its block sums and
        frame counts.
        """
        lookup: dict[str, Callable[..., Any]] = {
            'blockFrames': ResultsGroup.ndarray_sum}
        for species in self.groups:
            for direction in self._directions:
                lookup[f"{direction}HistT{species}"] = (
                    ResultsGroup.ndarray_vstack)
            for plane in _PLANES:
                lookup[f"{plane}HistBlock{species}"] = (
                    ResultsGroup.ndarray_sum)
        return ResultsGroup(lookup=lookup)

    def _prepare(self) -> None:
        """
        Allocate the per-frame rows of this run or worker and the block sums
        and frame counts of the whole run.
        """
        n_blocks = -(-self._run_n_frames // self.block_frames)
        self.results['blockFrames'] = np.zeros(n_blocks, dtype=np.int64)
        for species in self.groups:
            for direction in self._directions:
                self.results[f"{direction}HistT{species}"] = np.zeros(
                    (self.n_frames, self.edges[direction].size - 1),
                    dtype=np.int64)
            for plane in _PLANES:
                self.results[f"{plane}HistBlock{species}"] = np.zeros(
                    (n_blocks, self.edges[plane[0]].size - 1,
                     self.edges[plane[1]].size - 1), dtype=np.float64)

    def _shift(self) -> np.ndarray:
        """The vector subtracted from this frame's positions."""
        if self.recenter == 'none':
            return np.zeros(3)
        if self.recenter == 'half':
            return 0.5 * np.asarray(self._ts.dimensions[:3],
                                    dtype=np.float64)
        return np.asarray(self.recenter, dtype=np.float64)

    def _wrap(self, positions: np.ndarray) -> np.ndarray:
        """Wrap `positions` into the frame's box, :math:`[0, L)`, along the
        geometry's periodic directions, in place."""
        if self._ts.dimensions is None:
            raise ValueError(
                f"Frame {self._ts.frame} carries no box, which wrapping "
                f"along the periodic directions needs.")
        axes = list(_PERIODIC[self.geometry])
        lengths = np.asarray(self._ts.dimensions[:3],
                             dtype=np.float64)[axes]
        positions[:, axes] -= lengths * np.floor(positions[:, axes] / lengths)
        return positions

    def _single_frame(self) -> None:
        """
        Wrap, recenter, and bin every species on the current frame: the
        histograms into this frame's rows, the maps into this frame's block.
        """
        block = (self._run_offset + self._frame_index) // self.block_frames
        self.results['blockFrames'][block] += 1
        shift = self._shift()
        for species, group in self.groups.items():
            positions = self._wrap(
                group.positions.astype(np.float64)) - shift
            for direction in self._directions:
                self.results[f"{direction}HistT{species}"][
                    self._frame_index] = self._histograms[direction](
                        positions, self.edges[direction])
            for axis, key in enumerate(_CARTESIAN):
                edges = self.edges[key]
                outside = int(np.count_nonzero(
                    (positions[:, axis] < edges[0])
                    | (positions[:, axis] > edges[-1])))
                if outside:
                    raise ValueError(
                        f"{outside} of {len(positions)} {species} positions "
                        f"fall outside the {key} edges "
                        f"[{edges[0]}, {edges[-1]}] in frame "
                        f"{self._ts.frame}.")
            for plane in _PLANES:
                self.results[f"{plane}HistBlock{species}"][block] += _planar(
                    positions, plane, self.edges)

    def _conclude(self) -> None:
        """
        Turn each block's summed map into its mean and sum the counts over
        the frames.
        """
        counts = self.results['blockFrames']
        filled = counts > 0
        for species in self.groups:
            for direction in self._directions:
                self.results[f"{direction}Hist{species}"] = self.results[
                    f"{direction}HistT{species}"].sum(axis=0)
            for plane in _PLANES:
                key = f"{plane}HistBlock{species}"
                sums = self.results[key]
                self.results[f"{plane}Hist{species}"] = sums.sum(axis=0)
                means = np.zeros_like(sums)
                means[filled] = sums[filled] / counts[filled, None, None]
                self.results[key] = means

    def save(self, save_to: Pathish, name: str) -> list[Path]:
        """
        Save the edges, histograms, and maps under the research code's names.

        For each species ``S`` and direction ``d``: ``<name>-dEdgeS.npy``,
        ``<name>-dHistS.npy`` (summed over the frames), and
        ``<name>-dHistTS.npy`` (per frame); for each Cartesian axis ``c``:
        ``<name>-cEdge.npy``; for each plane ``p`` and species:
        ``<name>-pHistS.npy`` (summed) and ``<name>-pHistBlockS.npy`` (per
        block); and ``<name>-blockFrames.npy``.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        save_to : str or pathlib.Path
            An existing directory.
        name : str
            The artifact's name, ``<parser.name>-<parser.group>`` in the
            studies, such as ``'N200...ens1-all'``.

        Returns
        -------
        list of pathlib.Path
            The files written, in the order above.

        Raises
        ------
        RuntimeError
            If :meth:`run` has not been called.

        Examples
        --------
        >>> import tempfile
        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> universe = mda.Universe.empty(1, trajectory=True)
        >>> universe.load_new(np.array([[[1.0, 1.0, 1.0]]]), format='MEMORY',
        ...                   dimensions=np.array([[4.0, 4.0, 4.0, 90.0,
        ...                                         90.0, 90.0]]))
        <Universe with 1 atoms>
        >>> edges = {key: np.array([-2.0, 0.0, 2.0]) for key in 'xyz'}
        >>> edges['r'] = np.array([0.0, 2.0, 4.0])
        >>> analysis = SpatialHistograms(
        ...     {'Mon': universe.atoms}, edges, geometry='cubic',
        ...     block_frames=1).run()
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     paths = analysis.save(tmp, 'N1ens1-all')
        ...     len(paths), paths[0].name, paths[-1].name
        (13, 'N1ens1-all-rEdgeMon.npy', 'N1ens1-all-blockFrames.npy')
        """
        if 'blockFrames' not in self.results:
            raise RuntimeError("Call run() before save().")
        folder = Path(save_to)
        items: list[tuple[str, np.ndarray]] = []
        for species in self.groups:
            for direction in self._directions:
                items.append((f"{direction}Edge{species}",
                              self.edges[direction]))
                items.append((f"{direction}Hist{species}",
                              self.results[f"{direction}Hist{species}"]))
                items.append((f"{direction}HistT{species}",
                              self.results[f"{direction}HistT{species}"]))
        for key in _CARTESIAN:
            items.append((f"{key}Edge", self.edges[key]))
        for species in self.groups:
            for plane in _PLANES:
                items.append((f"{plane}Hist{species}",
                              self.results[f"{plane}Hist{species}"]))
                items.append((f"{plane}HistBlock{species}",
                              self.results[f"{plane}HistBlock{species}"]))
        items.append(('blockFrames', self.results['blockFrames']))
        paths = []
        for suffix, array in items:
            path = folder / f"{name}-{suffix}.npy"
            np.save(path, array)
            paths.append(path)
        return paths
