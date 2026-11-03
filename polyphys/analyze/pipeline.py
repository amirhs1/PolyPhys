"""
==========================================================
:mod:`polyphys.analyze.pipeline`
==========================================================

The :mod:`polyphys.analyze.pipeline` module reduces the probe artifacts of one
space down the lineage, ``segment`` to ``whole`` to ``ensemble`` to
``ensemble average``, for every property that a study's analysis table
declares. It holds no physics of its own: it selects the files of each
property, calls :mod:`polyphys.manage.organizer` to build the wholes, the
ensembles, and the ensemble averages, and lets a *post-hook* derive further
artifacts from the wholes, such as a time average, an autocorrelation
function, a density profile, or the pair-separation histograms of the foci.

A study describes its analysis as data: one :class:`AnalysisTable` per
particle group, listing a :class:`PropertySpec` per probe artifact. The
tables of the maintainer's studies live beside each study's
:class:`~polyphys.manage.projects.ProjectSpec` in the repository's
``projects/`` folder, so a data correction or a special case is code in the
study's folder, never an option of this module.

Notation
========

*property*
    The name of one probe artifact without its lineage name and group, such
    as ``gyrTMon`` or ``rHistCrd``: a *name* (``gyrT``, ``rHist``) followed
    by a *species* (``Mon``, ``Crd``). Its files end in ``-<property>.npy``
    (a whole) or ``.<property>.npy`` (a segment).
*relation*
    How the segments of a whole combine, as
    :func:`~polyphys.manage.organizer.whole_from_segments` defines it:
    ``'tseries'`` concatenates along the frames, ``'histogram'`` sums, and
    ``'bin_edge'`` keeps one copy of identical edges.
*whole type*
    The shape of one whole, as :func:`~polyphys.manage.organizer.ensemble`
    defines it: ``'vector'`` (a 1-D array), ``'matrix'`` (an N-D array,
    stacked along a new first axis), ``'bin_edge'`` (edges shared by the
    wholes), or ``'dataframe'`` (a table whose columns are elements, such as
    the pairs of foci). The ensemble average of a vector is a DataFrame of
    ``-mean``, ``-var``, and ``-sem`` columns; of a matrix, one array per
    statistic; of edges, the edges; of a DataFrame, a DataFrame with those
    three statistics per column.
*post-hook*
    A callable ``(wholes, context) -> tuple[Derived, ...]`` applied to a
    property's wholes before the ensembles are built. It returns every
    artifact to reduce, so a hook can pass the wholes through under the
    property's own name, add derived artifacts, or replace the wholes with a
    reduction of them.
*lineage of a table*
    ``'segment'`` when the probe phase wrote the group's artifacts in
    segments, ``'whole'`` when it wrote one file per whole. A property of
    the study and the group, not of this module.
:math:`r_k`
    The autocorrelation at lag :math:`k` frames, from
    :func:`~polyphys.analyze.correlations.acf`.
:math:`\\bar P(R)`
    The probability density of the separation vector of a pair of foci, in
    one over the length unit cubed, from
    :func:`~polyphys.analyze.clusters.pair_separation_density`.

Contents
========

Records:
:class:`PropertySpec` describes one probe artifact; :func:`tseries`,
:func:`histogram`, :func:`histogram_2d`, :func:`edge`, and
:func:`matrix_tseries` build the kinds the studies use. :class:`AnalysisTable`
lists the properties of one group. :class:`Derived` is one artifact a hook
returns, and :class:`HookContext` what a hook may read.

Hooks:
:class:`TimeAverage` averages a time series over its frames;
:class:`Autocorrelation` adds the autocorrelation function of a series;
:class:`LocalDensities` adds the number density and volume fraction profiles
of a histogram through an injected :class:`DensityKernel`;
:class:`FociPairDistances` turns the distance-matrix time series of the foci
into their pair histograms, separation densities, and separation time series.
:func:`half_contour` is the separation bound of a ring at full bond extension.

Drivers:
:func:`reduce_space` runs every property of a table over one space;
:func:`reduce_stamps` combines the stamps of a space;
:func:`equilibrium_wholes` reduces each whole's time series to its
equilibrium value with a blocked error, beside the whole's stamps.

Conventions
===========

- A property whose files are not found raises ``FileNotFoundError``; nothing
  is skipped in silence.
- The wholes of one ensemble must hold identical bin edges, or ``ValueError``
  names the ensemble and the two wholes; the organizer takes the first
  whole's edges without comparing them.
- Outputs go to ``<root>/analysis/<prefix>-<group>-<stage>`` for the stages
  ``wholeSim`` (the ``segment`` lineage only), ``ens``, and ``ensAvg``, as
  :func:`~polyphys.manage.organizer.make_database` builds them, with the
  names :func:`~polyphys.manage.organizer.save_artifact` gives.
- The ensemble of ``'dataframe'`` wholes is their long form, every whole's
  rows with a ``whole`` and a ``row`` column, so no whole is lost before the
  ensemble average; the average is taken row by row across the wholes.
- The error of a whole's time average is the blocked standard error of
  :func:`~polyphys.analyze.stats.sem_blocked`, since the frames of one
  run are correlated; across the wholes of an ensemble, independent
  runs, the organizer's plain standard error stays.
- Progress is logged at INFO level; no warning is raised for a renamed
  artifact, because every output name is declared.
- Complexities are stated in ``n_files`` (the files of one property),
  ``n_frames``, ``n_bins``, ``n_wholes``, and ``n_pairs``.

Dependencies
============

- Python standard library: `dataclasses`, `glob`, `logging`, `os`, `re`.
- `numpy`, `pandas`.
- :mod:`polyphys.analyze.clusters`: the foci pair kernels.
- :mod:`polyphys.analyze.correlations`:
  :func:`~polyphys.analyze.correlations.acf`.
- :mod:`polyphys.analyze.stats`:
  :func:`~polyphys.analyze.stats.sem_blocked`.
- :mod:`polyphys.manage.organizer`: the lineage reducers.
- :mod:`polyphys.manage.projects`:
  :class:`~polyphys.manage.projects.ProjectSpec`.

References
==========

- K. Kremer and G. S. Grest, "Dynamics of entangled linear polymer melts: A
  molecular-dynamics simulation", J. Chem. Phys. 92, 5057 (1990), for the
  FENE bond and its maximum extension :math:`R_0 = 1.5\\,\\sigma`.
  https://doi.org/10.1063/1.458541

Examples
========
A table that reads as data: the gyration radius time series with its
autocorrelation, the radial histogram of the crowders with their density
profiles, and the planar histograms with their edges:

>>> from polyphys.analyze.pipeline import (
...     AnalysisTable, Autocorrelation, LocalDensities, edge, histogram,
...     histogram_2d, tseries)
>>> table = AnalysisTable(
...     group='all', lineage='segment', stamps=False, properties=(
...         tseries('gyrT', 'Mon', post_hook=Autocorrelation(nlags=100)),
...         histogram('r', 'Crd', post_hook=LocalDensities()),
...         histogram_2d('xy', 'Crd'),
...         edge('x'),
...     ))
>>> [prop.prop for prop in table.properties]
['gyrTMon', 'rHistCrd', 'xyHistCrd', 'xEdge']

.. versionadded:: 0.6.0
"""
import glob
import logging
import os
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, Protocol, cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from polyphys.analyze import clusters, stats
from polyphys.analyze.correlations import acf
from polyphys.analyze.stats import PlateauMethod
from polyphys.manage import organizer
from polyphys.manage.parser import ParserBase
from polyphys.manage.projects import ProjectSpec
from polyphys.manage.types import (
    DirectionT,
    EntityT,
    GroupT,
    Pathish,
    PrimitiveLineageT,
    PropertyT,
    WholeRelationT,
)

__all__ = [
    'AnalysisTable',
    'Autocorrelation',
    'DensityKernel',
    'Derived',
    'FociPairDistances',
    'HookContext',
    'LocalDensities',
    'PropertySpec',
    'TimeAverage',
    'edge',
    'equilibrium_wholes',
    'half_contour',
    'histogram',
    'histogram_2d',
    'matrix_tseries',
    'reduce_space',
    'reduce_stamps',
    'tseries',
]

logger = logging.getLogger(__name__)

#: The shape of one whole; see "whole type" in the module notation.
WholeTypeT = Literal['vector', 'matrix', 'bin_edge', 'dataframe']

#: The organizer's ensemble-average type for each array whole type.
_ENS_TYPES: Mapping[str, Literal['dataframe', 'ndarray', 'bin_edge']] = {
    'vector': 'dataframe',
    'matrix': 'ndarray',
    'bin_edge': 'bin_edge',
}

#: The wholes of one property as loaded from its files, by whole name.
ArrayWholes = dict[str, np.ndarray]
#: The wholes a hook may return, by whole name.
Wholes = dict[str, np.ndarray] | dict[str, pd.DataFrame]
#: Bin edges by whole name.
Edges = Mapping[str, np.ndarray]

FloatArray = npt.NDArray[np.float64]


class DensityKernel(Protocol):
    """
    The volume-fraction kernel that :class:`LocalDensities` calls.

    A kernel turns the histogram of one species in one whole into its local
    number density and local volume fraction, bin by bin, from the
    geometry the whole's parser instance describes. The study's driver
    supplies one, so this module imports no density model.

    .. versionadded:: 0.6.0
    """

    def __call__(
        self,
        counts: np.ndarray,
        edges: np.ndarray,
        info: ParserBase,
        direction: str,
        diameter: float,
    ) -> tuple[npt.ArrayLike, npt.ArrayLike]:
        """
        Parameters
        ----------
        counts : numpy.ndarray
            The histogram of one whole, shape ``(n_bins,)``: the number of
            particle centres in each bin, summed over the frames.
        edges : numpy.ndarray
            Its bin edges, shape ``(n_bins + 1,)``, in the length unit.
        info : ParserBase
            The whole's parser instance, which holds the box sizes and the
            number of frames it needs.
        direction : str
            The direction of the histogram, such as ``'r'`` or ``'z'``.
        diameter : float
            The diameter of the species, in the length unit.

        Returns
        -------
        tuple of array_like
            The number density and the volume fraction, each of shape
            ``(n_bins,)``.
        """
        ...


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Derived:
    """
    One artifact a post-hook returns, ready for the ensemble stages.

    Attributes
    ----------
    name : str
        The output property, which names the files; for example
        ``'rRhoCrd'`` or ``'gyrTMon-acf'``.
    wholes : dict of str to numpy.ndarray or pandas.DataFrame
        The artifact of each whole, by whole name.
    whole_type : {'vector', 'matrix', 'bin_edge', 'dataframe'}
        The shape of each whole; see the module notation.
    edges : mapping of str to numpy.ndarray, optional
        The bin edges of each whole, for a ``'vector'`` histogram, so its
        ensembles carry a ``bin_center`` column.

    Raises
    ------
    ValueError
        If `whole_type` is unknown, or `edges` is given for a whole type
        other than ``'vector'``.

    Examples
    --------
    >>> import numpy as np
    >>> out = Derived('gyrMon', {'w1': np.array([1.0])}, 'vector')
    >>> out.name, out.whole_type
    ('gyrMon', 'vector')

    .. versionadded:: 0.6.0
    """

    name: str
    wholes: Wholes
    whole_type: WholeTypeT
    edges: Edges | None = None

    def __post_init__(self) -> None:
        _check_whole_type(self.whole_type)
        if self.edges is not None and self.whole_type != 'vector':
            raise ValueError(
                f"'{self.name}': edges are only valid with whole_type "
                f"'vector', not '{self.whole_type}'."
            )


@dataclass(frozen=True)
class PropertySpec:
    """
    One probe artifact of a particle group, and how it reduces.

    Attributes
    ----------
    name : str
        The property's name without its species, such as ``'gyrT'`` or
        ``'rHist'``.
    species : str
        The species the property is measured on, such as ``'Mon'``; the
        empty string for a property without one, such as an edge file.
    relation : {'tseries', 'histogram', 'bin_edge'}
        How the segments of a whole combine; see the module notation.
    whole_type : {'vector', 'matrix', 'bin_edge'}
        The shape of one whole as loaded; see the module notation.
    edge_prop : str, optional
        The property of the file that holds this histogram's bin edges, such
        as ``'rEdgeMon'`` for ``'rHistMon'``. Only for a ``'vector'``
        ``'histogram'``.
    post_hook : callable, optional
        A post-hook ``(wholes, context) -> tuple[Derived, ...]``; see the
        module notation. ``None`` passes the wholes through under `prop`.

    Raises
    ------
    ValueError
        If `name` is empty; if `relation` or `whole_type` is unknown; if
        exactly one of them is ``'bin_edge'``; or if `edge_prop` is given
        for anything but a ``'vector'`` ``'histogram'``.

    Notes
    -----
    The named constructors :func:`tseries`, :func:`histogram`,
    :func:`histogram_2d`, :func:`edge`, and :func:`matrix_tseries` build the
    kinds the studies use, so a table reads as data.

    Examples
    --------
    >>> spec = PropertySpec('rHist', 'Mon', 'histogram', 'vector',
    ...                     edge_prop='rEdgeMon')
    >>> spec.prop, spec.has_edges
    ('rHistMon', True)

    An edge file is a ``'bin_edge'`` in both its relation and its shape:

    >>> PropertySpec('xEdge', '', 'bin_edge', 'vector')
    Traceback (most recent call last):
    ...
    ValueError: 'xEdge': relation 'bin_edge' and whole_type 'vector' ...

    .. versionadded:: 0.6.0
    """

    name: str
    species: EntityT
    relation: WholeRelationT
    whole_type: Literal['vector', 'matrix', 'bin_edge']
    edge_prop: PropertyT | None = None
    post_hook: 'PostHook | None' = None

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("A property's name must not be empty.")
        if self.relation not in ('tseries', 'histogram', 'bin_edge'):
            raise ValueError(
                f"'{self.prop}': relation must be 'tseries', 'histogram', or "
                f"'bin_edge', not '{self.relation}'."
            )
        if self.whole_type not in ('vector', 'matrix', 'bin_edge'):
            raise ValueError(
                f"'{self.prop}': whole_type must be 'vector', 'matrix', or "
                f"'bin_edge', not '{self.whole_type}'."
            )
        if (self.relation == 'bin_edge') != (self.whole_type == 'bin_edge'):
            raise ValueError(
                f"'{self.prop}': relation '{self.relation}' and whole_type "
                f"'{self.whole_type}' do not match; an edge file is "
                "'bin_edge' in both."
            )
        if self.edge_prop is not None and not (
            self.relation == 'histogram' and self.whole_type == 'vector'
        ):
            raise ValueError(
                f"'{self.prop}': edge_prop is only valid for a 'vector' "
                f"'histogram', not a '{self.whole_type}' '{self.relation}'."
            )

    @property
    def prop(self) -> PropertyT:
        """
        The property as it appears in file names: `name` then `species`.
        """
        return self.name + self.species

    @property
    def has_edges(self) -> bool:
        """
        Whether the property has a bin-edge file.
        """
        return self.edge_prop is not None


@dataclass(frozen=True)
class HookContext:
    """
    What a post-hook may read besides the wholes.

    Attributes
    ----------
    prop : PropertySpec
        The property the wholes belong to.
    spec : ProjectSpec, optional
        The study's spec: its parser, topology, and species sizes.
    group : str
        The particle group.
    infos : mapping of str to ParserBase
        The parser instance of each whole, by whole name, which holds the
        whole's attributes such as sizes and counts.
    edges : mapping of str to numpy.ndarray, optional
        The bin edges of each whole, for a property with an edge file.
    density : DensityKernel, optional
        The injected volume-fraction kernel.

    .. versionadded:: 0.6.0
    """

    prop: PropertySpec
    spec: ProjectSpec | None = None
    group: GroupT = 'bug'
    infos: Mapping[str, ParserBase] = field(default_factory=dict)
    edges: Edges | None = None
    density: DensityKernel | None = None


#: A post-hook; see the module notation.
PostHook = Callable[[ArrayWholes, HookContext], tuple[Derived, ...]]


@dataclass(frozen=True)
class AnalysisTable:
    """
    The analysis of one particle group of a study, as data.

    Attributes
    ----------
    group : str
        The particle group whose artifacts the table reduces.
    lineage : {'segment', 'whole'}
        The lineage of the group's probe artifacts; see the module notation.
    stamps : bool
        Whether the group's stamp files are combined by
        :func:`reduce_stamps`.
    properties : tuple of PropertySpec
        The artifacts to reduce, each listed once.

    Raises
    ------
    ValueError
        If `lineage` is unknown, or two properties name the same artifact.

    Examples
    --------
    >>> table = AnalysisTable('bug', 'whole', True, (tseries('gyrT', 'Mon'),))
    >>> table.lineage, len(table.properties)
    ('whole', 1)

    Listing an artifact twice is an error, so a table cannot reduce one file
    two ways by accident:

    >>> AnalysisTable('bug', 'whole', True,
    ...               (tseries('gyrT', 'Mon'), tseries('gyrT', 'Mon')))
    Traceback (most recent call last):
    ...
    ValueError: Properties listed more than once: ['gyrTMon'].

    .. versionadded:: 0.6.0
    """

    group: GroupT
    lineage: PrimitiveLineageT
    stamps: bool
    properties: tuple[PropertySpec, ...]

    def __post_init__(self) -> None:
        if self.lineage not in ('segment', 'whole'):
            raise ValueError(
                f"lineage must be 'segment' or 'whole', not "
                f"'{self.lineage}'."
            )
        object.__setattr__(self, 'properties', tuple(self.properties))
        seen: set[str] = set()
        repeated: list[str] = []
        for prop in self.properties:
            if prop.prop in seen and prop.prop not in repeated:
                repeated.append(prop.prop)
            seen.add(prop.prop)
        if repeated:
            raise ValueError(f"Properties listed more than once: {repeated}.")


def tseries(
    name: str, species: EntityT, post_hook: 'PostHook | None' = None
) -> PropertySpec:
    """
    A scalar time series: one value per frame, concatenated over segments.

    Parameters
    ----------
    name : str
        The property's name, such as ``'gyrT'``.
    species : str
        The species, such as ``'Mon'``.
    post_hook : callable, optional
        A post-hook; see :class:`PropertySpec`.

    Returns
    -------
    PropertySpec
        A ``'tseries'`` ``'vector'`` property.

    Examples
    --------
    >>> tseries('gyrT', 'Mon').prop
    'gyrTMon'

    .. versionadded:: 0.6.0
    """
    return PropertySpec(name, species, 'tseries', 'vector',
                        post_hook=post_hook)


def histogram(
    direction: DirectionT, species: EntityT,
    post_hook: 'PostHook | None' = None,
) -> PropertySpec:
    """
    A one-dimensional histogram with its bin-edge file.

    Parameters
    ----------
    direction : str
        The direction, such as ``'r'``, ``'z'``, or ``'theta'``: the property
        is ``<direction>Hist<species>`` and its edges
        ``<direction>Edge<species>``.
    species : str
        The species, such as ``'Crd'``.
    post_hook : callable, optional
        A post-hook; see :class:`PropertySpec`.

    Returns
    -------
    PropertySpec
        A ``'histogram'`` ``'vector'`` property with `edge_prop` set.

    Examples
    --------
    >>> spec = histogram('r', 'Crd')
    >>> spec.prop, spec.edge_prop
    ('rHistCrd', 'rEdgeCrd')

    .. versionadded:: 0.6.0
    """
    return PropertySpec(f'{direction}Hist', species, 'histogram', 'vector',
                        edge_prop=f'{direction}Edge{species}',
                        post_hook=post_hook)


def histogram_2d(plane: str, species: EntityT) -> PropertySpec:
    """
    A two-dimensional histogram, summed over segments and stacked over wholes.

    Parameters
    ----------
    plane : str
        The plane as the probe phase names it, such as ``'xy'`` or
        ``'xz'``: the property is ``<plane>Hist<species>``. Its edges are
        the :func:`edge` properties of the plane's two directions.
    species : str
        The species, such as ``'Crd'``.

    Returns
    -------
    PropertySpec
        A ``'histogram'`` ``'matrix'`` property.

    Examples
    --------
    >>> histogram_2d('xy', 'Crd').prop
    'xyHistCrd'

    .. versionadded:: 0.6.0
    """
    return PropertySpec(f'{plane}Hist', species, 'histogram', 'matrix')


def edge(direction: DirectionT) -> PropertySpec:
    """
    A bin-edge file shared by the two-dimensional histograms of a direction.

    Parameters
    ----------
    direction : str
        The direction, such as ``'x'``: the property is
        ``<direction>Edge``, with no species.

    Returns
    -------
    PropertySpec
        A ``'bin_edge'`` property.

    Examples
    --------
    >>> edge('x').prop
    'xEdge'

    .. versionadded:: 0.6.0
    """
    return PropertySpec(f'{direction}Edge', '', 'bin_edge', 'bin_edge')


def matrix_tseries(
    name: str, species: EntityT, post_hook: 'PostHook | None' = None
) -> PropertySpec:
    """
    A time series of arrays: one array per frame, concatenated over segments.

    Parameters
    ----------
    name : str
        The property's name, such as ``'principalT'`` (one matrix per frame)
        or ``'clustersHistT'`` (one histogram per frame).
    species : str
        The species, such as ``'Mon'``.
    post_hook : callable, optional
        A post-hook; see :class:`PropertySpec`.

    Returns
    -------
    PropertySpec
        A ``'tseries'`` ``'matrix'`` property.

    Examples
    --------
    >>> matrix_tseries('principalT', 'Mon').whole_type
    'matrix'

    .. versionadded:: 0.6.0
    """
    return PropertySpec(name, species, 'tseries', 'matrix',
                        post_hook=post_hook)


# ---------------------------------------------------------------------------
# Hooks
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TimeAverage:
    """
    Average a time series over its frames; a post-hook.

    Attributes
    ----------
    output : str
        The name of the average, without the species, which is appended:
        ``TimeAverage('clustersHist')`` on ``'clustersHistTFoci'`` writes
        ``'clustersHistFoci'``.
    drop_first_column : bool, default False
        Drop column 0 of each frame before averaging. The size histograms of
        :func:`~polyphys.analyze.clusters.count_clusters` carry the count of
        empty clusters, always 0, in column 0.

    Notes
    -----
    The frames are axis 0; the average has the shape of one frame, so a
    ``(n_frames, n_sizes)`` histogram time series becomes an ``n_sizes``
    vector, reduced as a ``'vector'`` whole. Time and space are
    :math:`O(n_\\mathrm{frames} n_\\mathrm{sizes})` per whole.

    Examples
    --------
    >>> import numpy as np
    >>> hook = TimeAverage('clustersHist', drop_first_column=True)
    >>> context = HookContext(matrix_tseries('clustersHistT', 'Foci'))
    >>> wholes = {'w1': np.array([[0, 1, 2], [0, 3, 4], [0, 5, 6]])}
    >>> (out,) = hook(wholes, context)
    >>> out.name, out.wholes['w1']
    ('clustersHistFoci', array([3., 4.]))

    .. versionadded:: 0.6.0
    """

    output: str
    drop_first_column: bool = False

    def __call__(
        self, wholes: ArrayWholes, context: HookContext
    ) -> tuple[Derived, ...]:
        averaged: dict[str, np.ndarray] = {}
        for whole_name, whole in wholes.items():
            frames = np.asarray(whole, dtype=np.float64)
            if self.drop_first_column:
                frames = frames[:, 1:]
            averaged[whole_name] = frames.mean(axis=0)
        return (Derived(self.output + context.prop.species, averaged,
                        'vector'),)


@dataclass(frozen=True)
class Autocorrelation:
    """
    Pass a time series through and add its autocorrelation; a post-hook.

    Attributes
    ----------
    nlags : int
        The largest lag, in frames; see
        :func:`~polyphys.analyze.correlations.acf`.

    Notes
    -----
    Returns the series under its own name and :math:`r_k` for
    :math:`k = 0, \\dots,` `nlags` under ``<property>-acf``, both
    ``'vector'`` wholes. No confidence band is stored: Bartlett's band is a
    function of the stored :math:`r_k`, the number of frames, and the
    level, so :func:`~polyphys.analyze.correlations.bartlett_half_width`
    rebuilds it at any level when it is drawn. Time is
    :math:`O(n_\\mathrm{frames} \\log n_\\mathrm{frames})` per whole.

    Examples
    --------
    The alternating series has :math:`r_1 = -7/8` and :math:`r_2 = 6/8`:

    >>> import numpy as np
    >>> hook = Autocorrelation(nlags=2)
    >>> context = HookContext(tseries('gyrT', 'Mon'))
    >>> series, acf_out = hook({'w1': np.tile([1.0, -1.0], 4)}, context)
    >>> series.name, acf_out.name
    ('gyrTMon', 'gyrTMon-acf')
    >>> np.round(acf_out.wholes['w1'], 6)
    array([ 1.   , -0.875,  0.75 ])

    .. versionadded:: 0.6.0
    """

    nlags: int

    def __call__(
        self, wholes: ArrayWholes, context: HookContext
    ) -> tuple[Derived, ...]:
        prop = context.prop.prop
        acfs: dict[str, np.ndarray] = {
            whole_name: acf(whole, self.nlags, alpha=None).values
            for whole_name, whole in wholes.items()
        }
        return (Derived(prop, wholes, 'vector'),
                Derived(f'{prop}-acf', acfs, 'vector'))


@dataclass(frozen=True)
class LocalDensities:
    """
    Pass a histogram through and add its density profiles; a post-hook.

    On ``<direction>Hist<species>`` it adds ``<direction>Rho<species>``, the
    local number density, and ``<direction>Phi<species>``, the local volume
    fraction, from the :class:`DensityKernel` given to :func:`reduce_space`,
    with the species' diameter read from the attribute the study's
    :attr:`~polyphys.manage.projects.ProjectSpec.species_sizes` names.

    Raises
    ------
    ValueError
        When called: if the property is not a :func:`histogram`, if the
        context has no kernel, no spec, or no edges, or if the species is not
        in the spec's `species_sizes`.

    Notes
    -----
    The three outputs share the histogram's edges, so their ensembles carry
    the same ``bin_center`` column. Time is one kernel call per whole.

    Examples
    --------
    With a stub kernel that returns the counts and a constant, on one whole
    of the ``SumRuleCyl`` study:

    >>> import numpy as np
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> from polyphys.manage.projects import ProjectSpec
    >>> spec = ProjectSpec(
    ...     parser=SumRuleCyl, group='bug', divisor=0.025,
    ...     physical_attrs=('space',), species_sizes={'Crd': 'dcrowd'},
    ...     tseries_properties=())
    >>> whole = ('N20epsilon5.0r5.5lz25.0sig1.0nc100'
    ...          'dt0.005bdump1000adump5000ens1')
    >>> def kernel(counts, edges, info, direction, size):
    ...     return counts, np.full(counts.shape, size)
    >>> context = HookContext(
    ...     histogram('r', 'Crd'), spec=spec, group='all',
    ...     infos={whole: SumRuleCyl(whole, 'whole', 'all')},
    ...     edges={whole: np.array([0.0, 1.0, 2.0])}, density=kernel)
    >>> outputs = LocalDensities()({whole: np.array([3.0, 5.0])}, context)
    >>> [out.name for out in outputs]
    ['rHistCrd', 'rRhoCrd', 'rPhiCrd']
    >>> outputs[2].wholes[whole]
    array([1., 1.])

    .. versionadded:: 0.6.0
    """

    def __call__(
        self, wholes: ArrayWholes, context: HookContext
    ) -> tuple[Derived, ...]:
        prop = context.prop
        if not prop.name.endswith('Hist') or not prop.has_edges:
            raise ValueError(
                f"LocalDensities needs a histogram(...) property, not "
                f"'{prop.prop}'."
            )
        if context.density is None:
            raise ValueError(
                f"'{prop.prop}': LocalDensities needs the density kernel "
                "passed to reduce_space(..., density=...)."
            )
        if context.spec is None or context.edges is None:
            raise ValueError(
                f"'{prop.prop}': LocalDensities needs the study's spec and "
                "the histogram's edges in its context."
            )
        size_attr = context.spec.species_sizes.get(prop.species)
        if size_attr is None:
            raise ValueError(
                f"'{prop.prop}': species '{prop.species}' has no size in "
                f"{context.spec.name}.species_sizes "
                f"{dict(context.spec.species_sizes)}."
            )
        direction = prop.name[: -len('Hist')]
        rhos: dict[str, np.ndarray] = {}
        phis: dict[str, np.ndarray] = {}
        for whole_name, counts in wholes.items():
            info = context.infos[whole_name]
            rho, phi = context.density(
                np.asarray(counts), np.asarray(context.edges[whole_name]),
                info, direction, float(getattr(info, size_attr)),
            )
            rhos[whole_name] = np.asarray(rho, dtype=np.float64)
            phis[whole_name] = np.asarray(phi, dtype=np.float64)
        return (
            Derived(prop.prop, wholes, 'vector', context.edges),
            Derived(f'{direction}Rho{prop.species}', rhos, 'vector',
                    context.edges),
            Derived(f'{direction}Phi{prop.species}', phis, 'vector',
                    context.edges),
        )


#: A number, or a function of a whole's parser instance that gives one.
PerWhole = float | Callable[[ParserBase], float]


def _resolve(value: PerWhole, info: ParserBase | None, what: str) -> float:
    """Return `value`, calling it on `info` when it is a function."""
    if callable(value):
        if info is None:
            raise ValueError(
                f"{what} is a function of the whole, but the context has no "
                "parser instance for it."
            )
        return float(value(info))
    return float(value)


def half_contour(n_monomers: int, bond_max: float) -> float:
    """
    The largest separation two monomers of a ring can have.

    Parameters
    ----------
    n_monomers : int
        The number of monomers on the ring, at least 2.
    bond_max : float
        The largest length a bond can have, in the length unit; for the FENE
        bond, :math:`R_0 = 1.5\\,\\sigma_{ij}` with :math:`\\sigma_{ij}` the
        collision diameter of the bonded pair (Kremer and Grest, 1990).

    Returns
    -------
    float
        ``(n_monomers // 2) * bond_max``: the two monomers are at most
        ``n_monomers // 2`` bonds apart along the ring, and a path of bonds
        is longest when straight and fully extended.

    Raises
    ------
    ValueError
        If `n_monomers` is below 2 or `bond_max` is not positive.

    Notes
    -----
    A bound, not a typical value: it lets a pair-distance histogram cover
    every separation the chain allows, so that its counts conserve the
    number of frames and a separation above it marks a broken bond or an
    unwrapping error. Pass the largest bond's :math:`R_0` when the bonds
    differ; the bound is then loose by the ratio of the largest to the
    typical bond, which costs bins, not correctness.

    Examples
    --------
    A ring of 10 beads with bonds that stretch to 1.5:

    >>> half_contour(10, 1.5)
    7.5

    .. versionadded:: 0.6.0
    """
    if n_monomers < 2:
        raise ValueError(f'n_monomers must be at least 2; got {n_monomers}.')
    if not bond_max > 0:
        raise ValueError(f'bond_max must be positive; got {bond_max}.')
    return (n_monomers // 2) * bond_max


@dataclass(frozen=True)
class FociPairDistances:
    """
    Reduce the distance-matrix time series of the foci; a post-hook.

    On ``distMatT<species>``, whose whole is the ``(n_frames, n_foci,
    n_foci)`` matrix of centre-to-centre distances between the foci with the
    foci's backbone indices on the diagonal of frame 0, it writes three
    ``'dataframe'`` wholes with one column per pair, tagged
    ``loc<i>loc<j>genDist<d>`` with ``i < j`` the backbone indices and ``d``
    their genomic distance in bonds:

    - ``pairDistHist<species>``: the histogram of each pair's separation on
      one grid shared by every pair and whole, with a ``bin_center`` column;
    - ``pairDistRdf<species>``: :math:`\\bar P(R)`, the probability density
      of each pair's separation vector on that grid, with ``bin_center``;
    - ``pairDistT<species>``: each pair's separation over the frames.

    Attributes
    ----------
    bin_size : float or callable
        The bin width :math:`\\Delta r`, in the length unit, or a function of
        the whole's parser instance that gives it.
    r_max : float or callable
        The upper end of the grid, or a function of the parser instance that
        gives it, such as :func:`half_contour` of the ring; the grid runs to
        the first multiple of `bin_size` at or above it.
    n_monomers : int or callable
        The number of monomers on the backbone, or a function of the parser
        instance that gives it.

    Raises
    ------
    ValueError
        When called: if the context has no spec (the topology comes from
        it); if a whole is not a ``(n_frames, n_foci, n_foci)`` array; or if
        a separation lies above the grid, which names the whole and the
        separation, since a dropped frame would break the normalization.

    Notes
    -----
    Every pair has one separation per frame, so each histogram column sums
    to ``n_frames`` and :math:`\\sum_k \\bar P_k \\Delta V_k = 1` exactly.
    The research code set the grid from each whole's own largest separation,
    so wholes of one ensemble had different edges; a bound from the system
    gives one grid. Time is :math:`O(n_\\mathrm{pairs} n_\\mathrm{frames}
    \\log n_\\mathrm{bins})` and space :math:`O(n_\\mathrm{pairs}
    (n_\\mathrm{frames} + n_\\mathrm{bins}))` per whole.

    Examples
    --------
    Three foci at backbone indices 0, 4, and 8 of a ring of 10 monomers,
    over two frames, on a grid of unit bins to 3:

    >>> import numpy as np
    >>> from polyphys.manage.parser import SumRuleCubHeteroRing
    >>> from polyphys.manage.projects import ProjectSpec
    >>> spec = ProjectSpec(
    ...     parser=SumRuleCubHeteroRing, group='bug', divisor=0.025,
    ...     physical_attrs=('space',), species_sizes={},
    ...     tseries_properties=())
    >>> frame0 = np.array([[0.0, 1.5, 2.5], [1.5, 4.0, 0.5], [2.5, 0.5, 8.0]])
    >>> frame1 = np.array([[0.0, 1.5, 2.9], [1.5, 0.0, 2.5], [2.9, 2.5, 0.0]])
    >>> hook = FociPairDistances(bin_size=1.0, r_max=3.0, n_monomers=10)
    >>> context = HookContext(matrix_tseries('distMatT', 'Foci'), spec=spec)
    >>> hist, density, series = hook({'w1': np.stack([frame0, frame1])},
    ...                              context)
    >>> list(hist.wholes['w1'].columns)
    ['loc0loc4genDist4', 'loc0loc8genDist2', 'loc4loc8genDist4', 'bin_center']
    >>> hist.wholes['w1']['loc0loc8genDist2'].to_list()
    [0, 0, 2]
    >>> series.wholes['w1']['loc4loc8genDist4'].to_list()
    [0.5, 2.5]

    .. versionadded:: 0.6.0
    """

    bin_size: PerWhole
    r_max: PerWhole
    n_monomers: int | Callable[[ParserBase], int]

    def __call__(
        self, wholes: ArrayWholes, context: HookContext
    ) -> tuple[Derived, ...]:
        prop = context.prop
        if context.spec is None:
            raise ValueError(
                f"'{prop.prop}': FociPairDistances needs the study's spec in "
                "its context, for the topology."
            )
        # foci_pair_info checks the topology itself.
        topology = cast(clusters.Topology, context.spec.topology)
        hists: dict[str, pd.DataFrame] = {}
        densities: dict[str, pd.DataFrame] = {}
        series: dict[str, pd.DataFrame] = {}
        for whole_name, whole in wholes.items():
            matrix = np.asarray(whole, dtype=np.float64)
            if matrix.ndim != 3 or matrix.shape[1] != matrix.shape[2]:
                raise ValueError(
                    f"'{prop.prop}', whole '{whole_name}': expected a "
                    "(n_frames, n_foci, n_foci) distance matrix time series,"
                    f" got shape {matrix.shape}."
                )
            info = context.infos.get(whole_name)
            n_foci = matrix.shape[1]
            genomic_pos = np.rint(np.diagonal(matrix[0])).astype(np.int64)
            n_monomers = int(_resolve(self.n_monomers, info, 'n_monomers'))
            pairs = clusters.foci_pair_info(genomic_pos, n_monomers, topology)
            tags = [f'loc{i}loc{j}genDist{d}' for i, j, d in pairs]
            first, second = np.triu_indices(n_foci, 1)
            distances = matrix[:, first, second].T
            bin_size = _resolve(self.bin_size, info, 'bin_size')
            r_max = _resolve(self.r_max, info, 'r_max')
            upper = np.ceil(r_max / bin_size) * bin_size
            largest = float(distances.max())
            if largest > upper:
                raise ValueError(
                    f"'{prop.prop}', whole '{whole_name}': a pair separation "
                    f"of {largest:.6g} lies above r_max {upper:.6g}; the "
                    "chain cannot reach it, so check the whole."
                )
            counts, edges = clusters.pair_distance_histograms(
                distances, bin_size, r_max)
            counts = counts.reshape(len(tags), -1)
            density = clusters.pair_separation_density(counts, edges)
            centers = 0.5 * (edges[:-1] + edges[1:])
            hists[whole_name] = pd.DataFrame(
                dict(zip(tags, counts)) | {'bin_center': centers})
            densities[whole_name] = pd.DataFrame(
                dict(zip(tags, density)) | {'bin_center': centers})
            series[whole_name] = pd.DataFrame(dict(zip(tags, distances)))
        species = prop.species
        return (
            Derived(f'pairDistHist{species}', hists, 'dataframe'),
            Derived(f'pairDistRdf{species}', densities, 'dataframe'),
            Derived(f'pairDistT{species}', series, 'dataframe'),
        )


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------


SaveDirs = tuple[Pathish | None, Pathish | None, Pathish | None]


def reduce_stamps(
    stamp_files: Sequence[tuple[str]],
    lineage: PrimitiveLineageT,
    save_to: SaveDirs = (None, None, None),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Combine the stamps of one space into whole and ensemble stamps.

    Parameters
    ----------
    stamp_files : sequence of tuple of str
        The stamp CSV files, each path in a 1-tuple, as
        :func:`~polyphys.manage.utils.sort_filenames` returns them.
    lineage : {'segment', 'whole'}
        The lineage of the files.
    save_to : tuple of three paths or None, default (None, None, None)
        The ``wholeSim``, ``ens``, and ``ensAvg`` directories; ``None`` skips
        a stage's file. For ``'segment'`` stamps the combined segment
        stamps go to the first, the whole stamps to the second, and the
        ensemble stamps to the third; for ``'whole'`` stamps the whole
        stamps go to the second and the ensemble stamps to the third.

    Returns
    -------
    whole_stamps : pandas.DataFrame
        One row per whole.
    ensemble_stamps : pandas.DataFrame
        One row per ensemble.

    Raises
    ------
    ValueError
        As :func:`~polyphys.manage.organizer.children_stamps` and
        :func:`~polyphys.manage.organizer.parents_stamps`.

    Notes
    -----
    Segment stamps combine with
    :func:`~polyphys.manage.organizer.children_stamps`, then roll up once
    with :func:`~polyphys.manage.organizer.parents_stamps` into whole
    stamps; whole stamps combine directly. Either way the whole stamps roll
    up once more into ensemble stamps. Time and space are
    :math:`O(n_\\mathrm{files}\\, n_\\mathrm{columns})`.

    Examples
    --------
    .. code-block:: python

        whole_stamps, ensemble_stamps = reduce_stamps(
            [(path,) for path in sorted(stamp_paths)], 'whole')

    .. versionadded:: 0.6.0
    """
    whole_dir, ens_dir, ens_avg_dir = save_to
    if lineage == 'segment':
        segment_stamps = organizer.children_stamps(
            list(stamp_files), 'segment', save_to=whole_dir)
        whole_stamps = organizer.parents_stamps(
            segment_stamps, 'segment', save_to=ens_dir)
    else:
        whole_stamps = organizer.children_stamps(
            list(stamp_files), 'whole', save_to=ens_dir)
    ensemble_stamps = organizer.parents_stamps(
        whole_stamps, 'whole', save_to=ens_avg_dir)
    return whole_stamps, ensemble_stamps


def reduce_space(
    input_database: Pathish,
    spec: ProjectSpec,
    table: AnalysisTable,
    *,
    hierarchy: str | None = None,
    density: DensityKernel | None = None,
    save: bool = True,
) -> dict[str, dict[str, np.ndarray | pd.DataFrame]]:
    """
    Reduce every property of a table over one space of probe artifacts.

    Parameters
    ----------
    input_database : path-like
        The probe directory of the space,
        ``<root>/probes/<prefix>-<group>-probe``; it may hold a glob, as the
        research code's drivers pass it.
    spec : ProjectSpec
        The study's spec: its parser, topology, and species sizes.
    table : AnalysisTable
        The group's analysis table.
    hierarchy : str, optional
        The glob, relative to `input_database`, that lists the artifact
        files; the default is ``<spec.hierarchy>/<spec.hierarchy>``, the
        research code's ``dir/file`` layout.
    density : DensityKernel, optional
        The volume-fraction kernel for :class:`LocalDensities`.
    save : bool, default True
        Write the stages' files under ``<root>/analysis/``; ``False``
        returns the ensemble averages without writing anything.

    Returns
    -------
    dict of str to dict
        The ensemble averages by output property, then by ensemble name:
        a DataFrame for a ``'vector'`` or ``'dataframe'`` output, a dict of
        arrays for a ``'matrix'``, and the edges for a ``'bin_edge'``.

    Raises
    ------
    FileNotFoundError
        If the glob matches no file, or a property matches no file.
    ValueError
        If the wholes of one ensemble hold different bin edges, or if a
        histogram's edge files do not cover its wholes; as the organizer and
        the hooks raise.

    Notes
    -----
    For each property in turn: its files are selected by their last name
    token, wholes are built (for the ``'segment'`` lineage with
    :func:`~polyphys.manage.organizer.whole_from_segments`, written to
    ``wholeSim``; else with
    :func:`~polyphys.manage.organizer.whole_from_file`), the edges are
    loaded and checked, the post-hook is applied, and every returned
    artifact goes through :func:`~polyphys.manage.organizer.ensemble` and
    :func:`~polyphys.manage.organizer.ensemble_avg`, or their long-form
    counterparts for ``'dataframe'`` wholes. With `table.stamps`, the
    stamps are combined first by :func:`reduce_stamps`. One property's
    wholes are in memory at a time, so peak memory is the largest property's
    ``n_wholes`` wholes plus its ensembles; time is
    :math:`O(n_\\mathrm{files}\\, n_\\mathrm{frames})` per property plus the
    hook's.

    Examples
    --------
    Two replicates of one state point of the ``SumRuleCyl`` study, written
    as wholes and reduced without saving:

    >>> import tempfile
    >>> import numpy as np
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> from polyphys.manage.projects import ProjectSpec
    >>> spec = ProjectSpec(
    ...     parser=SumRuleCyl, group='bug', divisor=0.025,
    ...     physical_attrs=('space',), species_sizes={'Mon': 'dmon'},
    ...     tseries_properties=('gyrTMon',))
    >>> stem = 'N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000'
    >>> root = tempfile.mkdtemp()
    >>> np.save(f'{root}/{stem}ens1-bug-gyrTMon.npy', np.array([1.0, 2.0]))
    >>> np.save(f'{root}/{stem}ens2-bug-gyrTMon.npy', np.array([3.0, 6.0]))
    >>> table = AnalysisTable('bug', 'whole', False, (tseries('gyrT', 'Mon'),))
    >>> averages = reduce_space(root, spec, table, hierarchy='*', save=False)
    >>> averages['gyrTMon'][stem]
       gyrTMon-mean  gyrTMon-var  gyrTMon-sem
    0           2.0          2.0          1.0
    1           4.0          8.0          2.0

    .. versionadded:: 0.6.0
    """
    parser = spec.parser
    group = table.group
    if hierarchy is None:
        hierarchy = f'{spec.hierarchy}/{spec.hierarchy}'
    pattern = os.path.join(os.fspath(input_database), hierarchy)
    observations = sorted(glob.glob(pattern))
    if not observations:
        raise FileNotFoundError(f"No file matches '{pattern}'.")
    logger.info("Reducing %d files of group '%s' under '%s'.",
                len(observations), group, pattern)

    whole_dir: str | None = None
    ens_dir: str | None = None
    ens_avg_dir: str | None = None
    if save:
        if table.lineage == 'segment':
            whole_dir = organizer.make_database(
                input_database, 'analysis', 'wholeSim', group)
        ens_dir = organizer.make_database(
            input_database, 'analysis', 'ens', group)
        ens_avg_dir = organizer.make_database(
            input_database, 'analysis', 'ensAvg', group)

    if table.stamps:
        stamp_files = _select(observations, 'stamps', ext='csv')
        reduce_stamps(stamp_files, table.lineage,
                      save_to=(whole_dir, ens_dir, ens_avg_dir))

    results: dict[str, dict[str, np.ndarray | pd.DataFrame]] = {}
    for prop in table.properties:
        logger.info("Reducing property '%s'.", prop.prop)
        wholes = _load_wholes(prop.prop, prop.relation, observations, parser,
                              group, table.lineage, whole_dir)
        edges: dict[str, np.ndarray] | None = None
        if prop.edge_prop is not None:
            edges = _load_wholes(prop.edge_prop, 'bin_edge', observations,
                                 parser, group, table.lineage, whole_dir)
            _check_edges(prop, wholes, edges, parser, group)
        infos = {name: parser(name, 'whole', group) for name in wholes}
        context = HookContext(prop, spec=spec, group=group, infos=infos,
                              edges=edges, density=density)
        if prop.post_hook is None:
            outputs: tuple[Derived, ...] = (
                Derived(prop.prop, wholes, prop.whole_type, edges),)
        else:
            outputs = prop.post_hook(wholes, context)
        for output in outputs:
            if whole_dir is not None and output.name != prop.prop:
                for whole_name, data in output.wholes.items():
                    organizer.save_artifact(
                        organizer.create_fullname(whole_name, group,
                                                  output.name),
                        data, whole_dir)
            results[output.name] = _reduce_output(
                output, parser, group, ens_dir, ens_avg_dir)
    return results


def equilibrium_wholes(
    ens_database: Pathish,
    properties: Mapping[PropertyT, str],
    whole_stamps: pd.DataFrame,
    *,
    method: PlateauMethod = 'lee',
    topology: str | None = None,
) -> pd.DataFrame:
    """
    The equilibrium value of each whole's time series, with its error,
    beside the whole's stamps.

    For each property, every ensemble file of the ``ens`` stage is read
    once, and each of its columns, one whole's time series, is reduced to
    its time average, its variance, and the blocked standard error of the
    average. The results of all the properties are merged with
    `whole_stamps` on ``whole``, one row per whole: the table the
    research code's ``equilibrium_wholes`` built for the ``allInOne``
    phase.

    Parameters
    ----------
    ens_database : path-like
        The ``ens`` directory of one space and group,
        ``<root>/analysis/<prefix>-<group>-ens``, whose
        ``<ensemble_long>-<group>-<property>.csv`` files
        :func:`reduce_space` wrote.
    properties : mapping of str to str
        The time-series property of each ensemble file, mapped to its
        output name: ``{'gyrTMon': 'gyrMon'}`` reads the ``-gyrTMon.csv``
        files and writes ``gyrMon-mean``, ``gyrMon-var``, and
        ``gyrMon-sem``.
    whole_stamps : pandas.DataFrame
        The stamps of the space's wholes, one row per whole, with a
        ``whole`` column: the first result of :func:`reduce_stamps`.
    method : {'lee', 'sem_err'}, default 'lee'
        The plateau rule of :func:`~polyphys.analyze.stats.sem_blocked`;
        ``'lee'`` is the published rule.
    topology : str, optional
        Appended to each whole's name as ``.<topology>`` before the merge,
        for stamps that the probe phase keyed by the whole's file name,
        such as ``...ens1.ring``; the ensemble files name a whole without
        it.

    Returns
    -------
    pandas.DataFrame
        `whole_stamps` followed by, for each property in the order of
        `properties`, the columns ``<output>-mean``, ``<output>-var``, and
        ``<output>-sem``, with the rows in the order of `whole_stamps`.

    Raises
    ------
    FileNotFoundError
        If `ens_database` holds no CSV file, or a property matches none.
    ValueError
        If `whole_stamps` has no ``whole`` column or lists a whole twice;
        if two properties share an output name; if the wholes of a
        property's files differ from those of the stamps, which names the
        wholes on each side; or as
        :func:`~polyphys.analyze.stats.block_average` raises for a series
        that is too short or holds NaN.

    Warns
    -----
    RuntimeWarning
        As :func:`~polyphys.analyze.stats.select_plateau`: when a whole's
        series is too short for its correlation time, its ``-sem`` is NaN
        or unreliable.

    Notes
    -----
    ``-mean`` is the time average over the frames, ``-var`` the sample
    variance of the frames (``ddof = 1``), and ``-sem`` the standard error
    of the time average by block averaging, which accounts for the
    correlation between consecutive frames of one run (Flyvbjerg and
    Petersen, 1989; the plateau rule of Lee et al., 2011; both in
    :mod:`polyphys.analyze.stats`). A constant series has a ``-sem`` of 0.
    Blocking applies within a whole only: the wholes of an ensemble are
    independent runs, so :func:`~polyphys.manage.organizer.ensemble_avg`
    keeps the plain standard error across them. The research code took
    ``-sem`` as the per-frame standard error, which the blocked one
    exceeds by about :math:`\\sqrt{s}` for a statistical inefficiency
    :math:`s`, and its ``-var`` as the population variance.

    With :math:`n_\\mathrm{frames}` the frames of the longest series and
    :math:`n_\\mathrm{wholes}` the wholes of the space, time is
    :math:`O(n_\\mathrm{wholes}\\, n_\\mathrm{frames})` per property, and
    space is one ensemble file at a time plus a ``float64`` copy of it.

    Examples
    --------
    Two wholes of one ensemble whose series alternate between two values,
    so their time averages are exact and the blocked error vanishes at a
    block size of two frames:

    >>> import tempfile
    >>> import numpy as np
    >>> import pandas as pd
    >>> stem = 'N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000'
    >>> wholes = [f'{stem}ens1', f'{stem}ens2']
    >>> root = tempfile.mkdtemp()
    >>> pd.DataFrame({wholes[0]: np.tile([0.0, 2.0], 2048),
    ...               wholes[1]: np.tile([1.0, 3.0], 2048)}
    ...              ).to_csv(f'{root}/{stem}-bug-gyrTMon.csv', index=False)
    >>> stamps = pd.DataFrame({'whole': wholes, 'n_frames': [4096, 4096]})
    >>> table = equilibrium_wholes(root, {'gyrTMon': 'gyrMon'}, stamps)
    >>> table.drop(columns='whole').round(3)
       n_frames  gyrMon-mean  gyrMon-var  gyrMon-sem
    0      4096          1.0         1.0         0.0
    1      4096          2.0         1.0         0.0

    .. versionadded:: 0.6.0
    """
    if 'whole' not in whole_stamps.columns:
        raise ValueError("whole_stamps needs a 'whole' column.")
    expected = whole_stamps['whole'].astype(str).to_list()
    repeated_wholes = sorted({w for w in expected if expected.count(w) > 1})
    if repeated_wholes:
        raise ValueError(
            f"whole_stamps lists wholes more than once: {repeated_wholes}."
        )
    outputs = list(properties.values())
    repeated = sorted({out for out in outputs if outputs.count(out) > 1})
    if repeated:
        raise ValueError(f"Output names given more than once: {repeated}.")
    pattern = os.path.join(os.fspath(ens_database), '*.csv')
    observations = sorted(glob.glob(pattern))
    if not observations:
        raise FileNotFoundError(f"No file matches '{pattern}'.")
    suffix = '' if topology is None else f'.{topology}'

    frames: list[pd.DataFrame] = []
    for prop, output in properties.items():
        logger.info("Reducing the time series of '%s' to '%s'.", prop,
                    output)
        parts: list[pd.DataFrame] = []
        for (path,) in _select(observations, prop, ext='csv'):
            ensemble = pd.read_csv(path)
            values = ensemble.to_numpy(dtype=np.float64)
            parts.append(pd.DataFrame({
                f'{output}-mean': values.mean(axis=0),
                f'{output}-var': values.var(axis=0, ddof=1),
                f'{output}-sem': _blocked_errors(values, method),
            }, index=[f'{name}{suffix}' for name in ensemble.columns]))
        measured = pd.concat(parts)
        _check_wholes(prop, list(measured.index), expected)
        frames.append(measured.loc[expected])
    combined = pd.concat(frames, axis=1)
    combined.index = whole_stamps.index
    return pd.concat([whole_stamps, combined], axis=1)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _check_whole_type(whole_type: str) -> None:
    """Raise ``ValueError`` unless `whole_type` is a known whole type."""
    if whole_type not in ('vector', 'matrix', 'bin_edge', 'dataframe'):
        raise ValueError(
            "whole_type must be 'vector', 'matrix', 'bin_edge', or "
            f"'dataframe', not '{whole_type}'."
        )


def _natural_key(path: str) -> list[str | int]:
    """
    The sort key that orders a file name's digit runs as numbers, so
    ``j2`` comes before ``j10``.
    """
    return [int(token) if token.isdigit() else token
            for token in re.split(r'(\d+)', os.path.basename(path))]


def _select(
    observations: Sequence[str], prop: str, ext: str = 'npy'
) -> list[tuple[str]]:
    """
    The files of property `prop` among `observations`, as 1-tuples in
    natural order of their names.

    A file belongs to `prop` when its name ends in ``-<prop>.<ext>`` (a
    whole) or ``.<prop>.<ext>`` (a segment), so ``gyrTMon`` never matches
    ``xgyrTMon``. Raises ``FileNotFoundError`` when nothing matches.
    :func:`~polyphys.manage.utils.sort_filenames` is not used because it
    reads a name's stem up to its first dot, which a keyword value such as
    ``dt0.005`` holds.
    """
    suffixes = (f'-{prop}.{ext}', f'.{prop}.{ext}')
    matches = [path for path in observations if path.endswith(suffixes)]
    if not matches:
        raise FileNotFoundError(
            f"No file of property '{prop}': none of the {len(observations)}"
            f" files ends in '{suffixes[0]}' or '{suffixes[1]}'."
        )
    return [(path,) for path in sorted(matches, key=_natural_key)]


def _load_wholes(
    prop: str,
    relation: WholeRelationT,
    observations: Sequence[str],
    parser: type[ParserBase],
    group: GroupT,
    lineage: PrimitiveLineageT,
    whole_dir: str | None,
) -> dict[str, np.ndarray]:
    """Build the wholes of `prop` from its segment or whole files."""
    paths = _select(observations, prop)
    if lineage == 'segment':
        return organizer.whole_from_segments(
            prop, paths, parser, group, relation, save_to=whole_dir)
    return organizer.whole_from_file(paths, parser, group)


def _check_edges(
    prop: PropertySpec,
    wholes: Mapping[str, np.ndarray],
    edges: Mapping[str, np.ndarray],
    parser: type[ParserBase],
    group: GroupT,
) -> None:
    """
    Raise ``ValueError`` unless every whole has edges and the wholes of each
    ensemble hold identical edges.
    """
    missing = sorted(set(wholes) - set(edges))
    if missing:
        raise ValueError(
            f"'{prop.prop}': no '{prop.edge_prop}' file for the wholes "
            f"{missing}."
        )
    first_of: dict[str, str] = {}
    for whole_name in wholes:
        ens_name = str(getattr(parser(whole_name, 'whole', group),
                               'ensemble_long'))
        reference = first_of.setdefault(ens_name, whole_name)
        if not np.array_equal(edges[reference], edges[whole_name]):
            raise ValueError(
                f"'{prop.prop}': the wholes '{reference}' and '{whole_name}'"
                f" of ensemble '{ens_name}' hold different bin edges."
            )


def _reduce_output(
    output: Derived,
    parser: type[ParserBase],
    group: GroupT,
    ens_dir: str | None,
    ens_avg_dir: str | None,
) -> dict[str, np.ndarray | pd.DataFrame]:
    """Build and save the ensembles and ensemble averages of one output."""
    if output.whole_type == 'dataframe':
        frames = {name: data for name, data in output.wholes.items()
                  if isinstance(data, pd.DataFrame)}
        if len(frames) != len(output.wholes):
            raise ValueError(
                f"'{output.name}': every 'dataframe' whole must be a "
                "pandas.DataFrame."
            )
        long_ensembles = _frame_ensembles(output.name, frames, parser, group,
                                          ens_dir)
        return {
            ens_name: _frame_ensemble_avg(output.name, long_frame,
                                          ens_name, group, ens_avg_dir)
            for ens_name, long_frame in long_ensembles.items()
        }
    ensembles = organizer.ensemble(
        output.name, output.wholes, parser, group, output.whole_type,
        whole_edges=(dict(output.edges) if output.edges is not None
                     else None),
        save_to=ens_dir,
    )
    averages = organizer.ensemble_avg(
        output.name, ensembles, group, _ENS_TYPES[output.whole_type],
        save_to=ens_avg_dir,
    )
    return dict(averages)


def _frame_ensembles(
    prop: str,
    wholes: Mapping[str, pd.DataFrame],
    parser: type[ParserBase],
    group: GroupT,
    save_to: str | None,
) -> dict[str, pd.DataFrame]:
    """
    The long-form ensembles of ``'dataframe'`` wholes: every whole's rows,
    with a ``whole`` and a ``row`` column, saved as
    ``<ensemble_long>-<group>-<prop>.csv``.
    """
    members: dict[str, dict[str, pd.DataFrame]] = {}
    for whole_name, frame in wholes.items():
        ens_name = str(getattr(parser(whole_name, 'whole', group),
                               'ensemble_long'))
        members.setdefault(ens_name, {})[whole_name] = frame
    ensembles: dict[str, pd.DataFrame] = {}
    for ens_name, frames in members.items():
        long_frame = pd.concat(
            list(frames.values()), keys=list(frames),
            names=['whole', 'row'],
        ).reset_index()
        ensembles[ens_name] = long_frame
        if save_to is not None:
            organizer.save_artifact(
                organizer.create_fullname(ens_name, group, prop),
                long_frame, save_to)
    return ensembles


def _frame_ensemble_avg(
    prop: str,
    long_frame: pd.DataFrame,
    ens_name: str,
    group: GroupT,
    save_to: str | None,
    exclude: Sequence[str] = ('bin_center',),
) -> pd.DataFrame:
    """
    The row-by-row mean, variance, and standard error across the wholes of
    one long-form ensemble, as ``<column>-mean``, ``<column>-var``, and
    ``<column>-sem`` per column, after the `exclude` columns of the first
    whole; saved as ``<ensemble_long>-<group>-<prop>-ensAvg.csv``.
    """
    kept = [col for col in exclude if col in long_frame.columns]
    first_whole = long_frame['whole'].iloc[0]
    first = long_frame.loc[long_frame['whole'] == first_whole]
    average: pd.DataFrame = first.set_index('row').loc[:, kept]
    grouped = long_frame.drop(columns=['whole', *kept]).groupby('row')
    means = grouped.mean()
    variances = grouped.var(ddof=1)
    errors = grouped.sem(ddof=1)
    for column in means.columns:
        average[f'{column}-mean'] = means[column]
        average[f'{column}-var'] = variances[column]
        average[f'{column}-sem'] = errors[column]
    result: pd.DataFrame = average.reset_index(drop=True)
    if save_to is not None:
        organizer.save_artifact(
            organizer.create_fullname(ens_name, group, f'{prop}-ensAvg'),
            result, save_to)
    return result


def _blocked_errors(values: np.ndarray, method: PlateauMethod) -> np.ndarray:
    """
    The blocked standard error of each column of `values`, shape
    ``(n_frames, n_series)``, by :func:`~polyphys.analyze.stats.sem_blocked`;
    0 for a constant column, which
    :func:`~polyphys.analyze.stats.block_average` rejects. A column
    with NaN is passed on, so that it raises there.
    """
    errors = np.zeros(values.shape[1])
    varying = ~(values.var(axis=0) == 0.0)
    if varying.any():
        errors[varying] = np.asarray(
            stats.sem_blocked(values[:, varying], method=method))
    return errors


def _check_wholes(
    prop: str, found: Sequence[str], expected: Sequence[str]
) -> None:
    """
    Raise ``ValueError`` unless the wholes `found` in a property's ensemble
    files are exactly the `expected` wholes of the stamps.
    """
    missing = sorted(set(expected) - set(found))
    unknown = sorted(set(found) - set(expected))
    if missing or unknown:
        raise ValueError(
            f"'{prop}': the wholes of its ensemble files and of whole_stamps "
            f"differ. In the stamps only: {missing}. In the files only: "
            f"{unknown}. If the names differ by a suffix, pass topology."
        )
