"""\
==========================================================
:mod:`polyphys.manage.projects`
==========================================================

This :mod:`~polyphys.manage.projects` module defines :class:`ProjectSpec`,
the record of the per-project facts that the parser class does not hold, and
:class:`ProbeSpec`, the record of how the study's trajectories are probed. A
parser class declares how a study's filenames read; a :class:`ProjectSpec`
declares how the study's analysis reads its artifacts: which particle group
is the polymer, which parser attributes become columns, how the bulk volume
fraction of crowders is rounded, which time series exist, and which size
attribute scales each species' density profile. A :class:`ProbeSpec` declares
how :mod:`polyphys.probe` reads the study's trajectories: which MDAnalysis
selection is each species, which chain-size observables are measured, how
the chain is unwrapped, the bin sizes of the histograms, and, for a study
with binders, a :class:`BinderSpec` of how their binding is counted.

The facts the parser already holds, such as the geometry, the topology, and
the glob patterns of the ``space`` and ``ensemble_long`` names, are read from
the parser rather than stored twice.

PolyPhys holds no table of studies: a study builds its own spec where its
analysis lives. The specs of the maintainer's studies are in the repository's
top-level ``projects/`` folder, which is not part of the installed package.

Contents
========

- :class:`ProjectSpec`, a frozen record of one study.
- :class:`ProbeSpec`, a frozen record of how one study's trajectories are
  probed, built beside its :class:`ProjectSpec`.
- :data:`CHAIN_SIZE_OBSERVABLES`, the names a :class:`ProbeSpec` may list as
  chain-size observables.
- :class:`BinderSpec`, a frozen record of how the binding of a study's
  two-patch binders to its chain is counted, held by the :class:`ProbeSpec`
  of a study with binders.
- :data:`BINDER_SERIES`, the eight binding counts and the stems of their
  time series.

Dependencies
============
- Python standard library: `dataclasses`, `math`, `types`, `typing`.
- :mod:`polyphys.manage.parser`: :class:`~polyphys.manage.parser.ParserBase`.
- :mod:`polyphys.manage.types`: domain type aliases such as
  :py:data:`~polyphys.manage.types.GroupT`.

Notation
========
:math:`\\phi_c`
    The bulk volume fraction of crowders, the parser attribute
    ``phi_bulk_c``: the fraction of the volume available to crowder centres
    that the crowders fill (dimensionless).
:math:`\\Delta`
    The rounding step, `divisor`, of :math:`\\phi_c` in the
    ``phi_bulk_c_round`` column (dimensionless).

Examples
========
>>> from polyphys.manage.parser import SumRuleCyl
>>> spec = ProjectSpec(
...     parser=SumRuleCyl, group='bug', divisor=0.025,
...     physical_attrs=('space', 'nmon', 'dcyl', 'dcrowd'),
...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
...     tseries_properties=('gyrTMon',))
>>> spec.geometry, spec.space_pattern, spec.hierarchy
('cylindrical', 'N*D*ac*', 'N*')
"""
import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Callable, Literal, Mapping, get_args

from .parser import ParserBase
from .types import DirectionT, EntityT, GeometryT, GroupT, PropertyT, TopologyT

#: A quantity of a probe that is either the same number for every whole of a
#: study or a function of the whole's parser, such as a bin size that scales
#: with the crowder diameter. :meth:`ProbeSpec.resolve` turns either into the
#: number for one whole.
#:
#: .. versionadded:: 0.7.0
PerWhole = float | Callable[[ParserBase], float]

#: How a probe moves the positions of a frame before binning them. The
#: LAMMPS dump reader of MDAnalysis 2.4 and later translates every frame so
#: that the lower box bounds sit at the origin, which undoes the ``fix
#: recenter`` of the studies' LAMMPS scripts. ``'half'`` subtracts half the
#: box lengths of the frame, which restores the LAMMPS frame exactly;
#: ``'none'`` keeps the reader's frame; a triple is subtracted as given, so
#: the research code's shift by half the parser's ``(dcyl, dcyl, lcyl)`` or
#: ``lcube`` is reproducible.
#:
#: .. versionadded:: 0.7.0
RecentreT = Literal['none', 'half'] | tuple[float, float, float]

# The parser attribute that holds the bulk volume fraction of crowders, which
# the organizer adds and rounds into 'phi_bulk_c_round'.
_CROWDER_ATTR = 'phi_bulk_c'

# The directions of the density profiles in each geometry.
_DIRECTIONS: dict[GeometryT, tuple[DirectionT, ...]] = {
    'cubic': ('r',),
    'cylindrical': ('r', 'z'),
}

#: The chain-size observables that :class:`polyphys.probe.ChainSize`
#: measures per frame, in the research code's order: each is the name of a
#: time series without its species label, so ``'gyrT'`` with the chain label
#: ``'Mon'`` is the ``gyrTMon`` artifact. ``'principalT'`` is a per-frame
#: :math:`3 \\times 3` matrix; the others are scalars.
#:
#: .. versionadded:: 0.7.0
CHAIN_SIZE_OBSERVABLES: tuple[PropertyT, ...] = (
    'transSizeT', 'fsdT', 'gyrT', 'rfloryT', 'asphericityT', 'shapeT',
    'principalT')

# The chain-size observables that are per-frame matrices, which a project's
# `tseries_properties` does not list.
_MATRIX_OBSERVABLES: frozenset[PropertyT] = frozenset({'principalT'})

#: The eight per-frame binding counts that
#: :class:`polyphys.probe.BinderStatistics` measures, in the research code's
#: order: each attribute of
#: :class:`~polyphys.analyze.clusters.FrameBinding` mapped to the stem of its
#: time series and to which label completes the name, ``'patches'`` for the
#: three patch counts and ``'binders'`` for the five binder counts. With the
#: labels ``HnsPatch`` and ``HnsCore``, ``n_binders_bridging`` is the
#: ``nBridgeTHnsCore`` artifact.
#:
#: .. versionadded:: 0.7.0
BINDER_SERIES: Mapping[str, tuple[PropertyT, str]] = MappingProxyType({
    'n_patch_monomer_contacts': ('nBoundT', 'patches'),
    'n_patches_free': ('nFreeT', 'patches'),
    'n_patches_engaged': ('nEngagedT', 'patches'),
    'n_binders_free': ('nFreeT', 'binders'),
    'n_binders_bridging': ('nBridgeT', 'binders'),
    'n_binders_dangling': ('nDangleT', 'binders'),
    'n_pairs_cis': ('nCisT', 'binders'),
    'n_pairs_trans': ('nTransT', 'binders'),
})

# The number of patches per binder that binding_statistics counts.
_TWO_PATCHES = 2


@dataclass(frozen=True, slots=True)
class BinderSpec:
    """
    How the binding of a study's two-patch binders to its chain is counted:
    the facts :class:`polyphys.probe.BinderStatistics` needs beyond the
    selections of the :class:`ProbeSpec` that holds this record.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    patches : :py:data:`EntityT`
        The label in the probe spec's `selections` of the binders' patches,
        such as ``'HnsPatch'``, which names the three patch series
        ``nBoundT``, ``nFreeT``, and ``nEngagedT``. The patches of binder
        :math:`k` are the ones at positions :math:`2k` and :math:`2k + 1`
        of the selection, in atom-ID order.
    binders : :py:data:`EntityT`
        The label of the binders, such as ``'HnsCore'``, which names the
        five binder series ``nFreeT``, ``nBridgeT``, ``nDangleT``,
        ``nCisT``, and ``nTransT``. A label only: the binder cores are not
        selected.
    cutoff : float
        The binding threshold :math:`r_\\mathrm{thres}`, in units of the
        monomer diameter: a patch and a monomer closer than it are in
        contact. Finite and positive. A study constant that the spec
        states, never computed from the particle sizes; see the probe
        guide, "Choosing the cutoff".
    inclusive : bool, default True
        Whether a pair at exactly `cutoff` is a contact
        (:math:`r \\le r_\\mathrm{thres}`, the thesis criterion) or not
        (:math:`r < r_\\mathrm{thres}`, the research probe's).
    cis_threshold : int, default 4
        The largest genomic distance, in bonds, of a cis bridge;
        non-negative.
    patches_per_binder : int, default 2
        The number of patches per binder. Only ``2`` is accepted, which
        :func:`~polyphys.analyze.clusters.binding_statistics` counts; the
        field is kept so that binders of another valence can be declared
        once a kernel counts them.
    store_distances : bool, default False
        Whether the probe also keeps and saves the monomer-to-patch
        distance matrix of every frame, the input of the threshold sweep
        and of the binder-cluster post-step of the research code.

    Attributes
    ----------
    series : mapping of (str -> :py:data:`PropertyT`)
        Each attribute name of
        :class:`~polyphys.analyze.clusters.FrameBinding` in
        :data:`BINDER_SERIES` mapped to the full name of its time series,
        ``<stem><patches>`` or ``<stem><binders>``, in the research code's
        order.

    Raises
    ------
    ValueError
        If `patches` or `binders` is not a non-empty string or the two are
        equal; if `cutoff` is not a finite positive number; if
        `cis_threshold` is a negative integer or not an integer; or if
        `patches_per_binder` is not ``2``.

    Notes
    -----
    Construction and :attr:`series` are :math:`O(1)`.

    Examples
    --------
    >>> spec = BinderSpec(patches='HnsPatch', binders='HnsCore',
    ...                   cutoff=0.685)
    >>> spec.series['n_binders_bridging'], spec.series['n_patches_free']
    ('nBridgeTHnsCore', 'nFreeTHnsPatch')
    >>> BinderSpec('HnsPatch', 'HnsCore', cutoff=0.685, patches_per_binder=3)
    Traceback (most recent call last):
    ...
    ValueError: 'patches_per_binder' must be 2, not 3: ...
    """
    patches: EntityT
    binders: EntityT
    cutoff: float
    inclusive: bool = True
    cis_threshold: int = 4
    patches_per_binder: int = _TWO_PATCHES
    store_distances: bool = False

    def __post_init__(self) -> None:
        for name, label in (('patches', self.patches),
                            ('binders', self.binders)):
            if not (isinstance(label, str) and label):
                raise ValueError(
                    f"'{name}' must be a non-empty species label, not "
                    f"{label!r}.")
        if self.patches == self.binders:
            raise ValueError(
                f"'patches' and 'binders' must be different labels, not "
                f"both {self.patches!r}.")
        if not (isinstance(self.cutoff, (int, float))
                and not isinstance(self.cutoff, bool)
                and math.isfinite(self.cutoff) and self.cutoff > 0):
            raise ValueError(
                f"'cutoff' must be a finite positive number, not "
                f"{self.cutoff!r}.")
        if (isinstance(self.cis_threshold, bool)
                or not isinstance(self.cis_threshold, int)
                or self.cis_threshold < 0):
            raise ValueError(
                f"'cis_threshold' must be a non-negative integer, not "
                f"{self.cis_threshold!r}.")
        if (isinstance(self.patches_per_binder, bool)
                or self.patches_per_binder != _TWO_PATCHES):
            raise ValueError(
                f"'patches_per_binder' must be {_TWO_PATCHES}, not "
                f"{self.patches_per_binder!r}: binding_statistics counts "
                "two-patch binders only.")

    @property
    def series(self) -> Mapping[str, PropertyT]:
        """
        Each binding count's attribute name mapped to the full name of its
        time series, in the research code's order.
        """
        labels = {'patches': self.patches, 'binders': self.binders}
        return MappingProxyType({
            attr: f"{stem}{labels[which]}"
            for attr, (stem, which) in BINDER_SERIES.items()})


def _is_per_whole(value: object) -> bool:
    """Whether `value` is a finite positive number or a callable."""
    if callable(value):
        return True
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value > 0)


@dataclass(frozen=True, slots=True)
class ProjectSpec:
    """
    The per-project facts that the analysis of a study needs and its parser
    class does not hold.

    Every attribute a spec names is checked at construction against the
    attributes its parser sets at the ``ensemble_long`` lineage, the lineage
    whose names the organizer parses for these values.

    Parameters
    ----------
    parser : type of :class:`~polyphys.manage.parser.ParserBase`
        The study's concrete parser class.
    group : :py:data:`GroupT`
        The particle group of the polymer, such as ``'bug'`` or
        ``'nucleoid'``; one of the parser's groups, and not ``'all'``.
    divisor : float
        The step :math:`\\Delta > 0` to which the bulk volume fraction of
        crowders :math:`\\phi_c` (``phi_bulk_c``, dimensionless) is rounded
        in the ``phi_bulk_c_round`` column.
    physical_attrs : tuple of str
        The parser attributes added as columns to the ``space`` tables, in
        column order.
    species_sizes : mapping of (:py:data:`EntityT` -> str)
        Each species label of a density profile, such as ``'Mon'`` or
        ``'Crd'``, mapped to the parser attribute that holds its particle
        diameter, in units of the monomer diameter. Stored read-only, in the
        given order.
    tseries_properties : tuple of :py:data:`PropertyT`
        The per-frame time series of the study, named
        ``<property>T<Species>``, such as ``'gyrTMon'``.
    provisional : bool, default False
        Whether the entry is provisional: its data or values are not
        confirmed.

    Attributes
    ----------
    name : str
        The parser class name, which is the project name.
    geometry : :py:data:`GeometryT`
        The parser's geometry.
    topology : :py:data:`TopologyT`
        The parser's polymer topology.
    directions : tuple of :py:data:`DirectionT`
        The directions of the density profiles: ``('r',)`` in a cube and
        ``('r', 'z')`` in a cylinder.
    space_pattern : str
        A glob pattern that matches the study's ``space`` names: each
        keyword of the parser's ``space`` pattern followed by ``*``.
    space_hierarchy : str
        A glob prefix of the ``space`` names: their first keyword and ``*``.
    hierarchy : str
        A glob prefix of the ``ensemble_long`` names: their first keyword and
        ``*``.

    Raises
    ------
    TypeError
        If `parser` is not a concrete subclass of
        :class:`~polyphys.manage.parser.ParserBase`.
    ValueError
        If `group` is ``'all'`` or not one of the parser's groups; if
        `divisor` is not a finite positive number; or if `physical_attrs`
        or `species_sizes` names an attribute, or the parser lacks
        ``phi_bulk_c``, at ``ensemble_long``.

    Examples
    --------
    >>> from polyphys.manage.parser import HnsCub, SumRuleCyl
    >>> spec = ProjectSpec(
    ...     parser=HnsCub, group='nucleoid', divisor=0.01,
    ...     physical_attrs=('space', 'nmon', 'nhns', 'dcrowd'),
    ...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd', 'Hns': 'dhns'},
    ...     tseries_properties=('gyrTMon',))
    >>> spec.name, spec.geometry, spec.directions
    ('HnsCub', 'cubic', ('r',))
    >>> spec.space_pattern, spec.space_hierarchy, spec.hierarchy
    ('N*kbmm*nh*ac*epshc*', 'N*', 'N*')
    >>> spec.species_sizes['Hns']
    'dhns'

    A spec that names an attribute the parser does not set fails at once:

    >>> ProjectSpec(SumRuleCyl, 'bug', 0.025, ('nhns',), {}, ())
    Traceback (most recent call last):
    ...
    ValueError: SumRuleCyl does not set ['nhns'] at 'ensemble_long'.

    Notes
    -----
    Construction costs :math:`O(a + s + n)` for :math:`a` physical
    attributes, :math:`s` species, and the :math:`n` names that
    :meth:`~polyphys.manage.parser.ParserBase.attributes_at` collects; every
    derived attribute is :math:`O(1)`.
    """
    parser: type[ParserBase]
    group: GroupT
    divisor: float
    physical_attrs: tuple[str, ...]
    species_sizes: Mapping[EntityT, str] = field(hash=False)
    tseries_properties: tuple[PropertyT, ...]
    provisional: bool = False

    def __post_init__(self) -> None:
        # Copy the containers so that a caller's list or dict cannot change
        # the spec afterwards.
        object.__setattr__(self, 'physical_attrs', tuple(self.physical_attrs))
        object.__setattr__(self, 'tseries_properties',
                           tuple(self.tseries_properties))
        object.__setattr__(self, 'species_sizes',
                           MappingProxyType(dict(self.species_sizes)))
        parser = self.parser
        if not (isinstance(parser, type) and issubclass(parser, ParserBase)
                and parser._patterns is not None):
            raise TypeError(
                f"'parser' must be a concrete ParserBase subclass, not "
                f"{parser!r}.")
        name = parser.__name__
        if self.group == 'all' or self.group not in parser._groups:
            raise ValueError(
                f"{name}: 'group' must be the polymer group, one of "
                f"{[g for g in parser._groups if g != 'all']}, not "
                f"{self.group!r}.")
        if not (math.isfinite(self.divisor) and self.divisor > 0):
            raise ValueError(
                f"{name}: 'divisor' must be a finite positive number, not "
                f"{self.divisor!r}.")
        available = parser.attributes_at('ensemble_long')
        named = [*self.physical_attrs, *self.species_sizes.values(),
                 _CROWDER_ATTR]
        missing = sorted({attr for attr in named if attr not in available})
        if missing:
            raise ValueError(
                f"{name} does not set {missing} at 'ensemble_long'.")

    @property
    def name(self) -> str:
        """
        The parser class name, which is the project name.
        """
        return self.parser.__name__

    @property
    def geometry(self) -> GeometryT:
        """
        The parser's geometry.
        """
        return self.parser._geometry

    @property
    def topology(self) -> TopologyT:
        """
        The parser's polymer topology.
        """
        return self.parser._topology

    @property
    def directions(self) -> tuple[DirectionT, ...]:
        """
        The directions of the density profiles in the parser's geometry.
        """
        return _DIRECTIONS[self.geometry]

    @property
    def space_pattern(self) -> str:
        """
        A glob pattern that matches the study's ``space`` names.
        """
        return ''.join(f"{keyword}*"
                       for keyword in self.parser.pattern('space'))

    @property
    def space_hierarchy(self) -> str:
        """
        A glob prefix of the ``space`` names: their first keyword and ``*``.
        """
        return f"{self.parser.pattern('space')[0]}*"

    @property
    def hierarchy(self) -> str:
        """
        A glob prefix of the ``ensemble_long`` names: their first keyword
        and ``*``.
        """
        return f"{self.parser.pattern('ensemble_long')[0]}*"


@dataclass(frozen=True, slots=True)
class ProbeSpec:
    """
    How the trajectories of one study are probed: the facts that
    :mod:`polyphys.probe` needs beyond the study's :class:`ProjectSpec`.

    A probe spec is built beside the project spec, in the study's own
    folder, and reads the parser class and the geometry from it rather than
    storing them twice. Every observable it names is checked at construction
    against :data:`CHAIN_SIZE_OBSERVABLES` and, for the scalar ones, against
    the project's ``tseries_properties``.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    project : ProjectSpec
        The study's project spec.
    selections : mapping of (:py:data:`EntityT` -> str)
        Each species label, such as ``'Mon'`` or ``'Foci'``, mapped to the
        MDAnalysis selection string that picks its atoms, such as
        ``'resid 1'`` or ``'type 2'``. Stored read-only, in the given order.
    chain : :py:data:`EntityT`
        The label in `selections` of the polymer chain, whose chain-size
        time series are named ``<observable><chain>``.
    observables : tuple of :py:data:`PropertyT`
        The chain-size observables measured per frame, each a name in
        :data:`CHAIN_SIZE_OBSERVABLES`, in the order the probe stores them.
    unwrap : bool
        Whether the chain is unwrapped across the periodic boundaries before
        its gyration tensor is formed, the `unwrap` argument of
        ``AtomGroup.asphericity`` and ``AtomGroup.shape_parameter``:
        ``False`` for a cylinder, which is periodic along its axis only, and
        ``True`` for a periodic cube.
    bin_sizes : mapping of (:py:data:`DirectionT` -> :py:data:`PerWhole`), \
            default {}
        The bin size of the histogram along each direction, in units of the
        monomer diameter: a finite positive number, or a function of a
        whole's parser returning one, for a size that scales with a whole's
        particle diameters. Stored read-only. Left empty by a study whose
        histograms are not probed.
    histogram_species : tuple of :py:data:`EntityT`, default ()
        The labels in `selections` of the species whose spatial histograms
        :class:`polyphys.probe.SpatialHistograms` collects, each once. Left
        empty by a study whose histograms are not probed; otherwise
        `bin_sizes` and `block_frames` are required.
    block_frames : int or None, default None
        The block length :math:`b` of the two-dimensional maps: the number of
        consecutive frames whose maps are averaged into one stored map, a
        positive integer. Required with `histogram_species`.
    recentre : :py:data:`RecentreT`, default 'half'
        How the positions of a frame are moved before they are binned:
        ``'half'`` subtracts half the frame's box lengths, which restores the
        LAMMPS frame the dump reader translated; ``'none'`` keeps the reader's
        frame; a triple ``(dx, dy, dz)`` is subtracted as given. Stored as a
        tuple of floats when given as a triple.
    pair_species : :py:data:`EntityT` or None, default None
        The label in `selections` of the species whose pair distances and
        contacts :class:`polyphys.probe.PairDistances` collects, such as
        ``'Foci'``; ``None`` when the study has no such species. Requires
        `contact_cutoff`.
    contact_cutoff : :py:data:`PerWhole` or None, default None
        The centre-to-centre distance, in units of the monomer diameter, at
        or below which two particles of `pair_species` are in direct
        contact: a finite positive number or a function of a whole's parser,
        such as the sum of the large-monomer and crowder diameters.
    minimum_image : bool, default True
        Whether the pair distances use the minimum-image convention along
        the periodic directions of the geometry: the axis of a cylinder, all
        three directions of a cube.
    binders : BinderSpec, optional
        How the binding of the study's two-patch binders is counted, for a
        study with binders such as the H-NS ones; ``None`` for a study
        without. Its `patches` must be a label in `selections` other than
        `chain`, and each of its eight time series must be one of the
        project's ``tseries_properties``.

    Attributes
    ----------
    name : str
        The project name, from the project spec.
    parser : type of :class:`~polyphys.manage.parser.ParserBase`
        The study's parser class, from the project spec.
    geometry : :py:data:`GeometryT`
        The study's geometry, from the project spec.

    Raises
    ------
    TypeError
        If `project` is not a :class:`ProjectSpec`, or if `binders` is
        neither ``None`` nor a :class:`BinderSpec`.
    ValueError
        If `selections` is empty or has a label or a selection that is not
        a non-empty string; if `chain` is not a label in `selections`; if
        `observables` is empty, repeats a name, or names an observable that
        is not in :data:`CHAIN_SIZE_OBSERVABLES`; if a scalar observable
        with the chain label is not one of the project's
        ``tseries_properties``; if a bin size is keyed by an unknown
        direction or is neither a finite positive number nor callable; if a
        histogram species is not a label in `selections` or repeats, or is
        given without `bin_sizes` or `block_frames`; if `block_frames` is
        not a positive integer; if `recentre` is neither ``'none'``,
        ``'half'``, nor a triple of finite numbers; if `pair_species` is not
        a label in `selections` or is given without `contact_cutoff`; if
        `contact_cutoff` is neither a finite positive number nor callable;
        or if the patches of `binders` are not a label in `selections`
        other than `chain`, or one of its time series is not among the
        project's ``tseries_properties``.

    Examples
    --------
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> project = ProjectSpec(
    ...     parser=SumRuleCyl, group='bug', divisor=0.025,
    ...     physical_attrs=('space', 'nmon', 'dcyl', 'dcrowd'),
    ...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
    ...     tseries_properties=('fsdTMon', 'gyrTMon', 'rfloryTMon'))
    >>> probe = ProbeSpec(
    ...     project=project, selections={'Mon': 'resid 1'}, chain='Mon',
    ...     observables=('rfloryT', 'fsdT', 'gyrT'), unwrap=False)
    >>> probe.name, probe.geometry, probe.selections['Mon']
    ('SumRuleCyl', 'cylindrical', 'resid 1')

    An observable the project's time series do not list fails at once:

    >>> ProbeSpec(project, {'Mon': 'resid 1'}, 'Mon', ('shapeT',), False)
    Traceback (most recent call last):
    ...
    ValueError: SumRuleCyl: 'observables' names ['shapeTMon'], which ...

    A study whose histograms are probed names the species, the bin sizes,
    which may depend on the whole, and the block length; the bin sizes are
    resolved for one whole from its parser:

    >>> probe = ProbeSpec(
    ...     project=project, selections={'Mon': 'resid 1', 'Crd': 'resid 0'},
    ...     chain='Mon', observables=('gyrT',), unwrap=False,
    ...     bin_sizes={'r': lambda info: 0.1 * min(info.dmon, info.dcrowd),
    ...                'z': 0.5},
    ...     histogram_species=('Mon', 'Crd'), block_frames=100)
    >>> whole = SumRuleCyl(
    ...     'N200epsilon5r10.5lz50sig0.5nc100dt0.005bdump1adump1ens1',
    ...     'whole', 'all')
    >>> probe.resolve(probe.bin_sizes['r'], whole), probe.recentre
    (0.05, 'half')

    A study with binders adds a :class:`BinderSpec`:

    >>> hns = ProjectSpec(
    ...     parser=SumRuleCyl, group='bug', divisor=0.025,
    ...     physical_attrs=('space', 'nmon', 'dcyl', 'dcrowd'),
    ...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
    ...     tseries_properties=(
    ...         'gyrTMon', 'nBoundTHnsPatch', 'nFreeTHnsPatch',
    ...         'nEngagedTHnsPatch', 'nFreeTHnsCore', 'nBridgeTHnsCore',
    ...         'nDangleTHnsCore', 'nCisTHnsCore', 'nTransTHnsCore'))
    >>> probe = ProbeSpec(
    ...     project=hns, selections={'Mon': 'resid 1', 'HnsPatch': 'type 2'},
    ...     chain='Mon', observables=('gyrT',), unwrap=False,
    ...     binders=BinderSpec('HnsPatch', 'HnsCore', cutoff=0.685))
    >>> probe.binders.series['n_pairs_cis']
    'nCisTHnsCore'

    Notes
    -----
    Construction costs :math:`O(s + o + b + h)` for :math:`s` selections,
    :math:`o` observables, :math:`b` bin sizes, and :math:`h` histogram
    species; every derived attribute is :math:`O(1)`.
    """
    project: ProjectSpec
    selections: Mapping[EntityT, str] = field(hash=False)
    chain: EntityT
    observables: tuple[PropertyT, ...]
    unwrap: bool
    bin_sizes: Mapping[DirectionT, PerWhole] = field(
        default_factory=dict, hash=False)
    histogram_species: tuple[EntityT, ...] = ()
    block_frames: int | None = None
    recentre: RecentreT = 'half'
    pair_species: EntityT | None = None
    contact_cutoff: PerWhole | None = None
    minimum_image: bool = True
    binders: BinderSpec | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.project, ProjectSpec):
            raise TypeError(
                f"'project' must be a ProjectSpec, not {self.project!r}.")
        if self.binders is not None \
                and not isinstance(self.binders, BinderSpec):
            raise TypeError(
                f"'binders' must be a BinderSpec or None, not "
                f"{self.binders!r}.")
        name = self.project.name
        # Copy the containers so that a caller's dict or list cannot change
        # the spec afterwards.
        object.__setattr__(self, 'selections',
                           MappingProxyType(dict(self.selections)))
        object.__setattr__(self, 'observables', tuple(self.observables))
        object.__setattr__(self, 'bin_sizes',
                           MappingProxyType(dict(self.bin_sizes)))
        if not self.selections or any(
                not (isinstance(key, str) and key
                     and isinstance(value, str) and value)
                for key, value in self.selections.items()):
            raise ValueError(
                f"{name}: 'selections' must map each species label to a "
                f"non-empty selection string, not {dict(self.selections)!r}.")
        if self.chain not in self.selections:
            raise ValueError(
                f"{name}: 'chain' must be a label in 'selections', one of "
                f"{list(self.selections)}, not {self.chain!r}.")
        if not self.observables:
            raise ValueError(f"{name}: 'observables' must not be empty.")
        unknown = [obs for obs in self.observables
                   if obs not in CHAIN_SIZE_OBSERVABLES]
        if unknown:
            raise ValueError(
                f"{name}: 'observables' names {unknown}, which are not in "
                f"{list(CHAIN_SIZE_OBSERVABLES)}.")
        if len(set(self.observables)) != len(self.observables):
            raise ValueError(
                f"{name}: 'observables' repeats a name: "
                f"{list(self.observables)}.")
        unlisted = [f"{obs}{self.chain}" for obs in self.observables
                    if obs not in _MATRIX_OBSERVABLES
                    and f"{obs}{self.chain}"
                    not in self.project.tseries_properties]
        if unlisted:
            raise ValueError(
                f"{name}: 'observables' names {unlisted}, which are not "
                f"among the project's 'tseries_properties' "
                f"{list(self.project.tseries_properties)}.")
        directions = get_args(DirectionT)
        for direction, size in self.bin_sizes.items():
            if direction not in directions:
                raise ValueError(
                    f"{name}: 'bin_sizes' has the direction {direction!r}; "
                    f"each must be one of {list(directions)}.")
            if not _is_per_whole(size):
                raise ValueError(
                    f"{name}: the bin size along {direction!r} must be a "
                    f"finite positive number or a callable, not {size!r}.")
        object.__setattr__(self, 'histogram_species',
                           tuple(self.histogram_species))
        unselected = [species for species in self.histogram_species
                      if species not in self.selections]
        if unselected:
            raise ValueError(
                f"{name}: 'histogram_species' names {unselected}, which are "
                f"not labels in 'selections' {list(self.selections)}.")
        if len(set(self.histogram_species)) != len(self.histogram_species):
            raise ValueError(
                f"{name}: 'histogram_species' repeats a label: "
                f"{list(self.histogram_species)}.")
        if self.histogram_species and not self.bin_sizes:
            raise ValueError(
                f"{name}: 'histogram_species' needs 'bin_sizes'.")
        if self.histogram_species and self.block_frames is None:
            raise ValueError(
                f"{name}: 'histogram_species' needs 'block_frames'.")
        if self.block_frames is not None and not (
                isinstance(self.block_frames, int)
                and not isinstance(self.block_frames, bool)
                and self.block_frames > 0):
            raise ValueError(
                f"{name}: 'block_frames' must be a positive integer, not "
                f"{self.block_frames!r}.")
        if self.recentre not in ('none', 'half'):
            shift = self.recentre
            if not (isinstance(shift, (tuple, list)) and len(shift) == 3
                    and all(isinstance(value, (int, float))
                            and not isinstance(value, bool)
                            and math.isfinite(value) for value in shift)):
                raise ValueError(
                    f"{name}: 'recentre' must be 'none', 'half', or a triple "
                    f"of finite numbers, not {shift!r}.")
            object.__setattr__(self, 'recentre',
                               tuple(float(value) for value in shift))
        if self.pair_species is not None:
            if self.pair_species not in self.selections:
                raise ValueError(
                    f"{name}: 'pair_species' must be a label in "
                    f"'selections', one of {list(self.selections)}, not "
                    f"{self.pair_species!r}.")
            if self.contact_cutoff is None:
                raise ValueError(
                    f"{name}: 'pair_species' needs 'contact_cutoff'.")
        if self.contact_cutoff is not None and not _is_per_whole(
                self.contact_cutoff):
            raise ValueError(
                f"{name}: 'contact_cutoff' must be a finite positive number "
                f"or a callable, not {self.contact_cutoff!r}.")

        if self.binders is not None:
            if self.binders.patches not in self.selections \
                    or self.binders.patches == self.chain:
                raise ValueError(
                    f"{name}: the patches of 'binders' must be a label in "
                    f"'selections' other than the chain {self.chain!r}, "
                    f"not {self.binders.patches!r}.")
            unlisted = [
                series for series in self.binders.series.values()
                if series not in self.project.tseries_properties]
            if unlisted:
                raise ValueError(
                    f"{name}: 'binders' names the series {unlisted}, which "
                    f"are not among the project's 'tseries_properties' "
                    f"{list(self.project.tseries_properties)}.")

    @staticmethod
    def resolve(value: PerWhole, parser: ParserBase) -> float:
        """
        The number a :py:data:`PerWhole` quantity takes for one whole.

        Parameters
        ----------
        value : :py:data:`PerWhole`
            A number, returned as a float, or a function of `parser`.
        parser : ParserBase
            The whole's parser, passed to a callable `value`.

        Returns
        -------
        float
            The number for this whole.

        Raises
        ------
        ValueError
            If a callable `value` returns anything but a finite positive
            number.

        Examples
        --------
        >>> from polyphys.manage.parser import SumRuleCyl
        >>> whole = SumRuleCyl(
        ...     'N200epsilon5r10.5lz50sig0.5nc100dt0.005bdump1adump1ens1',
        ...     'whole', 'all')
        >>> ProbeSpec.resolve(0.5, whole)
        0.5
        >>> ProbeSpec.resolve(lambda info: info.dmon + info.dcrowd, whole)
        1.5

        Notes
        -----
        :math:`O(1)` beyond the cost of the callable.

        .. versionadded:: 0.7.0
        """
        if callable(value):
            result = value(parser)
            if not (isinstance(result, (int, float))
                    and not isinstance(result, bool)
                    and math.isfinite(result) and result > 0):
                raise ValueError(
                    f"A per-whole quantity of {parser.name} must resolve to "
                    f"a finite positive number, not {result!r}.")
            return float(result)
        return float(value)

    @property
    def name(self) -> str:
        """
        The project name, from the project spec.
        """
        return self.project.name

    @property
    def parser(self) -> type[ParserBase]:
        """
        The study's parser class, from the project spec.
        """
        return self.project.parser

    @property
    def geometry(self) -> GeometryT:
        """
        The study's geometry, from the project spec.
        """
        return self.project.geometry
