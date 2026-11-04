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
the chain is unwrapped, and the bin sizes of the histograms.

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
from typing import Mapping, get_args

from .parser import ParserBase
from .types import DirectionT, EntityT, GeometryT, GroupT, PropertyT, TopologyT

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
    bin_sizes : mapping of (:py:data:`DirectionT` -> float), default {}
        The bin size of the histogram along each direction, in units of the
        monomer diameter, each finite and positive. Stored read-only. Left
        empty by a study whose histograms are not yet probed.

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
        If `project` is not a :class:`ProjectSpec`.
    ValueError
        If `selections` is empty or has a label or a selection that is not
        a non-empty string; if `chain` is not a label in `selections`; if
        `observables` is empty, repeats a name, or names an observable that
        is not in :data:`CHAIN_SIZE_OBSERVABLES`; if a scalar observable
        with the chain label is not one of the project's
        ``tseries_properties``; or if a bin size is keyed by an unknown
        direction or is not a finite positive number.

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

    Notes
    -----
    Construction costs :math:`O(s + o + b)` for :math:`s` selections,
    :math:`o` observables, and :math:`b` bin sizes; every derived attribute
    is :math:`O(1)`.
    """
    project: ProjectSpec
    selections: Mapping[EntityT, str] = field(hash=False)
    chain: EntityT
    observables: tuple[PropertyT, ...]
    unwrap: bool
    bin_sizes: Mapping[DirectionT, float] = field(
        default_factory=dict, hash=False)

    def __post_init__(self) -> None:
        if not isinstance(self.project, ProjectSpec):
            raise TypeError(
                f"'project' must be a ProjectSpec, not {self.project!r}.")
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
            if not (isinstance(size, (int, float))
                    and math.isfinite(size) and size > 0):
                raise ValueError(
                    f"{name}: the bin size along {direction!r} must be a "
                    f"finite positive number, not {size!r}.")

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
