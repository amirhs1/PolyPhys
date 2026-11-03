"""\
==========================================================
:mod:`polyphys.manage.projects`
==========================================================

This :mod:`~polyphys.manage.projects` module defines :class:`ProjectSpec`,
the record of the per-project facts that the parser class does not hold. A
parser class declares how a study's filenames read; a :class:`ProjectSpec`
declares how the study's analysis reads its artifacts: which particle group
is the polymer, which parser attributes become columns, how the bulk volume
fraction of crowders is rounded, which time series exist, and which size
attribute scales each species' density profile.

The facts the parser already holds, such as the geometry, the topology, and
the glob patterns of the ``space`` and ``ensemble_long`` names, are read from
the parser rather than stored twice.

PolyPhys holds no table of studies: a study builds its own spec where its
analysis lives. The specs of the maintainer's studies are in the repository's
top-level ``projects/`` folder, which is not part of the installed package.

Contents
========

- :class:`ProjectSpec`, a frozen record of one study.

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
    ``phi_c_bulk_round`` column (dimensionless).

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
from typing import Mapping

from .parser import ParserBase
from .types import DirectionT, EntityT, GeometryT, GroupT, PropertyT, TopologyT

# The parser attribute that holds the bulk volume fraction of crowders, which
# the organizer adds and rounds into 'phi_c_bulk_round'.
_CROWDER_ATTR = 'phi_bulk_c'

# The directions of the density profiles in each geometry.
_DIRECTIONS: dict[GeometryT, tuple[DirectionT, ...]] = {
    'cubic': ('r',),
    'cylindrical': ('r', 'z'),
}


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
        in the ``phi_c_bulk_round`` column.
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
