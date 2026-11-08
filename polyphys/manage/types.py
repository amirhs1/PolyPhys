"""\
==========================================================
:mod:`polyphys.manage.types`
==========================================================

This :mod:`~polyphys.manage.types` module provides type definitions and type
aliases to support clarity, type safety, and documentation for the data
structures and classes used within the PolyPhys package.

Custom types here represent domain-specific concepts such as `PropertyT`,
`EntityT`, and `DirectionT`, which clarify the usage of fundamental
types like `str`, `int`, and `bool`. `PhaseT` and `StageT` are the phase and
stage vocabularies that :func:`~polyphys.manage.organizer.make_database`
checks its arguments against, and `NormalizerT` the ways a density profile
may be normalized.

Composite types such as tuples and dictionaries (`TimeSeriesT`,
`NonScalarHistT`, etc.) group related information for structured
data handling.

This module also defines type aliases for input sources and simulation
metadata.
"""
from typing import (
    TextIO, IO, Any, Callable, Literal, TYPE_CHECKING, TypeVar, TypeAlias,
    Iterable, Sequence
    )
from gzip import GzipFile
from os import PathLike
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .parser import ParserBase

__all__ = [
    'AnyPathIter',
    'AxisT',
    'BinT',
    'DirectionT',
    'EdgeDataT',
    'EdgeT',
    'EnsembleName',
    'EnsembleT',
    'EntityT',
    'FreqDataT',
    'GeometryT',
    'GroupT',
    'HasEdgeT',
    'HistogramT',
    'HnsStatdictT',
    'InputType',
    'LineageT',
    'NonScalarHistT',
    'NonScalarMatT',
    'NormalizerT',
    'ParserClassT',
    'ParserInstanceT',
    'PathGroup',
    'PathGroupIter',
    'PathIter',
    'Pathish',
    'PhaseT',
    'PlaneT',
    'PrimitiveLineageT',
    'PropertyT',
    'StageT',
    'TimeSeriesT',
    'TopologyT',
    'WholeName',
    'WholeRelationT',
    'WholeT',
]

#: A parser instance, an instance of a
#: :class:`~polyphys.manage.parser.ParserBase` subclass.
ParserInstanceT = TypeVar("ParserInstanceT", bound="ParserBase")
# A concrete alias, not type[ParserInstanceT]: an alias over a TypeVar is
# generic, and every bare use of it fails mypy --strict with [type-arg].
#: A parser class, a subclass of :class:`~polyphys.manage.parser.ParserBase`.
ParserClassT: TypeAlias = type["ParserBase"]


#: A file-system path, a ``str`` or an :class:`os.PathLike` of ``str``.
Pathish = str | PathLike[str]
#: An iterable of paths.
PathIter = Iterable[Pathish]
#: The paths of one group, such as the files that share a basename, e.g.
#: ``tuple[Path, ...]``.
PathGroup = Sequence[Pathish]
#: An iterable of path groups, e.g. ``list[tuple[Path, ...]]``.
PathGroupIter = Iterable[PathGroup]
#: An iterable of paths or of path groups, for user-facing functions.
AnyPathIter = PathIter | PathGroupIter

# --- Basic Type Aliases ---
#: A Cartesian axis, 0, 1, or 2 for x, y, or z.
AxisT = Literal[0, 1, 2]
#: How :func:`~polyphys.analyze.measurer.fixedsize_bins` fits a range to
#: whole bins, ``'ordinary'``, ``'nonnegative'``, or ``'periodic'``.
BinT = Literal['ordinary', 'nonnegative', 'periodic']
#: A direction of a histogram or profile, ``'x'``, ``'y'``, ``'z'``, the
#: radial ``'r'``, or the azimuthal ``'theta'``.
DirectionT = Literal['x', 'y', 'z', 'r', 'theta']
#: The name of an artifact at the ``ensemble`` level.
EnsembleName = str
#: A species label, such as ``'Mon'`` or ``'Crd'``.
EntityT = str
#: The geometry of a simulation box, a periodic ``'cubic'`` box or a
#: ``'cylindrical'`` one, periodic along its axis.
GeometryT = Literal['cubic', 'cylindrical']
#: A particle group, such as ``'bug'`` or ``'all'``.
GroupT = str
#: Whether a quantity has bin edges.
HasEdgeT = bool
#: A level of the artifact lineage.
LineageT = Literal['segment', 'whole', 'ensemble_long', 'ensemble', 'space']
#: How the segments of a whole combine;
#: :func:`~polyphys.manage.organizer.whole_from_segments` sums
#: ``'histogram'`` segments, joins ``'tseries'`` segments, and checks that
#: ``'bin_edge'`` segments agree.
WholeRelationT = Literal['histogram', 'tseries', 'bin_edge']
#: A phase directory of a project, in the order the data flows through them;
#: :func:`~polyphys.manage.organizer.make_database` accepts exactly these.
PhaseT = Literal['simsAll', 'simsCont', 'logs', 'trjs', 'probes', 'analysis',
                 'viz', 'galaxy', 'allInOne']
#: How a density profile is normalized, by the value of its first bin
#: (``r = 0``, the cylinder axis or the chain's center of mass), its last bin
#: (the far field along z or in a periodic cube), or its maximum; by a
#: number; or by a function of the scaled profile that returns the value.
NormalizerT = (Literal['first', 'last', 'max'] | float
               | Callable[[pd.DataFrame], float])
#: A Cartesian plane, ``'xy'``, ``'yz'``, or ``'zx'``.
PlaneT = Literal['xy', 'yz', 'zx']
#: A lineage level that the probe writes, ``'segment'`` or ``'whole'``.
PrimitiveLineageT = Literal['segment', 'whole']
#: The name of a physical property, such as ``'gyrT'``.
PropertyT = str
#: A stage of a phase directory;
#: :func:`~polyphys.manage.organizer.make_database` accepts exactly these.
StageT = Literal['segment', 'wholeSim', 'ens', 'ensAvg', 'space', 'galaxy']
#: The topology of a chain, such as ``'linear'`` or ``'ring'``.
TopologyT = str
#: The name of an artifact at the ``whole`` level.
WholeName = str

# --- Composite Data Structures ---
#: Bin edges keyed by name.
EdgeDataT = dict[str, np.ndarray]
#: A direction and a particle group.
EdgeT = tuple[DirectionT, GroupT]
#: An ensemble's name and its data.
EnsembleT = tuple[EnsembleName, np.ndarray | pd.DataFrame]
#: Frequencies keyed by name.
FreqDataT = dict[str, np.ndarray]
#: A direction, a species, and a particle group.
HistogramT = tuple[DirectionT, EntityT, GroupT]
#: Binding statistics of H-NS proteins keyed by name.
HnsStatdictT = dict[str, list[int] | np.ndarray]
#: A property, a species, a particle group, and an axis.
NonScalarHistT = tuple[PropertyT, EntityT, GroupT, AxisT]
#: A property, a species, and a particle group.
NonScalarMatT = tuple[PropertyT, EntityT, GroupT]
#: A property, a species, and a particle group.
TimeSeriesT = tuple[PropertyT, EntityT, GroupT]
#: The data of wholes keyed by whole name.
WholeT = dict[WholeName, np.ndarray | pd.DataFrame]

# --- IO Types ---
#: An open file, a :class:`gzip.GzipFile`, a text file, or any file object.
InputType = GzipFile | TextIO | IO[Any]
