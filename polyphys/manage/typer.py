"""\
==========================================================
:mod:`polyphys.manage.typer`
==========================================================

This :mod:`~polyphys.manage.typer` module provides type definitions and type
aliases to support clarity, type safety, and documentation for the data
structures and classes used within :mod:`polyphys` package.

Custom types here define domain-specific concepts such as `PropertyT`,
`EntityT`, `DirectionT`, and others, which clarify the usage of fundamental
types like `str`, `int`, and `bool` throughout the :mod:`polyphys` package.

The module defines composite types, such as tuples and dictionaries
(`TimeSeriesT`, `NonScalarHistT`, etc.), which group related information
together for structured data handling.

The module also provides type aliases for different classes, such as
`SumRuleCyl`. These types are used when passing the classes (not their
instances) to functions that need to distinguish between various parser types.
"""
from typing import Type, Union, Tuple, Dict, List, TextIO, IO, Any, Literal
from gzip import GzipFile
import numpy as np
import pandas as pd
from .parser import (
    ParserBase, SumRuleCyl, TransFociCyl, TransFociCub, HnsCyl, HnsCub,
    SumRuleCubHeteroLinear, SumRuleCubHeteroRing, TwoMonDepCub
)

# Type aliases for clarity
PropertyT = str
EntityT = str
GroupT = str
GeometryT = Literal['cubic', 'cylindrical']
PhaseT = str
# PhaseT = Literal['simAll', 'simCont', 'log', 'trj', 'probe', 'analysis',
#                 'viz', 'galaxy']
StageT = Literal['segment', 'wholeSim', 'ens', 'ensAvg', 'space', 'galaxy']
WholeName = str
EnsembleName = str
HasEdgeT = bool
TopologyT = str
DirectionT = Literal['x', 'y', 'z', 'r', 'theta']
LineageT = Literal['segment', 'whole', 'ensemble_long', 'ensemble', 'space']
PrimitiveLineageT = Literal['segment', 'whole']
PlaneT = Literal['xy', 'yz', 'zx']
BinT = Literal['ordinary', 'nonnegative', 'periodic']
AxisT = Literal[0, 1, 2]

# Type aliases for data structures
TimeSeriesT = Tuple[PropertyT, EntityT, GroupT]
NonScalarHistT = Tuple[PropertyT, EntityT, GroupT, AxisT]
NonScalarMatT = Tuple[PropertyT, EntityT, GroupT]
HistogramT = Tuple[DirectionT, EntityT, GroupT]
EdgeT = Tuple[DirectionT, GroupT]
WholeT = Dict[WholeName, Union[np.ndarray, pd.DataFrame]]
EnsembleT = Tuple[EnsembleName, Union[np.ndarray, pd.DataFrame]]
FreqDataT = Dict[str, np.ndarray]
EdgeDataT = Dict[str, np.ndarray]
HnsStatDictT = Dict[str, Union[List[int], np.ndarray]]

# Union types
InputType = Union[GzipFile, TextIO, IO[Any]]

# Type aliases for classes (used when passing the class itself)
ParserBaseT = Type[ParserBase]
TwoMonDepCubT = Type[TwoMonDepCub]
SumRuleCylT = Type[SumRuleCyl]
TransFociCylT = Type[TransFociCyl]
TransFociCubT = Type[TransFociCub]
SumRuleCubHeteroRingT = Type[SumRuleCubHeteroRing]
SumRuleCubHeteroLinearT = Type[SumRuleCubHeteroLinear]
HnsCubT = Type[HnsCub]
HnsCylT = Type[HnsCyl]

# Union types for classes
ClusterParserType = Union[
    TransFociCubT, TransFociCylT, SumRuleCubHeteroLinearT,
    SumRuleCubHeteroRingT
]
ParserType = Union[
    SumRuleCylT, TransFociCubT, TransFociCylT, HnsCubT, HnsCylT,
    SumRuleCubHeteroLinearT, SumRuleCubHeteroRingT, TwoMonDepCubT
]

# Type aliases for instances of the classes
ParserBaseInstance = ParserBase
TwoMonDepCubInstance = TwoMonDepCub
SumRuleCylInstance = SumRuleCyl
TransFociCylInstance = TransFociCyl
TransFociCubInstance = TransFociCub
SumRuleCubHeteroRingInstance = SumRuleCubHeteroRing
SumRuleCubHeteroLinearInstance = SumRuleCubHeteroLinear
HnsCubInstance = HnsCub
HnsCylInstance = HnsCyl

# Union types for class instances
ClusterParserInstance = Union[
    TransFociCubInstance, TransFociCylInstance, SumRuleCubHeteroLinearInstance,
    SumRuleCubHeteroRingInstance
]
ParserInstance = Union[
    SumRuleCylInstance, TransFociCubInstance, TransFociCylInstance,
    HnsCubInstance, HnsCylInstance, SumRuleCubHeteroLinearInstance,
    SumRuleCubHeteroRingInstance, TwoMonDepCubInstance
]
