"""\
==========================================================
:mod:`polyphys.manage.typer`
==========================================================

This :mod:`~polyphys.manage.typer` module provides type definitions and type
aliases to support clarity, type safety, and documentation for the data
structures and classes used within :mod:`polyphys` package.

Custom types here define domain-specific concepts such as `PropertyT`,
`SpeciesT`, `DirectionT`, and others, which clarify the usage of fundamental
types like `str`, `int`, and `bool` throughout the :mod:`polyphys` package.

The module defines composite types, such as tuples and dictionaries
(`TimeSeriesT`, `NonScalarHistT`, etc.), which group related information
together for structured data handling.

The module also provides type aliases for different classes, such as
`SumRuleCyl`. These types are used when passing the classes (not their
instances) to functions that need to distinguish between various parser types.

Dependencies
============
- `numpy`
- `pandas`
- `typing`
- `polyphys.manage.parser`: Import parser classes such as `SumRuleCyl`,
  `TransFociCub`, etc.

Usage
=====
Type aliases and class types are used throughout the :mode:`polyphys` package
to annotate function parameters, variables, and return types.

Notes
=====
Type hints defined in this module enable the use of static type checkers,
such as `mypy`, which verify the correctness of data types in the code
before runtime.
"""
from typing import Type, Union, Tuple, Dict, List, TextIO, IO, Any
from gzip import GzipFile
import numpy as np
import pandas as pd
from .parser import (
    ParserBase, SumRuleCyl, TransFociCyl, TransFociCub, HnsCyl, HnsCub,
    SumRuleCubHeteroLinear, SumRuleCubHeteroRing, TwoMonDep
)

# Type aliases for clarity
PropertyT = str
SpeciesT = str
GroupT = str
DirectionT = str
AxisT = int
PhaseT = str
WholeName = str
EnsembleName = str
HasEdgeT = bool

# Type aliases for data structures
TimeSeriesT = Tuple[PropertyT, SpeciesT, GroupT]
NonScalarHistT = Tuple[PropertyT, SpeciesT, GroupT, AxisT]
NonScalarMatT = Tuple[PropertyT, SpeciesT, GroupT]
HistogramT = Tuple[DirectionT, SpeciesT, GroupT]
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
TwoMonDepT = Type[TwoMonDep]
SumRuleCylT = Type[SumRuleCyl]
TransFociCylT = Type[TransFociCyl]
TransFociCubT = Type[TransFociCub]
SumRuleCubHeteroRingT = Type[SumRuleCubHeteroRing]
SumRuleCubHeteroLinearT = Type[SumRuleCubHeteroLinear]
HnsCubT = Type[HnsCub]
HnsCylT = Type[HnsCyl]

# Union types for clases
ParserType = Union[
    SumRuleCylT, TransFociCubT, TransFociCylT, HnsCubT, HnsCylT,
    SumRuleCubHeteroLinearT, SumRuleCubHeteroRingT, TwoMonDepT
]

# Type aliases for instances of the classes
ParserBaseInstance = ParserBase
TwoMonDepInstance = TwoMonDep
SumRuleCylInstance = SumRuleCyl
TransFociCylInstance = TransFociCyl
TransFociCubInstance = TransFociCub
SumRuleCubHeteroRingInstance = SumRuleCubHeteroRing
SumRuleCubHeteroLinearInstance = SumRuleCubHeteroLinear
HnsCubInstance = HnsCub
HnsCylInstance = HnsCyl

# Union types for class instances
ParserInstance = Union[
    SumRuleCylInstance, TransFociCubInstance, TransFociCylInstance,
    HnsCubInstance, HnsCylInstance, SumRuleCubHeteroLinearInstance,
    SumRuleCubHeteroRingInstance, TwoMonDepInstance
]
