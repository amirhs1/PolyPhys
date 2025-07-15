"""\
==========================================================
:mod:`polyphys.manage.types`
==========================================================

This :mod:`~polyphys.manage.types` module provides type definitions and type
aliases to support clarity, type safety, and documentation for the data
structures and classes used within :mod:`polyphys` package.

Custom types here represent domain-specific concepts such as `PropertyT`,
`EntityT`, and `DirectionT`, which clarify the usage of fundamental
types like `str`, `int`, and `bool`.

Composite types such as tuples and dictionaries (`TimeSeriesT`,
`NonScalarHistT`, etc.) group related information for structured
data handling.

This module also defines type aliases for input sources and simulation
metadata.
"""
from typing import (
    Union, Tuple, Dict, List, TextIO, IO, Any, Literal, TYPE_CHECKING,
    TypeAlias
    )
from gzip import GzipFile
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .parser import ParserBase

ParserType: TypeAlias = "ParserBase"


# --- Basic Type Aliases ---
AxisT = Literal[0, 1, 2]
BinT = Literal['ordinary', 'nonnegative', 'periodic']
DirectionT = Literal['x', 'y', 'z', 'r', 'theta']
EnsembleName = str
EntityT = str
GeometryT = Literal['cubic', 'cylindrical']
GroupT = str
HasEdgeT = bool
LineageT = Literal['segment', 'whole', 'ensemble_long', 'ensemble', 'space']
WholeRelationT = Literal['histogram', 'tseries', 'bin_edge']
PhaseT = str
# PhaseT = Literal['simAll', 'simCont', 'log', 'trj', 'probe', 'analysis',
#                 'viz', 'galaxy']
PlaneT = Literal['xy', 'yz', 'zx']
PrimitiveLineageT = Literal['segment', 'whole']
PropertyT = str
StageT = Literal['segment', 'wholeSim', 'ens', 'ensAvg', 'space', 'galaxy']
TopologyT = str
WholeName = str


# --- Composite Data Structures ---
EdgeDataT = Dict[str, np.ndarray]
EdgeT = Tuple[DirectionT, GroupT]
EnsembleT = Tuple[EnsembleName, Union[np.ndarray, pd.DataFrame]]
FreqDataT = Dict[str, np.ndarray]
HistogramT = Tuple[DirectionT, EntityT, GroupT]
HnsStatDictT = Dict[str, Union[List[int], np.ndarray]]
NonScalarHistT = Tuple[PropertyT, EntityT, GroupT, AxisT]
NonScalarMatT = Tuple[PropertyT, EntityT, GroupT]
TimeSeriesT = Tuple[PropertyT, EntityT, GroupT]
WholeT = Dict[WholeName, Union[np.ndarray, pd.DataFrame]]


# --- IO Types ---
InputType = Union[GzipFile, TextIO, IO[Any]]
