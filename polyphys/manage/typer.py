from typing import Type, Union, Tuple, Dict, NewType, TypeVar
import numpy as np
import pandas as pd
# pylint: disable-next=reportShadowedImports
from .parser import (
    ParserBase, TransFociCyl, TransFociCub, SumRuleCyl, HnsCub, Snapshot
)

# String types:

PropertyT = NewType('PropertyT', str)
SpeciesT = NewType('SpeciesT', str)
GroupT = NewType('GroupT', str)
DirectionT = NewType('DirectionT', str)
AxisT = NewType('AxisT', int)
WholeName = NewType('WholeName', str)
EnsembleName = NewType('EnsembleName', str)

# Data unit types
TimeSeriesT = Tuple[PropertyT, SpeciesT, GroupT]
NonScalarHistT = Tuple[PropertyT, SpeciesT, GroupT, AxisT]
NonScalarMatT = Tuple[PropertyT, SpeciesT, GroupT]
HistogramT = Tuple[DirectionT, SpeciesT, GroupT]
EdgeT = Union[Tuple[DirectionT, SpeciesT, GroupT],
              Tuple[DirectionT, GroupT]]
WholeT = Dict[WholeName, Union[np.ndarray, pd.DataFrame]]
EnsembleT = Tuple[EnsembleName, WholeT]
FreqDataT = Dict[str, np.ndarray]
EdgeDataT = Dict[str, np.ndarray]

# Parser types
ParserBase = TypeVar('ParserBase', bound=ParserBase)
TransFociCubT = Type[TransFociCub]
TransFociCylT = Type[TransFociCyl]
SumRuleCylT = Type[SumRuleCyl]
HnsCubT = Type[HnsCub]
ParserT = Union[TransFociCyl, TransFociCub, SumRuleCyl, HnsCub]

# Lammps Snapshot
SnapshotT = Type[Snapshot]
