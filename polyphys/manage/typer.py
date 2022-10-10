from .parser import TransFociCyl, SumRuleCyl, TransFociCubic
from typing import Type, Union, Tuple, Dict, NewType

import numpy as np
import pandas as pd

PropertyT = NewType('PropertyT', str)
SpeciesT = NewType('SpeciesT', str)
GroupT = NewType('GroupT', str)
DirectionT = NewType('DirectionT', str)
AxisT = NewType('AxisT', int)
WholeName = NewType('WholeName', str)
EnsembleName = NewType('EnsembleName', str)

TimeSeriesT = Tuple[PropertyT, SpeciesT, GroupT]
NonScalarHistT = Tuple[PropertyT, SpeciesT, GroupT, AxisT]
NonScalarMatT = Tuple[PropertyT, SpeciesT, GroupT]
HistogramT = Tuple[DirectionT, SpeciesT, GroupT]
EdgeT = Union[Tuple[DirectionT, SpeciesT, GroupT],
              Tuple[DirectionT, GroupT]]

TransFociCubicT = Type[TransFociCubic]
TransFociCylT = Type[TransFociCyl]
TransFociT = Union[TransFociCylT, TransFociCubicT]
SumRuleCylT = Type[SumRuleCyl]
ParserT = Union[TransFociCylT, SumRuleCylT, TransFociCubicT]

EnsWholes = Dict[WholeName, Union[np.ndarray, pd.DataFrame]]
EnsembleT = Tuple[EnsembleName, EnsWholes]

FreqDataT = Dict[str, np.ndarray]
EdgeDataT = Dict[str, np.ndarray]
