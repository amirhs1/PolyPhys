from .parser import TransFoci, SumRule, TransFociCubic
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

TransFociCubicParser = Type[TransFociCubic]
TransFociParser = Type[TransFoci]
SumRuleParser = Type[SumRule]
PolyPhysParser = Union[TransFociParser, SumRuleParser, TransFociCubicParser]

EnsWholes = Dict[WholeName, Union[np.ndarray, pd.DataFrame]]
EnsembleT = Tuple[EnsembleName, EnsWholes]
