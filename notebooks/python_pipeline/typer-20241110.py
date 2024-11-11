from typing import Type, Union, Tuple, Dict, NewType, List, TypeVar
import numpy as np
import pandas as pd
from .parser import (
    ParserBase, SumRuleCyl, TransFociCyl, TransFociCub, HnsCyl, HnsCub,
    SumRuleCubHeteroLinear, SumRuleCubHeteroRing, TwoMonDep
    )

# Define custom type aliases for clarity:

PropertyT = NewType('PropertyT', str)
SpeciesT = NewType('SpeciesT', str)
GroupT = NewType('GroupT', str)
DirectionT = NewType('DirectionT', str)
AxisT = NewType('AxisT', int)
WholeName = NewType('WholeName', str)
EnsembleName = NewType('EnsembleName', str)
HasEdgeT = NewType('HasEdgeT', bool)

# Define custom type aliases for various data structures:
TimeSeriesT = Tuple[PropertyT, SpeciesT, GroupT]
NonScalarHistT = Tuple[PropertyT, SpeciesT, GroupT, AxisT]
NonScalarMatT = Tuple[PropertyT, SpeciesT, GroupT]
HistogramT = Tuple[DirectionT, SpeciesT, GroupT]
EdgeT = Tuple[DirectionT, GroupT]
WholeT = Dict[WholeName, Union[np.ndarray, pd.DataFrame]]
EnsembleT = Tuple[EnsembleName, Union[np.ndarray, pd.DataFrame]]
FreqDataT = Dict[str, np.ndarray]
EdgeDataT = Dict[str, np.ndarray]
HnsStatDict = Dict[str, Union[List[int], np.ndarray]]

# Define custom type aliases for various parser classes
ParserBaseT = Type[ParserBase]
TwoMonDepT = Type[TwoMonDep]
TransFociCylT = Type[TransFociCyl]
TransFociCubT = Type[TransFociCub]
SumRuleCylT = Type[SumRuleCyl]
SumRuleCubHeteroRingT = Type[SumRuleCubHeteroRing]
SumRuleCubHeteroLinearT = Type[SumRuleCubHeteroLinear]
HnsCubT = Type[HnsCub]
HnsCylT = Type[HnsCyl]
TransFociT = Union[TransFociCubT, TransFociCylT]
HnsT = Union[HnsCubT, HnsCylT]
# ParserT = Union[SumRuleCylT, TransFociCubT, TransFociCylT, HnsCubT, HnsCylT]
ParserT = Union[SumRuleCyl, TransFociCub, TransFociCyl, HnsCub, HnsCyl,
                SumRuleCubHeteroLinear, SumRuleCubHeteroRing, TwoMonDep]
ParserT = Union[SumRuleCylT, TransFociCubT, TransFociCylT, HnsCubT, HnsCylT,
                SumRuleCubHeteroLinearT, SumRuleCubHeteroRingT, TwoMonDepT]
ParserCylT = Union[SumRuleCylT, TransFociCylT, HnsCylT]
ParserCubT = Union[TransFociCubT, HnsCubT]


TExcludedVolume = TypeVar("TExcludedVolume", bound="ExcludedVolume")
TFreeEnergyVirial = TypeVar("TFreeEnergyVirial", bound="FreeEnergyVirial")