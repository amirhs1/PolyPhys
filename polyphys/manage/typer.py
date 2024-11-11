"""\
==========================================================
:mod:`polyphys.manage.typer`
==========================================================

This :mod:`~polyphys.manage.typer` module provides type definitions and type
aliases to support clarity, type safety, and documentation for the data
structures and classes used within molecular simulation projects. These type
hints are used throughout the `polyphys.manage` package to specify expected
data types and class types, enabling more reliable code by supporting static
type-checking and improving code readability.

The module primarily defines type aliases for simulation properties, parser
classes, and custom types, facilitating a structured approach to representing
various domain-specific entities such as species, properties, and groups.
These aliases simplify the code by abstracting complex data structures and
parser class types, making the codebase easier to understand and maintain.

Custom Types
============
Custom types here define domain-specific concepts such as `PropertyT`,
`SpeciesT`, `DirectionT`, and others, which clarify the usage of fundamental
types like `str`, `int`, and `bool` within the simulation context.

The module also defines composite types, such as tuples and dictionaries
(`TimeSeriesT`, `NonScalarHistT`, etc.), which group related information
together for structured data handling.

Parser Classes
==============
The module provides type aliases for different parser classes used in
molecular simulation projects, such as `SumRuleCyl`, `TransFociCyl`, `HnsCyl`,
and others. These types are used when passing the classes (not their instances)
to functions that need to distinguish between various parser types.

Classes
-------
Each parser type alias is defined to match the parser classes used for
cylindrical and cubic geometries, supporting modular and type-safe handling of
different parser types.

Type Aliases
============
- PropertyT
- SpeciesT
- GroupT
- DirectionT
- AxisT
- WholeName
- EnsembleName
- HasEdgeT
- TimeSeriesT
- NonScalarHistT
- HistogramT
- ParserT, ParserCylT, ParserCubT

Dependencies
============
- `numpy`: For handling numerical data arrays.
- `pandas`: For data structures (e.g., `DataFrame`) used in simulations.
- `typing`: For type hinting support in type aliases.
- `polyphys.manage.parser`: Import parser classes such as `SumRuleCyl`,
  `TransFociCub`, etc.

Usage
=====
Type aliases and parser class types are used throughout the package to annotate
function parameters, variables, and return types. This setup enhances the
safety and readability of code by ensuring consistent types across different
modules and preventing type-related errors during runtime.

Examples
========
The following is an example of using `ParserT` as a function argument type,
indicating that it accepts any parser class (e.g., `SumRuleCyl`):

>>> def configure_parser(parser_class: ParserT):
>>>     parser = parser_class()
>>>     # Configure parser

In another example, composite types such as `WholeT` can be used to define data
structures that combine `numpy` arrays and `pandas` DataFrames, as in the
following:

>>> whole_data: WholeT = {
>>>     "density_data": np.array([1.0, 2.0, 3.0]),
>>>     "summary": pd.DataFrame({"A": [1, 2], "B": [3, 4]})
>>> }

Notes
=====
Type hints defined in this module enable the use of static type checkers,
such as `mypy`, which verify the correctness of data types in the code
before runtime. This improves code robustness, particularly for large projects
with complex data structures and diverse parser classes.
"""
from typing import Type, Union, Tuple, Dict, List
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

# Type aliases for parser classes
# Type aliases for individual parser classes
ParserBaseT = Type[ParserBase]
TwoMonDepT = Type[TwoMonDep]
SumRuleCylT = Type[SumRuleCyl]
TransFociCylT = Type[TransFociCyl]
TransFociCubT = Type[TransFociCub]
SumRuleCubHeteroRingT = Type[SumRuleCubHeteroRing]
SumRuleCubHeteroLinearT = Type[SumRuleCubHeteroLinear]
HnsCubT = Type[HnsCub]
HnsCylT = Type[HnsCyl]

# Union types for different parser categories
ParserType = Union[
    SumRuleCylT, TransFociCubT, TransFociCylT, HnsCubT, HnsCylT,
    SumRuleCubHeteroLinearT, SumRuleCubHeteroRingT, TwoMonDepT
]
ParserCylType = Union[SumRuleCylT, TransFociCylT, HnsCylT]
ParserCubType = Union[TransFociCubT, HnsCubT]
