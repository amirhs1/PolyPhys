"""
The `simulate` subpackage works with the files of the simulation engine,
LAMMPS: it reads logs and dumps and reads and writes data files, sharing the
box and the error types of :mod:`polyphys.simulate.core`.
"""
from polyphys.simulate import core
from polyphys.simulate import log
from polyphys.simulate import dump
from polyphys.simulate import data
from polyphys.simulate.core import Box, LammpsFormatError
from polyphys.simulate.data import (
    BrokenDataError,
    DataFile,
    read_data,
    with_frame,
    write_data,
)
from polyphys.simulate.dump import (
    BrokenDumpError,
    DumpFrame,
    iter_dump,
    read_dump_frame,
)
from polyphys.simulate.log import BrokenLogError, LogFile, parse_log, read_log

__all__ = [
    "core", "log", "dump", "data",
    "Box", "LammpsFormatError",
    "BrokenLogError", "LogFile", "parse_log", "read_log",
    "BrokenDumpError", "DumpFrame", "iter_dump", "read_dump_frame",
    "BrokenDataError", "DataFile", "read_data", "with_frame", "write_data",
]
