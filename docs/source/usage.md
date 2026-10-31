# Usage

To begin using **PolyPhys**, first install the package as described in the
[Installation](#installation_target) section below.

(installation_target)=
## Installation

PolyPhys is distributed from source and archived on Zenodo. It is not published
to PyPI, so install it directly from the repository:

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git"
```

To install a specific archived release, append the tag:

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git@v0.4.0"
```

PolyPhys requires Python 3.11 or newer. Its runtime dependencies are NumPy and
pandas.

### Development install

If you are working from a source checkout and want the documentation and
development tools:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[docs,dev]"
```

## Quickstart

`polyphys.manage.parser` decodes a simulation filename into typed physical
attributes and derives the artifact's position in the project hierarchy:

```python
>>> from polyphys.manage.parser import SumRuleCyl
>>> artifact = SumRuleCyl(
...     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1.j03',
...     lineage='segment',
...     group='all',
... )
>>> artifact.nmon
200
>>> float(artifact.dcyl)
20.0
>>> artifact.space
'N200D20ac2'

```

## Lineage names

A parser sets the artifact's name at its own lineage and at each parent
lineage, as the attributes `segment`, `whole`, `ensemble_long`, `ensemble`,
and `space`; the organizer uses them as join keys. Two rules fix their text.

- The long names, `segment`, `whole`, and `ensemble_long`, are cut from the
  parsed name, so they keep its text, and they carry no topology suffix. A
  segment `…ens1.j03.ring` and a whole `…ens1.ring` of one simulation both
  give the `whole` name `…ens1`.
- The short names, `ensemble` and `space`, are built from the parsed values.
  A float that is a whole number prints as an integer, and any other float
  prints as Python's `repr`, the shortest text that reads back as the same
  float. The `space` above is `N200D20ac2` although the name wrote `sig2.0`
  and `r10.5`, while a crowder size of `1.5` stays `ac1.5`.

```python
>>> from polyphys.manage.parser import HnsCub
>>> artifact = HnsCub('N200kbmm2nh12ac2l25epshc1.5nc0ens1.ring', 'whole',
...                   'nucleoid')
>>> artifact.whole, artifact.space
('N200kbmm2nh12ac2l25epshc1.5nc0ens1', 'N200kbmm2nh12ac2epshc1.5')

```

Names written by the research code that preceded PolyPhys can differ from
these in two ways: a `whole` it parsed from a whole file keeps `.ring` or
`.linear`, and a short name can print `ac2.0`. Parsing such a short name at
its lineage and reading the rebuilt name gives the current form.

## Declaring a parser for a new study

A parser class is a declaration: it lists the study's filename keywords, the
order in which each lineage writes them, any fixed values, the particle
species whose bulk number density and volume fraction to compute, and the
attribute that holds each group's dump period (here fixed in `_defaults`).
`ParserBase` does the parsing. Each keyword gets a `KeywordSpec` naming the
attribute that stores its value and the value's type; a `float` keyword may
also convert the value as `scale * value + offset`, here `lcube = 2 l` from
the half-side `l`:

```python
>>> from polyphys.manage.parser import KeywordSpec, ParserBase, Species
>>> class ChainCub(ParserBase):
...     _geometry = 'cubic'
...     _topology = 'linear'
...     _groups = ['bug', 'all']
...     _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
...     _keywords = {
...         'N': KeywordSpec('nmon', int),
...         'ac': KeywordSpec('dcrowd', float),
...         'nc': KeywordSpec('ncrowd', int),
...         'l': KeywordSpec('lcube', float, scale=2.0),
...         'ens': KeywordSpec('ensemble_id', int),
...         'j': KeywordSpec('segment_id', int),
...     }
...     _patterns = {
...         'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
...         'whole': ('N', 'ac', 'nc', 'l', 'ens'),
...         'ensemble_long': ('N', 'ac', 'nc', 'l'),
...         'ensemble': ('N', 'ac', 'nc'),
...         'space': ('N', 'ac'),
...     }
...     _defaults = {'dmon': 1.0, 'bdump': 1000, 'adump': 5000}
...     _species = (Species('m', 'nmon', 'dmon'),
...                 Species('c', 'ncrowd', 'dcrowd'))
...     _project_attributes = {
...         'segment': ['dmon', 'rho_bulk_m', 'rho_bulk_c'],
...         'whole': ['dmon', 'rho_bulk_m', 'rho_bulk_c'],
...         'ensemble_long': ['dmon', 'rho_bulk_m', 'rho_bulk_c'],
...         'ensemble': ['dmon'],
...         'space': ['dmon'],
...     }
>>> artifact = ChainCub('N100ac2.0nc50l10.0ens1.bug', 'whole', 'bug')
>>> artifact.lcube, artifact.ensemble, artifact.space
(20.0, 'N100ac2nc50', 'N100ac2')
>>> artifact.dump_keyword, artifact.bdump
('bdump', 1000)

```

The declaration is checked when the class is created: a pattern keyword
without a `KeywordSpec`, an attribute that appears twice in one pattern, a
species whose size or number no long-form name sets, or a group without a
dump keyword raises `TypeError`. A
name that lacks a keyword of its lineage's pattern, or whose value is not a
number, raises `ValueError` naming the keyword, the lineage, and the name.

A value that is a function of parsed values but fits none of the tables,
such as the complementary percentage `pct_two = 100 - pct_one` of
`SumRuleCubHeteroRingBiDisperse`, goes in an override of
`_derived_attributes`, which runs at every lineage after the bulk values.

See the [API reference](api.rst) for the full set of modules, and the
[README](https://github.com/amirhs1/PolyPhys#readme) for a longer introduction
to the artifact-lineage model.

## Reading LAMMPS log files

`read_log` reads a LAMMPS log into a `LogFile` record of three DataFrames:
`thermo` holds one row per thermo output step of every run, `run_stats` one
row per run from the summary LAMMPS prints after it, and `wall_time` the total
wall time. A log holds one run per `run` or `minimize` command, so a
production loop of `run` commands gives one run per iteration; the `run`
column, counted from 0 in log order, picks runs:

```python
from polyphys.simulate.log import read_log

log = read_log("log.lammps")
production = log.thermo[log.thermo["run"] >= 2]
timing = log.run_stats[["run", "loop_time", "timesteps_per_s"]]
```

The thermo styles `one`, `custom`, `multi`, and `yaml` are all read.
`WARNING` and `fix halt` lines inside a table are skipped and logged at level
INFO, which `logging.basicConfig(level=logging.INFO)` shows. A log cut off
before a run's `Loop time of` line, as after a walltime kill, or one that
stopped with an `ERROR` raises `BrokenLogError`, so a loop over many logs can
skip it:

```python
from polyphys.simulate.log import BrokenLogError, read_log

for path in paths:
    try:
        log = read_log(path)
    except BrokenLogError as err:
        print(f"skipped {path}: {err}")
```

`parse_log` does the same for text already in memory, such as a member of a
tarball.

## Reading LAMMPS dumps and writing data files

`iter_dump` reads a text dump of `dump atom` or `dump custom`, plain or
gzipped, one frame at a time, so a long trajectory is never in memory whole.
Each `DumpFrame` holds its `timestep`, its `box`, the column names, and the
per-atom `data`, with rows sorted by atom ID unless `sort_by_id=False` is
passed. `positions()` returns Cartesian coordinates from whichever columns
the dump has: `x y z`, or the scaled `xs ys zs` that `dump atom` writes by
default; with `unwrapped=True`, `xu yu zu`, `xsu ysu zsu`, or the wrapped
coordinates moved by the image flags `ix iy iz`. Orthogonal and restricted
triclinic boxes are both handled:

```python
from polyphys.simulate.dump import iter_dump

for frame in iter_dump("whole.bug.lammpstrj.gz"):
    unwrapped = frame.positions(unwrapped=True)
```

A frame cut off by a killed run raises `BrokenDumpError` once the complete
frames before it have been yielded. `read_dump_frame(path, index)` returns a
single frame.

`read_data` and `write_data` read and write LAMMPS data files in the
`atomic`, `bond`, `angle`, `molecular`, and `full` atom styles, as
`DataFile` records of numpy arrays with a `Box`. `with_frame` gives a dump the topology
it lacks: the atom types, masses, bonds, and angles come from a template, and
the positions, box, and image flags from a dump frame, matched by atom ID:

```python
from polyphys.simulate.data import read_data, with_frame, write_data
from polyphys.simulate.dump import read_dump_frame

template = read_data("initial.bug.data")
frame = read_dump_frame("whole.bug.lammpstrj", 0)
write_data(with_frame(template, frame), "whole.bug.data")
```

`write_data` also writes start configurations for LAMMPS. Build a
`DataFile` from arrays, leave `masses` and `velocities` as `None` for the
input script to set, and close a ring with one more bond. Sections of
force-field coefficients are skipped when read and never written, so the
input script sets the force field.

A dump frame and a data file carry the same `Box` from
`polyphys.simulate.core`, which converts scaled coordinates and applies image
flags; `DataFile.unwrapped_positions()` uses it too. The log, dump, and
data-file readers each read a `.gz` file through gzip, and each raises a
subclass of `LammpsFormatError`, itself a `ValueError`, for a broken file, so
one `except LammpsFormatError` skips a broken file of any kind.
