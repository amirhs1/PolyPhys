# Lineage, phases, and stages

This page explains how PolyPhys names and places every artifact of a
parameter sweep: the five levels of the *lineage*, the *phase* directories
the data passes through, and the *stage* directories inside a phase. The
README's "The core idea: artifact lineage" introduces the lineage; this page
is the reference.

## Notation

*artifact*
: One file of the sweep: a trajectory, a log, a probe's time series, an
  ensemble table, or a `space` table.

*lineage*
: The level of an artifact in the hierarchy
  `segment → whole → ensemble_long → ensemble → space`, the values of
  {data}`~polyphys.manage.types.LineageT`.

*state point*
: One set of input parameters of the sweep, such as one chain length, one
  cylinder diameter, and one number of crowders.

*replicate*
: One of the independent runs of a state point, started from its own random
  seed; the keyword `ens` numbers them in a filename.

*phase*
: A top-level directory of a sweep that holds one kind of data, such as
  trajectories or analysis tables, the values of
  {data}`~polyphys.manage.types.PhaseT`.

*stage*
: A directory inside a phase that holds the artifacts of one lineage level
  or one step of the reduction, the values of
  {data}`~polyphys.manage.types.StageT`.

*group*
: The particles a trajectory holds, such as `bug` (the chain) or `all`
  (every particle); a parser class lists its groups.

## The lineage

Each level aggregates the one before it:

| Level           | One artifact is                                         | Its name ends at        |
| --------------- | ------------------------------------------------------- | ----------------------- |
| `segment`       | a chunk of one run, as LAMMPS wrote it                  | the segment, `j03`      |
| `whole`         | one run, its segments joined                            | the replicate, `ens1`   |
| `ensemble_long` | one state point, every input parameter in the name      | the last input keyword  |
| `ensemble`      | one state point, its replicates collapsed               | the short state name    |
| `space`         | the parameter sweep: the state points of one study      | the space name          |

A `segment` is a chunk of a simulation run, not a segment of the polymer
chain, and an `ensemble` is a set of replicates of one state point, not a
statistical-mechanics ensemble such as NVT.

A parser, a subclass of {class}`~polyphys.manage.parser.ParserBase` that a
study declares, reads an artifact's name and
gives the name of the artifact at its own level and at every level above it,
so a file knows its replicate, its state point, and its sweep without
touching the file system:

```python
>>> from polyphys.manage.parser import KeywordSpec, ParserBase, Species
>>> long_keys = ('N', 'r', 'lz', 'ac', 'nc', 'dt', 'bdump', 'adump')
>>> bulk = ['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c']
>>> class ChainCyl(ParserBase):
...     """A linear chain among crowders in a cylinder, periodic along z."""
...     _geometry = 'cylindrical'
...     _topology = 'linear'
...     _groups = ['bug', 'all']
...     _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
...     _keywords = {
...         'N': KeywordSpec('nmon', int),
...         'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
...         'D': KeywordSpec('dcyl', float),
...         'lz': KeywordSpec('lcyl', float, scale=2.0),
...         'ac': KeywordSpec('dcrowd', float),
...         'nc': KeywordSpec('ncrowd', int),
...         'dt': KeywordSpec('dt', float),
...         'bdump': KeywordSpec('bdump', int),
...         'adump': KeywordSpec('adump', int),
...         'ens': KeywordSpec('ensemble_id', int),
...         'j': KeywordSpec('segment_id', int),
...     }
...     _patterns = {
...         'segment': (*long_keys, 'ens', 'j'),
...         'whole': (*long_keys, 'ens'),
...         'ensemble_long': long_keys,
...         'ensemble': ('N', 'D', 'ac', 'nc'),
...         'space': ('N', 'D', 'ac'),
...     }
...     _defaults = {'dmon': 1}
...     _species = (Species('m', 'nmon', 'dmon'),
...                 Species('c', 'ncrowd', 'dcrowd'))
...     _project_attributes = {
...         'segment': bulk, 'whole': bulk, 'ensemble_long': bulk,
...         'ensemble': ['dmon'], 'space': ['dmon']}
>>> artifact = ChainCyl(
...     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1.j03',
...     lineage='segment', group='all')
>>> for level in artifact.lineage_genealogy:
...     print(f'{level:13} {getattr(artifact, level)}')
segment       N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1.j03
whole         N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1
ensemble_long N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000
ensemble      N200D20ac2nc1000
space         N200D20ac2

```

The long names, `segment`, `whole`, and `ensemble_long`, are cut from the
parsed name; the short names, `ensemble` and `space`, are built from the
parsed values, so `r10.5` becomes the cylinder diameter `D20` and `ac2.0`
the crowder size `ac2`. "Lineage names" below gives the rules that fix their
text, and the [Studies](studies.md) page shows how a study declares which
keywords each level's name holds.

The organizer and the pipeline use these names as join keys: the segments of
one whole share its `whole` name, the replicates of one state point share
their `ensemble_long` name, and the state points of one sweep share their
`space` name.

### Lineage names

A parser sets the artifact's name at its own lineage and at each parent
lineage, as the attributes `segment`, `whole`, `ensemble_long`, `ensemble`,
and `space`; the organizer uses them as join keys. Two rules fix their text.

- The long names, `segment`, `whole`, and `ensemble_long`, are cut from the
  parsed name, so they keep its text, and they carry no topology suffix. A
  segment `…ens1.j03.linear` and a whole `…ens1.linear` of one simulation
  both give the `whole` name `…ens1`.
- The short names, `ensemble` and `space`, are built from the parsed values.
  A float that is a whole number prints as an integer, and any other float
  prints as Python's `repr`, the shortest text that reads back as the same
  float. The `space` above is `N200D20ac2` although the name wrote `ac2.0`
  and `r10.5`, while a crowder size of `1.5` stays `ac1.5`.

```python
>>> artifact = ChainCyl(
...     'N200r10.5lz25ac1.5nc0dt0.005bdump2000adump5000ens1.linear', 'whole',
...     'bug')
>>> artifact.whole, artifact.space
('N200r10.5lz25ac1.5nc0dt0.005bdump2000adump5000ens1', 'N200D20ac1.5')

```

Names written by the research code that preceded PolyPhys can differ from
these in two ways: a `whole` it parsed from a whole file keeps `.ring` or
`.linear`, and a short name can print `ac2.0`. Parsing such a short name at
its lineage and reading the rebuilt name gives the current form.

### Values a parser derives

Geometry is derived, not just read. The cylindrical wall is a layer of beads
of size `a_w = 1` whose centers lie on the cylinder of radius `r` set in
LAMMPS, so the cylinder size is `D = 2r - a_w`, stored as `dcyl`; `lz` is the
half-length, so the cylinder length is `Z = 2 lz`, stored as `lcyl`. For a
species `α` of size `a_α`, particle centers are confined to
`D_eff,α = D - a_α`, and its bulk number density and volume fraction both use
the effective volume `V_eff,α = π Z D_eff,α² / 4`:

```python
>>> artifact = ChainCyl(
...     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1.j03',
...     lineage='segment', group='all')
>>> float(artifact.dcyl)              # from 'r10.5': 2 * 10.5 - 1.0
20.0
>>> round(float(artifact.phi_bulk_c), 4)
0.3292

```

## Phases

A sweep's data lives under one root directory, in one subdirectory per
phase. {data}`~polyphys.manage.types.PhaseT` lists the phases in the order
the data flows through them:

| Phase      | Holds                                                         |
| ---------- | ------------------------------------------------------------- |
| `simsAll`, `simsCont` | simulation inputs                                  |
| `logs`     | LAMMPS log files                                              |
| `trjs`     | LAMMPS trajectories                                           |
| `probes`   | the per-frame observables and stamps that the probe writes    |
| `analysis` | the reduction of the probes down the lineage, and the `space` tables |
| `viz`      | data prepared for figures                                     |
| `galaxy`   | tables that combine spaces                                    |
| `allInOne` | a study's tables of one row per whole, such as the equilibrium values that {func}`~polyphys.analyze.pipeline.equilibrium_wholes` returns |

PolyPhys reads the files of the `logs` and `trjs` phases
({mod}`polyphys.simulate`, {mod}`polyphys.probe`) and writes the `analysis`
phase ({mod}`polyphys.analyze.pipeline`, {mod}`polyphys.manage.organizer`,
{mod}`polyphys.analyze.sumrule`). A study's driver decides where the probe
writes, and no function of PolyPhys writes the `viz`, `galaxy`, or
`allInOne` phase yet. The [Workflow](workflow.md) page follows the data from
`trjs` to `space`.

## Stages

A phase holds one directory per particle group and stage, named
`<prefix>-<group>-<stage>`; in the examples on these pages, `<prefix>` is
the space name. The probe's output directories are the exception: a study's
driver names them, and the examples on these pages end them in `-probe`,
which is not a stage. {func}`~polyphys.manage.organizer.make_database`
reads only their prefix.
{data}`~polyphys.manage.types.StageT` lists the stages:

| Stage      | Holds                                                          | Written by                                   |
| ---------- | -------------------------------------------------------------- | -------------------------------------------- |
| `segment`  | segment artifacts                                              | —                                            |
| `wholeSim` | whole artifacts joined from their segments                     | {func}`~polyphys.analyze.pipeline.reduce_space` |
| `ens`      | the replicates of one state point side by side, one file per property | {func}`~polyphys.analyze.pipeline.reduce_space` |
| `ensAvg`   | the mean, variance, and standard error across the replicates   | {func}`~polyphys.analyze.pipeline.reduce_space` |
| `space`    | one table per property for the whole sweep                     | {func}`~polyphys.manage.organizer.space_tseries`, {func}`~polyphys.manage.organizer.space_hists`, {func}`~polyphys.analyze.sumrule.space_profiles` |
| `galaxy`   | tables that combine spaces                                     | —                                            |

Three words belong to two vocabularies: `segment` and `space` are lineage
levels and stages, and `galaxy` is a phase and a stage, not a lineage level.
Where a word appears says which it is: as a `lineage` argument or the
lineage of a name it is a level; as a directory directly under the sweep's
root, a phase; as the last field of a directory inside a phase, a stage.

{func}`~polyphys.manage.organizer.make_database` builds a stage directory
from another directory of the same sweep: it keeps the parent of the old
phase directory, puts the new directory in the phase it is given, and takes
the prefix from the old directory's name, up to its first hyphen:

```python
>>> import tempfile
>>> from pathlib import Path
>>> from polyphys.manage.organizer import make_database
>>> tmp = tempfile.TemporaryDirectory()
>>> probe = Path(tmp.name) / 'probes' / 'N200D20ac2-bug-probe'
>>> ens_avg = make_database(probe, 'analysis', 'ensAvg', 'bug')
>>> Path(ens_avg).relative_to(tmp.name).as_posix()
'analysis/N200D20ac2-bug-ensAvg'
>>> tmp.cleanup()

```

The old directory must sit directly in a phase directory, so that the new
one lands at the same depth. A phase or stage outside the two vocabularies,
a group that contains `-`, or an old directory whose parent is not a phase
raises `ValueError`, so a misspelled or misplaced name fails before
anything is written.

## File names inside a stage

A file in a stage directory is named for its lineage name, its group, and
its property, joined by hyphens, as
{func}`~polyphys.manage.organizer.create_fullname` builds them:

| Stage      | File                                                          |
| ---------- | ------------------------------------------------------------- |
| `wholeSim` | `<whole>-<group>-<property>.npy`                              |
| `ens`      | `<ensemble_long>-<group>-<property>.csv` or `.npy`            |
| `ensAvg`   | `<ensemble_long>-<group>-<property>-ensAvg.csv`, or `-ensAvg-<statistic>.npy` per statistic |
| `space`    | `<space>-<group>-<property>-space.csv`                        |

A property is one observable of one species, such as `gyrTMon`, the radius
of gyration (`gyrT`, a time series) of the monomers (`Mon`). The
[Quickstart](quickstart.md) writes one file of each kind.

A file name carries the lineage name of its own level. It ends in a stage
word only where two stages hold files of one level, as `ens` and `ensAvg`
both hold `ensemble_long` files, or where the level is `space`.

{func}`~polyphys.analyze.pipeline.reduce_stamps` writes one stamps table
into each of the three reduction stages, named
`<space>-<group>-stamps-<level>.csv` for the lineage level of its rows:

| Stage      | Stamps table                                  | One row per    |
| ---------- | --------------------------------------------- | -------------- |
| `wholeSim` | `<space>-<group>-stamps-segment.csv`          | segment        |
| `ens`      | `<space>-<group>-stamps-whole.csv`            | whole          |
| `ensAvg`   | `<space>-<group>-stamps-ensemble_long.csv`    | state point    |

The `wholeSim` table exists only when the probe ran on segments.

### The grammar of a name

- Fields are joined by `-`: `<lineage name>-<group>-<property>`, then the
  stage word, if any, then the extension.
- A lineage name never contains `-`: a parser cuts a name at its first
  hyphen. A group cannot contain one either, which
  {func}`~polyphys.manage.organizer.make_database` checks.
- A property may contain `-`, which joins an observable to a measure of it,
  as in `gyrTMon-acf`. A property is therefore everything between the group
  and the stage word or the extension, not one hyphen-separated field.
- Inside a lineage name, a segment joins its `j<NN>` keyword with a dot, as
  in `...ens1.j03`, as the simulation writes it, and so do a topology suffix
  such as `.ring` and the file extension. A file is grouped on its full name
  with the matched extension removed, never on the part before its last dot,
  so `...epshc1.0nc0` and `...epshc1.0nc500` stay apart.
