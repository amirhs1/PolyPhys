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

A parser class of {mod}`polyphys.manage.parser` reads an artifact's name and
gives the name of the artifact at its own level and at every level above it,
so a file knows its replicate, its state point, and its sweep without
touching the file system:

```python
>>> from polyphys.manage.parser import SumRuleCyl
>>> artifact = SumRuleCyl(
...     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1.j03',
...     lineage='segment', group='all')
>>> for level in artifact.lineage_genealogy:
...     print(f'{level:13} {getattr(artifact, level)}')
segment       N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1.j03
whole         N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1
ensemble_long N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000
ensemble      N200D20ac2nc1000
space         N200D20ac2

```

The long names, `segment`, `whole`, and `ensemble_long`, are cut from the
parsed name; the short names, `ensemble` and `space`, are built from the
parsed values, so `r10.5` becomes the cylinder diameter `D20` and `sig2.0`
the crowder size `ac2`. The [Usage](usage.md) page, "Lineage names", gives
the rules that fix their text, and the [Studies](studies.md) page shows how a
study declares which keywords each level's name holds.

The organizer and the pipeline use these names as join keys: the segments of
one whole share its `whole` name, the replicates of one state point share
their `ensemble_long` name, and the state points of one sweep share their
`space` name.

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
the space name.
{data}`~polyphys.manage.types.StageT` lists the stages:

| Stage      | Holds                                                          | Written by                                   |
| ---------- | -------------------------------------------------------------- | -------------------------------------------- |
| `segment`  | segment artifacts                                              | —                                            |
| `wholeSim` | whole artifacts joined from their segments                     | {func}`~polyphys.analyze.pipeline.reduce_space` |
| `ens`      | the replicates of one state point side by side, one file per property | {func}`~polyphys.analyze.pipeline.reduce_space` |
| `ensAvg`   | the mean, variance, and standard error across the replicates   | {func}`~polyphys.analyze.pipeline.reduce_space` |
| `space`    | one table per property for the whole sweep                     | {func}`~polyphys.manage.organizer.space_tseries`, {func}`~polyphys.manage.organizer.space_hists`, {func}`~polyphys.analyze.sumrule.space_profiles` |
| `galaxy`   | tables that combine spaces                                     | —                                            |

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

A phase or stage outside the two vocabularies raises `ValueError`, so a
misspelt directory name fails before anything is written.

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
