# Workflow: from trajectories to a space

This page follows one study's data from its LAMMPS trajectories to the
`space` tables of its sweep: which function runs at each step, what it reads,
and what it writes. Each step has its own page; this one shows how they fit
together and describes the last step, the `space` stage, which has no page
of its own.

## Notation

*phase*, *stage*
: The directories of a sweep, such as `probes` and `analysis`, and the
  directories inside a phase, such as `ensAvg`; the [Lineage](lineage.md)
  page defines them.

*whole*, *ensemble*, *space*
: One run, the replicates of one state point, and the sweep; the
  [Lineage](lineage.md) page defines them with the other lineage levels.

*property*
: One observable of one species, such as `gyrTMon`, the radius of gyration
  time series (`gyrT`) of the monomers (`Mon`).

`SPEC`, `PROBE`, `ANALYSIS`
: A study's `ProjectSpec`, `ProbeSpec`, and analysis tables, which the
  [Studies](studies.md) page shows how to build.

$\phi_c$
: The bulk volume fraction of the crowders, the parser attribute
  `phi_bulk_c`.

$\Delta$
: The rounding step of $\phi_c$, the spec's `divisor`.

## The steps

| Step | Reads                                | Runs                                              | Writes                                   | Page |
| ---- | ------------------------------------ | ------------------------------------------------- | ---------------------------------------- | ---- |
| 1    | a whole or segment in `trjs`, with its data file | an analysis of {mod}`polyphys.probe`, built from `PROBE` | one `.npy` file per observable and a stamps file | [Probe](probe.md) |
| 2    | the probe files of one space         | {func}`~polyphys.analyze.pipeline.reduce_space` with one table of `ANALYSIS` | `analysis/<space>-<group>-ens` and `-ensAvg`, and `-wholeSim` for a `segment` table | [Pipeline](pipeline.md) |
| 3    | the `ensAvg` stage                   | {func}`~polyphys.manage.organizer.space_tseries`, {func}`~polyphys.manage.organizer.space_hists`, {func}`~polyphys.analyze.sumrule.space_profiles` | `analysis/<space>-<group>-space` | below |
| 4    | the `ens` stage and the stamps       | {func}`~polyphys.analyze.pipeline.equilibrium_wholes` | a table of one row per whole, returned | [Pipeline](pipeline.md) |

Step 1 is the only step that reads trajectories and the only one that needs
MDAnalysis, the `mda` extra; steps 2 to 4 read only the files the probe
wrote. The [Quickstart](quickstart.md) runs steps 2 and 3 on a synthetic
sweep.

### Step 1: probe

For each whole, a study's driver opens the trajectory and its data file as
one MDAnalysis universe with {func}`~polyphys.probe.core.make_universe`,
runs the analyses its `PROBE` asks for, saves their arrays, and writes the
whole's stamps, one row of its lineage names, its parser attributes, and
its number of frames, with {func}`~polyphys.probe.core.write_stamps`.

### Step 2: reduce down the lineage

{func}`~polyphys.analyze.pipeline.reduce_space` takes the probe directory of
one space, the study's `SPEC`, and one `AnalysisTable` of its `ANALYSIS`. For
each property of the table, it joins the segments of each whole, groups the
wholes of each state point into an ensemble, and averages each ensemble over
its replicates, writing each level to its stage. A study's folder is not an
installed package, so its driver puts `projects/` on the import path, as
[`projects/README.md`](https://github.com/amirhs1/PolyPhys/tree/main/projects)
shows. Each table reads the probe directory of its own group:

```python
import sys
sys.path.insert(0, "<path to PolyPhys>/projects")

from polyphys.analyze.pipeline import reduce_space
from trans_foci_cub import ANALYSIS, SPEC

for table in ANALYSIS:
    reduce_space(f'probes/ns400nl5al5ac1-{table.group}-probe', SPEC, table,
                 density=kernel)
```

`kernel` is the study's density kernel, which the `all` table's
`LocalDensities` hooks need; "Building the density kernel" on the [Pipeline](pipeline.md) page shows how to build it.

The ensemble average holds `<property>-mean`, `<property>-var`, and
`<property>-sem` per frame or bin. The wholes of an ensemble are independent
runs, so `-sem` is the plain standard error across them.

### Step 3: the space stage

Each `space_*` function reads every `ensAvg` file of one property in a
space, adds the state point's attributes as columns, and stacks the files
into one table, the input of a figure or a fit. With `is_save=True` it
writes the table to `analysis/<space>-<group>-space`. The three functions
differ in what a row is:

{func}`~polyphys.manage.organizer.space_tseries`
: A frame of a time series. It adds `t_index`, the frame index times the
  dump period of the group's trajectory, in time steps, and `time`,
  `t_index` times the time step `dt`.

{func}`~polyphys.manage.organizer.space_hists`
: A bin of a histogram. Its `bin_center` column is the one the `ensAvg`
  file carries, or the array passed as `bin_center`, which also adds
  `bin_center-norm`; with `normalize=True` it adds `<property>-norm`, the
  histogram divided by its sum.

{func}`~polyphys.analyze.sumrule.space_profiles`
: A bin of a density profile, scaled and normalized for the density sum
  rule, as the [Usage](usage.md) page, "The density sum rule", describes.

All three take the study's parser, its `hierarchy` glob, its
`physical_attrs`, its `divisor`, and its crowder attribute.
{func}`~polyphys.analyze.pipeline.space_tseries_from_spec`,
{func}`~polyphys.analyze.pipeline.space_hists_from_spec`, and
{func}`~polyphys.analyze.sumrule.space_profiles_from_spec` read them, and
the species sizes of the profiles, from the study's `SPEC`, so no study
choice comes from a library default:

```python
from polyphys.analyze.pipeline import space_tseries_from_spec
from polyphys.analyze.sumrule import space_profiles_from_spec
from trans_foci_cub import SPEC

table = space_tseries_from_spec(
    'analysis/ns400nl5al5ac1-bug-ensAvg', SPEC, 'gyrTMon',
    save=True)
profiles = space_profiles_from_spec(
    'analysis/ns400nl5al5ac1-all-ensAvg', SPEC, 'Phi', 'r', 'Crd',
    save=True)
```

The density profiles are in the `'all'` group, the default of
`space_profiles_from_spec`, not in the polymer group `SPEC.group`, which the
other two take by default.

Each also adds the column named by the spec's `crowder_attr`, such as
`phi_bulk_c`, and its rounded copy,
`phi_bulk_c_round`: $\phi_c$ rounded to the nearest multiple of $\Delta$ by
{func}`~polyphys.manage.utils.round_to_multiple`, then to `round_to`
decimal places. The unrounded fraction stays in the table beside it.

### Step 4: one value per whole

A time series keeps every frame through step 2. For tables of one number
per whole, {func}`~polyphys.analyze.pipeline.equilibrium_wholes` reduces each
whole's series in the `ens` stage to its time average, its variance, and the
blocked standard error of its time average, and merges them with the whole
stamps, the table `<space>-<group>-stamps-whole.csv` that step 2 writes into
the `ens` stage when the table has `stamps=True`.
{func}`~polyphys.analyze.pipeline.reduce_stamps` writes it, combining the
probe's stamps files with {func}`~polyphys.manage.organizer.children_stamps`
and rolling them up a level with
{func}`~polyphys.manage.organizer.parents_stamps`; the [Lineage](lineage.md)
page lists the stamps tables of each stage. The frames of one run are correlated, so this error comes from block
averaging within the whole, not from the per-frame standard error; the
[Statistics](statistics.md) page explains why and how.

## Which error bar at which level

| Level            | Average over                 | Error                                    |
| ---------------- | ---------------------------- | ---------------------------------------- |
| one whole        | its frames, correlated       | blocked standard error, {func}`~polyphys.analyze.stats.sem_blocked` |
| one ensemble     | its replicates, independent  | plain standard error, the `-sem` columns |
| a normalized profile | its replicates, each normalized | `<property>-norm-sem`, see [Usage](usage.md) |
