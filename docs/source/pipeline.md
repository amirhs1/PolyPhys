# Pipeline: reducing a space down the lineage

This page explains what {mod}`polyphys.analyze.pipeline` does with the
probe artifacts of one space, how a study describes its analysis as data,
what the post-hooks derive, where the outputs go, and what changed from the
research code it replaces.

## Notation

*segment*, *whole*, *ensemble*, *ensemble average*
: The levels of the lineage, from one chunk of one run up to the replicates
  of one state point collapsed into statistics. The README's "The core idea:
  artifact lineage" defines them; the pipeline walks them in that order.

*property*
: One probe artifact without its lineage name and group, such as `gyrTMon`:
  a *name* (`gyrT`) and a *species* (`Mon`). Its files end in
  `-gyrTMon.npy` (a whole) or `.gyrTMon.npy` (a segment).

*relation*
: How the segments of a whole combine: `'tseries'` concatenates along the
  frames, `'histogram'` sums, `'bin_edge'` keeps one copy of identical
  edges.

*whole type*
: The shape of one whole: `'vector'`, `'matrix'`, `'bin_edge'`, or
  `'dataframe'`. It decides how the wholes stack into an ensemble and what
  the ensemble average looks like.

*post-hook*
: A callable `(wholes, context) -> tuple[Derived, ...]` applied to a
  property's wholes before the ensembles are built. It returns every
  artifact to reduce.

$r_k$
: The autocorrelation at lag $k$ frames, from
  {func}`~polyphys.analyze.correlations.acf`.

$\bar P(R)$
: The probability density of the separation vector of a pair of foci, in
  one over the length unit cubed, from
  {func}`~polyphys.analyze.clusters.pair_separation_density`.

## A study's analysis as data

A study declares one {class}`~polyphys.analyze.pipeline.AnalysisTable` per
particle group: the group, the *lineage* of its probe artifacts (`'segment'`
when the probe phase wrote them in chunks, `'whole'` when it wrote one file
per run), whether its stamp files are combined, and a tuple of
{class}`~polyphys.analyze.pipeline.PropertySpec` records. Five constructors
build the kinds the studies use:

| Constructor                      | Property           | Relation      | Whole type   | Ensemble average                    |
| -------------------------------- | ------------------ | ------------- | ------------ | ----------------------------------- |
| `tseries('gyrT', 'Mon')`         | `gyrTMon`          | `'tseries'`   | `'vector'`   | `-mean`, `-var`, `-sem` per row     |
| `histogram('r', 'Crd')`          | `rHistCrd` + edges | `'histogram'` | `'vector'`   | the same, with `bin_center`         |
| `histogram_2d('xy', 'Crd')`      | `xyHistCrd`        | `'histogram'` | `'matrix'`   | one array per statistic             |
| `edge('x')`                      | `xEdge`            | `'bin_edge'`  | `'bin_edge'` | the edges                           |
| `matrix_tseries('principalT', 'Mon')` | `principalTMon` | `'tseries'` | `'matrix'`   | one array per statistic             |

A table lists each artifact once; listing one twice is an error, and so is
a property whose files are not found. The tables of the maintainer's
studies live beside each study's spec in the repository's `projects/`
folder, so a data correction or a special case is code in the study's
folder, never an option of the module.

```python
from polyphys.analyze.pipeline import (
    AnalysisTable, Autocorrelation, LocalDensities, edge, histogram,
    histogram_2d, tseries)

ANALYSIS = (
    AnalysisTable(group='bug', lineage='whole', stamps=True, properties=(
        tseries('gyrT', 'Mon', post_hook=Autocorrelation(nlags=7000)),
    )),
    AnalysisTable(group='all', lineage='segment', stamps=False, properties=(
        histogram('r', 'Crd', post_hook=LocalDensities()),
        histogram_2d('xy', 'Crd'),
        edge('x'),
    )),
)
```

## The post-hooks

A hook sees a property's wholes and a
{class}`~polyphys.analyze.pipeline.HookContext`: the study's spec, the
group, the parser instance of each whole (which holds its sizes and
counts), the bin edges, and the injected density kernel. It returns
{class}`~polyphys.analyze.pipeline.Derived` records, each an output
property with its wholes and whole type. Without a hook, the wholes pass
through under the property's own name.

{class}`~polyphys.analyze.pipeline.TimeAverage`
: Averages a per-frame array over the frames: a `(n_frames, n_sizes)`
  cluster-size histogram time series becomes one `n_sizes` vector per
  whole. The output name is given to the hook; `drop_first_column` removes
  the count of empty clusters that
  {func}`~polyphys.analyze.clusters.count_clusters` keeps in column 0.

{class}`~polyphys.analyze.pipeline.Autocorrelation`
: Passes a time series through and adds $r_k$ for $k = 0, \dots,$ `nlags`
  under `<property>-acf`. No confidence band is stored: Bartlett's band is
  a function of the stored $r_k$, the number of frames, and the level, so
  {func}`~polyphys.analyze.correlations.bartlett_half_width` rebuilds it at
  any level when it is drawn.

{class}`~polyphys.analyze.pipeline.LocalDensities`
: Passes a histogram through and adds `<direction>Rho<species>` and
  `<direction>Phi<species>`, the local number density and volume fraction,
  through the {class}`~polyphys.analyze.pipeline.DensityKernel` given to
  {func}`~polyphys.analyze.pipeline.reduce_space`. The species' diameter is
  read from the attribute the spec's `species_sizes` names. The module
  imports no density model; the study's driver supplies one.

{class}`~polyphys.analyze.pipeline.FociPairDistances`
: Turns the `(n_frames, n_foci, n_foci)` distance-matrix time series of
  the foci into three DataFrame wholes with one column per pair, tagged
  `loc<i>loc<j>genDist<d>` with $d$ the genomic distance in bonds: the
  histogram of each pair's separation, $\bar P(R)$, and the separation time
  series. Every pair has one separation per frame, so each histogram column
  sums to `n_frames` and $\sum_k \bar P_k \Delta V_k = 1$ exactly, provided
  the grid covers every separation. The grid's upper end comes from the
  system, not the data: {func}`~polyphys.analyze.pipeline.half_contour`
  gives the largest separation a ring allows, half its bonds times the
  largest bond extension, $R_0 = 1.5\,\sigma_{ij}$ for the FENE bond
  (Kremer and Grest, 1990). A separation above the grid raises and names
  the whole, since the chain cannot reach it.

## Running a space

```python
from polyphys.analyze.pipeline import reduce_space
from trans_foci_cub.spec import ANALYSIS, SPEC

for table in ANALYSIS:
    reduce_space('probes/al5.0nl5ml125ns400ac1.0-bug-probe', SPEC, table,
                 density=my_density_kernel)
```

`reduce_space` lists the files under `<input_database>/<hierarchy>` (the
default hierarchy is the spec's `dir/file` layout), combines the stamps if
the table asks for it, and then, one property at a time: selects the
property's files by their last name token, builds the wholes (segments are
combined and written to `wholeSim`), loads and checks the edges, applies
the hook, and builds the ensembles and ensemble averages. It returns the
ensemble averages by output property and ensemble name; `save=False` skips
the writing.

| Stage      | Directory                              | File                                            |
| ---------- | -------------------------------------- | ----------------------------------------------- |
| `wholeSim` | `<root>/analysis/<prefix>-<group>-wholeSim` | `<whole>-<group>-<property>.npy` (segment lineage only) |
| `ens`      | `<root>/analysis/<prefix>-<group>-ens`  | `<ensemble_long>-<group>-<property>.csv` or `.npy` |
| `ensAvg`   | `<root>/analysis/<prefix>-<group>-ensAvg` | `<ensemble_long>-<group>-<property>-ensAvg.csv`, or `-ensAvg-<statistic>.npy` per statistic |

The ensemble of DataFrame wholes is their long form, every whole's rows
with a `whole` and a `row` column, so no whole is lost before the ensemble
average; the average is taken row by row across the wholes and carries
`<column>-mean`, `<column>-var`, and `<column>-sem` per column.

## What changed from the research code

The module replaces the research code's `analyzer`, whose nine
near-identical blocks differed only in the six facts a
{class}`~polyphys.analyze.pipeline.PropertySpec` holds. The outputs match
it except where the change was deliberate:

- A property with no file raises; the research code wrote nothing and said
  nothing.
- The wholes of one ensemble must hold identical bin edges; the research
  code took the first whole's.
- The autocorrelation is stored alone; the `acfLowerCi` and `acfUpperCi`
  files are gone, since the band is a function of the stored values.
- The time average of a per-frame histogram is named by the hook, not by
  dropping every `T` from the property's name.
- The foci pair histograms share one grid up to a bound from the system;
  the research code set each whole's grid from its own largest separation,
  so the ensemble mean mixed different sets of wholes in different bins.
  `pairDistRdf` is now $\bar P(R)$ rather than a unit-sum curve, and the
  `genDist` tags count bonds, one more than the research code's monomers
  between.
- The ensemble of DataFrame wholes keeps every whole; the research code
  averaged them at the ensemble step, so the ensemble average had no
  variance or standard error.

## References

- K. Kremer and G. S. Grest, "Dynamics of entangled linear polymer melts:
  A molecular-dynamics simulation", J. Chem. Phys. 92, 5057 (1990).
  <https://doi.org/10.1063/1.458541>
- M. S. Bartlett, "On the theoretical specification and sampling
  properties of autocorrelated time-series", Supplement to the Journal of
  the Royal Statistical Society 8, 27-41 (1946).
  <https://doi.org/10.2307/2983611>
