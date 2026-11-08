# Usage

This page says where each task is explained. Install PolyPhys first, as
{ref}`Installation <installation_target>` shows; the [quickstart](quickstart.md)
then runs one complete synthetic sweep, from the files of its segments to its
`space` table.

## Tasks

| To | Read |
| --- | --- |
| Install PolyPhys, with pip, uv, or conda, or from a source checkout | [Installation](installation.md) |
| Run one sweep from its segments to its `space` table | [Quickstart](quickstart.md) |
| Decode a file name into its physical attributes and lineage names | [Lineage](lineage.md) |
| Go from trajectories to a `space` table, step by step | [Workflow](workflow.md) |
| Declare a study: its parser, its `ProjectSpec`, and its `ProbeSpec` | [Studies](studies.md) |
| Measure chain sizes, spatial histograms, pair distances, or binders from trajectories | [Probe](probe.md) |
| Reduce a space's probe artifacts down the lineage | [Pipeline](pipeline.md) |
| Estimate the error of a time average over correlated frames | [Statistics](statistics.md) |
| Count contacts, clusters, and binding modes | [Clusters](clusters.md) |
| Compute autocorrelations, decay fits, and persistence lengths | [Correlations](correlations.md) |
| Turn histograms into density profiles, and check the density sum rule | [Density](density.md) |
| Read LAMMPS logs, dumps, and data files, or build a start configuration | [Simulate](simulate.md) |
| Color a figure the way the studies' figures are colored | [Viz](viz.md) |

The [API reference](api.rst) documents every public name; the
[README](https://github.com/amirhs1/PolyPhys) introduces the
artifact-lineage model.

## Measures of one frame

{mod}`polyphys.analyze.measurer` holds the small functions the probes and
the analysis build on, each on the positions or values it is given. The
[API reference](api.rst) documents each one; this table says which to reach
for:

| Kind | Functions |
| --- | --- |
| Chain size, on one frame's positions | {func}`~polyphys.analyze.measurer.end_to_end`, {func}`~polyphys.analyze.measurer.transverse_size`, {func}`~polyphys.analyze.measurer.fsd`, {func}`~polyphys.analyze.measurer.max_distance` |
| Distances under periodic boundaries | {func}`~polyphys.analyze.measurer.apply_pbc_orthogonal`, {func}`~polyphys.analyze.measurer.pair_distance` |
| Bins | {func}`~polyphys.analyze.measurer.fixedsize_bins`, {func}`~polyphys.analyze.measurer.create_bin_edge_and_hist` |
| Histograms of one frame | {func}`~polyphys.analyze.measurer.radial_histogram` (spherical $r$), {func}`~polyphys.analyze.measurer.radial_cyl_histogram` (cylindrical $r$), {func}`~polyphys.analyze.measurer.axial_histogram`, {func}`~polyphys.analyze.measurer.azimuth_cyl_histogram`, {func}`~polyphys.analyze.measurer.planar_cartesian_histogram` |
| Volumes of a bead cut by planes or a sphere | {func}`~polyphys.analyze.measurer.spherical_segment`, {func}`~polyphys.analyze.measurer.sphere_sphere_intersection` |
| Statistics of independent samples | {func}`~polyphys.analyze.measurer.simple_stats`, {func}`~polyphys.analyze.measurer.sem` |

`simple_stats` and `sem` treat their samples as independent, as the
replicates of a state point are; for the correlated frames of one run, use
the blocked error of {mod}`polyphys.analyze.stats`, as the
[Statistics](statistics.md) page explains.

For example, the end-to-end distance and the span along `z` of three beads:

```python
>>> import numpy as np
>>> from polyphys.analyze import measurer
>>> positions = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [3.0, 4.0, 12.0]])
>>> float(measurer.end_to_end(positions))
13.0
>>> float(measurer.fsd(positions, 2))
12.0

```

## The density sum rule, in brief

{func}`~polyphys.analyze.sumrule.space_profiles` scales each species'
ensemble-averaged profile by its size and normalizes it, giving one row per
bin and state point with the columns `<prop>-scale`, `<prop>-normalizer`,
`<prop>-norm`, and, under the default normalizer, `<prop>-norm-sem`;
{func}`~polyphys.analyze.sumrule.sum_rule_profiles` sums the species into
`<prop>Sum-norm`. The inputs must be raw ensemble averages, and
`<prop>-norm-sem` is not generally an uncertainty estimate for the
normalized ensemble mean. The [Density](density.md) page, "The density sum
rule", gives the rule, the normalizers, and the error columns.
