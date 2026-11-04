# Probe: per-frame observables from LAMMPS trajectories

This page explains what {mod}`polyphys.probe` does with a study's LAMMPS
trajectories, how a study describes its probe as data, what the chain-size,
spatial-histogram, and pair-distance analyses measure, how they run in
serial or in parallel, where their outputs go, and what changed from the
research code they replace.

The probe phase is the one step of PolyPhys that reads trajectories. It is
built on [MDAnalysis](https://www.mdanalysis.org/), which the rest of the
package does not need, so the layer is an optional extra:

```bash
python -m pip install "polyphys[mda] @ git+https://github.com/amirhs1/PolyPhys.git"
```

Without MDAnalysis, `import polyphys.probe` raises `ImportError` with that
hint. The parsing and reduction layers work as before.

## Notation

*whole*, *segment*
: One run's trajectory, or one chunk of it, as the README's "The core
  idea: artifact lineage" defines them. A probe reads one whole or one
  segment and writes one file per observable plus a *stamps* file.

*selection*
: An MDAnalysis selection string, such as `'resid 1'` or `'type 2'`, that
  picks the atoms of one species from the universe.

*observable*
: One quantity measured per frame, named as its time series without the
  species label: `gyrT` for the file `gyrTMon.npy`.

$N$, $\mathbf{r}_i$, $m_i$
: The number of beads in the chain and the position and mass of bead $i$,
  in the dump's length and mass units.

$\mathbf{S}$, $\lambda_1 \ge \lambda_2 \ge \lambda_3$
: The mass-weighted gyration tensor of the chain about its centre of mass,
  and its eigenvalues.

$D$, $Z$, $L$, $a_w$
: The parser's `dcyl`, `lcyl`, and `lcube`, and the size $a_w = 1$ of the
  beads whose centres form a cylindrical wall, so that the LAMMPS box is
  $(D + a_w) \times (D + a_w) \times Z$ for a cylinder and $L^3$ for a cube.

$F$, $b$, $K = \lceil F / b \rceil$
: The frames of a run, the block length of the two-dimensional maps, and
  the number of blocks.

$d_{ij}$, $r_c$, $g_i$
: The distance between particles $i$ and $j$ of the pair species, the
  contact cutoff, and the 0-based backbone index of particle $i$.

## Reading a trajectory

{func}`~polyphys.probe.core.make_universe` opens a study's LAMMPS data file
and dump as one MDAnalysis universe with the four arguments the research
code used in each of its probe functions: the data file as a `DATA`
topology whose `Atoms` lines read `id resid type x y z`, so the molecule ID
of the `bond`, `angle`, or `molecular` atom styles becomes the residue ID
that `'resid 1'` selects; and the dump as a `LAMMPSDUMP` trajectory whose
coordinates are the wrapped, unscaled `x y z` columns. Any further keyword,
such as `dt`, passes through to the universe; the four fixed ones cannot be
overridden.

```python
from polyphys.probe import make_universe

universe = make_universe('N200...ens1.bug.data', 'N200...ens1.bug.lammpstrj')
chain = universe.select_atoms('resid 1')
```

Two facts about the files matter:

- The data file needs a `Masses` section. MDAnalysis weights the radius of
  gyration and the gyration tensor by mass, and without masses it cannot
  form the tensor. The studies' data files have one.
- MDAnalysis sorts the atoms of each dump frame by `id`, so `positions`
  are in atom-ID order, the order along the chain, which the end-to-end
  distance relies on.

This pull request does not convert frames to simulation time: the universe
is opened without `dt`, the dump reader warns that it has no time step, and
the `times` of an analysis are the dump's time steps. The research code
passed `dt` as the dump interval times the time step times the
Lennard-Jones time unit, which issue #56 lists as still to confirm.

## A study's probe as data

A study builds a {class}`~polyphys.manage.projects.ProbeSpec` beside its
{class}`~polyphys.manage.projects.ProjectSpec`, in its own folder, and
reads the parser class and the geometry from the project spec rather than
storing them twice:

```python
from polyphys.manage.projects import ProbeSpec

PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'Foci': 'type 2', 'Crd': 'resid 0',
                'Dna': 'type 1'},
    chain='Mon',
    observables=('asphericityT', 'fsdT', 'gyrT', 'shapeT', 'principalT'),
    unwrap=False,
    bin_sizes={'r': _r_bin, 'z': _z_bin, 'theta': np.pi / 36, 'x': _r_bin,
               'y': _r_bin},
    histogram_species=('Crd', 'Mon', 'Dna', 'Foci'),
    block_frames=100,
    recentre='half',
    pair_species='Foci',
    contact_cutoff=_contact_cutoff,
    minimum_image=True,
)
```

| Field               | Holds                                                                                                                                        |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `selections`        | Each species label mapped to its selection, read-only, in order                                                                              |
| `chain`             | The label in `selections` of the polymer chain; its series are named `<observable><chain>`                                                   |
| `observables`       | The chain-size observables to measure, each a name in `CHAIN_SIZE_OBSERVABLES`                                                              |
| `unwrap`            | Whether the chain is unwrapped across the periodic boundaries for its asphericity: `False` in a cylinder, `True` in a cube                   |
| `bin_sizes`         | The histogram bin size along each direction: a number, or a function of the whole's parser for a size that scales with its diameters (`PerWhole`) |
| `histogram_species` | The species whose spatial histograms are collected, each a label in `selections`; empty when the study's histograms are not probed           |
| `block_frames`      | The block length $b$ of the two-dimensional maps; required with `histogram_species`                                                          |
| `recentre`          | How a frame's positions are moved before binning: `'half'` (half the box), `'none'`, or a triple (`RecentreT`); see "Recentring a frame"     |
| `pair_species`      | The species whose pair distances and contacts are collected, such as the foci; `None` when the study has none                                |
| `contact_cutoff`    | The centre-to-centre distance at or below which two particles of `pair_species` touch, a number or a `PerWhole` function; required with it   |
| `minimum_image`     | Whether the pair distances use the minimum image along the periodic directions: the axis of a cylinder, every direction of a cube            |

The spec is checked at construction: the chain label must be a selection,
every observable must be known and listed once, every scalar observable
with the chain label must be one of the project's `tseries_properties`
(`principalT`, a per-frame matrix, is exempt), every bin size and the
contact cutoff must be a finite positive number or callable along a known
direction, every histogram species and the pair species must be a
selection, and the histogram species need `bin_sizes` and `block_frames`
as the pair species needs `contact_cutoff`. Every field after `unwrap` has
a default, so a record that probes only the chain size is unchanged.
{meth}`~polyphys.manage.projects.ProbeSpec.resolve` turns a `PerWhole`
value into the number for one whole, given its parser.

The `unwrap` value is the research code's: its probes read wrapped dumps,
and a cylinder is periodic along its axis only, so the chain was unwrapped
for the asphericity in the cubic studies and not in the cylindrical ones.
Whether the two conventions differ on a given study's frames is a check on
the data; switching a study is a spec edit, not a code change.

## The chain-size analysis

{class}`~polyphys.probe.chain_size.ChainSize` is an MDAnalysis
{class}`~MDAnalysis.analysis.base.AnalysisBase`. It measures the
observables of a `ProbeSpec`, or any sequence of names from
{data}`~polyphys.manage.projects.CHAIN_SIZE_OBSERVABLES`, on one chain in
every frame:

| Observable     | One frame's value                                                                                                 | From                                                       |
| -------------- | ----------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| `rfloryT`      | $\lvert \mathbf{r}_N - \mathbf{r}_1 \rvert$, the end-to-end distance                                               | {func}`~polyphys.analyze.measurer.end_to_end`              |
| `fsdT`         | $\max_i z_i - \min_i z_i$ along the cylinder axis, the span or Feret diameter                                     | {func}`~polyphys.analyze.measurer.fsd`                     |
| `transSizeT`   | Twice the mean distance from the centre of geometry in the plane across the axis                                  | {func}`~polyphys.analyze.measurer.transverse_size`         |
| `gyrT`         | $R_g = \sqrt{\lambda_1 + \lambda_2 + \lambda_3}$, the mass-weighted radius of gyration                            | `AtomGroup.radius_of_gyration`                             |
| `asphericityT` | $\tfrac{3}{2} \sum_k (\lambda_k - \bar\lambda)^2 / (\sum_k \lambda_k)^2$, from 0 (spherically symmetric) to 1 (a line) | `AtomGroup.asphericity(wrap=False, unwrap=unwrap)`        |
| `shapeT`       | $27 \prod_k (\lambda_k - \bar\lambda) / (\sum_k \lambda_k)^3$, from 0 to 2                                        | `AtomGroup.shape_parameter(wrap=False)`                    |
| `principalT`   | The $3 \times 3$ matrix of principal axes of the moment of inertia, one eigenvector per row, by decreasing eigenvalue | `AtomGroup.principal_axes(wrap=False)`                 |

The asphericity and shape parameter are those of Dima and Thirumalai
(2004), as MDAnalysis implements them. `shapeT` and `principalT` read the
positions as the dump wraps them whatever `unwrap` is, as the research
code did; `unwrap` reaches the asphericity only.

```python
from polyphys.manage.parser import TransFociCyl
from polyphys.probe import ChainSize, write_stamps
from trans_foci_cyl import PROBE

parser = TransFociCyl(trajectory, 'whole', 'bug', ispath=True)
analysis = ChainSize.from_spec(universe, PROBE).run()
analysis.save(save_to, f'{parser.name}-{parser.group}')
write_stamps(parser, analysis.n_frames, save_to)
```

After `run()`, `results` holds one float64 array per observable, keyed by
its name: shape `(n_frames,)` for a scalar series and `(n_frames, 3, 3)`
for `principalT`. `run(start, stop, step)` or `run(frames=...)` chooses
the frames as in any MDAnalysis analysis. `save` writes
`<name>-<observable><chain>.npy` per observable, the research code's
layout, and {func}`~polyphys.probe.core.write_stamps` writes
`<name>-<group>-stamps.csv`: one row of the artifact's lineage names, its
`group`, every parser attribute, and `n_frames`, which
{func}`~polyphys.analyze.pipeline.reduce_stamps` later combines across a
space.

### Parallel runs

The class is marked parallelizable and supports the `serial`,
`multiprocessing`, and `dask` backends of MDAnalysis 2.8 and later:

```python
analysis = ChainSize.from_spec(universe, PROBE).run(
    backend='multiprocessing', n_workers=4)
```

Each worker fills the arrays of its own frames, and the aggregator
concatenates them along the frame axis (`ResultsGroup.ndarray_hstack` for
a scalar series, `ResultsGroup.ndarray_vstack` for `principalT`), so a
parallel run gives the serial arrays.

### Cost

Every observable is a sum over the beads followed by at most a $3 \times 3$
eigenproblem, so each frame costs $O(N)$ time, plus the
$O(N + N_\mathrm{bonds})$ unwrap of the chain when `unwrap` is set, and
$O(N)$ extra space for the recentred positions. Over a run of $F$ frames
the series take $O(F)$ space each and $O(9F)$ for `principalT`. Positions
come from MDAnalysis as float32, so each value carries a relative round-off
of about $10^{-7}$.

## Recentring a frame

The studies' LAMMPS scripts set the box to $[-r, r]^2 \times [-l_z, l_z]$
for a cylinder and $[-l, l]^3$ for a cube, and `fix recenter` keeps the
polymer's centre of mass at $z = 0$ or at the origin. The LAMMPS dump
reader of MDAnalysis 2.4 and later translates every frame so that the lower
box bounds sit at the origin, which undoes that recentring, and no option
of the reader keeps the written frame; only the box lengths survive, as
`ts.dimensions`. A spatial histogram therefore has to move the positions
back, and `ProbeSpec.recentre` says how:

| `recentre`     | Subtracts, per frame                                     | Use                                                                 |
| -------------- | -------------------------------------------------------- | ------------------------------------------------------------------- |
| `'half'`       | Half the frame's box lengths, `ts.dimensions[:3] / 2`    | The default: restores the LAMMPS frame exactly, the axis at $(0, 0)$ |
| `'none'`       | Nothing                                                  | Keep the reader's frame, positions in $[0, L)$                      |
| `(dx, dy, dz)` | The triple as given                                      | Reproduce another convention, such as the research code's           |

The research code subtracted half the parser's $(D, D, Z)$ or $L$. In a
cube that is half the box. In a cylinder the box spans the wall-bead
centres, $D + a_w$ across, so the research shift of $D / 2$ left the
cylinder's axis at $(a_w / 2, a_w / 2) = (0.5, 0.5)$ rather than at the
origin, and the radial and azimuthal histograms of every cylindrical study
were measured about that displaced point. The maintainer chose `'half'` on
2026-11-04; the research numbers remain reproducible with the triple
`(D / 2, D / 2, Z / 2)`. The pair distances are translation invariant and
ignore `recentre`.

## The spatial histograms

{class}`~polyphys.probe.histograms.SpatialHistograms` bins every species of
a `ProbeSpec` in every frame, after recentring, along the directions of the
geometry and in the three Cartesian planes:

| Geometry      | Per frame, int64 `(n_frames, n_bins)` | Per block, float64 `(n_blocks, n_a, n_b)` |
| ------------- | ------------------------------------- | ----------------------------------------- |
| `cylindrical` | `rHistT`, `zHistT`, `thetaHistT`      | `xyHistBlock`, `xzHistBlock`, `yzHistBlock` |
| `cubic`       | `rHistT`                              | `xyHistBlock`, `xzHistBlock`, `yzHistBlock` |

each followed by the species label, with $r = \sqrt{x^2 + y^2}$ and
$\theta = \operatorname{atan2}(y, x)$ about the $z$ axis in a cylinder
and $r = |\mathbf{r}|$ from the box centre in a cube, through the helpers
of {mod}`polyphys.analyze.measurer`. The one-dimensional counts are kept
per frame, as the maintainer decided on 2026-11-03, so the time dependence
of a profile and its block statistics stay available. A map per frame
would not fit in memory, so the maps are averaged over blocks of
`block_frames` consecutive frames of the run, with `blockFrames` holding
the frames in each block; with $b = 100$ and $200 \times 200$ cells the
three planes cost about 1 MB per block per species. The `xz` map is
transposed to rows along $x$, as the research code stored it.

{func}`~polyphys.probe.histograms.histogram_edges` builds a whole's edges
from the spec's bin sizes, resolved for the whole, with
{func}`~polyphys.analyze.measurer.fixedsize_bins`:

| Geometry      | Direction          | Range                                    | Bins          |
| ------------- | ------------------ | ---------------------------------------- | ------------- |
| `cylindrical` | `r`                | $[0, (D + a_w) / 2]$, the wall surface   | `nonnegative` |
| `cylindrical` | `z`                | $[-Z / 2, Z / 2]$                        | `ordinary`    |
| `cylindrical` | `theta`            | $[-\pi, \pi]$                            | `periodic`    |
| `cylindrical` | `x`, `y`           | $[-(D + a_w) / 2, (D + a_w) / 2]$        | `ordinary`    |
| `cubic`       | `r`                | $[0, \sqrt{3} L / 2]$, the corner        | `nonnegative` |
| `cubic`       | `x`, `y`, `z`      | $[-L / 2, L / 2]$                        | `ordinary`    |

A position outside the edges of a direction, or outside the Cartesian
edges of the maps, raises and aborts the run, so every frame's counts sum
to the species' particle number. The research code's `np.histogram`
dropped such a position silently, and its ranges stopped at $D / 2$ in a
cylinder and $L / 2$ in a cube; the maintainer extended them to the wall
surface and the corner on 2026-11-04, since a repulsive wall lets a small
particle's centre approach the surface closer than half a monomer and a
cube's corners lie beyond $L / 2$. The shells of a cube beyond $L / 2$ are
only partly inside the box, which the density normalization downstream has
to account for.

```python
from polyphys.probe import SpatialHistograms

parser = TransFociCyl(trajectory, 'whole', 'all', ispath=True)
analysis = SpatialHistograms.from_spec(universe, PROBE, parser).run()
analysis.save(save_to, f'{parser.name}-{parser.group}')
```

`save` writes, per species and direction, `<name>-<dir>Edge<Species>.npy`,
`<name>-<dir>Hist<Species>.npy` (summed over the frames, the research
artifact), and `<name>-<dir>HistT<Species>.npy` (per frame); the Cartesian
edges `<name>-xEdge.npy`, `yEdge`, `zEdge` without a species; per plane and
species `<name>-<plane>Hist<Species>.npy` (summed) and
`<name>-<plane>HistBlock<Species>.npy` (per block); and
`<name>-blockFrames.npy`.

In parallel, the blocks are numbered over the frames of the whole run, so
a run split across workers gives the same blocks as a serial one: each
worker adds its frames' maps into the block sums and counts of the whole
run, the aggregator sums them across workers (`ResultsGroup.ndarray_sum`)
and concatenates the per-frame rows (`ResultsGroup.ndarray_vstack`), and
`_conclude` divides each block's sum by its frame count. Per frame the cost
is $O(S N \log B)$ time for $S$ species of $N$ particles and $B$ bins, from
the sorted-edge search of `np.histogram`, and $O(N)$ extra space for the
recentred positions; over a run $8 F B$ bytes per one-dimensional
histogram and $8 K B_a B_b$ per map.

## The pair distances

{class}`~polyphys.probe.pairs.PairDistances` measures, for the species
`pair_species` of a `ProbeSpec`, the foci of the *TransFoci* studies, four
quantities per frame through {mod}`polyphys.analyze.clusters`:

| Result                       | Shape, dtype                    | Holds                                                                                                   |
| ---------------------------- | ------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `distMatT<Species>`          | `(n_frames, n, n)`, float64     | $d_{ij}$ for $i < j$ above the diagonal, $g_i$ on it, zeros below                                        |
| `directContactsMatT<Species>` | `(n_frames, n, n)`, bool       | $d_{ij} \le r_c$, symmetric, `True` on the diagonal                                                     |
| `bondsHistT<Species>`        | `(n_frames, n)`, int64          | `[t, d]`: the particles with $d$ contacts in frame $t$                                                  |
| `clustersHistT<Species>`     | `(n_frames, n + 1)`, int64      | `[t, s]`: the clusters of size $s$ in frame $t$; `[t, 0]` is 0                                           |

The cutoff is `contact_cutoff` resolved for the whole, `dmon_large +
dcrowd` in the *TransFoci* studies: two foci touch when a crowder no longer
fits between them. With `minimum_image`, the distances use the minimum
image along the periodic directions of the geometry with the frame's box
lengths; the research code used none. `save` writes the four arrays as
`<name>-<quantity><Species>.npy`, the research names; the contact matrix
is bool here and was int64 there.

**Indexing.** LAMMPS numbers atoms from 1 (`read_data`: "a unique value
from 1 to Natoms for each atom"), and MDAnalysis keeps those numbers as
`AtomGroup.ids` while numbering its own atoms from 0 as
`AtomGroup.indices`. The research code wrote `foci.atoms.ids` on the
diagonal of the distance matrix, so the pair tags of its analysis,
`loc<i>loc<j>`, carry 1-based IDs. Here the diagonal carries each focus's
0-based position among the chain's atoms, which `from_spec` computes and
{func}`~polyphys.analyze.clusters.genomic_distance`, behind the pipeline's
{class}`~polyphys.analyze.pipeline.FociPairDistances`, requires: a 1-based
diagonal would raise there whenever a focus is the chain's last monomer. A
pair tag of this package is therefore one less than the research tag of
the same pair; the genomic distance itself is the same under either
convention. The maintainer settled this on 2026-11-04.

Per frame the cost is $O(n^2)$ time and space; over a run, $8 F n^2$ bytes
for the distance matrices, $F n^2$ for the contact matrices, and $8 F (2 n
+ 1)$ for the two histograms.

## What changed from the research code

- The research code shifted every position by half the box in every frame
  before measuring. Every chain-size observable is invariant under a
  translation of all positions, because each one recentres the chain on
  its centre of geometry or mass first, so the shift is dropped for this
  family. The spatial histograms recentre by `ProbeSpec.recentre`, whose
  default `'half'` restores the LAMMPS frame and corrects the research
  shift's half-bead offset in a cylinder (see "Recentring a frame").
- The spatial histograms keep the one-dimensional counts per frame and the
  two-dimensional maps per block of frames, where the research code summed
  both over the frames; the summed arrays are still written. Their ranges
  reach the wall surface and the cube's corner, and a position outside
  them raises instead of being dropped.
- The pair distances write 0-based backbone indices on the diagonal of the
  distance matrix, use the minimum image along the periodic directions,
  and store the contact matrix as bool (see "The pair distances").
- The HnsCyl radial bin size includes the H-NS core diameter in its
  minimum, as HnsCub's did; the maintainer confirmed the research code's
  omission as a mistake on 2026-11-04.
- One class, driven by a `ProbeSpec`, replaces the per-study loops; the
  collectors are preallocated per run or per worker rather than per
  function.
- `asphericity(wrap=False, unwrap=...)` takes `unwrap` from the spec.
- The stamps carry a `group` column, which
  {func}`~polyphys.manage.organizer.children_stamps` reads when it saves;
  the research stamps did not.
- Frames are not converted to simulation time (see "Reading a
  trajectory").

## References

- R. I. Dima and D. Thirumalai, "Asymmetry in the Shapes of Folded and
  Denatured States of Proteins", J. Phys. Chem. B 108, 6564-6570 (2004).
  <https://doi.org/10.1021/jp037128y>
- Y. Wang, I. Teraoka, F. Y. Hansen, G. H. Peters, and O. Hassager, "A
  Theoretical Study of the Separation Principle in Size Exclusion
  Chromatography", Macromolecules 43, 1651-1659 (2010).
  <https://doi.org/10.1021/ma902377g>
- MDAnalysis, "Analysis building blocks", version 2.11.0.
  <https://docs.mdanalysis.org/2.11.0/documentation_pages/analysis/base.html>
- MDAnalysis, "LAMMPS DCD trajectory and DATA I/O", `DumpReader`, version
  2.11.0: the reader translates the box to the origin since 2.4.0.
  <https://docs.mdanalysis.org/2.11.0/documentation_pages/coordinates/LAMMPS.html>
- LAMMPS, "read_data command": atom-IDs run from 1 to the number of atoms.
  <https://docs.lammps.org/read_data.html>
