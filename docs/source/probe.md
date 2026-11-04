# Probe: per-frame observables from LAMMPS trajectories

This page explains what {mod}`polyphys.probe` does with a study's LAMMPS
trajectories, how a study describes its probe as data, what the chain-size
analysis measures, how it runs in serial or in parallel, where its outputs
go, and what changed from the research code it replaces.

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
    selections={'Mon': 'resid 1', 'Foci': 'type 2'},
    chain='Mon',
    observables=('asphericityT', 'fsdT', 'gyrT', 'shapeT', 'principalT'),
    unwrap=False,
)
```

| Field         | Holds                                                                                                                  |
| ------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `selections`  | Each species label mapped to its selection, read-only, in order                                                        |
| `chain`       | The label in `selections` of the polymer chain; its series are named `<observable><chain>`                             |
| `observables` | The chain-size observables to measure, each a name in `CHAIN_SIZE_OBSERVABLES`                                        |
| `unwrap`      | Whether the chain is unwrapped across the periodic boundaries for its asphericity: `False` in a cylinder, `True` in a cube |
| `bin_sizes`   | The histogram bin size along each direction, for the histogram probes to come; empty until the sizes are confirmed     |

The spec is checked at construction: the chain label must be a selection,
every observable must be known and listed once, every scalar observable
with the chain label must be one of the project's `tseries_properties`
(`principalT`, a per-frame matrix, is exempt), and every bin size must be
finite and positive along a known direction.

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

## What changed from the research code

- The research code shifted every position by half the box in every frame
  before measuring. Every chain-size observable is invariant under a
  translation of all positions, because each one recentres the chain on
  its centre of geometry or mass first, so the shift is dropped for this
  family. The histogram probes, which will need it, follow in later pull
  requests.
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
