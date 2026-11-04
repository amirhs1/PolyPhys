# Probe: per-frame observables from LAMMPS trajectories

This page explains what {mod}`polyphys.probe` does with a study's LAMMPS
trajectories, how a study describes its probe as data, what the chain-size
analysis measures, what the binder analysis counts and how its cutoff was
chosen, how the analyses run in serial or in parallel, where their outputs
go, and what changed from the research code they replace.

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

$N_m$, $N_p$, $N_b = N_p / 2$
: The numbers of monomers, binder patches, and binders; the two patches of
  binder $k$ are the patches $2k$ and $2k + 1$ in atom-ID order.

$r_{ij}$, $r_\mathrm{thres}$
: The distance between monomer $i$ and patch $j$, with the minimum image
  along each periodic direction, and the binding threshold: the pair is
  in contact when $r_{ij} \le r_\mathrm{thres}$ (inclusive) or
  $r_{ij} < r_\mathrm{thres}$ (strict).

$d$
: The genomic distance of a bridged pair of monomers, in bonds along the
  backbone, the shorter way round a ring; a bridge is *cis* when
  $d \le$ `cis_threshold` and *trans* otherwise.

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
| `binders`     | A `BinderSpec` of how the study's two-patch binders are counted (see "The binder analysis"); `None` for a study without binders |

The spec is checked at construction: the chain label must be a selection,
every observable must be known and listed once, every scalar observable
with the chain label must be one of the project's `tseries_properties`
(`principalT`, a per-frame matrix, is exempt), every bin size must be
finite and positive along a known direction, and a `BinderSpec` must name
a selected label other than the chain as its patches and eight binding
series that the project's `tseries_properties` list.

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

## The binder analysis

{class}`~polyphys.probe.binders.BinderStatistics` counts, in every frame,
how the two-patch binders of a study bind its chain. The binder is the
bivalent patchy protein of Brackley et al. (2013): a repulsive core with a
sticky patch at each pole, which the H-NS studies follow (Jung, Sadeghi,
and Ha, 2024). The analysis takes the chain and the binders' patches, forms
the distance $r_{ij}$ between every monomer and every patch, marks the
pairs within the cutoff as contacts with
{func}`~polyphys.analyze.clusters.find_direct_contacts`, and counts them
with {func}`~polyphys.analyze.clusters.binding_statistics`, whose rules the
clusters guide states: a patch within the cutoff of two adjacent monomers
is one engaged patch, each bridging binder joins the lowest-index monomer
of each of its patches, and cis and trans count distinct monomer pairs. The
dangling, cis, and trans modes are those of Joyeux and Vreede (2013), whose
cis bridge joins sites "separated by at most two DNA beads", a genomic
distance of at most 3 bonds; `cis_threshold` of 4 reproduces that on the
studies' chains.

A study declares the facts the analysis needs in a
{class}`~polyphys.manage.projects.BinderSpec` inside its `ProbeSpec`:

```python
from polyphys.manage.projects import BinderSpec, ProbeSpec

PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'HnsPatch': 'type 2'},
    chain='Mon',
    observables=('asphericityT', 'gyrT', 'shapeT', 'principalT'),
    unwrap=True,
    binders=BinderSpec(patches='HnsPatch', binders='HnsCore', cutoff=0.685,
                       inclusive=True, cis_threshold=4,
                       store_distances=True),
)
```

| Field                | Holds                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `patches`            | The label in `selections` of the patches, in atom-ID order with the two patches of binder $k$ at $2k$ and $2k + 1$; names the three patch series |
| `binders`            | The label of the binders, a name only, since the cores are not selected; names the five binder series                              |
| `cutoff`             | The binding threshold $r_\mathrm{thres}$, in monomer diameters; a study constant, never computed from the particle sizes          |
| `inclusive`          | Whether a pair at exactly the cutoff is a contact: `True` by default, the thesis criterion; the research probe compared strictly   |
| `cis_threshold`      | The largest genomic distance, in bonds, of a cis bridge; 4 in the H-NS studies                                                     |
| `patches_per_binder` | 2, the only valence `binding_statistics` counts; kept as a field so that binders of another valence can be declared later          |
| `store_distances`    | Whether the monomer-to-patch distance matrices are kept and saved; `True` in the H-NS studies, which post-process them             |

```python
from polyphys.manage.parser import HnsCub
from polyphys.probe import BinderStatistics, write_stamps
from hns_cub import PROBE

parser = HnsCub(trajectory, 'whole', 'nucleoid', ispath=True)
analysis = BinderStatistics.from_spec(universe, PROBE).run()
analysis.save(save_to, f'{parser.name}-{parser.group}')
```

`from_spec` selects the chain and the patches, takes the topology from the
project, and takes the periodic directions from the geometry: every
direction in a periodic cube, the axis only in a cylinder, whose LAMMPS
boundary is `ff ff pp`. Along a periodic direction the distance takes the
minimum image with the frame's box length. The research probe called
`distance_array` without a box, so it took no image, and the `box`
argument of `distance_array` wraps every direction, which is wrong for a
cylinder, so the module forms the distances itself.

After `run()`, `results` holds one array per key, each with a leading frame
axis, which `save` writes as `<name>-<series>.npy` under the research
code's names:

| Key                                                                                                                                                                     | Shape, dtype                      | Saved as, with the H-NS labels                                                                                                                 |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `n_patch_monomer_contacts`, `n_patches_free`, `n_patches_engaged`                                                                                                     | `(n_frames,)`, int64              | `nBoundTHnsPatch`, `nFreeTHnsPatch`, `nEngagedTHnsPatch`                                                                                       |
| `n_binders_free`, `n_binders_bridging`, `n_binders_dangling`, `n_pairs_cis`, `n_pairs_trans`                                                                          | `(n_frames,)`, int64              | `nFreeTHnsCore`, `nBridgeTHnsCore`, `nDangleTHnsCore`, `nCisTHnsCore`, `nTransTHnsCore`                                                        |
| `loop_length_histogram`: row $t$, column $d$ is the number of distinct bridged pairs at $d$ bonds in frame $t$; $N_m // 2 + 1$ columns on a ring, $N_m$ on a linear chain | `(n_frames, n_lengths)`, int64    | `loopLengthHistTMon`                                                                                                                           |
| `distances`, with `store_distances`                                                                                                                                     | `(n_frames, N_m, N_p)`, float64   | `distMatTMonHnsPatch`                                                                                                                          |

The eight series are keyed by the attributes of
{class}`~polyphys.analyze.clusters.FrameBinding` and named by
{data}`~polyphys.manage.projects.BINDER_SERIES`, which maps each one to the
stem of its time series and to the label, patches or binders, that
completes the name. The loop-length rows are the maintainer's 2026-11-03
convention for a per-frame histogram, as the 1-D histogram probes keep
their int64 counts; the research probe saved their sum over the frames as
`loopLengthHistMon`, which is `rows.sum(axis=0)`. The H-NS analysis tables
in `projects/` still read that research artifact through
`tseries('loopLengthHist', 'Mon')`; reading the per-frame rows is an
analysis-side change for a later pull request.

The distance matrices are kept because two research post-steps read them
rather than the trajectory: the threshold sweep described next, and the
H-NS core clustering script, which recomputes the contacts and writes the
`bondsHistT` and `clustersHistT` families of `HnsCore` that the `HnsCub`
analysis table lists. They are float64, as `distance_array` returned them,
so a later sweep makes the same contact decisions as the research sweep
did; at $8 N_m N_p$ bytes per frame they are held in memory for the run,
as the research list was.

### Choosing the cutoff

The cutoff is a scientific choice that the spec states; the module never
computes it from the particle sizes, because three formulas give three
defensible numbers and the research probe picked the one that matches no
published table. The three sources of the H-NS model agree on the
patch-monomer geometry. Each starts from a patch sphere whose repulsive
range is $r_{\mathrm{cut},pp} = 0.2$, buried in a core of diameter $1$,
touching a monomer of diameter $1$, all in monomer diameters:

| Step                                                                     | Brackley et al. (2013), Table S2 and Fig. S4 | Jung, Sadeghi, and Ha (2024), Table 1 | Thesis, Table 5.1 and Fig. 5.2 |
| ------------------------------------------------------------------------ | -------------------------------------------- | ------------------------------------- | ------------------------------ |
| Patch-patch cutoff $r_{\mathrm{cut},pp}$                                 | 0.2                                          | 0.2                                   | 0.2                            |
| Patch size $\sigma_p = r_{\mathrm{cut},pp} / 2^{1/6}$                    | 0.178                                        | 0.178                                 | 0.178                          |
| Patch-to-core distance $(\sigma_h - r_{\mathrm{cut},pp}) / 2$            | 0.4 (the patch touches the pole)             | 0.4                                   | 0.4                            |
| Patch-monomer contact distance $d_{pm} = (\sigma_m + r_{\mathrm{cut},pp}) / 2$ | implied by $\sigma_{mp}$                | 0.6                                   | 0.6                            |
| Patch-monomer $\sigma_{mp} = d_{pm} / 2^{1/6}$                           | 0.535                                        | 0.535                                 | 0.535                          |
| Patch-monomer cutoff $r_{\mathrm{cut},mp} = d_{pm} + 0.085$              | 0.685                                        | 0.685                                 | 0.685                          |
| Well depth $\varepsilon_{mp}$                                            | 50                                           | 29                                    | 29                             |

Two checks tie the three together: the PNAS supporting information gives
the interaction range as 1.71 nm for 2.5 nm beads, and $0.685 \times 2.5$
is 1.7125; and the thesis bounds the attractive range $\delta_\mathrm{attr}
= 0.085$ by $\sigma_h (2^{1/6} - 1) \approx 0.122$, the reach of the core's
own repulsion past its contact distance, so that a monomer feels the patch
only where it would also feel the core. Outside $r_{\mathrm{cut},mp}$ the
force on the monomer is identically zero, which makes $0.685$ a natural
threshold for "bound": thesis Section 5.3 defines a bound H-NS as one
whose patch is within $r_\mathrm{thres} = r_{\mathrm{cut},mp}$ of a
monomer, with $r \le r_\mathrm{thres}$.

The research probe used a different number: $2^{1/6} (d_\mathrm{mon} +
d_\mathrm{hns\_patch}) / 2 = 2^{1/6} (1 + 0.178) / 2 = 0.66113$, rounded to
three decimals as $0.661$, with a strict comparison. That is the
Lennard-Jones minimum of an arithmetic-mean $\sigma$ of $0.589$, which no
source uses for this pair, since each sets $\sigma_{mp}$ from the contact
geometry instead. At $0.661$ the shifted potential is about $-3.1\,k_B T$
against $-8.7\,k_B T$ at the minimum, so a strict cut there drops the pairs
in the outer third of the well that the published criterion counts.

The research code also measured how much the counts depend on the
threshold. The sweep in SciOps `bash_pipeline/choice_of_threshold` and the
HnsCyl project notebook reloaded the saved distance matrices of the HnsCyl
wholes and recounted for every cutoff from $\sigma_{mp} = 0.535$ to
$r_{\mathrm{cut},mp} = 0.685$ in steps of $0.01$ and for a cis threshold of
2 and of 4. Its conclusions, recorded in the notebook: below a cutoff of
$0.575$ no loops of size 2 or more form; the loop-size distribution peaks
at a genomic distance of 4 for every cutoff; the normalized loop-size
distributions collapse for cutoffs of $0.645$ and above; and the cis
threshold of 2 against 4 changes nothing visible. Both $0.661$ and $0.685$
lie in that insensitive region, so the switch moves the ensemble averages
little while the per-frame counts differ by the pairs in the outer band of
the well.

The H-NS study specs therefore state the published criterion, `cutoff=0.685`
with `inclusive=True`, as the maintainer decided on 2026-11-04. A study
that must reproduce the research probe's artifacts exactly states
`cutoff=0.661` with `inclusive=False` instead; the change is one spec edit.

### Parallel runs and cost

The class is marked parallelizable with the same three backends as
`ChainSize`; each worker fills the arrays of its own frames, and the
aggregator concatenates them along the frame axis
(`ResultsGroup.ndarray_hstack` for the eight series,
`ResultsGroup.ndarray_vstack` for the histogram rows and the distances).
Per frame, the distances and the contacts cost $O(N_m N_p)$ time and
space, and making the bridged pairs distinct $O(N_b \log N_b)$. Over a run
of $F$ frames the series take $O(F)$ space each, the histogram
$O(F N_m / 2)$, and the distances $8 F N_m N_p$ bytes.

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
- The binder distances take the minimum image along the geometry's
  periodic directions; the research probe took none.
- The binding cutoff and its comparison are the spec's. The H-NS studies
  state the published criterion, $0.685$ inclusive, where the research
  probe used $0.661$ strict (see "Choosing the cutoff").
- The loop-length histogram is saved per frame as `loopLengthHistTMon`;
  the research probe summed it over the frames as `loopLengthHistMon`,
  which the H-NS analysis tables still read.
- The distance matrices are an option of the spec, on in the H-NS
  studies, rather than always written.

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
- C. A. Brackley, S. Taylor, A. Papantonis, P. R. Cook, and D. Marenduzzo,
  "Nonspecific bridging-induced attraction drives clustering of DNA-binding
  proteins and genome organization", Proc. Natl. Acad. Sci. U.S.A. 110,
  E3605-E3611 (2013). <https://doi.org/10.1073/pnas.1302950110>
- C. A. Brackley, "Polymer compaction and bridging-induced clustering of
  protein-inspired patchy particles", J. Phys.: Condens. Matter 32, 314002
  (2020). <https://doi.org/10.1088/1361-648X/ab7f6c>
- M. Joyeux and J. Vreede, "A Model of H-NS Mediated Compaction of
  Bacterial DNA", Biophys. J. 104, 1615-1622 (2013).
  <https://doi.org/10.1016/j.bpj.2013.02.043>
- Y. Jung, A. Sadeghi, and B.-Y. Ha, "Modeling the compaction of bacterial
  chromosomes by biomolecular crowding and the cross-linking protein H-NS",
  Sci. Rep. 14, 139 (2024). <https://doi.org/10.1038/s41598-023-50355-2>
- The maintainer's PhD thesis, "Computational modeling of bacterial
  chromosome organization: macromolecular crowding, chain heterogeneity,
  and chain cross-linking", University of Waterloo (2024), chapter 5.
