# Clusters: contacts, cluster sizes, and binding

This page explains what {mod}`polyphys.analyze.clusters` computes from one
frame of positions or contacts, how the cluster sizes are found, what the
pair-separation density is, and how the binding statistics of two-patch
binders are defined. Every routine works on one frame; reading frames and
accumulating over them is the caller's job.

## Notation

$n$
: The number of particles in the group, `n_particles`.

$r_c$
: The contact cutoff, `cutoff`, in the caller's length unit. Two particles
  are in *direct contact* when their center-to-center distance is at most
  $r_c$ (or below it, with `inclusive=False`).

$A$
: The *direct-contact matrix*, `direct_contacts`, boolean and
  $n \times n$: $A_{ij}$ is true when $i$ and $j$ are in direct contact. It
  is symmetric with a true diagonal. It is the adjacency matrix of the
  *contact graph*, whose vertices are the particles and whose edges are the
  direct contacts; $E$ is the number of edges, `n_contacts`.

*cluster*
: A connected component of the contact graph: a maximal set of particles
  joined by chains of direct contacts. Its *size* is its number of
  particles.

$C$
: The *contact matrix*, `contact_matrix`: $C_{ij}$ is true when $i$ and $j$
  are in the same cluster. It is the transitive closure of $A$.

*bonds*
: The number of direct contacts of one particle, its degree in the contact
  graph, from $0$ to $n - 1$.

*genomic distance*
: For two monomers of one chain, the number of backbone bonds between them:
  $|i - j|$ on a linear chain and $\min(|i - j|,\, N - |i - j|)$ on a ring
  of $N$ monomers.

$\vec R$, $R$
: The separation vector of one pair of monomers and its length,
  $R = |\vec R|$.

$P(\vec R)$, $\bar P(R)$
: The probability density of $\vec R$ over the chain's configurations,
  $\int P\, d^3R = 1$, and its average over the directions of $\vec R$.

*patch*, *binder*
: A binder is a protein with two binding sites, its patches; binder $k$ owns
  patches $2k$ and $2k + 1$. A binder is *free* when neither patch touches a
  monomer, *dangling* when one does, and *bridging* when both do. A bridge
  joins two monomers; it is *cis* when their genomic distance is at most
  `cis_threshold` bonds and *trans* when it is larger.

## From positions to clusters

A frame goes through four steps:

```python
import numpy as np
from polyphys.analyze import clusters

dist = clusters.pair_distances(positions, box=[0.0, 0.0, 2 * lz])
direct = clusters.find_direct_contacts(dist, cutoff=d_large + d_crowd)
bonds_hist = clusters.count_bonds(direct)        # [d]: particles with d bonds
cluster_hist = clusters.count_clusters(direct)   # [s]: clusters of size s
```

{func}`~polyphys.analyze.clusters.pair_distances` applies the minimum image
in every direction whose box length is positive and finite, so a pair that
straddles a periodic boundary of wrapped coordinates gets its true
separation. In the cylindrical projects only the axial direction is
periodic.

{func}`~polyphys.analyze.clusters.find_direct_contacts` is the contact
criterion of Sevick, Monson, and Ottino (1988). The cutoff is a project
constant. In the foci projects it is the sum of a large monomer's and a
crowder's diameter: the largest separation at which no crowder fits between
the two, so that depletion holds them together.

The clusters are the connected components of the contact graph, found by
{func}`scipy.sparse.csgraph.connected_components`
({func}`~polyphys.analyze.clusters.cluster_labels`). Finding the components
by depth-first search costs time and space proportional to $n + E$
(Hopcroft and Tarjan, 1973); converting the dense matrix to a sparse graph
first costs $O(n^2)$, so a frame costs $O(n^2)$ in all. The research code
merged pairs of matrix columns in Python, $O(n^3)$, and its single pass
could stop before the closure was reached: on the 7-node graph with edges
(0, 6), (1, 6), (2, 3), (3, 4), (4, 6) it reported clusters of sizes 1, 1, 5,
and 7, fourteen particles out of seven. That graph is a regression test.

{func}`~polyphys.analyze.clusters.count_clusters` histograms the component
sizes with two calls to {func}`numpy.bincount`. The research code read the
sizes off the eigenvalues of $C$, which is $O(n^3)$ and correct only when
$C$ is exactly block all-ones. Two conventions to note:

- The histogram is indexed by the size itself and has length $n + 1$;
  index 0 is always 0. The research code's index 0 held the number of zero
  eigenvalues, which the analysis stage dropped before averaging.
- Labels from `cluster_labels` are `0, ..., n_clusters - 1`, in an order
  SciPy does not specify. Compare `labels[i] == labels[j]`, never the label
  values.

## The pair-separation density

For a tagged pair of monomers, the probe saves the time series of $R$.
{func}`~polyphys.analyze.clusters.pair_distance_histograms` bins the series
of every pair on one grid, from 0 to the first multiple of `bin_size` at or
above `r_max`. Pass the same `r_max` for every whole of a project, chosen
from the system (the largest separation the box allows), not from the data:
a grid set by each whole's own maximum gives wholes different numbers of
bins, and an average over wholes is then taken over different sets of
wholes in different bins. A separation beyond the grid raises rather than
being dropped, since a dropped frame would break the normalization below.

{func}`~polyphys.analyze.clusters.pair_separation_density` turns the
counts $H_k$ of one pair, out of $N_f$ frames, into an estimate of
$\bar P(R)$. In spherical coordinates the scalar separation has density
$P(R) = 4\pi R^2 \bar P(R)$; the factor $4\pi R^2$ is the Jacobian and
nothing else. Integrating $P(\vec R)$ over the shell between the edges
$R_k$ and $R_{k+1}$, of volume $\Delta V_k = \tfrac{4\pi}{3}(R_{k+1}^3 -
R_k^3)$, gives the probability of the bin, so

$$
\bar P_k = \frac{H_k}{N_f\, \Delta V_k},
\qquad
\sum_k \bar P_k\, \Delta V_k = 1 ,
\qquad
\sum_k \bar P_k\, \Delta V_k\, R_k^2 \approx \langle R^2 \rangle .
$$

$\bar P$ has units of one over the length unit cubed, and
$-k_B T \ln \bar P(R)$ is the free energy of the pair at separation $R$ up
to a constant. For an ideal chain of $n$ Kuhn segments of length $b$,
$P(\vec R) = (3 / 2\pi n b^2)^{3/2} \exp(-3R^2 / 2 n b^2)$ (Rubinstein and
Colby, 2003, Chapter 2); because the estimate carries no free constant, it
can be laid over that reference directly. The estimate itself assumes no
chain model and no geometry: the accessible volume of a confinement shapes
$\bar P$, not its definition. The exact shell volume is used; the midpoint
shell $4\pi R_c^2 \Delta r$ is 25% too small in the first bin.

This is a probability density over the chain's configurations, not the
radial distribution function $g(r)$ of liquid-state theory, whose reference
state, two independent uniform points in the container, is not defined
here. For spatial densities of particles in the confinement, see
{mod}`polyphys.analyze.density` and the [Density](density.md) page.

## Binding statistics

The probe saves a monomer-patch contact matrix per frame, with the monomers
in backbone order along the rows and the two patches of binder $k$ in
columns $2k$ and $2k + 1$. {func}`~polyphys.analyze.clusters.binding_statistics`
returns a {class}`~polyphys.analyze.clusters.FrameBinding` with eight
counts and the bridged pairs; {class}`~polyphys.analyze.clusters.BindingStats`
accumulates frames and the loop-length histogram:

```python
stats = clusters.BindingStats(n_monomers=n_mon, topology='ring',
                              cis_threshold=4)
for frame in trajectory:
    stats.update(find_direct_contacts(dist_monomer_patch, cutoff))
series = stats.as_dict()                    # eight (n_frames,) arrays
loops = stats.loop_length_histogram         # [d]: bridged pairs at d bonds
```

The counting rules, kept from the research code so that its numbers are
reproduced:

- The patch and binder counts use the contacts as given, so a patch marked
  on two adjacent monomers (a cutoff artifact) is one engaged patch.
- The bridged pairs use the contacts after
  {func}`~polyphys.analyze.clusters.enforce_single_patch_contact`, which
  keeps the lowest-index monomer of each patch, so each bridging binder
  joins exactly two monomers.
- Cis and trans count *distinct monomer pairs*: two binders bridging the
  same pair count once, and a binder whose two patches hold one monomer is
  bridging but bridges no pair.
- In {func}`~polyphys.analyze.clusters.binder_contacts`, which links binders
  on the same or neighboring monomers for cluster counting, every binder is
  a cluster of at least size 1, bound or not.

Names changed from the research code's `hns_binding` dictionary:

| Research code     | {class}`~polyphys.analyze.clusters.FrameBinding` |
| ----------------- | ------------------------------------------------ |
| `n_m_hpatch_bound` | `n_patch_monomer_contacts`                      |
| `n_hpatch_engaged` | `n_patches_engaged`                             |
| `n_hpatch_free`    | `n_patches_free`                                |
| `n_hcore_free`     | `n_binders_free`                                |
| `n_hcore_dangle`   | `n_binders_dangling`                            |
| `n_hcore_bridge`   | `n_binders_bridging`                            |
| `n_hcore_cis`      | `n_pairs_cis`                                   |
| `n_hcore_trans`    | `n_pairs_trans`                                 |

Genomic distances are in bonds everywhere in the module. The research
code's `foci_info` returned the number of monomers between a pair, one less;
{func}`~polyphys.analyze.clusters.foci_pair_info` returns bonds and takes
the topology, so a linear chain no longer gets the ring's shorter way round.

### The threshold sweep

The binding criterion has two parameters, the cutoff and the cis threshold,
and the research code chose them with a sweep that PolyPhys reproduces:
{func}`~polyphys.analyze.clusters.threshold_sweep` takes the float64
monomer-to-patch distance matrices that
{class}`~polyphys.probe.binders.BinderStatistics` stores when the study's
`BinderSpec.store_distances` is on, shape `(n_frames, n_monomers,
n_patches)`, and recounts every frame at every `(cutoff, cis_threshold)`
pair of a grid:

```python
cutoffs = np.round(0.535 + 0.01 * np.arange(16), 3)   # 0.535, ..., 0.685
counts, loops = clusters.threshold_sweep(
    np.load('N200...ens1-nucleoid-distMatTMonHnsPatch.npy'),
    cutoffs, [2, 4], 'ring', inclusive=True)
counts            # one row per (cutoff, cis_threshold): eight frame means
loops             # (n_rows, n_lengths): the frame-averaged loop-length histogram
pd.concat([counts[['cutoff', 'cis_threshold']], pd.DataFrame(loops)], axis=1)
```

`counts` has the columns `cutoff`, `cis_threshold`, and the eight
{class}`~polyphys.analyze.clusters.FrameBinding` attribute names, each a
mean over the frames; its rows follow the cutoffs as given, with the cis
thresholds inside each cutoff, and no row is repeated. `loops` has one row
per row of `counts`, so the last line rebuilds the research
`loopLengthHistMon` table; the cis threshold does not enter the histogram,
so the rows of one cutoff are equal. Each frame is counted once per cutoff
and its bridged pairs are split once per cis threshold, since only the
cis-trans split depends on the threshold.

Three things to know before comparing a sweep with other numbers:

- `inclusive` is a required keyword. The research sweep compared strictly,
  `dist < cutoff`, while the H-NS study specs count with `inclusive=True`
  at 0.685 (the probe guide, "Choosing the cutoff"). A strict sweep's row
  at 0.685 omits the pairs at exactly the cutoff, so only an inclusive sweep
  reproduces the series that `BinderStatistics` saved, and only at the
  spec's cutoff and cis threshold.
- The research grid, `np.arange(0.535, 0.685 + 0.01, 0.01)`, ran one step
  past its end to 0.695 because of float round-off, so its tables hold
  seventeen cutoffs; the example above builds the grid from a count and
  stops at 0.685.
- The frame means treat the frames as exchangeable counts and are not
  uncertainties. An error bar on a swept count needs the per-frame series
  of {class}`~polyphys.analyze.clusters.BindingStats` and
  {func}`~polyphys.analyze.stats.block_average`.

The research script wrote two CSVs per whole,
`*-choice_of_threshold-bindingStat.csv` and
`*-choice_of_threshold-loopLengthHistMon.csv`, each rebuilt from lists that
were never reset between cutoffs, so every table repeated the rows of the
earlier cutoffs. The sweep returns each row once; the metadata columns the
script added (`dcrowd`, `phi_c_bulk`, `nhns`, `whole`, `ensemble_id`) are
the whole's parser attributes, which the caller attaches.

## References

- E. M. Sevick, P. A. Monson, and J. M. Ottino, "Monte Carlo calculations of
  cluster statistics in continuum models of composite morphology",
  J. Chem. Phys. 88, 1198-1206 (1988). <https://doi.org/10.1063/1.454720>
- J. Hopcroft and R. Tarjan, "Algorithm 447: efficient algorithms for graph
  manipulation", Commun. ACM 16, 372-378 (1973).
  <https://doi.org/10.1145/362248.362272>
- SciPy, {func}`scipy.sparse.csgraph.connected_components`.
- M. Rubinstein and R. H. Colby, Polymer Physics (Oxford University Press,
  2003), Chapter 2.
