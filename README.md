# PolyPhys

[![CI](https://github.com/amirhs1/PolyPhys/actions/workflows/ci.yaml/badge.svg)](https://github.com/amirhs1/PolyPhys/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/amirhs1/PolyPhys/branch/main/graph/badge.svg)](https://app.codecov.io/gh/amirhs1/PolyPhys)
[![Docs](https://img.shields.io/readthedocs/polyphys/latest)](https://polyphys.readthedocs.io/en/latest/)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](https://github.com/amirhs1/PolyPhys)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15858407.svg)](https://doi.org/10.5281/zenodo.15858407)

**Research data management and analysis for large-scale polymer
molecular-dynamics simulations**, built around the study of bacterial
chromosome organization under macromolecular crowding.

## The problem

A single parameter sweep in a coarse-grained LAMMPS study produces tens of
thousands of files. The physics of each run — chain length, confinement
geometry, crowder size and count, timestep, dump frequency — is encoded in its
_filename_, and each physical state point is repeated across independent
replicates that must eventually collapse into one number with an honest error
bar.

That leaves two recurring problems, and PolyPhys is built around them:

1. **The filesystem is the database, but nothing can query it.** Filenames
   carry the physics, yet they are opaque strings.
2. **Replicates must be reduced correctly.** Averaging correlated
   molecular-dynamics frames with a naive standard error understates the
   uncertainty.

## The core idea: artifact lineage

Every simulation artifact sits at one of five levels of a hierarchy, and each
level is a well-defined aggregation of the one before it:

```text
segment  →  whole  →  ensemble_long  →  ensemble  →  space
 a chunk    a full     one state         replicates  the whole
 of a run   trajectory point, verbose    collapsed   parameter sweep
```

`polyphys.manage.parser` decodes a filename into typed physical attributes
_and_ derives that artifact's ancestry, so a file knows which state point and
which sweep it belongs to. Once artifacts know their lineage,
`polyphys.manage.organizer` can walk the hierarchy and reduce it level by
level.

## Installation

PolyPhys is distributed from source and archived on Zenodo; it is not on PyPI.

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git"
```

For a development checkout:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[dev]"
```

Requires Python 3.11+. The runtime dependencies are NumPy, pandas, and SciPy.
The probe layer, which reads LAMMPS trajectories, needs MDAnalysis as well;
the `mda` extra installs it:

```bash
python -m pip install "polyphys[mda] @ git+https://github.com/amirhs1/PolyPhys.git"
```

## Quickstart

Decode a simulation filename into physics, and read off where it sits in the
hierarchy:

```python
>>> from polyphys.manage.parser import SumRuleCyl
>>> artifact = SumRuleCyl(
...     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1.j03',
...     lineage='segment',
...     group='all',
... )

>>> artifact.nmon                     # monomers in the chain
200
>>> artifact.ncrowd                   # crowder particles
1000
>>> float(artifact.dcrowd)            # crowder diameter
2.0

```

Geometry is derived, not just read. The cylindrical wall is a layer of beads
of size `a_w = 1` whose centres lie on the cylinder of radius `r` set in
LAMMPS, so the cylinder size is `D = 2r - a_w`, stored as `dcyl`; `lz` is the
half-length, so the cylinder length is `Z = 2 lz`, stored as `lcyl`. For a
species `α` of size `a_α`, particle centres are confined to
`D_eff,α = D - a_α`, and its bulk number density and volume fraction both use
the effective volume `V_eff,α = π Z D_eff,α² / 4`:

```python
>>> float(artifact.dcyl)              # from 'r10.5': 2 * 10.5 - 1.0
20.0
>>> round(float(artifact.phi_bulk_c), 4)
0.3292

```

Each artifact also knows its ancestors, so a segment can be grouped with its
replicates and its sweep without touching the filesystem:

```python
>>> artifact.lineage_genealogy
['segment', 'whole', 'ensemble_long', 'ensemble', 'space']
>>> artifact.whole
'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1'
>>> artifact.ensemble
'N200D20ac2nc1000'
>>> artifact.space
'N200D20ac2'

```

`polyphys.analyze.measurer` provides the structural observables computed per
frame:

```python
>>> import numpy as np
>>> from polyphys.analyze import measurer
>>> positions = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [3.0, 4.0, 12.0]])
>>> float(measurer.end_to_end(positions))
13.0
>>> float(measurer.fsd(positions, 2))   # farthest-site distance along z
12.0

```

## What's in the package

| Module             | Purpose                                                                                                                                                                                        |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `manage.parser`    | Filename → typed physical attributes and lineage. One `ParserBase` subclass per project, declaring its own geometry, topology, and attribute schema.                                           |
| `manage.organizer` | Walks the lineage, combining segments into wholes and reducing replicate ensembles to averaged measurements.                                                                                   |
| `manage.projects`  | `ProjectSpec`, what a study's analysis needs beyond its parser: the polymer group, the rounding step of the crowder volume fraction, the attributes that become columns, the species-to-size map, and the time series; and `ProbeSpec`, what its probe needs: the selection of each species, the chain-size observables, the unwrap convention, the histogram species, bin sizes, and block length, how a frame is recentred, the pair species and its contact cutoff, and, for a study with binders, a `BinderSpec` of the binding cutoff and cis threshold. Each study's own specs live in `projects/`, outside the package. |
| `manage.utils`     | Filename sorting, safe IO, and number-density / volume-fraction conversions for cubic and cylindrical geometries.                                                                              |
| `manage.types`     | Domain type aliases (`LineageT`, `GeometryT`, `PhaseT`, …) that keep the vocabulary explicit and checkable.                                                                                    |
| `simulate.core`    | `Box`, the LAMMPS simulation box, and `LammpsFormatError`, the base class of the errors of the log, dump, and data-file readers.                                                               |
| `simulate.log`     | Reads a LAMMPS log file into pandas DataFrames: the thermodynamic output of every run, each run's summary, and the total wall time.                                                            |
| `simulate.dump`    | Reads the text dumps of the LAMMPS `dump atom` and `dump custom` commands into numpy arrays, one frame at a time, without MDAnalysis.                                                           |
| `simulate.data`    | Reads and writes LAMMPS data files of atoms, bonds, and angles, both as start configurations and as topologies for analysis tools such as MDAnalysis.                                          |
| `simulate.build`   | Builds LAMMPS start configurations of a bead-spring polymer, a helix or a closed plectoneme, and checks that LAMMPS can start from them.                                                        |
| `probe.core`       | Opens a study's LAMMPS data file and dump as one MDAnalysis universe, and writes the stamps of a probed artifact. Needs the `mda` extra.                                                      |
| `probe.chain_size` | `ChainSize`, an MDAnalysis analysis that measures a chain's end-to-end distance, span, transverse size, radius of gyration, asphericity, shape parameter, and principal axes in every frame, driven by a study's `ProbeSpec`. |
| `probe.binders`    | `BinderStatistics`, an MDAnalysis analysis that counts, in every frame, the patches of two-patch binders in contact with the chain, the free, dangling, and bridging binders, the cis and trans bridged pairs, and the loop-length histogram, with the minimum image along the geometry's periodic directions, driven by a study's `BinderSpec`. |
| `probe.histograms` | `SpatialHistograms`, an MDAnalysis analysis that bins every species of a whole along the directions of its geometry per frame and in the three Cartesian planes per block of frames, after recentring each frame, and `histogram_edges`, the whole's bin edges from its `ProbeSpec` and parser. |
| `probe.pairs`      | `PairDistances`, an MDAnalysis analysis that measures the pair distances, direct contacts, bonds, and clusters of one species, such as the foci, in every frame, with the minimum image along the periodic directions. |
| `analyze.measurer` | Per-frame structural observables — end-to-end distance, farthest-site distance, transverse size, maximum extent — plus fixed-size binning and radial, axial, azimuthal, and planar histograms. |
| `analyze.stats`    | The uncertainty of a time average over correlated frames by the blocking method of Flyvbjerg and Petersen: the blocked standard error, its plateau, and the statistical inefficiency.          |
| `analyze.clusters` | Contact and cluster statistics of a particle group (direct contacts under a cutoff, connected-component clusters, size and bond histograms), genomic distances, the probability density of a pair's separation, the free, dangling, bridging, cis, and trans statistics of two-patch binders, and the binding-threshold sweep on saved distance matrices. |
| `analyze.correlations` | The autocorrelation function of a time series by the FFT with Bartlett's confidence band, exponential and stretched-exponential decay fits with a recorded starting point, and the bond cosine correlation of a linear or ring chain, from which a persistence length follows. |
| `analyze.sumrule`  | The density sum rule of Jeon et al. (2016) on a space's ensemble-averaged density profiles: each species' profile scaled by its size, normalized, folded, and summed per state point.          |
| `analyze.density`  | Local number density and volume fraction profiles from a histogram of bead centres along one coordinate, with the share of each finite-size bead in the bins it reaches.                       |
| `analyze.pipeline` | Reduces one space's probe artifacts down the lineage for every property of a study's `AnalysisTable`, with post-hooks that derive time averages, autocorrelations, density profiles through an injected kernel, and the pair distances of the foci; and reduces each whole's time series to its equilibrium value with a blocked error. |

Nine parser subclasses ship with the package, covering cylindrical and cubic
confinement: linear and ring chains, some with two monomer sizes, among
crowders of one or two sizes, and a pair of monomers held at a fixed distance
(`TwoMonDepCub`).

## Documentation

The full documentation, including the API reference, is at
[polyphys.readthedocs.io](https://polyphys.readthedocs.io/en/latest/). Its
quickstart runs a synthetic sweep from segments to a `space` table, and its
lineage, workflow, and studies pages explain the phases, the steps from
trajectories to a space, and how to add a study. To build it locally:

```bash
python -m pip install -e ".[docs]"
python -m sphinx -W -b html docs/source docs/_build
```

## Project status

PolyPhys is under active development toward a stable 1.0 API. It is used to
produce published research results, and the parsing and measurement layers are
covered by tests that run against Python 3.11, 3.12, 3.13, and 3.14 on every
pull request and every push to `main`; the probe layer's tests run on Python
3.14 with MDAnalysis.

The examples above are executed as part of the test suite, so anything shown in
this README is guaranteed to run against the current code.

Current focus, in order:

1. Expand `analyze` with the polymer- and crowding-specific observables that
   general MD toolkits do not provide, and finish the trajectory probes of
   `probe`.
2. Check `analyze.density` against the research code's density profiles.
3. Move the study parsers out of `manage.parser` into their study folders,
   so a new study _uses_ PolyPhys from outside the package rather than
   extending it; the study specs have already moved to `projects/`.

## Citation

If you use PolyPhys in published work, please cite the archived release. The
DOI above resolves to the latest version; see [`CITATION.cff`](CITATION.cff)
for machine-readable metadata, or use GitHub's _Cite this repository_ button.

## AI assistance

PolyPhys is developed with the help of AI coding assistants. Every change is
reviewed and merged by the maintainer, who is accountable for the result; no AI
system is an author of this software or of any release. Numerical fixtures,
formulas, and citations are derived or verified rather than accepted from a
model — see [`AI-POLICY.md`](AI-POLICY.md) for the full policy, including what
is expected of AI-assisted contributions.

Each part has a tier, set by whether I can evaluate AI output there. The table
records the checks or human review applied to each part:

- **Instrumented** — I could write it myself. AI is used for review,
  refactoring, and alternative implementations, not first drafts of core logic.
- **Supervised** — AI drafts; I read every line and set the acceptance criteria
  and test values.
- **Delegated** — AI generates; I can't fully evaluate it. It is covered by
  tests, kept isolated and low-risk, and not presented as my work.

| Part                                                      | Tier         | Checks or human review                                                                                                   |
| --------------------------------------------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------ |
| `polyphys/analyze/`                                       | Instrumented | CI: flake8, mypy, pytest with doctests on Python 3.11–3.14; maintainer review                                            |
| `polyphys/manage/`                                        | Instrumented | CI: flake8, mypy, pytest with doctests on Python 3.11–3.14; parser tests for filename-parsing changes; maintainer review |
| `polyphys/simulate/`                                      | Instrumented | CI: flake8, mypy, pytest with doctests on Python 3.11–3.14; tests against LAMMPS-written files; maintainer review        |
| `polyphys/probe/`                                         | Instrumented | CI: flake8, mypy; pytest with doctests on Python 3.14 with MDAnalysis installed; maintainer review                        |
| `polyphys/tests/`, `.github/workflows/`, `pyproject.toml` | Instrumented | CI on every pull request: lint, type-check, tests, docs build, package build; maintainer review                          |
| `projects/`                                               | Supervised   | CI: flake8, mypy, pytest of each spec against its parser and the research code's values; maintainer review of every line |
| `docs/`                                                   | Instrumented | CI: Sphinx build with warnings as errors for `docs/source/`; maintainer review                                           |

Nothing that handles security, credentials, private data, or published results,
or that can block a merge, is Delegated. Where a tier is unclear, I treat the
part as Supervised. Tiers last reviewed: 2026-11-06.

Instructions for AI agents: [`AGENTS.md`](AGENTS.md).

## Contributing

Bug reports and questions are welcome via
[issues](https://github.com/amirhs1/PolyPhys/issues). Please see
[`SECURITY.md`](SECURITY.md) for reporting suspected vulnerabilities privately.
See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to contribute.

## License

MIT — see [`LICENSE`](LICENSE).
