# PolyPhys

[![CI](https://github.com/amirhs1/PolyPhys/actions/workflows/ci.yaml/badge.svg)](https://github.com/amirhs1/PolyPhys/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/amirhs1/PolyPhys/branch/main/graph/badge.svg)](https://app.codecov.io/gh/amirhs1/PolyPhys)
[![Docs](https://img.shields.io/readthedocs/polyphys/latest)](https://polyphys.readthedocs.io/en/latest/)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](https://github.com/amirhs1/PolyPhys)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](https://github.com/amirhs1/PolyPhys/blob/main/LICENSE)
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

From v1.0.0, PolyPhys is on PyPI, and each release is archived on Zenodo:

```bash
python -m pip install polyphys
```

With [uv](https://docs.astral.sh/uv/), install it into the active environment,
or add it as a dependency of a uv project:

```bash
uv pip install polyphys
uv add polyphys
```

PolyPhys is not on conda-forge yet; a conda-forge package is planned after
v1.0.0 is on PyPI. Until it is listed, create a conda environment with the
dependencies from conda-forge, then install PolyPhys into it with pip, as
[conda's guide](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
advises: conda packages first, then pip. `mdanalysis` is needed only by the
probe layer and `matplotlib-base` only by the viz layer's figure helpers; leave
out the ones you do not need:

```bash
conda create -n polyphys -c conda-forge "python>=3.11" "numpy>=2.0" "pandas>=3.0" "scipy>=1.14.1" "mdanalysis>=2.8" "matplotlib-base>=3.8.4"
conda activate polyphys
python -m pip install polyphys
```

Releases before v1.0.0 are not on PyPI; install one from its tag:

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git@v0.9.0"
```

For a development checkout:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[dev,viz]"
```

The `viz` extra installs matplotlib, which the doctests of `polyphys/viz/`
import; without it the documented test command stops with an `ImportError`.

Requires Python 3.11+. The runtime dependencies are NumPy 2.0 or newer, pandas
3.0 or newer, and SciPy 1.14.1 or newer. The probe layer, which reads LAMMPS
trajectories, needs MDAnalysis as well; the `mda` extra installs it:

```bash
python -m pip install "polyphys[mda]"
```

The viz layer's figure helpers need matplotlib 3.8.4 or newer; the `viz` extra
installs it (the color lists in `viz.colors` need nothing extra):

```bash
python -m pip install "polyphys[viz]"
```

## Quickstart

Decode a simulation filename into physics, and read off where it sits in the
hierarchy. A study declares its parser, a subclass of `ParserBase` that lists
the keywords of its names; this one reads a linear chain among crowders in a
cylinder, and the [Studies](https://polyphys.readthedocs.io/en/latest/studies.html) page shows how to declare
one:

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
of size `a_w = 1` whose centers lie on the cylinder of radius `r` set in
LAMMPS, so the cylinder size is `D = 2r - a_w`, stored as `dcyl`; `lz` is the
half-length, so the cylinder length is `Z = 2 lz`, stored as `lcyl`. For a
species `α` of size `a_α`, particle centers are confined to
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
'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1'
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
| `spec`             | `ProjectSpec`, what a study's analysis needs beyond its parser: the polymer group, the rounding step of the crowder volume fraction, the attributes that become columns, the species-to-size map, and the time series; and `ProbeSpec`, what its probe needs: the selection of each species, the chain-size observables, the unwrap convention, the histogram species, bin sizes, and block length, how a frame is recentered, the pair species and its contact cutoff, and, for a study with binders, a `BinderSpec` of the binding cutoff and cis threshold; and `AnalysisTable`, the artifacts of one particle group to reduce. It imports no MDAnalysis. Each study's own specs live in `projects/`, outside the package. |
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
| `probe.histograms` | `SpatialHistograms`, an MDAnalysis analysis that bins every species of a whole along the directions of its geometry per frame and in the three Cartesian planes per block of frames, after recentering each frame, and `histogram_edges`, the whole's bin edges from its `ProbeSpec` and parser. |
| `probe.pairs`      | `PairDistances`, an MDAnalysis analysis that measures the pair distances, direct contacts, bonds, and clusters of one species, such as the foci, in every frame, with the minimum image along the periodic directions. |
| `analyze.measurer` | Per-frame structural observables — end-to-end distance, farthest-site distance, transverse size, maximum extent — plus fixed-size binning and radial, axial, azimuthal, and planar histograms. |
| `analyze.stats`    | The uncertainty of a time average over correlated frames by the blocking method of Flyvbjerg and Petersen: the blocked standard error, its plateau, and the statistical inefficiency.          |
| `analyze.clusters` | Contact and cluster statistics of a particle group (direct contacts under a cutoff, connected-component clusters, size and bond histograms), genomic distances, the probability density of a pair's separation, the free, dangling, bridging, cis, and trans statistics of two-patch binders, and the binding-threshold sweep on saved distance matrices. |
| `analyze.correlations` | The autocorrelation function of a time series by the FFT with Bartlett's confidence band, exponential and stretched-exponential decay fits with a recorded starting point, and the bond cosine correlation of a linear or ring chain, from which a persistence length follows. |
| `analyze.sumrule`  | The density sum rule of Jeon et al. (2016) on a space's ensemble-averaged density profiles: each species' profile scaled by its size, normalized, folded, and summed per state point.          |
| `analyze.density`  | Local number density and volume fraction profiles from a histogram of bead centers along one coordinate, with the share of each finite-size bead in the bins it reaches.                       |
| `analyze.pipeline` | Reduces one space's probe artifacts down the lineage for every property of a study's `AnalysisTable`, with post-hooks that derive time averages, autocorrelations, density profiles through an injected kernel, and the pair distances of the foci; and reduces each whole's time series to its equilibrium value with a blocked error. |
| `viz.colors`       | The qualitative color lists of the studies' figures as `#RRGGBB` tuples (Paul Tol's bright, vibrant, and muted schemes, the Dutch Field palette, and the maintainer's own list), which need no plotting library, and `truncate_colormap`, a matplotlib color map cut down to a sub-range of its domain. |
| `viz.helpers`      | `set_xticks` and `set_yticks`, which place an axis's major ticks every `step` from `lower` to `upper` inside a small frame; proxy legend handles keyed by color or by line style and marker; a legend moved from a figure into chosen panels; and a band shaded under each line in its own color. Needs the `viz` extra. |

The package ships no study's parser: each study declares its own
`ParserBase` subclass in its own folder. The parsers of the maintainer's
published studies are in `projects/`, outside the installed package.

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

A study _uses_ PolyPhys from outside the package rather than extending it:
its parser and its specs live in its own folder, as `projects/` shows.

## Citation

If you use PolyPhys in published work, please cite the archived release. The
DOI above resolves to the latest version; see [`CITATION.cff`](https://github.com/amirhs1/PolyPhys/blob/main/CITATION.cff)
for machine-readable metadata, or use GitHub's _Cite this repository_ button.

## AI assistance

PolyPhys is developed with the help of AI coding assistants. Every change is
reviewed and merged by the maintainer, who is accountable for the result; no AI
system is an author of this software or of any release. Numerical fixtures,
formulas, and citations are derived or verified rather than accepted from a
model — see [`AI-POLICY.md`](https://github.com/amirhs1/PolyPhys/blob/main/AI-POLICY.md) for the full policy, including what
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
| `polyphys/viz/`                                           | Instrumented | CI: flake8, pytest with doctests on Python 3.11–3.14 with matplotlib installed; maintainer review                         |
| `polyphys/tests/`, `.github/workflows/`, `pyproject.toml` | Instrumented | CI on every pull request: lint, type-check, tests, docs build, package build; maintainer review                          |
| `projects/`                                               | Supervised   | CI: flake8, mypy, pytest of each spec against its parser and the research code's values; maintainer review of every line |
| `docs/`                                                   | Instrumented | CI: Sphinx build with warnings as errors for `docs/source/`; maintainer review                                           |

Nothing that handles security, credentials, private data, or published results,
or that can block a merge, is Delegated. Where a tier is unclear, I treat the
part as Supervised. Tiers last reviewed: 2026-11-06.

Instructions for AI agents: [`AGENTS.md`](https://github.com/amirhs1/PolyPhys/blob/main/AGENTS.md).

## Contributing

Bug reports and questions are welcome via
[issues](https://github.com/amirhs1/PolyPhys/issues). Please see
[`SECURITY.md`](https://github.com/amirhs1/PolyPhys/blob/main/SECURITY.md) for reporting suspected vulnerabilities privately.
See [`CONTRIBUTING.md`](https://github.com/amirhs1/PolyPhys/blob/main/CONTRIBUTING.md) for how to contribute.

## License

MIT — see [`LICENSE`](https://github.com/amirhs1/PolyPhys/blob/main/LICENSE).
