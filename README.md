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
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git@v1.0.0"
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

A study declares a parser, a subclass of `ParserBase` that lists the
keywords of its file names. This short one reads a linear chain among
crowders in a cylinder; the
[Studies](https://polyphys.readthedocs.io/en/latest/studies.html) page shows
the full declaration:

```python
>>> from polyphys.manage.parser import KeywordSpec, ParserBase
>>> keys = ('N', 'D', 'ac', 'nc', 'bdump')
>>> class ChainCyl(ParserBase):
...     """A linear chain among crowders in a cylinder."""
...     _geometry = 'cylindrical'
...     _topology = 'linear'
...     _groups = ['bug']
...     _dump_keywords = {'bug': 'bdump'}
...     _keywords = {
...         'N': KeywordSpec('nmon', int),
...         'D': KeywordSpec('dcyl', float),
...         'ac': KeywordSpec('dcrowd', float),
...         'nc': KeywordSpec('ncrowd', int),
...         'bdump': KeywordSpec('bdump', int),
...         'ens': KeywordSpec('ensemble_id', int),
...         'j': KeywordSpec('segment_id', int),
...     }
...     _patterns = {
...         'segment': (*keys, 'ens', 'j'), 'whole': (*keys, 'ens'),
...         'ensemble_long': keys, 'ensemble': ('N', 'D', 'ac', 'nc'),
...         'space': ('N', 'D', 'ac'),
...     }
...     _project_attributes = dict.fromkeys(
...         ('segment', 'whole', 'ensemble_long', 'ensemble', 'space'), [])

```

A file name then decodes into typed physical attributes, and into the name of
the artifact at every level of the lineage above it, without touching the
file system:

```python
>>> artifact = ChainCyl('N200D20ac2nc1000bdump2000ens1.j03',
...                     lineage='segment', group='bug')
>>> artifact.nmon, artifact.ncrowd, float(artifact.dcrowd)
(200, 1000, 2.0)
>>> for level in artifact.lineage_genealogy:
...     print(f'{level:13} {getattr(artifact, level)}')
segment       N200D20ac2nc1000bdump2000ens1.j03
whole         N200D20ac2nc1000bdump2000ens1
ensemble_long N200D20ac2nc1000bdump2000
ensemble      N200D20ac2nc1000
space         N200D20ac2

```

The segments of one run share their `whole` name, the replicates of one
state point their `ensemble_long` name, and the state points of one sweep
their `space` name, so the organizer groups and averages them by name. The
[quickstart](https://polyphys.readthedocs.io/en/latest/quickstart.html) of
the documentation runs a small synthetic sweep from its segments to its
`space` table, and the
[Lineage](https://polyphys.readthedocs.io/en/latest/lineage.html) page shows
what else a parser derives, such as the cylinder size and the bulk volume
fractions.

## What's in the package

One short purpose per module; each guide page and the API reference give the
details. `probe` needs the `mda` extra and `viz.helpers` the `viz` extra.

| Module | Purpose | Guide |
| --- | --- | --- |
| `spec` | Study, probe, and analysis records shared across the workflow; no MDAnalysis. | [Studies](https://polyphys.readthedocs.io/en/latest/studies.html) |
| `manage.parser` | File name to typed physical attributes and lineage names. | [Lineage](https://polyphys.readthedocs.io/en/latest/lineage.html) |
| `manage.organizer` | Join segments into wholes and reduce replicates to ensemble averages. | [Workflow](https://polyphys.readthedocs.io/en/latest/workflow.html) |
| `manage.utils` | File sorting, file opening, and number-density and volume-fraction conversions. | [API](https://polyphys.readthedocs.io/en/latest/api.html) |
| `manage.types` | Type aliases of the vocabulary, such as `LineageT` and `PhaseT`. | [API](https://polyphys.readthedocs.io/en/latest/api.html) |
| `simulate.core` | The LAMMPS box and the readers' base error. | [Simulate](https://polyphys.readthedocs.io/en/latest/simulate.html) |
| `simulate.log` | Read LAMMPS log files into DataFrames. | [Simulate](https://polyphys.readthedocs.io/en/latest/simulate.html) |
| `simulate.dump` | Read LAMMPS text dumps, one frame at a time. | [Simulate](https://polyphys.readthedocs.io/en/latest/simulate.html) |
| `simulate.data` | Read and write LAMMPS data files. | [Simulate](https://polyphys.readthedocs.io/en/latest/simulate.html) |
| `simulate.build` | Build bead-spring start configurations for LAMMPS. | [Simulate](https://polyphys.readthedocs.io/en/latest/simulate.html) |
| `probe.core` | Open a whole's data file and dump as an MDAnalysis universe. | [Probe](https://polyphys.readthedocs.io/en/latest/probe.html) |
| `probe.chain_size` | Chain size and shape in every frame. | [Probe](https://polyphys.readthedocs.io/en/latest/probe.html) |
| `probe.binders` | Binding of two-patch binders to the chain in every frame. | [Probe](https://polyphys.readthedocs.io/en/latest/probe.html) |
| `probe.histograms` | Spatial histograms of every species. | [Probe](https://polyphys.readthedocs.io/en/latest/probe.html) |
| `probe.pairs` | Pair distances, contacts, and clusters of one species. | [Probe](https://polyphys.readthedocs.io/en/latest/probe.html) |
| `analyze.measurer` | Per-frame structural measures and histograms of one frame. | [Usage](https://polyphys.readthedocs.io/en/latest/usage.html) |
| `analyze.stats` | Errors of time averages over correlated frames, by blocking. | [Statistics](https://polyphys.readthedocs.io/en/latest/statistics.html) |
| `analyze.clusters` | Contacts, clusters, and binding statistics of one frame. | [Clusters](https://polyphys.readthedocs.io/en/latest/clusters.html) |
| `analyze.correlations` | Autocorrelations, decay fits, and bond correlations. | [Correlations](https://polyphys.readthedocs.io/en/latest/correlations.html) |
| `analyze.sumrule` | The density sum rule on a space's density profiles. | [Density](https://polyphys.readthedocs.io/en/latest/density.html) |
| `analyze.density` | Local number-density and volume-fraction profiles. | [Density](https://polyphys.readthedocs.io/en/latest/density.html) |
| `analyze.pipeline` | Reduce a space's probe artifacts down the lineage. | [Pipeline](https://polyphys.readthedocs.io/en/latest/pipeline.html) |
| `viz.colors` | The color lists of the studies' figures. | [Viz](https://polyphys.readthedocs.io/en/latest/viz.html) |
| `viz.helpers` | Matplotlib helpers for ticks, legends, and bands. | [Viz](https://polyphys.readthedocs.io/en/latest/viz.html) |

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

From 1.0.0, PolyPhys follows semantic versioning: a backward-incompatible
change to a stable module comes only in a major release, while `polyphys.probe`
and `polyphys.viz` stay provisional, as the
[API stability](https://polyphys.readthedocs.io/en/latest/stability.html) page
sets out. It is used to produce published research results, and the parsing
and measurement layers are covered by tests that run against Python 3.11,
3.12, 3.13, and 3.14 on every pull request and every push to `main`; the probe
layer's tests run on Python 3.14 with MDAnalysis.

The Python examples above run as doctests in the test suite.

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

| Part                                                      | Tier         | Checks or human review                                                                                                                 |
| --------------------------------------------------------- | ------------ | -------------------------------------------------------------------------------------------------------------------------------------- |
| `polyphys/analyze/`                                       | Supervised   | CI: flake8, mypy, pytest with doctests on Python 3.11–3.14; maintainer review of every line                                            |
| `polyphys/manage/`                                        | Supervised   | CI: flake8, mypy, pytest with doctests on Python 3.11–3.14; parser tests for filename-parsing changes; maintainer review of every line |
| `polyphys/simulate/`                                      | Supervised   | CI: flake8, mypy, pytest with doctests on Python 3.11–3.14; tests against LAMMPS-written files; maintainer review of every line        |
| `polyphys/probe/`                                         | Supervised   | CI: flake8, mypy; pytest with doctests on Python 3.14 with MDAnalysis installed; maintainer review of every line                       |
| `polyphys/viz/`                                           | Instrumented | CI: flake8, pytest with doctests on Python 3.11–3.14 with matplotlib installed; maintainer review                                      |
| `polyphys/tests/`, `.github/workflows/`, `pyproject.toml` | Instrumented | CI on every pull request: lint, type-check, tests, docs build, package build; maintainer review                                        |
| `projects/`                                               | Supervised   | CI: flake8, mypy, pytest of each spec against its parser and the research code's values; maintainer review of every line               |
| `docs/`                                                   | Instrumented | CI: Sphinx build with warnings as errors for `docs/source/`; maintainer review                                                         |

Nothing that handles security, credentials, private data, or published results,
or that can block a merge, is Delegated. Where a tier is unclear, I treat the
part as Supervised. Tiers last reviewed: 2026-11-08.

Instructions for AI agents: [`AGENTS.md`](https://github.com/amirhs1/PolyPhys/blob/main/AGENTS.md).

## Contributing

Bug reports and questions are welcome via
[issues](https://github.com/amirhs1/PolyPhys/issues). Please see
[`SECURITY.md`](https://github.com/amirhs1/PolyPhys/blob/main/SECURITY.md) for reporting suspected vulnerabilities privately.
See [`CONTRIBUTING.md`](https://github.com/amirhs1/PolyPhys/blob/main/CONTRIBUTING.md) for how to contribute.

## License

MIT — see [`LICENSE`](https://github.com/amirhs1/PolyPhys/blob/main/LICENSE).
