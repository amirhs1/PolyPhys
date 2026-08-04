# PolyPhys

[![CI](https://github.com/amirhs1/poly-phys/actions/workflows/ci.yaml/badge.svg)](https://github.com/amirhs1/poly-phys/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/amirhs1/poly-phys/branch/main/graph/badge.svg)](https://app.codecov.io/gh/amirhs1/poly-phys)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](https://github.com/amirhs1/poly-phys)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15858407.svg)](https://doi.org/10.5281/zenodo.15858407)

**Research data management and analysis for large-scale polymer molecular-dynamics
simulations**, built around the study of bacterial chromosome organization under
macromolecular crowding.

## The problem

A single parameter sweep in a coarse-grained LAMMPS study produces tens of thousands
of files. The physics of each run — chain length, confinement geometry, crowder size
and count, timestep, dump frequency — is encoded in its *filename*, and each physical
state point is repeated across independent replicates that must eventually collapse
into one number with an honest error bar.

That leaves two recurring problems, and PolyPhys is built around them:

1. **The filesystem is the database, but nothing can query it.** Filenames carry the
   physics, yet they are opaque strings.
2. **Replicates must be reduced correctly.** Averaging correlated molecular-dynamics
   frames with a naive standard error understates the uncertainty.

## The core idea: artifact lineage

Every simulation artifact sits at one of five levels of a hierarchy, and each level
is a well-defined aggregation of the one before it:

```
segment  →  whole  →  ensemble_long  →  ensemble  →  space
 a chunk    a full     one state         replicates  the whole
 of a run   trajectory point, verbose    collapsed   parameter sweep
```

`polyphys.manage.parser` decodes a filename into typed physical attributes *and*
derives that artifact's ancestry, so a file knows which state point and which sweep
it belongs to. Once artifacts know their lineage, `polyphys.manage.organizer` can
walk the hierarchy and reduce it level by level.

## Installation

PolyPhys is distributed from source and archived on Zenodo; it is not on PyPI.

```bash
python -m pip install "git+https://github.com/amirhs1/poly-phys.git"
```

For a development checkout:

```bash
git clone https://github.com/amirhs1/poly-phys.git
cd poly-phys
python -m pip install -e ".[dev]"
```

Requires Python 3.11+. The runtime dependencies are NumPy and pandas.

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

Geometry is derived, not just read. The cylindrical wall is built from particles of
size 1.0, so the usable confinement diameter is `D = 2r - 1.0`, and the bulk crowder
volume fraction follows from the confining volume:

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
'N200epsilon5.0r20.0lz50.0sig2.0nc1000dt0.005bdump2000adump5000ens1'
>>> artifact.ensemble
'N200D20.0ac2.0nc1000'
>>> artifact.space
'N200D20.0ac2.0'

```

`polyphys.analyze.measurer` provides the structural observables computed per frame:

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

| Module | Purpose |
|---|---|
| `manage.parser` | Filename → typed physical attributes and lineage. One `ParserBase` subclass per project, declaring its own geometry, topology, and attribute schema. |
| `manage.organizer` | Walks the lineage, combining segments into wholes and reducing replicate ensembles to averaged measurements. |
| `manage.utils` | Filename sorting, safe IO, and number-density / volume-fraction conversions for cubic and cylindrical geometries. |
| `manage.types` | Domain type aliases (`LineageT`, `GeometryT`, `PhaseT`, …) that keep the vocabulary explicit and checkable. |
| `analyze.measurer` | Per-frame structural observables — end-to-end distance, farthest-site distance, transverse size, maximum extent — plus fixed-size binning and radial, axial, azimuthal, and planar histograms. |

Eight parser subclasses ship with the package, covering cylindrical and cubic
confinement with linear, ring, and bidisperse chain topologies.

## Project status

PolyPhys is under active development toward a stable 1.0 API. It is used to produce
published research results, and the parsing and measurement layers are covered by
tests that run against Python 3.11, 3.12, and 3.13 on every push.

The examples above are executed as part of the test suite, so anything shown in this
README is guaranteed to run against the current code.

Current focus, in order:

1. Broaden test coverage and executable examples for `manage.organizer`, including
   autocorrelation-aware error estimation for correlated MD frames.
2. Expand `analyze` with the polymer- and crowding-specific observables that general
   MD toolkits do not provide.
3. Reduce the per-project coupling in the `manage` layer so a new study can *use*
   PolyPhys from outside the package rather than extend it.

## Citation

If you use PolyPhys in published work, please cite the archived release. The DOI
above resolves to the latest version; see [`CITATION.cff`](CITATION.cff) for
machine-readable metadata, or use GitHub's *Cite this repository* button.

## AI assistance

PolyPhys is developed with the help of AI coding assistants. Every change is
reviewed and merged by the maintainer, who is accountable for the result; no AI
system is an author of this software or of any release. Numerical fixtures,
formulas, and citations are derived or verified rather than accepted from a
model — see [`AI-POLICY.md`](AI-POLICY.md) for the full policy, including what
is expected of AI-assisted contributions.

## Contributing

Bug reports and questions are welcome via
[issues](https://github.com/amirhs1/poly-phys/issues). Please see
[`SECURITY.md`](SECURITY.md) for reporting suspected vulnerabilities privately.

## License

MIT — see [`LICENSE`](LICENSE).
