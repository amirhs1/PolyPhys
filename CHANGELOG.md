# Changelog

All notable changes to PolyPhys are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the major version is `0`, the public API is not yet stable and a minor
version increase may introduce incompatible changes.

## [Unreleased]

Targeting `0.4.0`.

### Fixed

- Removed the `polyphys = polyphys.cli:main` console-script entry point. The
  `polyphys.cli` module does not exist, so every install shipped a command that
  raised `ModuleNotFoundError` when run. A CLI will be reintroduced when there is
  a module behind it.
- Removed the `[tool.setuptools.package-data]` entry for `polyphys.test_data`,
  which is not a package in this tree.
- Reduced the runtime dependencies from ten to the two the package actually
  imports, `numpy` and `pandas`. The removed entries — `scipy`, `matplotlib`,
  `seaborn`, `MDAnalysis`, `statsmodels`, `pyarrow`, `packaging`, and `click` —
  were never imported anywhere in the package.
- Corrected the project URLs and README badges, which pointed at the
  pre-rename `amirhs1/PolyPhys` path. The coverage badge additionally tracked a
  stale `master` branch and reported a coverage figure that no longer reflected
  the default branch.
- Documentation no longer instructs users to `pip install polyphys`; the package
  is distributed from source and archived on Zenodo, not published to PyPI.
- `docs/source/conf.py` now uses the standard-library `tomllib` instead of the
  undeclared third-party `tomli` dependency.

### Added

- `polyphys/tests/test_packaging.py`, validating the packaging contract against
  installed distribution metadata: declared console scripts must resolve to a
  callable, and `__version__` must be valid semantic versioning.
- `CITATION.cff`, providing machine-readable citation metadata and enabling
  GitHub's *Cite this repository* button.
- This changelog.
- A README that states the problem the package solves, documents installation
  and the artifact-lineage model, and carries a quickstart whose examples are
  executed by the test suite.

## [0.3.0] - 2025-07-10

First public pre-release, archived on Zenodo.

- Test suite introduced.

## Earlier history

Versions `0.1.0` (2022-01-28, initial project layout) and `0.2.0` (2024-11-18,
first public publication on GitHub and major structural changes) were recorded in
the previous `HISTORY.rst` as development milestones. Neither was tagged or
released, so `0.3.0` is the only release preceding this changelog.

[Unreleased]: https://github.com/amirhs1/poly-phys/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/amirhs1/poly-phys/releases/tag/v0.3.0
