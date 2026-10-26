# Changelog

All notable changes to PolyPhys are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the major version is `0`, the public API is not yet stable and a minor
version increase may introduce incompatible changes.

## [Unreleased]

### Added

- Shared test fixtures in `polyphys/tests/conftest.py`: `make_lineage_tree`
  writes a synthetic lineage tree of tiny `.npy` or `.csv` files named for any
  parser class, and `make_lammps_data` and `make_lammps_dump` write tiny LAMMPS
  data files (`atomic`, `bond`, and `angle` atom styles) and dump files. Their
  tests are in `polyphys/tests/test_conftest.py`.
- A `private_data` pytest marker for tests that read the maintainer's local
  data: such a test is skipped when `POLYPHYS_DATA_DIR` is unset and errors
  when it names something that is not a directory.
- A Testing page in the documentation, describing the three test-data tiers,
  the kinds of test, the shared fixtures, and the pytest settings.

### Changed

- pytest now runs with `--strict-markers`, `--strict-config`,
  `xfail_strict = true`, `-ra`, and `--import-mode=importlib`, and loads the
  `pytester` plugin. A misspelt marker is an error rather than a warning.

## [0.4.0] - 2026-10-26

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
- `space_tseries` and `space_hists` in `polyphys.manage.organizer` find their
  `*-<prop>-ensAvg.csv` input files again. They passed that suffix to
  `sort_filenames` as a file extension, so no file matched and both raised
  `ValueError: No objects to concatenate`. They also pass `ispath=True` to the
  parser, which otherwise split the full path at its first hyphen. They read
  the crowder volume fraction from a new keyword parameter, `crowder_attr`,
  defaulting to the parsers' `phi_bulk_c`, instead of the research-code name
  `phi_c_bulk`, and add that column themselves, so it no longer has to be
  listed in `physical_attrs`
  ([#44](https://github.com/amirhs1/PolyPhys/issues/44)).
- `ens_from_bin_edge` and `ens_avg_from_bin_edge` raise `ValueError` when the
  wholes of an ensemble have different bin edges, instead of merging them with
  `np.unique` into a longer array that misaligned the paired histogram.
  `ens_avg_from_df` no longer modifies the caller's DataFrame; it returns a
  new one ([#45](https://github.com/amirhs1/PolyPhys/issues/45)).
- `fixedsize_bins` now allocates `collector` as `int64` and `collector_std` as
  `float64`, and `create_bin_edge_and_hist` allocates its histogram as `int64`.
  They were `int16`, so per-frame counts added in place wrapped silently past
  32,767 counts per bin ([#43](https://github.com/amirhs1/PolyPhys/issues/43)).
- `HnsCub` now parses `bend_mm` (`kbmm`) and `eps_hc` (`epshc`) as floats
  instead of truncating them to integers. Two filenames that differed only in a
  non-integer value, such as `epshc1.5` and `epshc1.0`, produced identical
  `whole`, `ensemble_long`, `ensemble`, and `space` names, so distinct state
  points could be averaged together. Integer-valued inputs now give names such
  as `kbmm2.0` and `epshc1.0`, matching `HnsCyl`
  ([#41](https://github.com/amirhs1/PolyPhys/issues/41)).
- Corrected the `attributes` tables of the four hetero parser classes, which
  listed `phi_bulk_m_large` and `rho_bulk_m` twice and omitted
  `rho_bulk_m_large`. `HnsCub` and `HnsCyl` now set the `ndump` attribute
  their tables list, in place of `bdump`
  ([#42](https://github.com/amirhs1/PolyPhys/issues/42)).
- `whole_from_segments`, `whole_from_file`, and `whole_from_dist_mat_t` in
  `polyphys.manage.organizer` pass `ispath=True` to the parser. Without it the
  parser took the full path as the name and split it at its first hyphen, so
  any file under a hyphenated directory, as every lineage directory is, raised
  `AttributeError` ([#66](https://github.com/amirhs1/PolyPhys/issues/66)).
- `ensemble` in `polyphys.manage.organizer` again adds a `bin_center` column,
  the midpoints of the shared bin edges, to each histogram ensemble when
  `whole_edges` is given. It computed the centres and then discarded them, so
  `whole_edges` had no effect and ensemble-averaged histograms had no x-axis.
  The centres come from the first whole of each ensemble, whose wholes share
  their edges by design. Passing `whole_edges` with a `whole_type` other than
  `"vector"` now raises `ValueError`
  ([#70](https://github.com/amirhs1/PolyPhys/issues/70)).

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
- `AI-POLICY.md`, stating how AI tools may be used here: whoever submits a
  change is responsible for it, test values and citations are verified
  independently, and AI-assisted commits carry an `Assisted-by:` trailer. The
  README records each part's AI tier and checks; `AGENTS.md` and `CLAUDE.md`
  carry the matching agent rules.

## [0.3.0] - 2025-07-10

First public pre-release, archived on Zenodo.

- Test suite introduced.

## Earlier history

Versions `0.1.0` (2022-01-28, initial project layout) and `0.2.0` (2024-11-18,
first public publication on GitHub and major structural changes) were recorded in
the previous `HISTORY.rst` as development milestones. Neither was tagged or
released, so `0.3.0` is the only release preceding this changelog.

[Unreleased]: https://github.com/amirhs1/PolyPhys/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/amirhs1/PolyPhys/releases/tag/v0.4.0
[0.3.0]: https://github.com/amirhs1/PolyPhys/releases/tag/v0.3.0
