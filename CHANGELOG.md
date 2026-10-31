# Changelog

All notable changes to PolyPhys are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the major version is `0`, the public API is not yet stable and a minor
version increase may introduce incompatible changes.

## [0.5.0](https://github.com/amirhs1/PolyPhys/releases/tag/v0.5.0) - 2026-10-31

### Changed

- **Breaking:** Rename `polyphys.manage.utils.round_up_nearest` to `round_to_multiple`, which rounds to the nearest multiple of `divider`, ties to even; values are unchanged ([#96](https://github.com/amirhs1/PolyPhys/pull/96), `98fe6bb`)
- Raise `ValueError` in `whole_from_segments` when a whole's segments are not in strictly increasing `segment_id` order, or `relation` is unknown ([#93](https://github.com/amirhs1/PolyPhys/pull/93), `144e691`)
- Run every `polyphys.manage.organizer` example as a doctest, and state each reducer's time and memory cost ([#93](https://github.com/amirhs1/PolyPhys/pull/93), `144e691`)
- Cover every public organizer function in `test_organizer.py`, raising its line coverage from 65% to 98% ([#93](https://github.com/amirhs1/PolyPhys/pull/93), `144e691`)
- Run pytest with `--strict-markers`, `--strict-config`, `xfail_strict = true`, `-ra`, and `--import-mode=importlib` ([#87](https://github.com/amirhs1/PolyPhys/pull/87), `1a08b44`)
- Type-check `polyphys/analyze` and `polyphys/manage` under `mypy --strict`, set in `pyproject.toml` ([#97](https://github.com/amirhs1/PolyPhys/pull/97), `af5266c`, `7f3a08b`)
- Fail the CI test job when total coverage falls below 95% ([#98](https://github.com/amirhs1/PolyPhys/pull/98), `6e9f3d9`)
- Take `html_baseurl` from Read the Docs' `READTHEDOCS_CANONICAL_URL` ([#91](https://github.com/amirhs1/PolyPhys/pull/91), `017fcb5`)
- Stop ignoring `*.png` in `.gitignore`, so documentation figures can be committed ([#88](https://github.com/amirhs1/PolyPhys/pull/88), `60bb1a7`)

### Added

- Add `polyphys.analyze.stats`, Flyvbjerg-Petersen block averaging for the uncertainty of a time average over correlated frames, warning when a run fails the acceptance criterion of Lee et al. (2011) ([#88](https://github.com/amirhs1/PolyPhys/pull/88), [#90](https://github.com/amirhs1/PolyPhys/pull/90), `60bb1a7`, `db3f43b`)
- Add a Statistics documentation page, with figures from `docs/scripts/make_stats_figures.py` and a demo in `docs/scripts/stats_demo.py` ([#88](https://github.com/amirhs1/PolyPhys/pull/88), `60bb1a7`)
- Add the shared test fixtures `make_lineage_tree`, `make_lammps_data`, and `make_lammps_dump` in `polyphys/tests/conftest.py` ([#87](https://github.com/amirhs1/PolyPhys/pull/87), `1a08b44`, `ad65894`)
- Add a `private_data` pytest marker, skipped when `POLYPHYS_DATA_DIR` is unset ([#87](https://github.com/amirhs1/PolyPhys/pull/87), `1a08b44`)
- Add a Testing documentation page ([#87](https://github.com/amirhs1/PolyPhys/pull/87), `1a08b44`)
- Publish the documentation on Read the Docs at https://polyphys.readthedocs.io/ ([#91](https://github.com/amirhs1/PolyPhys/pull/91), `017fcb5`)

### Removed

- **Breaking:** Remove the `bin_range` and `bin_ranges` arguments of the `polyphys.analyze.measurer` histogram functions, which NumPy ignored; counts are unchanged ([#96](https://github.com/amirhs1/PolyPhys/pull/96), `4f01879`)

### Fixed

- **Breaking:** Raise `ValueError` in the `measurer` histogram functions when positions fall outside the edges, instead of dropping them silently ([#96](https://github.com/amirhs1/PolyPhys/pull/96), `4f01879`)
- Make the end edges of `create_bin_edge_and_hist` and `fixedsize_bins` enclose `lmin` and `lmax` exactly, with no extra edge ([#96](https://github.com/amirhs1/PolyPhys/pull/96), `4f01879`)
- **Breaking:** Raise `ValueError` in the `periodic` branch of `fixedsize_bins` when the period is not a whole number of bins ([#96](https://github.com/amirhs1/PolyPhys/pull/96), `4f01879`)
- Keep bins of width `bin_size` in the `nonnegative` branch of `fixedsize_bins` when it clips at zero, and reject a negative `lmin` ([#96](https://github.com/amirhs1/PolyPhys/pull/96), `4f01879`)
- Render the Statistics page's mathematics by enabling MyST `dollarmath` ([#95](https://github.com/amirhs1/PolyPhys/pull/95), `090e17f`)

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
  points could be averaged together. Integer-valued inputs now give `ensemble`
  and `space` names such as `kbmm2.0` and `epshc1.0`, matching `HnsCyl`
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
- The parsers rebuild the `segment`, `whole`, and `ensemble_long` names by
  cutting the artifact's own name after the parent's last keyword, instead of
  writing converted attributes back under the filename keywords. Every parser
  converts some keywords (`r` to `dcyl = 2r - 1`, `lz` to `lcyl = 2 lz`, and
  `l` or `hl` to `lcube = 2 l`), so each re-parse converted them again: a
  `SumRuleCyl` segment with `r10.5lz25.0` reached `ensemble_long` with
  `dcyl = 77.0` and `lcyl = 200.0`, which makes the volume fraction of
  crowders of size 1 64 times too small, and `ensemble`, `space_tseries`, and `space_hists` carried the drift
  into their output. The rebuilt names now keep the filename's tokens and
  number formatting, such as `r10.5` and `.j01`. The `ensemble` and `space`
  names are built from the attributes as before
  ([#85](https://github.com/amirhs1/PolyPhys/issues/85)).

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

[0.4.0]: https://github.com/amirhs1/PolyPhys/releases/tag/v0.4.0
[0.3.0]: https://github.com/amirhs1/PolyPhys/releases/tag/v0.3.0
