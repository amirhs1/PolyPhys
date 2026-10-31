# Changelog

All notable changes to PolyPhys are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the major version is `0`, the public API is not yet stable and a minor
version increase may introduce incompatible changes.

## [Unreleased]

### Added

- `polyphys.manage.lammps_log`: `read_log` and `parse_log` read a LAMMPS log
  into a `LammpsLog` record of three DataFrames: `thermo`, one row per thermo
  output step of every run in the `one`, `custom`, `multi`, or `yaml` style,
  with a `run` column; `run_stats`, one row per run with its loop time,
  performance, CPU use, MPI timing breakdown, and dangerous builds; and
  `wall_time`. `WARNING` and `fix halt` lines inside a table are skipped and
  logged, and a cut-off or malformed log raises `BrokenLogError` (#55). The
  test logs in `polyphys/tests/data/lammps_log/` are written by LAMMPS from
  the input scripts beside them. A Usage section describes the module.
- `KeywordSpec` and `Species` in `polyphys.manage.parser`, with which a
  parser class is a declaration: one `KeywordSpec` per filename keyword (the
  attribute, `int` or `float`, and an optional `scale * value + offset`
  conversion), the keyword pattern of each lineage, fixed values, and the
  species whose bulk number density and volume fraction to compute.
  `ParserBase` parses the name and computes the bulk values. A declaration
  is checked when its class is created, and an incomplete or inconsistent
  one raises `TypeError`. A new study declares a parser instead of copying
  one; the Usage page shows how
  ([#49](https://github.com/amirhs1/PolyPhys/issues/49)).
- `SumRuleCubHeteroRingBiDisperse` in `polyphys.manage.parser`, the parser
  of the project with a ring of small and large monomers among crowders of
  two sizes in a periodic cube, ported from the research code. `pctc` is the
  percentage of the total crowder volume fraction that the `ac` crowders
  carry, and `pct_two = 100 - pct_one`. For scripts that migrate from the
  research code: the bulk values are named `rho_bulk_<species>` and
  `phi_bulk_<species>` (`phi_c_bulk_one` → `phi_bulk_c_one`, `rho_c_bulk` →
  `rho_bulk_c`, `phi_m_bulk` → `phi_bulk_m`, `rho_m_bulk` → `rho_bulk_m`),
  with per-species monomer values added; the constructor is
  `(artifact, lineage, group, ispath=False)`, without `geometry` and
  `topology`; `nmon`, `ncrowd`, and `eps_others` are gone (the docstring
  states that every LJ pair has strength 1); and `mmon_large`, `mcrowd_one`,
  and `mcrowd_two` exist only at `segment`, `whole`, and `ensemble_long`
  ([#50](https://github.com/amirhs1/PolyPhys/issues/50)).
- `ParserBase._derived_attributes`, an empty method that a parser class
  overrides to set a value derived from parsed values outside the
  declaration tables; it runs at every lineage after the bulk values
  ([#50](https://github.com/amirhs1/PolyPhys/issues/50)).
- Hypothesis in the `dev` extra, for property tests that generate their
  inputs. The first one checks that every position in a histogram's range is
  counted exactly once over the edges `create_bin_edge_and_hist` builds.
  Hypothesis's built-in `ci` profile, active whenever `CI` is set, makes CI
  runs draw the same inputs every time. The Testing page describes both
  ([#83](https://github.com/amirhs1/PolyPhys/issues/83)).

- `ParserBase.dump_keywords` in `polyphys.manage.parser` maps each particle
  group to the attribute that holds the dump period, in time steps, of its
  trajectory, and `dump_keyword` gives the one for a parser's own group:
  `bdump` for `bug`, `ndump` for `nucleoid`, and `adump` for `all` in every
  project parser ([#18](https://github.com/amirhs1/PolyPhys/issues/18)).
- Python 3.14: CI runs the test suite on it, the lint, type-check, docs, and
  build jobs move to it from 3.13, and `pyproject.toml` lists its classifier
  ([#118](https://github.com/amirhs1/PolyPhys/issues/118)).

### Changed

- **Breaking:** a parser in `polyphys.manage.parser` raises `ValueError` when
  a name lacks a keyword of its lineage's pattern, naming the keyword, its
  attribute, the lineage, and the name, and when a keyword's value is not a
  number, naming the keyword and the value. Before, it printed a message and
  failed later with an `AttributeError` that named neither the keyword nor
  the file ([#47](https://github.com/amirhs1/PolyPhys/issues/47)).
- **Breaking:** a non-integral value of an integer keyword, such as `nc1.5`,
  raises `ValueError` instead of being truncated to an integer.
- **Breaking:** a parser in `polyphys.manage.parser` built from a name, with
  `ispath=False` or through `from_name`, sets `filepath` and `ext` to `None`
  instead of the string `'N/A'`, which read as a real directory or extension.
  Built from a file path, `filepath` is now a `pathlib.Path` instead of a
  `str`; `ext` stays a `str`, and is `''` for a file with no extension
  ([#100](https://github.com/amirhs1/PolyPhys/issues/100)).
- The eight parser classes of `polyphys.manage.parser` are declarations;
  `ParserBase` holds the one copy of the parsing and bulk-value code that
  each class repeated, so the per-class lists of float attributes, which had
  drifted out of step with the keyword tables, are gone. Subclasses no longer
  define `__init__`, `_initiate_attributes`, `_parse_name`, or
  `_dependent_attributes`. `ParserBase` is no longer an `abc.ABC`;
  instantiating it, or a subclass that declares no patterns, raises
  `TypeError`. For every name the previous parsers accepted, the attribute
  values and their types are unchanged. Parsing splits a name once, in time
  linear in its length
  ([#49](https://github.com/amirhs1/PolyPhys/issues/49)).
- The parser docstrings now match the code: `SumRuleCyl.dmon` and the
  `HnsCub`/`HnsCyl` time step and dump periods are fixed values, not
  keywords; `HnsCub` and `HnsCyl` read `nmon` from `N` and describe the
  `nucleoid` group; the `TransFociCyl` space pattern includes `D` and its
  wall energies are `epsilon_small` and `epsilon_large`; the hetero classes
  do not set `nmon`; `from_name` sets `filepath` and `ext` to `'N/A'`.
- `polyphys/tests/manage/test_parser.py` follows the test-writing rules:
  shared cases in fixtures, one behaviour per test, and a round-trip test
  that builds a name from every class's declared pattern.
- `polyphys.manage` no longer prints. `make_database` logs an INFO record
  when it reuses an existing directory, and `children_stamps` and
  `parents_stamps` log DEBUG records for an absent column, to the
  `polyphys.manage.organizer` logger; with no logging configured, nothing is
  shown. `openany_context` no longer catches an error raised inside its
  `with` block to print "Error opening ..." before re-raising it; the error
  propagates unchanged
  ([#47](https://github.com/amirhs1/PolyPhys/issues/47)).
- **Breaking:** the `ensemble` and `space` names that a parser in
  `polyphys.manage.parser` rebuilds print a float keyword value that is a
  whole number as an integer, and any other float as Python's `repr`, the
  shortest text that reads back as the same float. One rule holds for every
  float keyword in every class, whatever the parsed name wrote: `ac2.0`
  becomes `ac2`, `epshc0.0` becomes `epshc0`, `D20.0` becomes `D20`, and
  `pctc40.0` becomes `pctc40`; `epshc1.5` and `ac1.000001` are unchanged.
  The long names, `segment`, `whole`, and `ensemble_long`, are cut from the
  parsed name and keep its text. The new names match the SciOps generator
  scripts and the maintainer's local `analysis` directories; a short key in
  existing stamps, such as `N200kbmm2nh12ac2.0epshc0` in `HnsCub`, maps to
  the new name by parsing it at its lineage and reading the rebuilt name
  ([#102](https://github.com/amirhs1/PolyPhys/issues/102)).
- The rebuilt `segment` and `whole` names carry no topology suffix, as
  since v0.3.0, and this is now the documented rule: a segment
  `…ens1.j07.ring` and a whole `…ens1.ring` of one simulation both give the
  `whole` `…ens1`. Stamps that SciOps wrote for an artifact parsed from a
  whole file key the whole as `…ens1.ring`; append `.<topology>` to the new
  name, or drop it from the old one, to match them. Stamps that SciOps
  derived from segments already match
  ([#102](https://github.com/amirhs1/PolyPhys/issues/102)).
- A test marked `private_data` compares the rebuilt `HnsCub` names with the
  SciOps stamps under `POLYPHYS_DATA_DIR`, in any directory whose name starts
  with `HnsCub`, under the two mappings above; the comparison itself is
  tested on synthetic stamp files, and the data-reading test is excluded
  from coverage because CI never runs it
  ([#102](https://github.com/amirhs1/PolyPhys/issues/102)).

- **Breaking:** every concrete parser in `polyphys.manage.parser` declares
  `_dump_keywords`, one entry per group, each naming an attribute that the
  `segment`, `whole`, and `ensemble_long` patterns or `_defaults` set; a
  subclass that lacks it, or whose entries do not match `_groups`, raises
  `TypeError` when it is created. To migrate a parser declared outside
  PolyPhys, add `_dump_keywords`, for example
  `{'bug': 'bdump', 'all': 'adump'}`
  ([#18](https://github.com/amirhs1/PolyPhys/issues/18)).

### Removed

- **Breaking:** the `dump_key` argument of `space_tseries` in
  `polyphys.manage.organizer`. The frame index is now scaled by the dump
  period of the `group` trajectory, which the parser names through
  `dump_keyword`, and the arguments after `group` are keyword-only, so a call
  that still passes a dump keyword by position raises `TypeError`. To
  migrate, drop the dump keyword and pass `divisor`, `round_to`, `is_save`,
  and `crowder_attr` by name. A call that paired a group with another
  group's dump keyword, such as `'bug'` with `'adump'`, now gets the group's
  own period, so its `t_index` and `time` change; a call that passed the
  group's own keyword gets the same values
  ([#18](https://github.com/amirhs1/PolyPhys/issues/18)).

### Fixed

- `parents_stamps` in `polyphys.manage.organizer` drops the child lineage
  column (`segment` or `whole`) from the parent stamps when the stamps have
  no `lineage_name` column, as the stamps the probe phase writes do not.
  It removed the two columns together, so a missing `lineage_name` left the
  child column in, holding the last child's name: whole stamps built from
  segments carried a stale `segment` column, and ensemble-average stamps a
  stale `whole` column. Tables produced before this fix keep those columns
  ([#47](https://github.com/amirhs1/PolyPhys/issues/47)).
- `children_stamps(..., "whole")` drops whichever of `segment` and
  `segment_id` is present; it dropped neither when only one was there.
- The two-size parsers `SumRuleCubHeteroRing`, `SumRuleCubHeteroLinear`,
  `TransFociCyl`, and `TransFociCub` list `phi_bulk_m`, the bulk volume
  fraction of all monomers, in `attributes` at the `segment`, `whole`, and
  `ensemble_long` lineages, just before `rho_bulk_m`, as the other parsers
  do. They computed it but left it out, so their stamps carried
  `rho_bulk_m` and not `phi_bulk_m`. Stamps written from now on gain a
  `phi_bulk_m` column; no value changes
  ([#103](https://github.com/amirhs1/PolyPhys/issues/103)).

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
- A written naming convention for branches, commits, pull requests, issues,
  labels, milestones, and tags, built on two axes, `type` and `area`, and set
  out in `CONTRIBUTING.md`, "Names". The flat label set is replaced by
  `type:*` and `area:*` labels, renamed in place so existing assignments are
  kept. The issue templates prefill an `[area]` title and a `type:*` label,
  and Dependabot declares its labels explicitly, so it no longer recreates
  its default `dependencies` label
  ([#13](https://github.com/amirhs1/PolyPhys/pull/13)).

## [0.3.0] - 2025-07-10

First public pre-release, archived on Zenodo.

- Test suite introduced.

## Earlier history

Versions `0.1.0` (2022-01-28, initial project layout) and `0.2.0` (2024-11-18,
first public publication on GitHub and major structural changes) were recorded in
the previous `HISTORY.rst` as development milestones. Neither was tagged or
released, so `0.3.0` is the only release preceding this changelog.

[Unreleased]: https://github.com/amirhs1/PolyPhys/compare/v0.5.0...HEAD
[0.4.0]: https://github.com/amirhs1/PolyPhys/releases/tag/v0.4.0
[0.3.0]: https://github.com/amirhs1/PolyPhys/releases/tag/v0.3.0
