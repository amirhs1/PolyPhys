# Changelog

All notable changes to PolyPhys are documented here. The format is
[Common Changelog](https://common-changelog.org/), with one departure: an
`[Unreleased]` section collects changes between releases, and each release
adds the commit references its entries lack. This project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

From 1.0.0, the public API and its deprecation policy are those of the
"API stability" page of the documentation (`docs/source/stability.md`). A
deprecation is listed under "Changed" and starts with `**Deprecated:**`.

## [Unreleased]

### Added

- Declare the public API, the `__all__` of every public module, with each module stable or provisional, and a deprecation policy under Semantic Versioning on the "API stability" page ([#81](https://github.com/amirhs1/PolyPhys/issues/81))

### Removed

- **Breaking:** Remove `polyphys.manage.projects` and the table records and builders that `polyphys.analyze.pipeline` served, the import paths that 0.9.0 deprecated; import them from `polyphys.spec` ([#81](https://github.com/amirhs1/PolyPhys/issues/81))
- **Breaking:** Remove `whole_from_dist_mat_t` and `find_unique_properties` from `polyphys.manage.organizer`, which nothing in PolyPhys called ([#81](https://github.com/amirhs1/PolyPhys/issues/81))

### Migrating to 1.0.0

Each breaking change since 0.4.0, oldest last, with its migration path.

| Release | Change | Migration |
| ------- | ------ | --------- |
| 1.0.0 | `polyphys.manage.projects` and the table names of `polyphys.analyze.pipeline` removed | Import them from `polyphys.spec`. |
| 1.0.0 | `whole_from_dist_mat_t` and `find_unique_properties` removed | None; nothing in PolyPhys called them. |
| 0.9.0 | `ProjectSpec` requires `crowder_attr`, keyword-only, with no default | Pass `crowder_attr='phi_bulk_c'` to keep the old behavior. |
| 0.9.0 | `polyphys.manage.utils.sort_filenames_str` removed | Use `sort_filenames`. |
| 0.9.0 | Stamps tables are named `<space>-<group>-stamps-<level>.csv` | New runs write the new names; tables already on disk keep their old names, so code that finds them by name matches both. |
| 0.9.0 | The `lineage` argument of `children_stamps` and `parents_stamps` is `child_lineage` | Rename the keyword. |
| 0.9.0 | `make_database` raises `ValueError` for a `group` with `-` and for an old directory outside a phase directory | Use a group name without `-`, and pass an old directory whose parent is a phase directory. |
| 0.9.0 | `whole_from_segments` raises `ValueError` when a whole skips a `segment_id` | Supply the missing segment, or leave the incomplete whole out. |
| 0.9.0 | `ens_from_vec`, `ens_from_df`, and `reduce_space` raise `ValueError` for wholes of unequal lengths | Trim the wholes to a common length, or leave the short ones out, before reducing. |
| 0.9.0 | `parents_stamps` and `sum_rule_profiles` raise `ValueError` for NaN in a group key, and `parents_stamps` for a `properties` entry it aggregates itself | Fill or drop the NaN rows before the call, and leave those columns out of `properties`. |
| 0.9.0 | `whole_from_file` raises `ValueError` when two files parse to one whole | Pass one file per whole. |
| 0.9.0 | `volume_fraction_cube` and `volume_fraction_cylinder` return a fraction above one instead of raising; `number_density_cube` and `number_density_cylinder` raise `ValueError` when a particle leaves its center no room | Check `phi > 1` where the old error was relied on, and catch the new `ValueError` of the number densities. |
| 0.8.0 | `RecentreT`, `recentre`, and `'centre'` are `RecenterT`, `recenter`, and `'center'` | Rename them; no file on disk stores these names. |
| 0.7.0 | The rounded crowder volume fraction column is `phi_bulk_c_round` | Rename the column; `rename_bulk_columns` maps the research-code bulk names, with `BULK_NAME_MAP`. |
| 0.6.0 | A parser raises `ValueError` for a name that lacks a keyword or has a non-numeric value | Catch `ValueError` instead of `AttributeError`. |
| 0.6.0 | A parser raises `ValueError` for a non-integral integer keyword, such as `nc1.5` | Fix the name; the value is no longer truncated. |
| 0.6.0 | A parser built from a name has `filepath` and `ext` of `None`, not `'N/A'`; `filepath` is a `pathlib.Path` | Test for `None`, and call `str()` where a string path is needed. |
| 0.6.0 | Rebuilt `ensemble` and `space` names print a whole-number float as an integer, such as `ac2` for `ac2.0` | Map an old short key by parsing it at its lineage and reading the rebuilt name. |
| 0.6.0 | Every concrete parser needs `_dump_keywords` | Add one entry per group, such as `{'bug': 'bdump', 'all': 'adump'}`. |
| 0.6.0 | `PhaseT` holds the phase vocabulary; `make_database` checks `phase` against `PhaseT` and `stage` against `StageT` | Pass a phase and a stage from those vocabularies. |
| 0.6.0 | `space_tseries` has no `dump_key`, and its arguments after `group` are keyword-only | Drop the dump keyword and pass `divisor`, `round_to`, `is_save`, and `crowder_attr` by name. |
| 0.5.0 | `round_up_nearest` is `round_to_multiple` | Rename the call; values are unchanged. |
| 0.5.0 | The `measurer` histogram functions have no `bin_range` or `bin_ranges` | Drop the argument; counts are unchanged. |
| 0.5.0 | The `measurer` histogram functions raise `ValueError` for positions outside the edges | Widen the edges, or filter the positions, before binning. |
| 0.5.0 | `fixedsize_bins` raises `ValueError` in its `periodic` branch when the period is not a whole number of bins | Choose a bin size that divides the period. |

## [0.9.0](https://github.com/amirhs1/PolyPhys/releases/tag/v0.9.0) - 2026-11-07

### Changed

- **Breaking:** Require `ProjectSpec` to name `crowder_attr`, the parser attribute of the bulk volume fraction of crowders, as a keyword-only field with no default, where it was always `phi_bulk_c` before ([#195](https://github.com/amirhs1/PolyPhys/issues/195), [#204](https://github.com/amirhs1/PolyPhys/pull/204), `f40d710`)
- **Breaking:** Remove `polyphys.manage.utils.sort_filenames_str`, which nothing called and which zipped separately sorted lists, so one missing file shifted every later pair, and the unused private alias `_normalize_exts`; use `sort_filenames` ([#177](https://github.com/amirhs1/PolyPhys/issues/177), [#206](https://github.com/amirhs1/PolyPhys/pull/206), `6fdeb17`)
- **Breaking:** Name every stamps table `<space>-<group>-stamps-<level>.csv`, with `<level>` the lineage of its rows: `children_stamps` and `parents_stamps` write `-stamps-segment.csv`, `-stamps-whole.csv`, and `-stamps-ensemble_long.csv` instead of `-segment-stamps.csv`, `-whole-stamps.csv` or `-stamps-whole.csv` for the same `ens` table, and `-stamps-ensAvg.csv`; existing analysis directories keep their old names ([#92](https://github.com/amirhs1/PolyPhys/issues/92), [#206](https://github.com/amirhs1/PolyPhys/pull/206), `179d0c4`)
- **Breaking:** Rename the `lineage` argument of `children_stamps` and `parents_stamps` to `child_lineage`, with no alias, since it is the level of the child stamps ([#92](https://github.com/amirhs1/PolyPhys/issues/92), [#206](https://github.com/amirhs1/PolyPhys/pull/206), `179d0c4`)
- **Breaking:** Raise `ValueError` in `make_database` for a `group` that contains `-` and for an old directory whose parent is not a phase directory, before any directory is made ([#92](https://github.com/amirhs1/PolyPhys/issues/92), [#206](https://github.com/amirhs1/PolyPhys/pull/206), `179d0c4`)
- **Breaking:** Raise `ValueError` in `whole_from_segments` when the segments of a whole skip a `segment_id`, such as `j01` then `j03`, instead of joining them as if adjacent ([#180](https://github.com/amirhs1/PolyPhys/issues/180), [#207](https://github.com/amirhs1/PolyPhys/pull/207), `10155d9`)
- **Breaking:** Raise `ValueError` naming the ensemble and each whole's length when the wholes of an ensemble have unequal lengths, in `ens_from_vec`, `ens_from_df`, and the `'dataframe'` ensembles of `reduce_space`, instead of averaging the tail over fewer replicates or taking the first whole's rows ([#180](https://github.com/amirhs1/PolyPhys/issues/180), [#207](https://github.com/amirhs1/PolyPhys/pull/207), `10155d9`)
- **Breaking:** Raise `ValueError` in `parents_stamps` and `sum_rule_profiles` for a row with NaN in its group key instead of dropping it, and in `parents_stamps` for a `properties` entry that it aggregates itself instead of overwriting it ([#180](https://github.com/amirhs1/PolyPhys/issues/180), [#207](https://github.com/amirhs1/PolyPhys/pull/207), `10155d9`)
- **Breaking:** Raise `ValueError` naming both files in `whole_from_file` when two files parse to one whole name, instead of keeping the last ([#180](https://github.com/amirhs1/PolyPhys/issues/180), [#207](https://github.com/amirhs1/PolyPhys/pull/207), `10155d9`)
- **Breaking:** Return a bulk volume fraction above one from `volume_fraction_cube` and `volume_fraction_cylinder`, and so from the parsers, as the nominal fraction it is, instead of raising `ValueError`; `number_density_cube` and `number_density_cylinder` raise `ValueError` when a particle leaves its center no room, and a parser names the species, the artifact, and the lineage in it, instead of dividing by zero ([#181](https://github.com/amirhs1/PolyPhys/issues/181), [#201](https://github.com/amirhs1/PolyPhys/pull/201), `75cd943`)
- Raise `ValueError` in `ChainSize.run` on a frame in which a chain bond spans half the box or more, since every observable except an unwrapped `asphericityT` needs the dump to hold the chain whole; the check is skipped without a box or bonds ([#183](https://github.com/amirhs1/PolyPhys/issues/183), [#203](https://github.com/amirhs1/PolyPhys/pull/203), `9cf6dd4`)
- Move `ProjectSpec`, `ProbeSpec`, `BinderSpec`, `PerWhole`, `RecenterT`, `CHAIN_SIZE_OBSERVABLES`, and `BINDER_SERIES` from `polyphys.manage.projects`, and `AnalysisTable`, `PropertySpec`, `PostHook`, `tseries`, `histogram`, `histogram_2d`, `edge`, and `matrix_tseries` from `polyphys.analyze.pipeline`, to the new module `polyphys.spec`, which imports no MDAnalysis, and define `PerWhole` once; the old import paths warn with `DeprecationWarning` until 1.0.0 ([#196](https://github.com/amirhs1/PolyPhys/issues/196), [#198](https://github.com/amirhs1/PolyPhys/pull/198), `818aca9`)
- Read the tier C test files from one folder per study under `POLYPHYS_DATA_DIR`, which now names the maintainer's `polyphys_test_files` folder, each with the `make_database` phase folders: the probe real frames from `<project>/trjs/<whole>/` instead of `polyphys_test_files/probe_real_frames/<whole>/`, and the density golden inputs from `TransFociCyl/golden/density/<direction>/` instead of `density/<direction>/` ([#213](https://github.com/amirhs1/PolyPhys/pull/213), `cb55905`)
- Exclude the tier C tests and their fixtures from the coverage measured without `POLYPHYS_DATA_DIR`, as CI runs, with `# pragma: no cover - needs data` on each `private_data` marker, so tests CI cannot run no longer count as uncovered ([#216](https://github.com/amirhs1/PolyPhys/pull/216), `d78e50b`)

### Added

- Add `space_tseries_from_spec` and `space_hists_from_spec` to `polyphys.analyze.pipeline`, which build the `space` tables of `organizer.space_tseries` and `organizer.space_hists` with the parser, the glob prefix, the columns, the rounding step, and the crowder attribute of a `ProjectSpec`, so no study choice comes from the organizer's defaults ([#195](https://github.com/amirhs1/PolyPhys/issues/195), [#204](https://github.com/amirhs1/PolyPhys/pull/204), `f40d710`)
- Add `space_profiles_from_spec` to `polyphys.analyze.sumrule`, which builds the sum-rule profiles of `space_profiles` with the parser, the glob prefix, the columns, the species sizes, the rounding step, and the crowder attribute of a `ProjectSpec`, so no study choice comes from the defaults of `space_profiles` ([#209](https://github.com/amirhs1/PolyPhys/issues/209), [#211](https://github.com/amirhs1/PolyPhys/pull/211), `10ff1a3`)
- Add a tier A test that runs `ChainSize` on chains written in the test and reduces its files with `reduce_space`, checking the stage files' names and the ensemble mean; `make_chain_files` moves to the shared test `conftest.py` ([#191](https://github.com/amirhs1/PolyPhys/issues/191), [#208](https://github.com/amirhs1/PolyPhys/pull/208), `907ad4c`)
- Add a release step: the maintainer runs the full test command with `POLYPHYS_DATA_DIR` set and LAMMPS available, since CI runs neither the `private_data` nor the `lammps` tests, and the release commit records the result in a `Checks-run:` line; in `CONTRIBUTING.md`, "Releases", and the `draft-release` skill ([#190](https://github.com/amirhs1/PolyPhys/issues/190), [#200](https://github.com/amirhs1/PolyPhys/pull/200), `f413c1b`)
- Add tier C tests that compare each stage of the TransFociCub `bug` pipeline with the research code's output on the four tier B wholes: the chain-size observables, the principal axes up to sign, the foci bonds and clusters histograms, and the stamps against the research probe; `reduce_space` on the research probe, and on the package's own probe of the same dumps, against the research `ensAvg`; `space_tseries_from_spec` against the research space tables; and a figure of a space table with the viz helpers ([#215](https://github.com/amirhs1/PolyPhys/pull/215), `0d9acc2`)

### Fixed

- Return the rows of `sum_rule_profiles` in the input order, as documented, when a state point's rows are not contiguous ([#180](https://github.com/amirhs1/PolyPhys/issues/180), [#207](https://github.com/amirhs1/PolyPhys/pull/207), `10155d9`)
- Convert a `KeywordSpec` value with a scale or an offset in exact decimal arithmetic, so `r4.1` gives the cylinder size `7.2` and the short name `D7.2` instead of `D7.199999999999999` ([#181](https://github.com/amirhs1/PolyPhys/issues/181), [#201](https://github.com/amirhs1/PolyPhys/pull/201), `75cd943`)
- Strip an extension with `ispath=True` only when the suffix holds a letter, so a name ending with a value, such as `ac1.5`, or with a dot keeps it ([#181](https://github.com/amirhs1/PolyPhys/issues/181), [#201](https://github.com/amirhs1/PolyPhys/pull/201), `75cd943`)
- Reduce the research code's probe stamps, which have no `group` column: `children_stamps` adds the column from the file names, `<name>-<group>-stamps.csv`, instead of raising `KeyError`, so `reduce_space` runs on an archived `-probe` folder; `parents_stamps` raises `ValueError` naming the missing `space` or `group` column when asked to save without one ([#214](https://github.com/amirhs1/PolyPhys/pull/214), `eafda92`)
- Key each file in `sort_filenames` on its full name without the matched extension instead of the part before its last dot, so the HnsCyl names `...epshc1.0nc0` and `...epshc1.0nc500` no longer collide and a `.data` file no longer pairs with the `.lammpstrj` file of another state point; `space_tseries` and `space_hists` join their glob with `os.path.join` ([#177](https://github.com/amirhs1/PolyPhys/issues/177), [#206](https://github.com/amirhs1/PolyPhys/pull/206), `6fdeb17`)
- Describe the lineage, phases, stages, and naming patterns of the `organizer` module docstring and of `make_database` as the code builds them, and state in `lineage.md` which vocabulary a shared word belongs to, the grammar of a name, and the stamps tables ([#92](https://github.com/amirhs1/PolyPhys/issues/92), [#206](https://github.com/amirhs1/PolyPhys/pull/206), `179d0c4`)
- Report a `fit_decay` fit as failed, with `converged=False` and NaN parameters, when the model is not finite at the start or the fit does not move from its start and has no finite covariance, instead of reporting it converged under the default method and raising `ValueError` under `method='trf'` ([#178](https://github.com/amirhs1/PolyPhys/issues/178), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `a470534`)
- Count the bins of `pair_distance_histograms` and `FociPairDistances` with the round-off tolerance of `create_bin_edge_and_hist`, so a range of a whole number of bins, such as `r_max=2.1` with `bin_size=0.3`, ends at `r_max` instead of one bin past it or an ulp below it ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Treat an `inf` box length as non-periodic and accept integer positions in `measurer.pair_distance`, as `clusters.pair_distances` does ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Format the axis and shape in the `IndexError` of `transverse_size` and `fsd` ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Raise `ValueError` for a single sample in `simple_stats` and `sem` instead of returning NaN ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Raise `ValueError` for a cis threshold that is not a whole number in `threshold_sweep` instead of truncating it ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Correct the memory note of `pair_distances` to the four `(n, n)` arrays it holds at once ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Reject bin centers that are not symmetric about `z = 0` when `normalize_z` folds, and fold `<prop>-scale-normalized_curve` by summing mirror bins, so it is the distribution of `|z|` and sums to one ([#184](https://github.com/amirhs1/PolyPhys/issues/184), [#202](https://github.com/amirhs1/PolyPhys/pull/202), `9b4dc51`)
- Reject an `AnalysisTable` whose hooks declare one output name twice, or an output that another entry also yields, through a new `outputs(prop)` method of `TimeAverage`, `Autocorrelation`, `LocalDensities`, and `FociPairDistances`, and raise in `reduce_space` when a hook that declares none repeats a name, instead of overwriting the earlier output ([#185](https://github.com/amirhs1/PolyPhys/issues/185), [#205](https://github.com/amirhs1/PolyPhys/pull/205), `d45ec7a`)
- Select every property's files in `reduce_space` before any analysis directory is made, so a typo in a table raises `FileNotFoundError` with nothing written ([#185](https://github.com/amirhs1/PolyPhys/issues/185), [#205](https://github.com/amirhs1/PolyPhys/pull/205), `d45ec7a`)
- Raise `BrokenLogError` from `parse_log` for a `Loop time of` line of another form than LAMMPS writes and for a thermo header that repeats a keyword, and `BrokenDumpError` from the dump reader for an atom line of the wrong length even when the frame's total is right; `with_frame` documents that it needs wrapped coordinates and says so when a frame has unwrapped ones only ([#179](https://github.com/amirhs1/PolyPhys/issues/179), [#199](https://github.com/amirhs1/PolyPhys/pull/199), `817f5a1`)
- Take one color, a name or an RGB(A) tuple, for every label in `color_handles`, and one NumPy integer or dash tuple for every label in `line_handles`, instead of reading it as one entry per label ([#186](https://github.com/amirhs1/PolyPhys/issues/186), [#199](https://github.com/amirhs1/PolyPhys/pull/199), `2bfe9de`)
- Compute each part of a split probe run on its own copy of the analysis, so that `ChainSize`, `BinderStatistics`, `SpatialHistograms`, and `PairDistances` with the `serial` backend and `n_parts > 1` give the results, frames, and times of the other backends instead of repeating the last part ([#182](https://github.com/amirhs1/PolyPhys/issues/182), [#203](https://github.com/amirhs1/PolyPhys/pull/203), `35780aa`)
- Raise `ValueError` instead of `TypeError` in `PairDistances` on a frame with no box when a direction is periodic, and accept NumPy booleans as its periodic flags, as `BinderStatistics` does ([#182](https://github.com/amirhs1/PolyPhys/issues/182), [#203](https://github.com/amirhs1/PolyPhys/pull/203), `35780aa`)
- Correct the numbers the documentation states: the error of the sphere kernel (about 0.047 a/R, not 0.09), the conservation error of the centered form (up to 2.65 % cylindrical and 7.27 % spherical at h = a/4), the outermost radial bins where a uniform fluid falls short of ρv, the two-level gap in the AR(1) figure, the `shapeT` range of −1/4 to 2, the 313 bytes per frame of `PairDistances`, the n/50 rule of the `select_plateau` warning, and the cost of the cylindrical share matrix; the half-space test of `sphere_cylinder_intersection` tightens to the corrected error ([#187](https://github.com/amirhs1/PolyPhys/issues/187), [#212](https://github.com/amirhs1/PolyPhys/pull/212), `db7da9d`)
- Install the `viz` extra in the documented development installs, which the doctests of `polyphys/viz/` need; reduce each analysis table from the probe directory of its own group, under the space name the parser gives, with `projects/` on the import path, in the workflow and pipeline examples; state that the recentering triple reproduces the research shift but not its histogram ranges; pin `usage.md` to `v0.8.0`; and drop stale text: the finished items of the README's current focus, mypy for `polyphys/viz/` in the README's tier table, pull-request wording on the probe page, and a check of `bond_cosine_correlation` against `PersistenceLength` that no test makes ([#188](https://github.com/amirhs1/PolyPhys/issues/188), [#212](https://github.com/amirhs1/PolyPhys/pull/212), `d81d56d`)
- Make the pages agree with the code and fill their gaps: the study parsers that still ship in `manage.parser`, the `density` module for spatial densities, `wholeSim` for `segment` tables only, where `space_hists` takes `bin_center`, the tier B real-frame checks, the one fixture that is not a factory, the doctest options and `test_docs.py`, the `-norm-sem` column of density profiles, the stamps tables and functions of workflow step 4, `provisional` and `space_hierarchy` of `ProjectSpec`, a table of the `measurer` functions, and how to build the `DensityKernel` that `LocalDensities` needs, whose docstring no longer says a parser holds the number of frames ([#189](https://github.com/amirhs1/PolyPhys/issues/189), [#212](https://github.com/amirhs1/PolyPhys/pull/212), `40315cb`)

## [0.8.0](https://github.com/amirhs1/PolyPhys/releases/tag/v0.8.0) - 2026-11-06

### Changed

- **Breaking:** Rename `polyphys.manage.projects.RecentreT` to `RecenterT`, the `recentre` field of `ProbeSpec` and the `recentre` keyword and attribute of `SpatialHistograms` to `recenter`, and the `'centre'` value of `polyphys.analyze.density.ShareMode`, which `share_matrix` and `density_profiles` take as `shares`, to `'center'`, with no alias; no file the probes or the organizer write stores these names ([#168](https://github.com/amirhs1/PolyPhys/issues/168), [#172](https://github.com/amirhs1/PolyPhys/pull/172), `04115c5`)
- Adopt American English as the spelling of PolyPhys in `AGENTS.md` ([#168](https://github.com/amirhs1/PolyPhys/issues/168), [#167](https://github.com/amirhs1/PolyPhys/pull/167), `17c918e`)
- Convert the prose, comments, docstrings, test docstrings, and private and test names to American English spelling ([#168](https://github.com/amirhs1/PolyPhys/issues/168), [#170](https://github.com/amirhs1/PolyPhys/pull/170), `9c3a17f`)
- Take the commit template's `Why:` line from the issue, the maintainer, or an outside report, mark a breaking change with the changelog's `**Breaking:**` instead of `!`, and allow `Refs #n` in the pull request template ([#166](https://github.com/amirhs1/PolyPhys/issues/166), [#167](https://github.com/amirhs1/PolyPhys/pull/167), `3e7fe0b`)
- Read the tier C probe test files from `<POLYPHYS_DATA_DIR>/polyphys_test_files/probe_real_frames/<whole>/` instead of `<POLYPHYS_DATA_DIR>/fixtures/<whole>/`; `POLYPHYS_DATA_DIR` stays the datasets root ([#174](https://github.com/amirhs1/PolyPhys/pull/174), `5092443`)

### Added

- Add probe tests on real frames of the four live studies: 50 tier B files from 22 wholes in `polyphys/tests/data/probe/`, polymer dumps opened with the data file of the uncrowded whole of their family, and a tier C half on the crowded windows under `POLYPHYS_DATA_DIR`; `testing.md` describes the layout and raises the tier B size guide to about 300 KB ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#173](https://github.com/amirhs1/PolyPhys/pull/173), `958f0bb`)
- Add `polyphys.viz`, the viz layer behind the new `viz` optional extra, `matplotlib>=3.7`: `viz.colors` holds the Tol bright, vibrant, and muted schemes, the Dutch Field palette, and the maintainer's `STUDY_COLORS` as `#RRGGBB` tuples, importable without a plotting library, and `truncate_colormap`, which cuts a color map to a sub-range; `viz.helpers` holds `set_xticks` and `set_yticks`, which place an axis's ticks every `step` between two limits inside a rounded-down-tenth-of-a-step frame, the proxy-handle makers `color_handles` and `line_handles`, `move_legend_to_axes`, which moves a legend into chosen panels, and `band_under_lines`, which shades a band under each line. These port the generic helpers and palettes of the research code's `visualize` package, with its removed `Legend.legendHandles` use fixed; its study tables, figure scripts, and unused or uncalled helpers are not ported. CI installs the extra in every test job, and `CONTRIBUTING.md` gains the area `viz` ([#80](https://github.com/amirhs1/PolyPhys/issues/80), [#169](https://github.com/amirhs1/PolyPhys/pull/169), `a175e39`)
- Add `benchmarks/organizer_timing.py`, a script outside the package that times `whole_from_segments`, `ensemble`, and `ensemble_avg` on one property of one space and prints timings, peak memory, and counts only, for step 1 of [#69](https://github.com/amirhs1/PolyPhys/issues/69) ([#171](https://github.com/amirhs1/PolyPhys/pull/171), `c8287f7`)

### Fixed

- Correct the parser count, the module lists of the README and `AGENTS.md`, the contents of a study folder, the CI and coverage descriptions, and four spellings, now American, in the prose files ([#166](https://github.com/amirhs1/PolyPhys/issues/166), [#167](https://github.com/amirhs1/PolyPhys/pull/167), `6a9a186`)
- Wrap each position into the box along the periodic directions in `SpatialHistograms` before binning, since a LAMMPS dump can place an atom slightly outside the box between neighbor-list rebuilds; real HnsCub frames raised `ValueError` before ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#173](https://github.com/amirhs1/PolyPhys/pull/173), `89f8d72`)
- Update the README's current focus and tier-review date, the release tag of the install example in `usage.md`, the stage-directory wording for the probe's outputs in `lineage.md`, and a spelling in `CITATION.cff` ([#166](https://github.com/amirhs1/PolyPhys/issues/166), [#167](https://github.com/amirhs1/PolyPhys/pull/167), `3e7fe0b`)
- Compare the density golden test's axial profiles with the research code's times `(dcyl / lcyl)**2`, since the research code's axial slab volume takes `lcyl` where its integrand expects `dcyl`, and read one folder per direction under `POLYPHYS_DATA_DIR/density/` ([#144](https://github.com/amirhs1/PolyPhys/issues/144), [#175](https://github.com/amirhs1/PolyPhys/pull/175), `5b6b780`)

## [0.7.0](https://github.com/amirhs1/PolyPhys/releases/tag/v0.7.0) - 2026-11-05

### Changed

- **Breaking:** Rename the column that `space_tseries`, `space_hists`, and `space_profiles` round the crowder volume fraction into from `phi_c_bulk_round` to `phi_bulk_c_round`, so that every bulk name in PolyPhys reads `<quantity>_bulk_<species>` ([#105](https://github.com/amirhs1/PolyPhys/issues/105), [#150](https://github.com/amirhs1/PolyPhys/pull/150), `fdea860`). The research-code names map to the PolyPhys names as `phi_m_bulk` → `phi_bulk_m`, `rho_m_bulk` → `rho_bulk_m`, `phi_c_bulk` → `phi_bulk_c`, `rho_c_bulk` → `rho_bulk_c`, `phi_c_bulk_one` → `phi_bulk_c_one`, `phi_c_bulk_two` → `phi_bulk_c_two`, `rho_c_bulk_one` → `rho_bulk_c_one`, `rho_c_bulk_two` → `rho_bulk_c_two`, `phi_hns_bulk` → `phi_bulk_hns`, `rho_hns_bulk` → `rho_bulk_hns`, and `phi_c_bulk_round` → `phi_bulk_c_round`.
- Add area `probe` in `CONTRIBUTING.md`, "Names", for `polyphys/probe/` and `polyphys/tests/probe/`, and give `polyphys/probe/` the Instrumented tier in `AGENTS.md` and the README ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#153](https://github.com/amirhs1/PolyPhys/pull/153), `82401b6`)
- Measure the CI matrix jobs' coverage floor without the probe layer and its tests, through `.coveragerc-no-mda`, since those jobs run without MDAnalysis and cannot execute them; the MDAnalysis job keeps measuring the probe layer, and Codecov merges the reports ([#156](https://github.com/amirhs1/PolyPhys/issues/156), [#157](https://github.com/amirhs1/PolyPhys/pull/157), `6deb071`)
- Warn with a `UserWarning`, naming the file, in `polyphys.analyze.sumrule.space_profiles` when a profile's `<prop>-mean` sums to one within `1e-6`, the signature of a profile normalized by its sum, and state in `sumrule`, `usage.md`, and `density.md` that the sum rule takes the raw ensemble-averaged profiles, since a species normalized on its own enters the species sum with an arbitrary weight ([#140](https://github.com/amirhs1/PolyPhys/issues/140), [#151](https://github.com/amirhs1/PolyPhys/pull/151), `90d163f`)
- Restructure `test_measurer.py` and `test_utils.py` into one behaviour per test with stated tolerances, add the cases that bring `polyphys.analyze.measurer` to full line coverage, and correct the walled-cube volume-fraction expectation to its closed form, `1000 / 9**3 * pi / 6` ([#82](https://github.com/amirhs1/PolyPhys/issues/82), [#161](https://github.com/amirhs1/PolyPhys/pull/161), `6dfefd2`)

### Added

- Add the `mda` optional extra, `MDAnalysis>=2.8`, and a CI job on Python 3.14 that installs `.[dev,mda]` and runs the test suite with the probe tests, which the matrix jobs skip; the matrix jobs are unchanged ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#153](https://github.com/amirhs1/PolyPhys/pull/153), `42bf61d`, `a1bccf5`)
- Add `ProbeSpec` and `CHAIN_SIZE_OBSERVABLES` to `polyphys.manage.projects`: the record, built beside a study's `ProjectSpec`, of the MDAnalysis selection of each species, the chain's label, the chain-size observables to measure, the `unwrap` convention of the asphericity, and the histogram bin sizes, checked against the project's time series at construction ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#153](https://github.com/amirhs1/PolyPhys/pull/153), `f04285b`)
- Add `polyphys.probe`, the probe layer on MDAnalysis: `make_universe` opens a LAMMPS data file and dump in the studies' convention, `write_stamps` writes the stamps of a probed artifact, and `ChainSize`, an `AnalysisBase` subclass driven by a `ProbeSpec`, measures a chain's end-to-end distance, span, transverse size, radius of gyration, asphericity, shape parameter, and principal axes per frame, in serial or across workers, and saves each series under the research code's name ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#153](https://github.com/amirhs1/PolyPhys/pull/153), `636beac`, `6a58280`). For scripts that migrate from the research code: no per-frame shift of the positions is applied, since every chain-size observable is translation invariant; `asphericity` takes `unwrap` from the spec; and the frame times are not converted to simulation time.
- Add the `PROBE` record of each study in `projects/`, with the research code's selections, chain-size observables, and unwrap convention, and add `transSizeTMon` to the `TransFociCyl` spec's time series, which the research registry omitted although its probe wrote the series and its analysis table reads it, so the pipeline reduces one more time series for that study than the research code did ([#56](https://github.com/amirhs1/PolyPhys/issues/56), [#153](https://github.com/amirhs1/PolyPhys/pull/153), `a8eeadc`; [#154](https://github.com/amirhs1/PolyPhys/pull/154), `118095b`)
- Add `BinderSpec` and `BINDER_SERIES` to `polyphys.manage.projects`, and the optional `binders` field of `ProbeSpec`: the record of how a study's two-patch binders are counted, its patch and binder labels, binding cutoff, inclusive or strict comparison, cis threshold, patches per binder, and whether the distance matrices are kept, checked against the selections and the project's time series at construction ([#156](https://github.com/amirhs1/PolyPhys/issues/156), [#157](https://github.com/amirhs1/PolyPhys/pull/157), `5196e8f`)
- Add `BinderStatistics` to `polyphys.probe`, an `AnalysisBase` subclass driven by a `BinderSpec` that counts per frame the patch-monomer contacts, the free and engaged patches, the free, dangling, and bridging binders, and the cis and trans bridged pairs with `binding_statistics`, writes the loop-length histogram as per-frame int64 rows, `loopLengthHistT<chain>`, and optionally the float64 monomer-to-patch distance matrices, `distMatT<chain><patches>`, in serial or across workers, under the research code's file names ([#156](https://github.com/amirhs1/PolyPhys/issues/156), [#157](https://github.com/amirhs1/PolyPhys/pull/157), `08d83a1`). For scripts that migrate from the research code: the distances take the minimum image along the geometry's periodic directions, which the research probe did not; the comparison with the cutoff is inclusive by default where the research probe's was strict; and the loop-length histogram is saved per frame where the research probe saved its sum over the frames.
- Add the `binders` record to the `HnsCub` and `HnsCyl` probe specs in `projects/`: the `HnsPatch` and `HnsCore` labels, the research cis threshold of 4, the distance matrices kept, and the published binding criterion of the thesis, a distance of at most 0.685 monomer diameters, in place of the research probe's 0.661 with a strict comparison; the probe guide, "Choosing the cutoff", derives both from the model's Lennard-Jones parameters and records the research threshold sweep ([#156](https://github.com/amirhs1/PolyPhys/issues/156), [#157](https://github.com/amirhs1/PolyPhys/pull/157), `e8ce7e2`, `c217a2f`)
- Add to `ProbeSpec` the fields that drive the spatial histograms and the pair distances of a study: `histogram_species`, `block_frames`, `recentre`, `pair_species`, `contact_cutoff`, and `minimum_image`, each with a default so that an existing record is unchanged; let a bin size or the contact cutoff be a number or a function of the whole's parser (`PerWhole`), resolved for one whole by `ProbeSpec.resolve`; and add `RecentreT`, the three ways a probe moves a frame's positions before binning them ([#155](https://github.com/amirhs1/PolyPhys/issues/155), [#159](https://github.com/amirhs1/PolyPhys/pull/159), `3b5b269`)
- Add `polyphys.probe.histograms`: `SpatialHistograms`, an `AnalysisBase` subclass driven by a `ProbeSpec`, bins every species of a whole along the directions of its geometry (`r`, `z`, `theta` in a cylinder, `r` in a cube) with int64 counts kept per frame, and in the `xy`, `xz`, and `yz` planes with the mean map of each block of `block_frames` frames, after recentring the frame as the spec says; it raises on a position outside the edges, where the research code dropped it silently, so every frame conserves the particle number; and `histogram_edges` builds a whole's edges from the spec's bin sizes over ranges that reach the cylinder's wall surface, `(dcyl + 1) / 2`, and the cube's corner, `sqrt(3) lcube / 2`, which the maintainer extended from the research code's `dcyl / 2` and `lcube / 2` on 2026-11-04. The saved files are the research code's `<dir>Edge<Species>`, `<dir>Hist<Species>`, `<x|y|z>Edge`, and `<plane>Hist<Species>`, plus `<dir>HistT<Species>`, `<plane>HistBlock<Species>`, and `blockFrames` ([#155](https://github.com/amirhs1/PolyPhys/issues/155), [#159](https://github.com/amirhs1/PolyPhys/pull/159), `c6fb60d`, `aca6368`). For scripts that migrate from the research code: the research shift by half the parser's `(dcyl, dcyl, lcyl)` put a cylinder's axis at `(0.5, 0.5)` rather than the origin, since the LAMMPS box spans `dcyl + 1`; `recentre='half'` restores the LAMMPS frame exactly, and a triple reproduces the research shift.
- Add `polyphys.probe.pairs`: `PairDistances`, an `AnalysisBase` subclass driven by a `ProbeSpec`, measures the pair distances, direct contacts, bonds, and clusters of one species per frame, with the minimum image along the periodic directions when the spec asks, and saves them as the research code's `distMatT<Species>`, `directContactsMatT<Species>` (bool here, int64 there), `bondsHistT<Species>`, and `clustersHistT<Species>` ([#155](https://github.com/amirhs1/PolyPhys/issues/155), [#159](https://github.com/amirhs1/PolyPhys/pull/159), `88550d2`). For scripts that migrate from the research code: the diagonal of the distance matrix carries each particle's 0-based backbone index, which `genomic_distance` requires, where the research code wrote the 1-based LAMMPS atom ID, so a pipeline pair tag `loc<i>loc<j>` is one less than the research tag of the same pair.
- Extend the `PROBE` record of each study in `projects/` with the all-group selections of the research probe (`'Crd': 'resid 0'`, `'Dna': 'type 1'` for the TransFoci studies, `'Hns': 'type 3'` for the Hns studies), its bin sizes as functions of the whole's diameters, the histogram species the analysis tables read, `block_frames=100`, `recentre='half'`, and, for the TransFoci studies, the foci as the pair species with the research cutoff `dmon_large + dcrowd`; the HnsCyl radial bin size now includes the H-NS core diameter in its minimum as HnsCub's did, which the maintainer confirmed as the research code's mistake on 2026-11-04 ([#155](https://github.com/amirhs1/PolyPhys/issues/155), [#159](https://github.com/amirhs1/PolyPhys/pull/159), `1acc7d5`)
- Add `threshold_sweep` to `polyphys.analyze.clusters`, the research code's choice-of-threshold sweep: it recounts the monomer-to-patch distance matrices that `BinderStatistics` stores at every cutoff and cis threshold of a grid and returns one row per pair with the frame means of the eight binding counts, plus the frame-averaged loop-length histogram, with the comparison (`inclusive`) a required keyword ([#158](https://github.com/amirhs1/PolyPhys/issues/158), [#160](https://github.com/amirhs1/PolyPhys/pull/160), `9347e01`). For scripts that migrate from the research code: the research sweep compared strictly where the H-NS specs count inclusively, so its row at 0.685 does not match the stored series; its `arange` grid ran one step past 0.685 to 0.695; and its two CSVs repeated the rows of earlier cutoffs, which the sweep does not.
- Add `BULK_NAME_MAP` and `rename_bulk_columns` to `polyphys.manage.organizer`, which rename the bulk-quantity columns of a stamp, ensemble-average, or `space` table the research code wrote to the PolyPhys names and leave every other column unchanged ([#105](https://github.com/amirhs1/PolyPhys/issues/105), [#150](https://github.com/amirhs1/PolyPhys/pull/150), `fdea860`)
- Write `<prop>-norm-sem`, the standard error over the wholes of each whole's density profile divided by its own value in the bin that the sum rule's default normalizer names, beside the ensAvg columns in `polyphys.analyze.pipeline`, through a `normalizer` field of `Derived` that `LocalDensities` sets; carry it through `polyphys.analyze.sumrule.space_profiles`, folded along `z` and dropped with a warning under an explicit normalizer; and add the species' columns in quadrature as `<prop>Sum-norm-sem` in `sum_rule_profiles` ([#143](https://github.com/amirhs1/PolyPhys/issues/143), [#151](https://github.com/amirhs1/PolyPhys/pull/151), `fb5bc51`, `f195b79`)
- Add Quickstart, Lineage, Workflow, and Studies documentation pages, move the study sections of the Usage page to Studies, and run the examples of the Quickstart, Lineage, Studies, and Usage pages as tests ([#79](https://github.com/amirhs1/PolyPhys/issues/79), [#162](https://github.com/amirhs1/PolyPhys/pull/162), `4118a0c`)
- Add the "Simulation design" guide page, `docs/source/design.md`, which lists for each study with a parser the design values its LAMMPS inputs set and its names do not carry (interaction strengths and cutoffs, wall strengths, bonds, angles, masses, thermostat, and time steps), and point each parser's docstring to it ([#104](https://github.com/amirhs1/PolyPhys/issues/104), [#163](https://github.com/amirhs1/PolyPhys/pull/163), `b11bbd2`)

### Fixed

- Describe `SumRuleCyl.epsilon` as the wall-monomer LJ strength, not the wall-particle one, since the wall-crowder strength is fixed at 1, and describe `HnsCub.eps_hc` and `HnsCyl.eps_hc` as a name keyword that the project's LAMMPS inputs do not read, not as an attractive core-crowder strength ([#104](https://github.com/amirhs1/PolyPhys/issues/104), [#163](https://github.com/amirhs1/PolyPhys/pull/163), `b11bbd2`)
- Add Python 3.14 to the README's Python badge, which listed 3.11 to 3.13 while the package has supported 3.14 since it joined the CI matrix ([#149](https://github.com/amirhs1/PolyPhys/issues/149), [#152](https://github.com/amirhs1/PolyPhys/pull/152), `7f08516`, `a76728b`)
- Serve the README's Docs badge from shields.io, since GitHub's image proxy returned an error for the Read the Docs badge URL and the README showed a broken image ([#164](https://github.com/amirhs1/PolyPhys/pull/164), `2b20e19`, `f89c7bf`)

## [0.6.0](https://github.com/amirhs1/PolyPhys/releases/tag/v0.6.0) - 2026-11-03

### Changed

- **Breaking:** Raise `ValueError` in the `polyphys.manage.parser` parsers when a name lacks a keyword of its lineage's pattern, or a keyword's value is not a number, naming the keyword; before, a parser printed a message and failed later with an `AttributeError` ([#99](https://github.com/amirhs1/PolyPhys/pull/99), `a18845a`)
- **Breaking:** Raise `ValueError` for a non-integral value of an integer keyword, such as `nc1.5`, instead of truncating it ([#99](https://github.com/amirhs1/PolyPhys/pull/99), `a18845a`)
- **Breaking:** Set a parser's `filepath` and `ext` to `None` when it is built from a name, instead of `'N/A'`, and make `filepath` a `pathlib.Path` when it is built from a file path ([#101](https://github.com/amirhs1/PolyPhys/pull/101), `5786ca0`)
- **Breaking:** Print a whole-number float keyword as an integer, and any other float as its `repr`, in the `ensemble` and `space` names a parser rebuilds ([#108](https://github.com/amirhs1/PolyPhys/pull/108), `f7f44dc`). `ac2.0` becomes `ac2` and `pctc40.0` becomes `pctc40`; a short key in existing stamps maps to the new name by parsing it at its lineage and reading the rebuilt name.
- **Breaking:** Require `_dump_keywords`, one entry per group, in every concrete parser, and raise `TypeError` when a class lacks it or its entries do not match `_groups` ([#113](https://github.com/amirhs1/PolyPhys/pull/113), `4eb9973`). To migrate a parser declared outside PolyPhys, add it, for example `{'bug': 'bdump', 'all': 'adump'}`.
- Make the eight parser classes declarations over `ParserBase`, which holds the one copy of the parsing and bulk-value code, and make `ParserBase` a plain class instead of an `abc.ABC`; attribute values and types are unchanged for every name the previous parsers accepted ([#99](https://github.com/amirhs1/PolyPhys/pull/99), `a18845a`)
- Correct the parser docstrings to match the code ([#99](https://github.com/amirhs1/PolyPhys/pull/99), `a18845a`)
- Restructure `test_parser.py` into fixtures and one behaviour per test, with a round-trip test over every class's declared patterns ([#99](https://github.com/amirhs1/PolyPhys/pull/99), `a18845a`)
- Log instead of printing in `polyphys.manage`, and let an error raised inside `openany_context` propagate unchanged ([#107](https://github.com/amirhs1/PolyPhys/pull/107), `496ca9b`)
- Document that rebuilt `segment` and `whole` names carry no topology suffix, as since v0.3.0 ([#108](https://github.com/amirhs1/PolyPhys/pull/108), `f7f44dc`). SciOps stamps of an artifact parsed from a whole file key the whole as `…ens1.ring`; append `.<topology>` to the new name to match them.
- Define each commit type in `CONTRIBUTING.md`, "Names", drop type `ci`, and rename area `docs` to `guide` ([#123](https://github.com/amirhs1/PolyPhys/pull/123), `03f8478`)
- Add area `simulate` in `CONTRIBUTING.md`, "Names", for `polyphys/simulate/` and `polyphys/tests/simulate/`, and give `polyphys/simulate/` the Instrumented tier in `AGENTS.md` and the README ([#125](https://github.com/amirhs1/PolyPhys/pull/125), `05384bb`)
- Type-check `polyphys/simulate` in CI, in `pyproject.toml`, and in the pull request template ([#126](https://github.com/amirhs1/PolyPhys/pull/126), `e3a1689`)
- Add area `projects` in `CONTRIBUTING.md`, "Names", for `projects/`, and list `projects/` in the tier table of `AGENTS.md` and the README ([#76](https://github.com/amirhs1/PolyPhys/issues/76), [#132](https://github.com/amirhs1/PolyPhys/pull/132), `45be4dc`)
- Lint, type-check, and test `projects/` in CI, in `pyproject.toml`, and in the pull request template ([#76](https://github.com/amirhs1/PolyPhys/issues/76), [#132](https://github.com/amirhs1/PolyPhys/pull/132), `45be4dc`)
- State in `CONTRIBUTING.md`, "Changelog", which pull requests add a changelog entry and when, with a matching checkbox in the pull request template and a check in the `open-pull-request` skill ([#130](https://github.com/amirhs1/PolyPhys/issues/130), [#141](https://github.com/amirhs1/PolyPhys/pull/141), `51355fa`)
- Add a Releases section to `CONTRIBUTING.md`: a milestone is released from `main` as soon as its last issue closes and never from an older commit, and the changelog's format is Common Changelog, with the `[Unreleased]` section and the release-time commit references as the one departure ([#124](https://github.com/amirhs1/PolyPhys/pull/124), `b1d317d`, `9143524`)
- Ask for `git config fetch.prune true` once per clone, in the `CONTRIBUTING.md` checklist, so a fetch drops the remote-tracking refs of deleted branches ([#136](https://github.com/amirhs1/PolyPhys/pull/136), `811d9cd`)
- Bump `actions/setup-python` from 7.0.0 to 7.1.0 in the CI workflow ([#127](https://github.com/amirhs1/PolyPhys/pull/127), `284b2b8`)
- **Breaking:** Declare `PhaseT` in `polyphys.manage.types` as the phase vocabulary of `AGENTS.md`, `simsAll`, `simsCont`, `logs`, `trjs`, `probes`, `analysis`, `viz`, `galaxy`, and `allInOne`, which `make_database` already accepted, instead of the singular names it never accepted; `make_database` checks `phase` against `PhaseT` and, new, `stage` against `StageT`, and takes neither as a plain `str` ([#34](https://github.com/amirhs1/PolyPhys/issues/34), [#138](https://github.com/amirhs1/PolyPhys/pull/138), `1c6df51`)

### Added

- Add `polyphys.analyze.sumrule`, whose `space_profiles`, `normalize_r`, and `normalize_z` are the space reducers of the density sum rule ported from the research code: a species' ensemble-averaged `Phi` or `Rho` profiles scaled by its diameter, normalized by a reference value, folded about `z = 0` along `z`, and recentred on the wall for `r` in a cylinder, one function for one or several crowder species, with the species-to-size map of the study's `ProjectSpec` ([#54](https://github.com/amirhs1/PolyPhys/issues/54), [#138](https://github.com/amirhs1/PolyPhys/pull/138), `0a35d00`). For scripts that migrate from the research code: `space_sum_rule` and `space_sum_rule_bi_disperse` are one function, whose `normalizer` argument takes `'first'`, `'last'`, `'max'`, a number, or a function and defaults to the research rule in `PROFILE_NORMALIZERS`; the bulk volume fraction column is `phi_bulk_c`; the recentred columns are added for `r` only; and the saved file's name carries the direction, so `r` and `z` no longer collide.
- Add `sum_rule_profiles` to `polyphys.analyze.sumrule`, the sum of several species' scaled profiles state point by state point, normalized by its own reference value from `SUM_NORMALIZERS` or the `normalizer` argument, which the research code's notebooks built by hand ([#54](https://github.com/amirhs1/PolyPhys/issues/54), [#138](https://github.com/amirhs1/PolyPhys/pull/138), `0a35d00`)
- Add `NormalizerT` to `polyphys.manage.types` ([#54](https://github.com/amirhs1/PolyPhys/issues/54), [#138](https://github.com/amirhs1/PolyPhys/pull/138), `0a35d00`)
- Add `polyphys.analyze.density`, whose `bin_volumes`, `share_matrix`, `number_density`, `volume_fraction`, and `density_profiles` turn a histogram of bead centres into the local number density and volume fraction along `r` in a cube or a cylinder and along the periodic `z` of a cylinder, ported from the research code's `distributions.py` as numpy and scipy kernels with no parser dependency ([#52](https://github.com/amirhs1/PolyPhys/issues/52), [#142](https://github.com/amirhs1/PolyPhys/pull/142), `b9836bb`). The bin volumes are closed forms instead of quadratures; the volume shares are a sparse band matrix built once per bins and bead size; by default the test bead is averaged over its bin (`shares='average'`), which conserves the total bead volume to round-off where the research code's bin-centre placement (`shares='centre'`) was off by about 2 % of a bead volume at a bin width of a quarter diameter, and the radial shells of a cylinder use the exact ball-in-cylinder volume, `sphere_cylinder_intersection` (`intersection='exact'`), instead of the sphere-sphere formula (`intersection='sphere'`), whose error is about `0.09 a / R` of a bead volume. For scripts that migrate from the research code: `rho` is divided by `n_frames`; the slabs along `z` take the effective diameter `D - a` of the parsers' bulk values, so a crowder's far-field value is `(D / (D - a))**2` times the old one before normalization; `normalized=True` is `normalize='sum'`; and `n_particles` makes a lost bead an error.
- Add `polyphys.analyze.pipeline`, whose `reduce_space` reduces one space's probe artifacts from segments or wholes to ensembles and ensemble averages for every property of a study's `AnalysisTable`, and `reduce_stamps` combines its stamps, ported from the research code's `analyzer` ([#78](https://github.com/amirhs1/PolyPhys/issues/78), [#145](https://github.com/amirhs1/PolyPhys/pull/145), `fdef05b`, `617d579`). A property is a `PropertySpec` built by `tseries`, `histogram`, `histogram_2d`, `edge`, or `matrix_tseries`; a property with no file raises instead of being skipped; the wholes of an ensemble must share their bin edges; and the ensemble of DataFrame wholes keeps every whole, with the mean, variance, and standard error of each column in the ensemble average.
- Add the post-hooks `TimeAverage`, `Autocorrelation`, `LocalDensities`, and `FociPairDistances` to `polyphys.analyze.pipeline`, with `half_contour` as the separation bound of a ring ([#78](https://github.com/amirhs1/PolyPhys/issues/78), [#145](https://github.com/amirhs1/PolyPhys/pull/145), `fdef05b`, `617d579`). For scripts that migrate from the research code: the time average of a per-frame histogram is named by the hook, not by dropping `T`; the autocorrelation is stored alone, without the `acfLowerCi` and `acfUpperCi` bands; the density profiles come from a `DensityKernel` the caller injects; and the foci pair histograms share one grid up to a bound from the system, so `pairDistRdf` is the separation density of `clusters.pair_separation_density` and the `genDist` tags count bonds.
- Add `bartlett_half_width` to `polyphys.analyze.correlations`, the confidence band of a stored autocorrelation at any level, which `acf` now calls ([#78](https://github.com/amirhs1/PolyPhys/issues/78), [#145](https://github.com/amirhs1/PolyPhys/pull/145), `fdef05b`, `617d579`)
- Add an `ANALYSIS` tuple of `AnalysisTable` per particle group to the `TransFociCub`, `TransFociCyl`, `HnsCub`, and `HnsCyl` specs in `projects/`, with the values of the research code's analysis tables and each artifact listed once ([#78](https://github.com/amirhs1/PolyPhys/issues/78), [#145](https://github.com/amirhs1/PolyPhys/pull/145), `fdef05b`, `617d579`)
- Add `equilibrium_wholes` to `polyphys.analyze.pipeline`, which reads each ensemble file of the `ens` stage once and reduces every whole's time series to its time average, its sample variance, and the blocked standard error of `polyphys.analyze.stats.sem_blocked`, merged with the whole stamps, ported from the research code's `equilibrium_wholes` ([#78](https://github.com/amirhs1/PolyPhys/issues/78), [#146](https://github.com/amirhs1/PolyPhys/pull/146), `eb04b81`). For scripts that migrate from the research code: `-sem` is the blocked error with the `'lee'` plateau rule by default, larger than the per-frame `measurer.sem` by about the square root of the statistical inefficiency; `-var` has `ddof = 1`, where `np.var` on a pandas column gave the population variance; the output name is given per property, not made by dropping the last `T`; and a whole in the stamps or in the files only raises instead of being dropped.
- Add `polyphys.simulate.build`, whose `helix`, `plectoneme_ring`, and `resample_curve` place beads an exact bond length apart along a helix, along the closed plectoneme of Neukirch and Starostin, or along any curve given as points, and whose `chain_bonds`, `chain_angles`, and `chain_data` give a linear or ring chain its topology as a `DataFile` for `write_data` ([#114](https://github.com/amirhs1/PolyPhys/issues/114), [#131](https://github.com/amirhs1/PolyPhys/pull/131), `94b2e8e`, `831b9ae`). The default bond length is 0.97, the Kremer-Grest start value; the research builders' 0.95 R0 = 1.425 is above the length at which LAMMPS warns.
- Add `check_conformation`, `ConformationReport`, and `ConformationError` to `polyphys.simulate.build`, which test a start configuration against the conditions under which LAMMPS warns or stops: the FENE warning length, the ghost cutoff, coincident atoms, the non-bonded distance, the box on non-periodic axes, and a cylinder's effective radius ([#114](https://github.com/amirhs1/PolyPhys/issues/114), [#131](https://github.com/amirhs1/PolyPhys/pull/131), `94b2e8e`, `831b9ae`)
- Add the `lammps` pytest marker, the `lammps_executable` fixture, and `polyphys/tests/data/lammps_build/`, whose `in.kg` runs a built start under a chosen integrator; the marked tests run when `LMP` or `lmp_serial` names an executable and are skipped otherwise, as in CI ([#114](https://github.com/amirhs1/PolyPhys/issues/114), [#131](https://github.com/amirhs1/PolyPhys/pull/131), `94b2e8e`, `831b9ae`)
- Add `polyphys.analyze.clusters`, the contact and cluster statistics of a particle group, ported from the research code: `pair_distances`, `find_direct_contacts`, `cluster_labels`, `count_clusters`, `contact_matrix`, and `count_bonds`, which take the clusters as the connected components of the contact graph ([#51](https://github.com/amirhs1/PolyPhys/issues/51), [#128](https://github.com/amirhs1/PolyPhys/pull/128), `a0a9a56`). For scripts that migrate from the research code: `count_clusters` is indexed by size and its index 0 is always 0, where `count_foci_clusters` put the number of zero eigenvalues. The research `generate_contact_matrix` could stop short of the transitive closure, which reached the foci cluster-size histograms but not the H-NS cis, trans, and loop-length statistics; no published number came through it, as the maintainer confirmed on 2026-11-02 ([#129](https://github.com/amirhs1/PolyPhys/pull/129), `491081e`).
- Add `genomic_distance` and `foci_pair_info` to `polyphys.analyze.clusters`, in bonds on linear and ring chains ([#128](https://github.com/amirhs1/PolyPhys/pull/128), `a0a9a56`). `foci_pair_info` takes the topology and returns bonds, one more than the research `foci_info`'s monomers between.
- Add `pair_distance_histograms` and `pair_separation_density` to `polyphys.analyze.clusters`: the separations of every pair binned on one shared grid, and the probability density of a pair's separation vector, `H_k / (N_f dV_k)` with the exact shell volume ([#128](https://github.com/amirhs1/PolyPhys/pull/128), `a0a9a56`). It replaces the research `foci_rdf`, whose midpoint shell and unit-sum rescaling it drops: every stored `pairDistRdfFoci` value changes by one constant per curve plus the first-bin shell correction, and free-energy differences within a curve do not.
- Add `binding_statistics`, `FrameBinding`, and `BindingStats` to `polyphys.analyze.clusters`, the free, dangling, bridging, cis, and trans statistics of two-patch binders and their loop-length histogram, with `enforce_single_patch_contact`, `monomer_binder_contacts`, `split_binders`, and `binder_contacts` ([#128](https://github.com/amirhs1/PolyPhys/pull/128), `a0a9a56`). The counting rules are those of the research `hns_binding`; the Clusters page maps its keys to the `FrameBinding` fields.
- Add a Clusters documentation page ([#128](https://github.com/amirhs1/PolyPhys/pull/128), `a0a9a56`)
- Add a Correlations documentation page and the `analyze.correlations` row of the README's module table ([#53](https://github.com/amirhs1/PolyPhys/issues/53), [#134](https://github.com/amirhs1/PolyPhys/pull/134), `2a7871f`)
- Add `polyphys.analyze.correlations`, ported from the research code with numpy and scipy only: `acf` estimates the autocorrelation of a time series by the FFT with Bartlett's confidence band as an `AcfResult`, replacing `statsmodels.tsa.stattools.acf` with the same `1/N` estimator and band convention ([#53](https://github.com/amirhs1/PolyPhys/issues/53), [#133](https://github.com/amirhs1/PolyPhys/pull/133), `6ad5846`)
- Add `fit_decay`, `initial_guess`, `DecayFit`, and the `DECAY_MODELS` table to `polyphys.analyze.correlations`: one `scipy.optimize.curve_fit` entry point over the pure, stretched, scaled, and offset exponentials, written in the decay time `tau` rather than the rate `omega` of the research `mono_*` functions, with a data-driven starting point recorded in the result ([#53](https://github.com/amirhs1/PolyPhys/issues/53), [#133](https://github.com/amirhs1/PolyPhys/pull/133), `6ad5846`). A failed fit returns `converged=False` and NaN parameters of the model's length, where the research `fit_wholes` and `fit_exp_wholes` padded four zeros whatever the model; and `stderr` is the square root of the covariance diagonal, where those drivers stored the variance under a `-std` name.
- Add `bond_vectors` and `bond_cosine_correlation` to `polyphys.analyze.correlations`, the bond vectors of a linear or ring chain and the mean cosine between bonds at every lag by the FFT, `O(N log N)` per frame against the research `bond_info`'s `O(N^2)`, which raised on every linear chain ([#53](https://github.com/amirhs1/PolyPhys/issues/53), [#133](https://github.com/amirhs1/PolyPhys/pull/133), `6ad5846`). On a ring every lag now averages all `N` pairs going round the ring, so `C(k) = C(N - k)` and the curve does not depend on which monomer is labeled first; `bond_info` summed the `N - k` open pairs, so stored `bondCosineCorrVecMon` values differ beyond lag 0 by the pairs that cross the label origin, with the same expectation for an equilibrated ring.
- Add `scipy>=1.14.1` as a runtime dependency, the first SciPy with wheels for every Python in the CI matrix, and `scipy-stubs` to the `dev` extra ([#128](https://github.com/amirhs1/PolyPhys/pull/128), `a0a9a56`)
- Add `polyphys.simulate`, a subpackage for the files of LAMMPS runs, whose `core` module holds `Box`, the simulation box that converts scaled coordinates and applies image flags, and `LammpsFormatError`, the base of `BrokenLogError`, `BrokenDumpError`, and `BrokenDataError` ([#126](https://github.com/amirhs1/PolyPhys/pull/126), `e3a1689`)
- Add `polyphys.simulate.log`, whose `read_log` and `parse_log` read a LAMMPS log, plain or gzipped, into a `LogFile` of its thermo output, per-run statistics, and wall time, and raise `BrokenLogError` for a cut-off or malformed log ([#111](https://github.com/amirhs1/PolyPhys/pull/111), `944801f`, `71eb401`; moved from `polyphys.manage.lammps_log` in [#126](https://github.com/amirhs1/PolyPhys/pull/126))
- Add `polyphys.simulate.dump`, whose `iter_dump` and `read_dump_frame` read a LAMMPS text dump, plain or gzipped, one frame at a time into `DumpFrame` records of numpy arrays, with Cartesian positions from scaled, wrapped, or unwrapped columns in orthogonal and restricted triclinic boxes, and raise `BrokenDumpError` for a cut-off or malformed dump; it needs numpy only ([#126](https://github.com/amirhs1/PolyPhys/pull/126), `7baf0d4`)
- Add `polyphys.simulate.data`, whose `read_data` and `write_data` read and write LAMMPS data files in the `atomic`, `bond`, `angle`, `molecular`, and `full` atom styles as `DataFile` records, and whose `with_frame` gives a dump frame its topology from a template, matching atoms by ID ([#126](https://github.com/amirhs1/PolyPhys/pull/126), `7baf0d4`)
- Add `KeywordSpec` and `Species` to `polyphys.manage.parser`, with which a parser class declares its filename keywords, lineage patterns, fixed values, and bulk species ([#99](https://github.com/amirhs1/PolyPhys/pull/99), `a18845a`)
- Add the `SumRuleCubHeteroRingBiDisperse` parser, for a ring of small and large monomers among crowders of two sizes in a periodic cube ([#106](https://github.com/amirhs1/PolyPhys/pull/106), `bd1c257`). For scripts that migrate from the research code: the bulk values are named `rho_bulk_<species>` and `phi_bulk_<species>` (`phi_c_bulk_one` → `phi_bulk_c_one`, `rho_c_bulk` → `rho_bulk_c`, `phi_m_bulk` → `phi_bulk_m`, `rho_m_bulk` → `rho_bulk_m`); the constructor is `(artifact, lineage, group, ispath=False)`; `nmon`, `ncrowd`, and `eps_others` are gone; and `mmon_large`, `mcrowd_one`, and `mcrowd_two` exist only at `segment`, `whole`, and `ensemble_long`.
- Add `ParserBase._derived_attributes`, which a parser class overrides to set values derived from parsed ones ([#106](https://github.com/amirhs1/PolyPhys/pull/106), `bd1c257`)
- Add `ParserBase.dump_keywords` and `dump_keyword`, which name the attribute that holds each group's dump period ([#113](https://github.com/amirhs1/PolyPhys/pull/113), `4eb9973`)
- Add `polyphys.manage.projects`, whose `ProjectSpec` records what a study's analysis needs beyond its parser: the polymer group, the rounding step `divisor` of the bulk crowder volume fraction, the attributes added as columns, the species-to-size map of the density profiles, and the time series; it reads the geometry, topology, profile directions, and the `space` and `ensemble_long` glob patterns from the parser, and raises `ValueError` at construction for an attribute the parser does not set at `ensemble_long` ([#76](https://github.com/amirhs1/PolyPhys/issues/76), [#132](https://github.com/amirhs1/PolyPhys/pull/132), `45be4dc`)
- Add a top-level `projects/` folder, outside the installed package, with one folder per study holding its spec: `TransFociCyl`, `TransFociCub`, `HnsCub`, and the provisional `HnsCyl`; CI lints, type-checks, and tests it ([#76](https://github.com/amirhs1/PolyPhys/issues/76), [#132](https://github.com/amirhs1/PolyPhys/pull/132), `45be4dc`)
- Add `ParserBase.pattern` and `ParserBase.attributes_at` to `polyphys.manage.parser`, class methods that give a lineage's filename keywords and the names of the attributes a parser of that lineage sets, without parsing a name ([#76](https://github.com/amirhs1/PolyPhys/issues/76), [#132](https://github.com/amirhs1/PolyPhys/pull/132), `45be4dc`)
- Add a Design decisions page to the documentation, recording that `manage` builds on no workflow framework such as signac ([#19](https://github.com/amirhs1/PolyPhys/issues/19), [#132](https://github.com/amirhs1/PolyPhys/pull/132), `45be4dc`), and that each study's code lives in a folder outside the package ([#76](https://github.com/amirhs1/PolyPhys/issues/76))
- Add Hypothesis to the `dev` extra, with a property test that a histogram counts every position in its range exactly once ([#109](https://github.com/amirhs1/PolyPhys/pull/109), `4a6a968`)
- Add a `private_data` test that compares the rebuilt `HnsCub` names with the SciOps stamps under `POLYPHYS_DATA_DIR` ([#108](https://github.com/amirhs1/PolyPhys/pull/108), `f7f44dc`, `0e212ae`)
- Add Python 3.14 to the CI test matrix and the `pyproject.toml` classifiers, and run the lint, type-check, docs, and build jobs on it instead of 3.13 ([#119](https://github.com/amirhs1/PolyPhys/pull/119), `758c5ff`)

### Removed

- **Breaking:** Remove the `dump_key` argument of `space_tseries` in `polyphys.manage.organizer`, which now scales the frame index by the `group`'s own dump period, and make its arguments after `group` keyword-only ([#113](https://github.com/amirhs1/PolyPhys/pull/113), `4eb9973`). To migrate, drop the dump keyword and pass `divisor`, `round_to`, `is_save`, and `crowder_attr` by name. A call that paired a group with another group's dump keyword, such as `'bug'` with `'adump'`, gets new `t_index` and `time` values.

### Fixed

- Drop the child lineage column in `parents_stamps` when the stamps have no `lineage_name` column ([#107](https://github.com/amirhs1/PolyPhys/pull/107), `496ca9b`). Tables produced before this fix keep a stale `segment` or `whole` column.
- Drop whichever of `segment` and `segment_id` is present in `children_stamps(..., "whole")` ([#107](https://github.com/amirhs1/PolyPhys/pull/107), `496ca9b`)
- List `phi_bulk_m` in the `attributes` of `SumRuleCubHeteroRing`, `SumRuleCubHeteroLinear`, `TransFociCyl`, and `TransFociCub`, so their stamps gain a `phi_bulk_m` column; no value changes ([#110](https://github.com/amirhs1/PolyPhys/pull/110), `2d0b394`)
- Replace the eight `inplace=True` pandas calls in `polyphys.manage.organizer`, which pandas 3.1 deprecates with `Pandas4Warning`; results are unchanged ([#115](https://github.com/amirhs1/PolyPhys/issues/115), [#120](https://github.com/amirhs1/PolyPhys/pull/120), `3fea0e8`)

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
