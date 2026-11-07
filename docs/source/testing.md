# Testing

This page describes how PolyPhys is tested: which data a test may use, the
kinds of test the project writes, the shared fixtures, and the pytest settings.

## Running the tests

From the repository root, in an environment with the `dev` and `viz` extras
installed; `--doctest-modules` imports `polyphys/viz/`, which needs
matplotlib, the `viz` extra:

```bash
python -m pytest polyphys README.md projects --cov=polyphys --cov-report=term-missing --doctest-modules --doctest-glob='README.md'
```

The README's examples run as doctests, and so do the `Examples` sections of the
docstrings in `polyphys/`.

CI's Python-matrix jobs, which have no MDAnalysis, run the same command
with `--cov-config=.coveragerc-no-mda`, which leaves the probe layer out of
the measured set, and `--cov-fail-under=95`, so a pull request fails when
total coverage drops below 95 percent; the MDAnalysis job measures the probe
layer and has no floor. coverage.py rounds the total to a whole percent
before it compares, so 94.6 percent passes. Add
`--cov-fail-under=95` to the command above to run the same check locally.

A CI job on Python 3.11 tests the lowest versions that `pyproject.toml`
allows of NumPy, pandas, SciPy, and the `viz` extra's matplotlib, with
[uv](https://docs.astral.sh/uv/)'s `--resolution lowest-direct`, and runs the
suite without MDAnalysis or coverage. To run it locally, with uv installed:

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python --resolution lowest-direct -e ".[viz]"
uv pip install --python .venv/bin/python pytest hypothesis
.venv/bin/python -m pytest polyphys README.md projects --doctest-modules --doctest-glob='README.md'
```

## Test data: three tiers

Every test uses data from one of three tiers.

| Tier | Data                                                                                                                                                   | Where it lives                                   | Who runs it                    |
| ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------ | ------------------------------ |
| A    | Synthetic inputs a test generates itself: lineage trees of tiny files, short LAMMPS files, random or analytic arrays                                   | the test module or `polyphys/tests/conftest.py`  | everyone, including CI         |
| B    | Small trimmed real files, such as a few dump frames or a filename listing, used where a file format is under test and, for the probes, to check results on real frames (see "Real frames of the probes"); the maintainer approves each one | `polyphys/tests/data/`, each about 300 KB or less | everyone, including CI         |
| C    | Full-size processed artifacts, such as a whole space                                                                                                   | outside the repository, under `POLYPHYS_DATA_DIR` | the maintainer, locally; never CI |

Tier A is the default. The algorithms are what is under test, and properties
such as conservation, symmetry, invariance, or agreement with an analytic result
do not need real data. Tier C tests compare a result against the maintainer's
processed data, for example when a module replaces code that produced published
numbers.

## Kinds of test

| Kind                       | Catches                                               | Notes                                                                                                                   |
| -------------------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Unit test                  | a function breaking its own contract                  | tier A data; one behavior per test                                                                                     |
| Known-answer test          | code that runs but computes the wrong physics         | the expected value comes from an analytic result or a cited reference                                                   |
| Property test              | edge cases nobody wrote down                          | checks an invariant, such as conservation of histogram counts, on hand-picked or Hypothesis-generated inputs            |
| Integration test           | lineage levels that pass alone but disagree at a seam | runs a lineage tree from `segment` up through the organizer                                                             |
| Regression test            | a fixed bug returning                                 | written before the fix and seen to fail against the unfixed code                                                        |
| Comparison against old code | a port that changed a number                         | tier C; the old and new implementations run on the same data, and any number that changes is explained in the pull request |
| Doctest                    | examples drifting from behavior                      | every public name has an `Examples` section                                                                             |
| Static check               | type and style errors                                 | `flake8` and `mypy`                                                                                                     |

Expected values never come from the code under test. Floating-point results are
compared with a stated tolerance, using `numpy.testing.assert_allclose` for
arrays and `pytest.approx` for scalars. Uncertainties of quantities measured over
correlated molecular-dynamics frames use block averaging or another
autocorrelation-aware estimate, not a per-frame standard error.

## Property tests with Hypothesis

A property test can let [Hypothesis](https://hypothesis.readthedocs.io/)
generate its inputs: the `@given` decorator draws them from strategies, and
when an input breaks the property, Hypothesis shrinks it to a minimal failing
example and reports it. For example,
`test_axial_histogram_counts_every_position_in_range_once` in
`polyphys/tests/analyze/test_measurer.py` draws a bin size, a range, and
positions within that range, and checks that each position is counted exactly
once.

Hypothesis is in the `dev` extra. Its built-in `ci` settings profile is active
whenever the `CI` environment variable is set, as it is on GitHub Actions: CI
runs then draw the same inputs every time (`derandomize=True`) and keep no
example database. Local runs draw new inputs each time and store failing
examples under `.hypothesis/`, which git ignores, so a failure found once is
tried first on the next run. To run a test as CI does, set `CI=true` for that
run. See
[Built-in profiles](https://hypothesis.readthedocs.io/en/latest/reference/api.html#built-in-profiles)
in the Hypothesis documentation.

## Shared fixtures

`polyphys/tests/conftest.py` provides tier A fixtures to every test module.

`make_lineage_tree(parser, ...)`
: Writes a lineage tree of tiny `.npy` or `.csv` files, named for a parser
  class such as those of `polyphys.tests.synthetic`, under the test's
  `tmp_path`. The
  tree has `n_ensembles` state points that differ in one filename keyword,
  `n_wholes` replicates per state point, and, at the `segment` level,
  `n_segments` segments per replicate. Each file's values encode its indices
  (`1000 * ensemble + 100 * whole + 10 * segment + row`), so expected sums and
  concatenations can be written down by hand.

`make_lammps_data(atom_style=..., ...)`
: Writes a LAMMPS data file of one linear chain in the `atomic`, `bond`, or
  `angle` atom style, following the formats of the LAMMPS
  [read_data](https://docs.lammps.org/read_data.html) command.

`make_lammps_dump(...)`
: Writes a LAMMPS dump file in the `atom` style with unscaled coordinates
  (`ITEM: ATOMS id type x y z`) and chosen boundary flags, following the
  LAMMPS [dump](https://docs.lammps.org/dump.html) command.

`make_temp_file(filename)`
: Creates an empty file with a given name.

`lammps_executable`
: The LAMMPS executable, for a test marked `lammps` (see below).

The fixtures other than `lammps_executable`, which gives the executable as a
string, are factories, called inside a test with the arguments that test
needs; see [Factories as fixtures](https://docs.pytest.org/en/stable/how-to/fixtures.html#factories-as-fixtures)
in the pytest documentation. Their tests are in `polyphys/tests/test_conftest.py`.

## Real frames of the probes

The probe tests in `polyphys/tests/probe/` also run on real frames of the
four live studies, 22 wholes in all: per study, two replicates without
crowders and two at a bulk crowder volume fraction of about 0.36 to 0.38;
for the Hns studies, also the wholes without H-NS; and for HnsCyl, crowders
of size 1 as well as 2. The maintainer trimmed them from the studies' runs.

- Tier B, in `polyphys/tests/data/probe/<whole>/`: every whole's last 50
  frames of the polymer group dump,
  `<whole>.ring.<group>.last50.lammpstrj.gz`; and, for the wholes without
  crowders and the HnsCyl wholes at `ac2 nc386`, the `all` frames of the
  same time window, `<whole>.j<NN>.ring.all.window.lammpstrj.gz`, and the
  whole's data file, `<whole>.ring.all.data.gz`.
- Tier C, in `<POLYPHYS_DATA_DIR>/<project>/trjs/<whole>/`, with
  `<project>` the study's project name, such as `HnsCub`: the `all` windows
  and data files of the other crowded wholes, under the runs' own names.

The HnsCyl runs name their wholes without `epshc`, which the `HnsCyl` parser
reads; the tier B copies insert `epshc1`, the value of their H-NS
core-crowder pair coefficients.

A polymer dump holds only the polymer group, so it is opened with the data
file of the whole without crowders of its family and replicate, whose atoms
are that group; the polymer atoms carry the IDs 1 to n in every system. The
`RealWhole` cases and the fixtures that open the files are in
`polyphys/tests/probe/conftest.py`.

The four TransFociCub wholes also have a tier C reference: the research
code's probe of the same 50 frames of their polymer dumps, under
`<POLYPHYS_DATA_DIR>/TransFociCub/probes/<space>-bug-probe/<whole>.ring/`,
and its reduction of that probe, under
`<POLYPHYS_DATA_DIR>/TransFociCub/analysis/<space>-bug-{ens,ensAvg,space}/`.
Tier C tests compare each stage of the `bug` pipeline with it: the probe
tests (`test_chain_size.py`, `test_pairs.py`, `test_core.py`) compare the
chain-size observables, the principal axes up to sign, the foci bonds and
clusters histograms, and the stamps with the research probe;
`test_pipeline.py` reduces the research probe, and the package's own probe
of the same dumps, to the research `ensAvg`, and builds the research space
tables from the research `ensAvg`; and `viz/test_helpers.py` draws a space
table with the helpers. The foci pair distances and the autocorrelation
functions are left out on purpose: PolyPhys takes the minimum image, which
the research code did not, and `acf` rejects lags beyond a 50-frame
series, which the research code's `statsmodels` call clamped. The
`research_probe_file` fixture finds the probe files.

## Tests that need LAMMPS

A test that runs LAMMPS is marked `lammps` and takes the `lammps_executable`
fixture:

```python
import pytest


@pytest.mark.lammps
def test_lammps_runs_from_a_built_start(lammps_executable, tmp_path):
    subprocess.run([lammps_executable, "-in", "in.kg", ...], cwd=tmp_path)
```

- The executable is the one `LMP` names, else `lmp_serial` on `PATH`.
- When neither exists, the test is skipped. CI has no LAMMPS.
- When `LMP` names something that is not an executable, the test errors.

The input scripts live in `polyphys/tests/data/lammps_build/`; the tests
build their start files in `tmp_path` and read the log LAMMPS writes there.
Each such test takes LAMMPS itself as its oracle: the run finishes with no
error and no unexpected warning.

## Tests that need local data

A tier C test is marked `private_data` and takes the `private_data_dir`
fixture:

```python
import pytest


@pytest.mark.private_data
def test_profiles_match_the_research_code(private_data_dir):
    folder = private_data_dir / "TransFociCyl" / "golden" / "density" / "r"
    ...
```

`POLYPHYS_DATA_DIR` names the maintainer's tier C folder, which holds one
folder per study, named by its project (`HnsCub`, `HnsCyl`, `TransFociCub`,
`TransFociCyl`). Each study folder holds the phase folders of
`make_database`, such as `trjs/`, `probes/`, and `analysis/`, and `golden/`
for test inputs derived from them by the maintainer.

- When `POLYPHYS_DATA_DIR` is unset, the test is skipped. CI never sets it.
- When it names a directory, the test runs and `private_data_dir` is that
  directory.
- When it names something that is not a directory, the test errors, so a
  misconfigured machine is not mistaken for one without data.

Keep the path out of tracked files; set the variable in your shell.

## pytest settings

`pyproject.toml` sets these options under `[tool.pytest.ini_options]`:

| Setting                  | Effect                                                                                       |
| ------------------------ | -------------------------------------------------------------------------------------------- |
| `--strict-markers`       | An unregistered marker, such as a misspelled `private_data`, is an error rather than a warning |
| `--strict-config`        | A warning while reading the pytest configuration is an error                                 |
| `xfail_strict = true`    | A test marked `xfail` that passes fails the run                                              |
| `-ra`                    | The summary lists every test that did not pass, with the reason for each skip                |
| `--import-mode=importlib` | Test modules are imported without changing `sys.path`                                       |
| `-p pytester`            | Loads the `pytester` fixture, which `test_conftest.py` uses to run nested pytest sessions    |
| `--doctest-modules`      | Runs the `Examples` sections of the docstrings in the modules pytest collects               |
| `--doctest-glob=README.md` | Runs the examples of `README.md` as doctests                                              |
| `testpaths = ["polyphys/tests"]` | The directory pytest collects when no path is given                                 |
| `filterwarnings = ["error"]` | Every warning raised during a test is an error; a test that expects one wraps the call in `pytest.warns`, and a test module that must ignore one names it in a `filterwarnings` mark |

The documented command passes `polyphys README.md projects`, so the
modules' docstrings, the README, and the studies are collected as well as
`polyphys/tests`. The examples of five pages, `quickstart.md`,
`lineage.md`, `studies.md`, `usage.md`, and `viz.md`, run through
`polyphys/tests/test_docs.py`.

See [Working with custom markers](https://docs.pytest.org/en/stable/how-to/mark.html)
and the [configuration reference](https://docs.pytest.org/en/stable/reference/reference.html)
in the pytest documentation.

`filterwarnings` reaches only warnings that pytest's own process raises
while a test runs. The resource tracker of `multiprocessing` reports
leaked semaphores from a process of its own after pytest has exited, so
`test_multiprocessing_run_leaves_no_leaked_semaphore` in
`polyphys/tests/probe/test_core.py` runs a parallel probe in a subprocess
and checks its standard error instead.
