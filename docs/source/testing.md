# Testing

This page describes how PolyPhys is tested: which data a test may use, the
kinds of test the project writes, the shared fixtures, and the pytest settings.

## Running the tests

From the repository root, in an environment with the `dev` extra installed:

```bash
python -m pytest polyphys README.md --cov=polyphys --cov-report=term-missing --doctest-modules --doctest-glob='README.md'
```

The README's examples run as doctests, and so do the `Examples` sections of the
docstrings in `polyphys/`.

CI runs the same command with `--cov-fail-under=95`, so a pull request fails
when total coverage drops below 95 percent. coverage.py rounds the total to a
whole percent before it compares, so 94.6 percent passes. Add
`--cov-fail-under=95` to the command above to run the same check locally.

## Test data: three tiers

Every test uses data from one of three tiers.

| Tier | Data                                                                                                                                                   | Where it lives                                   | Who runs it                    |
| ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------ | ------------------------------ |
| A    | Synthetic inputs a test generates itself: lineage trees of tiny files, short LAMMPS files, random or analytic arrays                                   | the test module or `polyphys/tests/conftest.py`  | everyone, including CI         |
| B    | Small trimmed real files, used only where a file format is under test, such as a few dump frames or a filename listing; the maintainer approves each one | `polyphys/tests/data/`, each about 200 KB or less | everyone, including CI         |
| C    | Full-size processed artifacts, such as a whole space                                                                                                   | outside the repository, under `POLYPHYS_DATA_DIR` | the maintainer, locally; never CI |

Tier A is the default. The algorithms are what is under test, and properties
such as conservation, symmetry, invariance, or agreement with an analytic result
do not need real data. Tier C tests compare a result against the maintainer's
processed data, for example when a module replaces code that produced published
numbers.

## Kinds of test

| Kind                       | Catches                                               | Notes                                                                                                                   |
| -------------------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Unit test                  | a function breaking its own contract                  | tier A data; one behaviour per test                                                                                     |
| Known-answer test          | code that runs but computes the wrong physics         | the expected value comes from an analytic result or a cited reference                                                   |
| Property test              | edge cases nobody wrote down                          | checks an invariant, such as conservation of histogram counts, on hand-picked inputs                                    |
| Integration test           | lineage levels that pass alone but disagree at a seam | runs a lineage tree from `segment` up through the organizer                                                             |
| Regression test            | a fixed bug returning                                 | written before the fix and seen to fail against the unfixed code                                                        |
| Comparison against old code | a port that changed a number                         | tier C; the old and new implementations run on the same data, and any number that changes is explained in the pull request |
| Doctest                    | examples drifting from behaviour                      | every public name has an `Examples` section                                                                             |
| Static check               | type and style errors                                 | `flake8` and `mypy`                                                                                                     |

Expected values never come from the code under test. Floating-point results are
compared with a stated tolerance, using `numpy.testing.assert_allclose` for
arrays and `pytest.approx` for scalars. Uncertainties of quantities measured over
correlated molecular-dynamics frames use block averaging or another
autocorrelation-aware estimate, not a per-frame standard error.

## Shared fixtures

`polyphys/tests/conftest.py` provides tier A fixtures to every test module.

`make_lineage_tree(parser, ...)`
: Writes a lineage tree of tiny `.npy` or `.csv` files, named for a parser
  class from {mod}`polyphys.manage.parser`, under the test's `tmp_path`. The
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

The fixtures are factories, called inside a test with the arguments that test
needs; see [Factories as fixtures](https://docs.pytest.org/en/stable/how-to/fixtures.html#factories-as-fixtures)
in the pytest documentation. Their tests are in `polyphys/tests/test_conftest.py`.

## Tests that need local data

A tier C test is marked `private_data` and takes the `private_data_dir`
fixture:

```python
import pytest


@pytest.mark.private_data
def test_space_matches_processed_data(private_data_dir):
    space_dir = private_data_dir / "SumRuleCyl-allInOne"
    ...
```

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
| `--strict-markers`       | An unregistered marker, such as a misspelt `private_data`, is an error rather than a warning |
| `--strict-config`        | A warning while reading the pytest configuration is an error                                 |
| `xfail_strict = true`    | A test marked `xfail` that passes fails the run                                              |
| `-ra`                    | The summary lists every test that did not pass, with the reason for each skip                |
| `--import-mode=importlib` | Test modules are imported without changing `sys.path`                                       |
| `-p pytester`            | Loads the `pytester` fixture, which `test_conftest.py` uses to run nested pytest sessions    |

See [Working with custom markers](https://docs.pytest.org/en/stable/how-to/mark.html)
and the [configuration reference](https://docs.pytest.org/en/stable/reference/reference.html)
in the pytest documentation. Warnings are not yet turned into errors.
