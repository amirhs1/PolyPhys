"""Tests for :mod:`polyphys.simulate.log`.

Tier A tests parse synthetic logs assembled in the test; each expected value
is a number the test writes into the log, so the parser must read it back
exactly, or, for the wall time, an analytic conversion.

Tier B tests read the logs in ``polyphys/tests/data/lammps_log/``, written
by LAMMPS (22 Jul 2025 - Update 4) from the inputs beside them. Their
expected values were read from the log text with ``grep`` and ``awk``, not
from the parser.
"""

import gzip
import logging
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd
import pytest

from polyphys.simulate.log import (
    BrokenLogError,
    _lammps_version,
    parse_log,
    read_log,
)


@pytest.fixture
def make_log() -> Callable[..., str]:
    """Return a factory that assembles a synthetic LAMMPS log.

    Each table, a sequence of lines, becomes one run: a run-start line, the
    table, and that run's summary lines (by default a ``Loop time of`` line
    for 100 steps with 20 atoms on 1 process). The log starts with
    `first_line` and ends with `wall`.
    """
    def _make(
        *tables: Sequence[str],
        summaries: Sequence[Sequence[str]] | None = None,
        first_line: str = "LAMMPS (22 Jul 2025 - Update 4)",
        run_start: str = ("Per MPI rank memory allocation (min/avg/max) = "
                          "4.9 | 4.9 | 4.9 Mbytes"),
        wall: str = "Total wall time: 0:00:01",
    ) -> str:
        if summaries is None:
            summaries = [
                ["Loop time of 0.5 on 1 procs for 100 steps with 20 atoms"]
            ] * len(tables)
        lines = [first_line]
        for table, summary in zip(tables, summaries, strict=True):
            lines += [run_start, *table, *summary]
        lines.append(wall)
        return "\n".join(lines) + "\n"
    return _make


@pytest.fixture
def log_dir() -> Path:
    """Return the directory of the LAMMPS log fixtures."""
    return Path(__file__).resolve().parents[1] / "data" / "lammps_log"


# Thermo tables --------------------------------------------------------------


def test_parse_log_one_style_reads_keywords_and_values(make_log):
    """A one-style table becomes one column per keyword, Step as int64."""
    text = make_log(["   Step   Temp   Press",
                     "      0    1.0     0.5",
                     "    100    0.9     0.4"])

    log = parse_log(text)

    expected = pd.DataFrame({"run": [0, 0], "Step": [0, 100],
                             "Temp": [1.0, 0.9], "Press": [0.5, 0.4]})
    pd.testing.assert_frame_equal(
        log.thermo, expected.astype({"run": "int64", "Step": "int64"}))


def test_parse_log_custom_keywords_kept_verbatim(make_log):
    """Custom keywords with digits and brackets keep their exact names."""
    text = make_log(["Step c_gyr[1] c_gyr[2] v_x1 c_1 f_3[2]",
                     "0 1.5 2.5 3.5 4.5 5.5"])

    log = parse_log(text)

    assert list(log.thermo.columns) == [
        "run", "Step", "c_gyr[1]", "c_gyr[2]", "v_x1", "c_1", "f_3[2]"]


def test_parse_log_multi_style_reads_sections(make_log):
    """A multi-style section gives Step, CPU, and each Keyword = value."""
    text = make_log([
        "------------ Step              0 ----- CPU =            0 (sec) "
        "-------------",
        "TotEng   =        -1.5000 KinEng   =         1.5000 Temp     = "
        "        1.0000 ",
        "PotEng   =        -3.0000 Press    =         0.2500",
        "------------ Step            100 ----- CPU =      0.01234 (sec) "
        "-------------",
        "TotEng   =        -1.2500 KinEng   =         1.2500 Temp     = "
        "        0.8000 ",
        "PotEng   =        -2.5000 Press    =         0.1250",
    ])

    log = parse_log(text)

    expected = pd.DataFrame({
        "run": [0, 0], "Step": [0, 100], "CPU": [0.0, 0.01234],
        "TotEng": [-1.5, -1.25], "KinEng": [1.5, 1.25], "Temp": [1.0, 0.8],
        "PotEng": [-3.0, -2.5], "Press": [0.25, 0.125],
    })
    pd.testing.assert_frame_equal(
        log.thermo, expected.astype({"run": "int64", "Step": "int64"}))


def test_parse_log_yaml_style_reads_rows(make_log):
    """A yaml-style table takes its keywords from the keywords line."""
    text = make_log(["---",
                     "keywords: ['Step', 'Temp', 'c_rg[1]', ]",
                     "data:",
                     "  - [0, 1, 2.5, ]",
                     "  - [100, 0.75, 2.25, ]",
                     "..."])

    log = parse_log(text)

    expected = pd.DataFrame({"run": [0, 0], "Step": [0, 100],
                             "Temp": [1.0, 0.75], "c_rg[1]": [2.5, 2.25]})
    pd.testing.assert_frame_equal(
        log.thermo, expected.astype({"run": "int64", "Step": "int64"}))


@pytest.mark.parametrize(
    ("table", "expected_steps"),
    [
        pytest.param(
            ["Step Temp",
             "0 1.0",
             "WARNING: Lost atoms: original 20 current 19 "
             "(src/thermo.cpp:532)",
             "For more information see https://docs.lammps.org/err0008 "
             "(src/thermo.cpp:532)",
             "100 0.9"],
            [0, 100], id="warning-in-one-table"),
        pytest.param(
            ["WARNING: Communication cutoff is shorter than a bond length "
             "(src/comm.cpp:743)",
             "Step Temp",
             "0 1.0"],
            [0], id="warning-before-header"),
        pytest.param(
            ["------------ Step              0 ----- CPU =            0 "
             "(sec) -------------",
             "Temp     =         1.0000 ",
             "WARNING: Lost atoms: original 20 current 19 "
             "(src/thermo.cpp:532)",
             "------------ Step            100 ----- CPU =         0.01 "
             "(sec) -------------",
             "Temp     =         0.9000 "],
            [0, 100], id="warning-in-multi-table"),
        pytest.param(
            ["Step Temp",
             "0 1.0",
             "Fix halt condition for fix-id F_halt met on step 50 with "
             "value 50"],
            [0], id="fix-halt-message"),
    ],
)
def test_parse_log_skips_messages_in_table(make_log, table, expected_steps):
    """WARNING and fix halt lines inside a table are not read as data."""
    text = make_log(table)

    log = parse_log(text)

    assert log.thermo["Step"].tolist() == expected_steps


def test_parse_log_logs_skipped_message_with_line_number(make_log, caplog):
    """A skipped WARNING line is logged at INFO with its line number."""
    text = make_log(["Step Temp", "0 1.0",
                     "WARNING: Lost atoms: original 20 current 19"])

    with caplog.at_level(logging.INFO, logger="polyphys.simulate.log"):
        parse_log(text)

    assert "line 5" in caplog.text
    assert "WARNING: Lost atoms" in caplog.text


def test_parse_log_reads_nan_and_inf(make_log):
    """nan, -nan, inf, and -inf values are read as such, not dropped."""
    text = make_log(["Step A B C D",
                     "0 nan -nan inf -inf"])

    log = parse_log(text)

    np.testing.assert_array_equal(
        log.thermo.loc[0, ["A", "B", "C", "D"]].to_numpy(dtype=float),
        [np.nan, np.nan, np.inf, -np.inf])


def test_parse_log_joins_runs_with_run_column(make_log):
    """Runs are joined in log order, keeping each run's repeated first step."""
    text = make_log(["Step Temp", "0 1.0", "100 0.9"],
                    ["Step Temp", "100 0.9", "200 0.8"])

    log = parse_log(text)

    assert log.thermo["run"].tolist() == [0, 0, 1, 1]
    assert log.thermo["Step"].tolist() == [0, 100, 100, 200]


def test_parse_log_keyword_absent_from_run_is_nan(make_log):
    """A keyword missing from one run is NaN in that run's rows."""
    text = make_log(["Step Temp Press", "0 1.0 0.5"],
                    ["Step Temp", "100 0.9"])

    log = parse_log(text)

    np.testing.assert_array_equal(log.thermo["Press"].to_numpy(),
                                  [0.5, np.nan])


@pytest.mark.parametrize(
    ("first_line", "expected"),
    [
        pytest.param("LAMMPS (22 Jul 2025 - Update 4)",
                     "22 Jul 2025 - Update 4", id="recent"),
        pytest.param("LAMMPS (29 Oct 2020)", "29 Oct 2020", id="older"),
        pytest.param("units lj", None, id="opened-by-log-command"),
    ],
)
def test_parse_log_reads_version_from_first_line(make_log, first_line,
                                                 expected):
    """The version comes from a first line LAMMPS (<version>), else None."""
    text = make_log(["Step", "0"], first_line=first_line)

    log = parse_log(text)

    assert log.lammps_version == expected


def test_lammps_version_of_empty_log_is_none():
    """An empty list of lines has no version line."""
    assert _lammps_version([]) is None


def test_parse_log_accepts_old_run_start_line(make_log):
    """The run-start line of old versions also opens a run."""
    text = make_log(["Step Temp", "0 1.0"],
                    run_start="Memory usage per processor = 2.1 Mbytes")

    log = parse_log(text)

    assert log.thermo["Step"].tolist() == [0]


# Broken logs -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("tables", "options", "match"),
    [
        pytest.param(
            [["Step Temp", "0 1.0"]], {"summaries": [[]]},
            r"Run 0 \(line 2\) has no 'Loop time of' line", id="no-footer"),
        pytest.param(
            [["Step Temp", "0 1.0"]], {"wall": ""},
            r"no 'Total wall time:' line", id="no-wall-time"),
        pytest.param(
            [["Step Temp", "0 1.0"]],
            {"wall": "# run ends with Total wall time: 0:00:01"},
            r"no 'Total wall time:' line", id="wall-time-only-in-comment"),
        pytest.param(
            [], {}, r"The log has no run", id="no-run"),
        pytest.param(
            [["Step Temp", "0 1.0",
              "ERROR: Lost atoms: original 20 current 19"]], {},
            r"Line 5 reports an error", id="error-line"),
        pytest.param(
            [["Step Temp", "0 1.0", "fix print: Rg 5.2"]], {},
            r"Line 5 is inside a thermo table but is not data",
            id="stray-line"),
        pytest.param(
            [["Step Temp", "0 1.0", "100"]], {},
            r"Line 5 is inside a thermo table but is not data",
            id="short-row"),
        pytest.param(
            [["Step Temp", "0 abc"]], {},
            r"Line 4 is inside a thermo table but is not data",
            id="non-numeric-value"),
        pytest.param(
            [["WARNING: Lost atoms: original 20 current 19"]], {},
            r"Run 0 \(line 2\) has no thermo header", id="no-header"),
        pytest.param(
            [["------------ Step              0 ----- CPU =            0 "
              "(sec) -------------",
              "fix print: Rg 5.2"]], {},
            r"Line 4 is inside a thermo table but is not data",
            id="multi-stray-line"),
        pytest.param(
            [["---", "keywords: ['Step', ]", "data:", "fix print: Rg 5.2"]],
            {}, r"Line 6 is inside a thermo table but is not data",
            id="yaml-stray-line"),
        pytest.param(
            [["---", "data:", "..."]], {},
            r"Run 0 \(line 2\) has no 'keywords:' line",
            id="yaml-without-keywords"),
    ],
)
def test_parse_log_raises_broken_log_error(make_log, tables, options, match):
    """A cut-off or malformed log raises BrokenLogError naming the problem."""
    text = make_log(*tables, **options)

    with pytest.raises(BrokenLogError, match=match):
        parse_log(text)


# Run summary -----------------------------------------------------------------


def test_parse_log_reads_run_summary(make_log):
    """Every field of a full run summary lands in its run_stats column."""
    summary = [
        "Loop time of 2.5 on 4 procs for 1000 steps with 2000 atoms",
        "",
        "Performance: 172800.000 tau/day, 400.000 timesteps/s, "
        "800.000 katom-step/s",
        "100.0% CPU use with 4 MPI tasks x 1 OpenMP threads",
        "",
        "MPI task timing breakdown:",
        "Section |  min time  |  avg time  |  max time  |%varavg| %total",
        "---------------------------------------------------------------",
        "Pair    | 1.0        | 1.2        | 1.4        |   5.0 | 48.00",
        "Other   |            | 0.05       |            |       |  2.00",
        "",
        "Dangerous builds = 3",
    ]
    text = make_log(["Step", "0"], summaries=[summary])

    log = parse_log(text)

    assert log.run_stats.iloc[0].to_dict() == {
        "run": 0, "loop_time": 2.5, "n_cores": 4, "loop_timesteps": 1000,
        "n_atoms": 2000, "perf_per_day": 172800.0, "perf_unit": "tau/day",
        "timesteps_per_s": 400.0, "cpu_use": 100.0, "mpi_task": 4,
        "openmp_threads": 1, "Pair_avg_time": 1.2, "Pair_total_pct": 48.0,
        "Other_avg_time": 0.05, "Other_total_pct": 2.0,
        "dangerous_builds": 3,
    }


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        pytest.param(
            "Performance: 172800.000 tau/day, 400.000 timesteps/s, "
            "800.000 katom-step/s",
            (172800.0, "tau/day", 400.0), id="lj"),
        pytest.param(
            "Performance: 172800.000 tau/day, 400.000 timesteps/s",
            (172800.0, "tau/day", 400.0), id="lj-without-atom-steps"),
        pytest.param(
            "Performance: 0.512 ns/day, 46.875 hours/ns, 5.926 timesteps/s, "
            "11.852 katom-step/s",
            (0.512, "ns/day", 5.926), id="real"),
        pytest.param(
            "Performance: 3.200 fs/day, 7.500 hours/fs, 37.037 timesteps/s, "
            "74.074 katom-step/s",
            (3.2, "fs/day", 37.037), id="electron"),
    ],
)
def test_parse_log_reads_performance_in_any_unit(make_log, line, expected):
    """The performance line gives its rate, unit, and timesteps per second."""
    summary = ["Loop time of 2.5 on 4 procs for 1000 steps with 2000 atoms",
               "", line]
    text = make_log(["Step", "0"], summaries=[summary])

    log = parse_log(text)

    row = log.run_stats.iloc[0]
    assert (row["perf_per_day"], row["perf_unit"],
            row["timesteps_per_s"]) == expected


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        pytest.param("100.0% CPU use with 4 MPI tasks x 1 OpenMP threads",
                     (100.0, 4, 1), id="full-use"),
        pytest.param("9.5% CPU use with 2 MPI tasks x no OpenMP threads",
                     (9.5, 2, 0), id="no-openmp"),
    ],
)
def test_parse_log_reads_cpu_line(make_log, line, expected):
    """The CPU line gives use, MPI tasks, and threads (0 without OpenMP)."""
    summary = ["Loop time of 2.5 on 4 procs for 1000 steps with 2000 atoms",
               "", line]
    text = make_log(["Step", "0"], summaries=[summary])

    log = parse_log(text)

    row = log.run_stats.iloc[0]
    assert (row["cpu_use"], row["mpi_task"],
            row["openmp_threads"]) == expected


@pytest.mark.parametrize(
    ("header", "row"),
    [
        pytest.param(
            "Section |  min time  |  avg time  |  max time  |%varavg| %total",
            "Pair    | 1.0        | 1.2        | 1.4        |   5.0 | 48.00",
            id="without-cpu-column"),
        pytest.param(
            "Section |  min time  |  avg time  |  max time  |%varavg|  %CPU "
            "| %total",
            "Pair    | 1.0        | 1.2        | 1.4        |   5.0 |  99.0 "
            "| 48.00",
            id="with-cpu-column"),
    ],
)
def test_parse_log_finds_timing_columns_by_name(make_log, header, row):
    """Average time and percent of total are found whichever columns exist."""
    summary = ["Loop time of 2.5 on 4 procs for 1000 steps with 2000 atoms",
               "", "MPI task timing breakdown:", header, "-" * 63, row, ""]
    text = make_log(["Step", "0"], summaries=[summary])

    log = parse_log(text)

    row_stats = log.run_stats.iloc[0]
    assert (row_stats["Pair_avg_time"], row_stats["Pair_total_pct"]) == (
        1.2, 48.0)


def test_parse_log_skips_malformed_timing_rows(make_log):
    """Timing rows without a section or with too few cells are ignored."""
    summary = [
        "Loop time of 2.5 on 4 procs for 1000 steps with 2000 atoms", "",
        "MPI task timing breakdown:",
        "Section |  min time  |  avg time  |  max time  |%varavg| %total",
        "-" * 63,
        "Pair    | 1.0        | 1.2",
        "        | 1.0        | 1.2        | 1.4        |   5.0 | 48.00",
        "Neigh   | 1.0        | n/a        | 1.4        |   5.0 | 12.00",
        "",
    ]
    text = make_log(["Step", "0"], summaries=[summary])

    log = parse_log(text)

    row = log.run_stats.iloc[0]
    timing = [col for col in log.run_stats.columns
              if col.endswith(("_avg_time", "_total_pct"))]
    assert timing == ["Neigh_avg_time", "Neigh_total_pct"]
    assert np.isnan(row["Neigh_avg_time"])
    assert row["Neigh_total_pct"] == 12.0


def test_parse_log_missing_summary_fields_are_nan(make_log):
    """A run without a performance or dangerous-builds count gets NaN."""
    minimize = ["Loop time of 0.1 on 1 procs for 5 steps with 20 atoms",
                "Dangerous builds not checked"]
    run = ["Loop time of 0.5 on 1 procs for 100 steps with 20 atoms",
           "Performance: 172800.000 tau/day, 400.000 timesteps/s",
           "Dangerous builds = 2"]
    text = make_log(["Step", "0"], ["Step", "5"], summaries=[minimize, run])

    log = parse_log(text)

    np.testing.assert_array_equal(log.run_stats["perf_per_day"].to_numpy(),
                                  [np.nan, 172800.0])
    np.testing.assert_array_equal(
        log.run_stats["dangerous_builds"].to_numpy(), [np.nan, 2.0])


def test_parse_log_wall_time_in_seconds_and_hours(make_log):
    """1:02:03 is 3723 s; cores and atoms come from the last run."""
    summaries = [["Loop time of 0.5 on 4 procs for 100 steps with 20 atoms"],
                 ["Loop time of 0.5 on 8 procs for 100 steps with 19 atoms"]]
    text = make_log(["Step", "0"], ["Step", "100"], summaries=summaries,
                    wall="Total wall time: 1:02:03")

    log = parse_log(text)

    row = log.wall_time.iloc[0]
    assert row["wall_time_s"] == 3723
    # 3723 / 3600 h; tolerance covers float64 round-off only.
    assert row["wall_time_hr"] == pytest.approx(3723 / 3600, rel=1e-12)
    assert (row["n_cores"], row["n_atoms"]) == (8, 19)


# Reading files ---------------------------------------------------------------


@pytest.mark.parametrize("as_type", [pytest.param(str, id="str"),
                                     pytest.param(Path, id="path")])
def test_read_log_accepts_str_and_path(make_log, tmp_path, as_type):
    """read_log takes a str or a Path and stores the file as a Path."""
    path = tmp_path / "log.lammps"
    path.write_text(make_log(["Step", "0"]))

    log = read_log(as_type(path))

    assert log.path == path


def test_read_log_reads_a_gzipped_log(make_log, tmp_path):
    """A log whose name ends in .gz is read through gzip."""
    path = tmp_path / "log.lammps.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(make_log(["Step Temp", "0 1.5", "100 0.5"]))

    log = read_log(path)

    assert log.thermo["Temp"].tolist() == [1.5, 0.5]


# LAMMPS-written fixtures (tier B) -------------------------------------------


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        pytest.param("log.one", [2, 7, 7, 7, 7, 12], id="one"),
        pytest.param("log.multi", [6, 6, 6, 6, 6], id="multi"),
        pytest.param("log.custom", [6, 6, 3, 7], id="custom"),
        pytest.param("log.yaml", [6, 6], id="yaml"),
    ],
)
def test_read_log_fixture_rows_per_run(log_dir, name, expected):
    """Each run of a LAMMPS-written log keeps all its thermo rows."""
    log = read_log(log_dir / name)

    assert log.thermo.groupby("run").size().tolist() == expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        pytest.param("log.one", [0, 5, 505, 1005, 1505, 2005], id="one"),
        pytest.param("log.multi", [0, 0, 500, 1000, 1500],
                     id="multi-reset-timestep"),
        pytest.param("log.custom", [0, 500, 1000, 1250], id="custom"),
        pytest.param("log.yaml", [0, 500], id="yaml"),
    ],
)
def test_read_log_fixture_first_step_per_run(log_dir, name, expected):
    """Each run starts at the step LAMMPS printed first for it."""
    log = read_log(log_dir / name)

    assert log.thermo.groupby("run")["Step"].first().tolist() == expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        pytest.param(
            "log.one",
            ["run", "Step", "Temp", "E_pair", "E_mol", "TotEng", "Press"],
            id="one"),
        pytest.param(
            "log.multi",
            ["run", "Step", "CPU", "TotEng", "KinEng", "Temp", "PotEng",
             "E_bond", "E_angle", "E_dihed", "E_impro", "E_vdwl", "E_coul",
             "E_long", "Press"],
            id="multi"),
        pytest.param(
            "log.custom",
            ["run", "Step", "KinEng", "PotEng", "E_tail", "Ecouple",
             "Econserve", "E_vdwl", "E_bond", "c_rg", "c_rg[1]", "c_rg[3]",
             "v_rg2"],
            id="custom"),
        pytest.param(
            "log.yaml",
            ["run", "Step", "Temp", "KinEng", "PotEng", "E_bond", "E_angle",
             "E_dihed", "E_impro", "E_vdwl", "E_coul", "E_long", "Press"],
            id="yaml"),
    ],
)
def test_read_log_fixture_keywords(log_dir, name, expected):
    """The thermo columns are the keywords LAMMPS printed, in order."""
    log = read_log(log_dir / name)

    assert list(log.thermo.columns) == expected


@pytest.mark.parametrize(
    ("name", "steps", "atoms"),
    [
        pytest.param("log.one", [5, 500, 500, 500, 500, 1000],
                     [170, 170, 170, 170, 170, 123], id="one"),
        pytest.param("log.multi", [500] * 5, [320] * 5, id="multi"),
        pytest.param("log.custom", [500, 500, 250, 500], [320] * 4,
                     id="custom-fix-halt"),
        pytest.param("log.yaml", [500, 500], [320, 320], id="yaml"),
    ],
)
def test_read_log_fixture_loop_lines(log_dir, name, steps, atoms):
    """Each run's steps and atoms come from its Loop time line."""
    log = read_log(log_dir / name)

    assert log.run_stats["loop_timesteps"].tolist() == steps
    assert log.run_stats["n_atoms"].tolist() == atoms


def test_read_log_fixture_first_row_values(log_dir):
    """The first row of log.one equals the numbers LAMMPS printed."""
    log = read_log(log_dir / "log.one")

    assert log.thermo.iloc[0].tolist() == [
        0, 0, 0.0, 0.020390867, 2.2622954, 2.2826862, -0.0049223331]


def test_read_log_fixture_unchecked_builds_are_nan(log_dir):
    """The run with neighbour checks off has no dangerous-builds count."""
    log = read_log(log_dir / "log.custom")

    np.testing.assert_array_equal(
        log.run_stats["dangerous_builds"].to_numpy(), [0, np.nan, 0, 0])


@pytest.mark.parametrize(
    ("name", "match"),
    [
        pytest.param("log.error",
                     r"Line \d+ reports an error: ERROR: Lost atoms",
                     id="error"),
        pytest.param("log.killed",
                     r"Run 0 \(line \d+\) has no 'Loop time of' line",
                     id="killed"),
    ],
)
def test_read_log_fixture_broken_log_raises(log_dir, name, match):
    """A LAMMPS run that stopped with an error or was killed raises."""
    with pytest.raises(BrokenLogError, match=match):
        read_log(log_dir / name)
