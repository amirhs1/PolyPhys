"""\
=====================================================================
LAMMPS log files --- :mod:`polyphys.manage.lammps_log`
=====================================================================

The :mod:`~polyphys.manage.lammps_log` module reads a LAMMPS log file into
pandas DataFrames: the thermodynamic output of every run, the summary LAMMPS
prints at the end of each run, and the total wall time.

A log holds one *run* per ``run`` or ``minimize`` command. Each run starts
with the line ``Per MPI rank memory allocation ...`` (``Memory usage per
processor ...`` in old versions), continues with the thermo table, and ends
with the line ``Loop time of ...``, after which comes the run summary. The
log ends with ``Total wall time: H:MM:SS``. Every marker is matched at the
start of a line only, because the log also echoes the input script, whose
comments may contain the same words.

Contents
========

:func:`read_log` reads a file and :func:`parse_log` reads text; both return
a :class:`LammpsLog`. :class:`BrokenLogError` reports a log that is cut off
or does not follow the format.

Thermo styles
=============

- ``one`` and ``custom``: a header line of keywords, then one row of values
  per output step, separated by spaces.
- ``multi``: one section per output step, opened by a line
  ``------------ Step N ----- CPU = t (sec) -------------`` and followed by
  ``Keyword = value`` pairs, three per line.
- ``yaml``: ``---``, ``keywords: [...]``, ``data:``, one ``  - [...]`` line
  per output step, and ``...``.

``WARNING`` lines and the ``Fix halt condition ...`` message may appear
inside a table, and recent versions follow each with a ``For more
information see ...`` line; all three are skipped, and the first two are
logged. Blank lines are skipped. Any other line inside a table that is not
data raises :class:`BrokenLogError`.

Complexity
==========

Time and memory are :math:`O(L)` in the length :math:`L` of the log: the
file is read once and each line is visited a bounded number of times.

References
==========

The formats follow the LAMMPS documentation of ``thermo_style``,
``thermo_modify``, and ``run``, and these output statements in the LAMMPS
source, release ``stable_30Sep2026``:

- Run start: `output.cpp, line 1214
  <https://github.com/lammps/lammps/blob/8de817dd79bfe4525d5d39246a212d833e6dee07/src/output.cpp#L1214>`_.
- Thermo header and rows: `thermo.cpp, lines 91, 256 to 325, and
  387 to 419
  <https://github.com/lammps/lammps/blob/8de817dd79bfe4525d5d39246a212d833e6dee07/src/thermo.cpp#L387-L419>`_.
- Run summary: `finish.cpp, lines 128 to 185, 321 to 356, and 543 to 544
  <https://github.com/lammps/lammps/blob/8de817dd79bfe4525d5d39246a212d833e6dee07/src/finish.cpp#L128-L185>`_.
- Total wall time: `lammps.cpp, line 796
  <https://github.com/lammps/lammps/blob/8de817dd79bfe4525d5d39246a212d833e6dee07/src/lammps.cpp#L796>`_.
- Fix halt message: `fix_halt.cpp, lines 313 to 320
  <https://github.com/lammps/lammps/blob/8de817dd79bfe4525d5d39246a212d833e6dee07/src/fix_halt.cpp#L313-L320>`_.

Dependencies
============
- Python standard library: `dataclasses`, `enum`, `logging`, `pathlib`,
  `re`
- Third-party: `pandas`
"""

import logging
import math
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterator, Sequence

import pandas as pd

__all__ = ["BrokenLogError", "LammpsLog", "parse_log", "read_log"]

logger = logging.getLogger(__name__)

# Line-start markers, from the LAMMPS source lines cited in the module
# docstring.
_RUN_START = ("Per MPI rank memory allocation", "Memory usage per processor")
_RUN_END = "Loop time of"
_WALL_TIME = "Total wall time:"
_ERROR = "ERROR"
_SKIPPED_MESSAGES = ("WARNING", "Fix halt condition")
_MORE_INFO = "For more information see"
_TIMING_TABLE = "MPI task timing breakdown:"

_RUN_STATS_LEADING = (
    "run", "loop_time", "n_cores", "loop_timesteps", "n_atoms",
    "perf_per_day", "perf_unit", "timesteps_per_s",
    "cpu_use", "mpi_task", "openmp_threads",
)

_VERSION = re.compile(r"LAMMPS \((?P<version>.+)\)\s*")
_MULTI_HEADER = re.compile(
    r"-+ Step\s+(?P<step>\S+) -+ CPU =\s+(?P<cpu>\S+) \(sec\) -+\s*"
)
_MULTI_PAIRS = re.compile(r"\s*(?:\S+\s+=\s+\S+\s*)+")
_MULTI_PAIR = re.compile(r"(?P<key>\S+)\s+=\s+(?P<value>\S+)")
_LOOP = re.compile(
    r"Loop time of (?P<time>\S+) on (?P<procs>\d+) procs for "
    r"(?P<steps>\d+) steps with (?P<atoms>\d+) atoms"
)
_PERFORMANCE = re.compile(
    r"Performance: (?P<per_day>\S+) (?P<unit>\S+/day), (?:.*?, )?"
    r"(?P<steps_per_s>\S+) timesteps/s"
)
_CPU = re.compile(
    r"(?P<cpu>\S+)% CPU use with (?P<tasks>\d+) MPI tasks x "
    r"(?P<threads>\d+|no) OpenMP threads"
)
_DANGEROUS = re.compile(r"Dangerous builds = (?P<count>\d+)")
_WALL = re.compile(
    r"Total wall time: (?P<hours>\d+):(?P<minutes>\d{2}):(?P<seconds>\d{2})"
)


class BrokenLogError(Exception):
    """A LAMMPS log is cut off or does not follow the log format.

    Examples
    --------
    A log killed during its first run has no ``Loop time of`` line:

    >>> text = '''LAMMPS (22 Jul 2025 - Update 4)
    ... Per MPI rank memory allocation (min/avg/max) = 4.9 | 4.9 | 4.9 Mbytes
    ...    Step          Temp
    ...          0   1
    ...        100   0.5
    ... '''
    >>> try:
    ...     parse_log(text)
    ... except BrokenLogError as err:
    ...     print(err)
    Run 0 (line 2) has no 'Loop time of' line; the log may be cut off.
    """


class _ThermoStyle(Enum):
    """Layout of the thermo table of one run."""

    ONE = "one"  # also custom: a header line, then one row per output step
    MULTI = "multi"
    YAML = "yaml"


@dataclass(frozen=True)
class _RunBlock:
    """The lines of one run.

    Attributes
    ----------
    index : int
        Position of the run in the log, from 0.
    first_line : int
        Line number, from 1, of the run-start line in the log; error
        messages give line numbers counted from it.
    thermo : tuple of str
        The lines after the run-start line, up to the ``Loop time of`` line.
    summary : tuple of str
        The lines from the ``Loop time of`` line up to the next run-start
        line or the end of the log.
    """

    index: int
    first_line: int
    thermo: tuple[str, ...]
    summary: tuple[str, ...]

    def numbered_thermo(self) -> Iterator[tuple[int, str]]:
        """Yield each thermo line with its line number in the log.

        Skipped messages and blank lines are left out; ``WARNING`` and
        ``Fix halt condition`` lines are logged at level INFO.
        """
        for offset, line in enumerate(self.thermo, start=1):
            if not line.strip() or line.startswith(_MORE_INFO):
                continue
            if line.startswith(_SKIPPED_MESSAGES):
                logger.info(
                    "Run %d, line %d: skipped %r",
                    self.index, self.first_line + offset, line,
                )
                continue
            yield self.first_line + offset, line


@dataclass(frozen=True)
class LammpsLog:
    """The parsed contents of one LAMMPS log.

    The DataFrames can still be changed in place; copy one before editing it.

    Attributes
    ----------
    path : pathlib.Path or None
        The log file, or None when the log was parsed from text.
    lammps_version : str or None
        The version on the first line of the log, such as
        ``'22 Jul 2025 - Update 4'``, or None if the first line is not
        ``LAMMPS (<version>)``.
    thermo : pandas.DataFrame
        One row per thermo output step of every run, in log order. Column
        ``run`` (int64) is the position of the run in the log, from 0; the
        other columns are the thermo keywords exactly as LAMMPS prints them,
        such as ``Step``, ``Temp``, or ``c_gyr[1]``. ``Step`` is int64 and
        every other keyword float64. A keyword absent from a run is NaN in
        that run's rows. The first row of each run repeats the step at which
        it starts, so a run that follows another without ``reset_timestep``
        repeats the previous run's last step.
    run_stats : pandas.DataFrame
        One row per run, from its summary:

        - ``run``: position of the run in the log, from 0.
        - ``loop_time``: wall time of the run, in seconds.
        - ``n_cores``: processes times OpenMP threads.
        - ``loop_timesteps``: number of timesteps of the run.
        - ``n_atoms``: number of atoms at the end of the run.
        - ``perf_per_day``, ``perf_unit``, ``timesteps_per_s``: simulated
          time per day of wall time, its unit (``'tau/day'``, ``'fs/day'``,
          or ``'ns/day'``), and timesteps per second of wall time.
        - ``cpu_use``: CPU use, in percent.
        - ``mpi_task``: number of MPI tasks.
        - ``openmp_threads``: OpenMP threads per task, 0 for a build without
          OpenMP.
        - ``<section>_avg_time`` and ``<section>_total_pct``: average time,
          in seconds, and percent of the loop time spent in each section of
          the MPI task timing breakdown, such as ``Pair_avg_time``.
        - ``dangerous_builds``: dangerous neighbour-list builds.

        A field the run does not report is NaN: LAMMPS prints no
        performance line after ``minimize`` or in ``cgs`` and ``si`` units,
        and no count of dangerous builds when neighbour checks are off.
    wall_time : pandas.DataFrame
        One row: ``wall_time_s`` (int64, seconds) and ``wall_time_hr``
        (float64, hours) from the ``Total wall time`` line, and ``n_cores``
        and ``n_atoms`` from the summary of the last run.

    Examples
    --------
    >>> text = '''LAMMPS (22 Jul 2025 - Update 4)
    ... Per MPI rank memory allocation (min/avg/max) = 4.9 | 4.9 | 4.9 Mbytes
    ...    Step          Temp
    ...          0   1
    ...        100   0.5
    ... Loop time of 0.01 on 1 procs for 100 steps with 20 atoms
    ...
    ... Per MPI rank memory allocation (min/avg/max) = 4.9 | 4.9 | 4.9 Mbytes
    ...    Step          Temp
    ...        100   0.5
    ...        200   0.25
    ... Loop time of 0.02 on 1 procs for 100 steps with 20 atoms
    ...
    ... Total wall time: 0:00:01
    ... '''
    >>> log = parse_log(text)
    >>> log.thermo
       run  Step  Temp
    0    0     0  1.00
    1    0   100  0.50
    2    1   100  0.50
    3    1   200  0.25

    Keep the second run only:

    >>> log.thermo[log.thermo["run"] == 1]["Temp"].tolist()
    [0.5, 0.25]
    >>> log.run_stats[["run", "loop_time", "loop_timesteps"]]
       run  loop_time  loop_timesteps
    0    0       0.01             100
    1    1       0.02             100
    """

    path: Path | None
    lammps_version: str | None
    thermo: pd.DataFrame
    run_stats: pd.DataFrame
    wall_time: pd.DataFrame


def read_log(path: str | Path) -> LammpsLog:
    """Read a LAMMPS log file.

    Parameters
    ----------
    path : str or pathlib.Path
        The log file, read as UTF-8; an undecodable byte becomes U+FFFD.

    Returns
    -------
    LammpsLog
        The parsed log, with ``path`` set to `path` as a
        :class:`pathlib.Path`.

    Raises
    ------
    BrokenLogError
        If the log is cut off or does not follow the format; see
        :func:`parse_log`.

    Notes
    -----
    Time and memory are :math:`O(L)` in the length :math:`L` of the file.

    Examples
    --------
    >>> import tempfile
    >>> text = '''LAMMPS (22 Jul 2025 - Update 4)
    ... Per MPI rank memory allocation (min/avg/max) = 4.9 | 4.9 | 4.9 Mbytes
    ...    Step          Temp
    ...          0   1
    ...        100   0.5
    ... Loop time of 0.01 on 1 procs for 100 steps with 20 atoms
    ...
    ... Total wall time: 0:01:05
    ... '''
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "log.lammps"
    ...     _ = path.write_text(text)
    ...     log = read_log(path)
    >>> log.path.name
    'log.lammps'
    >>> log.wall_time[["wall_time_s", "n_atoms"]]
       wall_time_s  n_atoms
    0           65       20
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    return parse_log(text, path)


def parse_log(text: str, path: Path | None = None) -> LammpsLog:
    """Parse the text of a LAMMPS log.

    Parameters
    ----------
    text : str
        The whole log.
    path : pathlib.Path or None, default None
        The file the text came from, stored in the result.

    Returns
    -------
    LammpsLog
        The parsed log.

    Raises
    ------
    BrokenLogError
        If the log has no run; a run has no ``Loop time of`` line, as after a
        walltime kill; an ``ERROR`` line appears; a thermo table holds a line
        that is neither data nor a skipped message; or the log has no
        ``Total wall time`` line.

    Notes
    -----
    Time and memory are :math:`O(L)` in the length :math:`L` of `text`.

    Examples
    --------
    >>> text = '''LAMMPS (22 Jul 2025 - Update 4)
    ... Per MPI rank memory allocation (min/avg/max) = 4.9 | 4.9 | 4.9 Mbytes
    ...    Step          Temp          c_rg[1]
    ...          0   1                 nan
    ... WARNING: Lost atoms: original 20 current 19 (src/thermo.cpp:532)
    ... For more information see https://docs.lammps.org/err0008
    ...        100   0.5               2.5
    ... Loop time of 0.01 on 1 procs for 100 steps with 19 atoms
    ...
    ... Total wall time: 0:00:01
    ... '''
    >>> log = parse_log(text)
    >>> log.lammps_version
    '22 Jul 2025 - Update 4'
    >>> log.thermo
       run  Step  Temp  c_rg[1]
    0    0     0   1.0      NaN
    1    0   100   0.5      2.5
    """
    lines = text.splitlines()
    blocks = _split_runs(lines)
    wall_time_s = _parse_wall_time(lines)
    thermo_tables = []
    stats_rows: list[dict[str, float | str]] = []
    for block in blocks:
        table = _parse_thermo(block)
        table.insert(0, "run", block.index)
        thermo_tables.append(table)
        stats_rows.append({"run": block.index, **_parse_run_summary(block)})
    run_stats = pd.DataFrame(stats_rows)
    # Fixed columns first, whichever run reported them first; then the
    # timing sections in log order; then dangerous_builds.
    leading = [col for col in _RUN_STATS_LEADING if col in run_stats]
    trailing = [col for col in ("dangerous_builds",) if col in run_stats]
    sections = [col for col in run_stats
                if col not in leading and col not in trailing]
    run_stats = run_stats[leading + sections + trailing]
    last_run = stats_rows[-1]
    wall_time = pd.DataFrame(
        {
            "wall_time_s": [wall_time_s],
            "wall_time_hr": [wall_time_s / 3600],
            "n_cores": [last_run["n_cores"]],
            "n_atoms": [last_run["n_atoms"]],
        }
    )
    return LammpsLog(
        path=path,
        lammps_version=_lammps_version(lines),
        thermo=pd.concat(thermo_tables, ignore_index=True),
        run_stats=run_stats,
        wall_time=wall_time,
    )


def _lammps_version(lines: Sequence[str]) -> str | None:
    """Return the version from the first line, ``LAMMPS (<version>)``.

    Parameters
    ----------
    lines : sequence of str
        The lines of the log.

    Returns
    -------
    str or None
        The text between the parentheses, or None if the log is empty or
        its first line has another form, as in a log opened by the ``log``
        command.
    """
    if not lines:
        return None
    match = _VERSION.fullmatch(lines[0])
    return match["version"] if match else None


def _split_runs(lines: Sequence[str]) -> list[_RunBlock]:
    """Split the lines of a log into runs.

    Parameters
    ----------
    lines : sequence of str
        The lines of the log, without line endings.

    Returns
    -------
    list of _RunBlock
        The runs in log order.

    Raises
    ------
    BrokenLogError
        If there is no run-start line, a run has no ``Loop time of`` line,
        or a line starts with ``ERROR``.

    Notes
    -----
    Time is :math:`O(L)` in the number :math:`L` of lines.
    """
    starts = []
    for number, line in enumerate(lines):
        if line.startswith(_ERROR):
            raise BrokenLogError(f"Line {number + 1} reports an error: {line}")
        if line.startswith(_RUN_START):
            starts.append(number)
    if not starts:
        raise BrokenLogError(
            "The log has no run: no line starts with "
            f"{' or '.join(map(repr, _RUN_START))}."
        )
    blocks = []
    bounds = starts + [len(lines)]
    for index, (start, stop) in enumerate(zip(bounds[:-1], bounds[1:])):
        end = next(
            (n for n in range(start + 1, stop)
             if lines[n].startswith(_RUN_END)),
            None,
        )
        if end is None:
            raise BrokenLogError(
                f"Run {index} (line {start + 1}) has no {_RUN_END!r} line; "
                "the log may be cut off."
            )
        blocks.append(
            _RunBlock(
                index=index,
                first_line=start + 1,
                thermo=tuple(lines[start + 1:end]),
                summary=tuple(lines[end:stop]),
            )
        )
    return blocks


def _detect_style(lines: Sequence[str]) -> _ThermoStyle:
    """Return the thermo style of a run's thermo lines.

    Parameters
    ----------
    lines : sequence of str
        The thermo lines of one run, as in :attr:`_RunBlock.thermo`.

    Returns
    -------
    _ThermoStyle
        ``MULTI`` if the first line that is neither blank nor a skipped
        message is a ``------------ Step`` line, ``YAML`` if it is ``---``,
        and ``ONE`` otherwise, including when there is no such line.
    """
    for line in lines:
        if (not line.strip() or line.startswith(_SKIPPED_MESSAGES)
                or line.startswith(_MORE_INFO)):
            continue
        if _MULTI_HEADER.fullmatch(line):
            return _ThermoStyle.MULTI
        if line.strip() == "---":
            return _ThermoStyle.YAML
        return _ThermoStyle.ONE
    return _ThermoStyle.ONE


def _parse_thermo(block: _RunBlock) -> pd.DataFrame:
    """Parse the thermo table of one run.

    Parameters
    ----------
    block : _RunBlock
        One run.

    Returns
    -------
    pandas.DataFrame
        One row per output step and one column per keyword, named as LAMMPS
        prints it; ``Step`` is int64 and the other columns float64. ``nan``
        and ``inf`` values are read as such. No ``run`` column.

    Raises
    ------
    BrokenLogError
        If the table has no header, or a line is neither data, blank, nor a
        skipped message; the message names its line number in the log.

    Notes
    -----
    Time is :math:`O(n)` in the number :math:`n` of thermo lines.
    """
    style = _detect_style(block.thermo)
    if style is _ThermoStyle.MULTI:
        table = _parse_multi(block)
    elif style is _ThermoStyle.YAML:
        table = _parse_yaml(block)
    else:
        table = _parse_one(block)
    if "Step" in table.columns:
        table["Step"] = table["Step"].astype("int64")
    return table


def _parse_one(block: _RunBlock) -> pd.DataFrame:
    """Parse a thermo table in the ``one`` or ``custom`` style."""
    numbered = block.numbered_thermo()
    header = next(numbered, None)
    if header is None:
        raise BrokenLogError(
            f"Run {block.index} (line {block.first_line}) has no thermo "
            "header."
        )
    keywords = header[1].split()
    rows = [_floats(number, line.split(), len(keywords), line)
            for number, line in numbered]
    return pd.DataFrame(rows, columns=keywords, dtype="float64")


def _parse_multi(block: _RunBlock) -> pd.DataFrame:
    """Parse a thermo table in the ``multi`` style."""
    records: list[dict[str, float]] = []
    for number, line in block.numbered_thermo():
        header = _MULTI_HEADER.fullmatch(line)
        if header:
            records.append({
                "Step": _float(number, header["step"], line),
                "CPU": _float(number, header["cpu"], line),
            })
        elif records and _MULTI_PAIRS.fullmatch(line):
            for pair in _MULTI_PAIR.finditer(line):
                records[-1][pair["key"]] = _float(number, pair["value"], line)
        else:
            raise _not_data(number, line)
    return pd.DataFrame.from_records(records).astype("float64")


def _parse_yaml(block: _RunBlock) -> pd.DataFrame:
    """Parse a thermo table in the ``yaml`` style."""
    keywords: list[str] | None = None
    rows: list[list[float]] = []
    for number, line in block.numbered_thermo():
        stripped = line.strip()
        if stripped in ("---", "data:", "..."):
            continue
        if stripped.startswith("keywords: [") and stripped.endswith("]"):
            keywords = [word.strip().strip("'\"")
                        for word in _yaml_items(stripped[len("keywords:"):])]
        elif (keywords is not None and stripped.startswith("- [")
                and stripped.endswith("]")):
            rows.append(_floats(number, _yaml_items(stripped[1:]),
                                len(keywords), line))
        else:
            raise _not_data(number, line)
    if keywords is None:
        raise BrokenLogError(
            f"Run {block.index} (line {block.first_line}) has no "
            "'keywords:' line."
        )
    return pd.DataFrame(rows, columns=keywords, dtype="float64")


def _yaml_items(flow_list: str) -> list[str]:
    """Return the items of a YAML flow list such as ``[0, 1.5, ]``."""
    inner = flow_list.strip()[1:-1]
    return [item.strip() for item in inner.split(",") if item.strip()]


def _floats(number: int, fields: Sequence[str], expected: int,
            line: str) -> list[float]:
    """Convert the fields of one data row, checking their count."""
    if len(fields) != expected:
        raise _not_data(number, line)
    return [_float(number, field, line) for field in fields]


def _float(number: int, field: str, line: str) -> float:
    """Convert one field, reporting its line on failure."""
    try:
        return float(field)
    except ValueError:
        raise _not_data(number, line) from None


def _not_data(number: int, line: str) -> BrokenLogError:
    """Return the error for a line of a thermo table that is not data."""
    return BrokenLogError(
        f"Line {number} is inside a thermo table but is not data: {line!r}"
    )


def _parse_run_summary(block: _RunBlock) -> dict[str, float | str]:
    """Parse the summary that ends one run.

    Parameters
    ----------
    block : _RunBlock
        One run.

    Returns
    -------
    dict of str to float or str
        The fields of one row of :attr:`LammpsLog.run_stats`, without
        ``run``; a field the summary lacks is absent from the dict. When a
        line occurs more than once, the first occurrence counts.

    Notes
    -----
    Time is :math:`O(n)` in the number :math:`n` of summary lines.
    """
    fields: dict[str, float | str] = {}
    lines = block.summary
    number = 0
    while number < len(lines):
        line = lines[number]
        number += 1
        if (match := _LOOP.match(line)) and "loop_time" not in fields:
            fields.update(
                loop_time=float(match["time"]),
                n_cores=int(match["procs"]),
                loop_timesteps=int(match["steps"]),
                n_atoms=int(match["atoms"]),
            )
        elif (match := _PERFORMANCE.match(line)) and \
                "perf_per_day" not in fields:
            fields.update(
                perf_per_day=float(match["per_day"]),
                perf_unit=match["unit"],
                timesteps_per_s=float(match["steps_per_s"]),
            )
        elif (match := _CPU.match(line)) and "cpu_use" not in fields:
            threads = match["threads"]
            fields.update(
                cpu_use=float(match["cpu"]),
                mpi_task=int(match["tasks"]),
                openmp_threads=0 if threads == "no" else int(threads),
            )
        elif (match := _DANGEROUS.match(line)) and \
                "dangerous_builds" not in fields:
            fields["dangerous_builds"] = int(match["count"])
        elif line.startswith(_TIMING_TABLE) and number + 2 <= len(lines):
            # The column header, a line of dashes, then one row per section
            # until a blank line.
            columns = [cell.strip() for cell in lines[number].split("|")]
            number += 2
            while number < len(lines) and lines[number].strip():
                cells = [cell.strip() for cell in lines[number].split("|")]
                number += 1
                if len(cells) != len(columns):
                    continue
                row = dict(zip(columns, cells))
                section = row.get("Section")
                if not section:
                    continue
                fields.setdefault(f"{section}_avg_time",
                                  _float_or_nan(row.get("avg time", "")))
                fields.setdefault(f"{section}_total_pct",
                                  _float_or_nan(row.get("%total", "")))
    return fields


def _float_or_nan(cell: str) -> float:
    """Return a timing-table cell as a float, NaN if blank or not numeric."""
    try:
        return float(cell)
    except ValueError:
        return math.nan


def _parse_wall_time(lines: Sequence[str]) -> int:
    """Return the total wall time, in seconds.

    Parameters
    ----------
    lines : sequence of str
        The lines of the log.

    Returns
    -------
    int
        Seconds, from the line ``Total wall time: H:MM:SS``.

    Raises
    ------
    BrokenLogError
        If no line starts with ``Total wall time:``.
    """
    for line in reversed(lines):
        if line.startswith(_WALL_TIME) and (match := _WALL.match(line)):
            return (int(match["hours"]) * 3600 + int(match["minutes"]) * 60
                    + int(match["seconds"]))
    raise BrokenLogError(f"The log has no {_WALL_TIME!r} line.")
