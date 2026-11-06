"""
Time the organizer stages on one space, for amirhs1/PolyPhys#69, step 1.

Runs ``whole_from_segments`` -> ``ensemble`` -> ``ensemble_avg`` serially on
the probe files of one property in one space, as
:func:`polyphys.analyze.pipeline.reduce_space` calls them, and prints a
summary of timings, peak memory, and counts. The summary holds no file
names, no lineage names, and no array contents, so it can be pasted into
the issue.

Nothing is written: every stage runs with ``save_to=None``, so the time to
save the artifacts is not measured.

What is timed
=============
Each stage is timed end to end through its public function. Where the
stage interleaves its steps in one loop, the steps are also timed by
replaying the same calls outside the function:

``whole_from_segments``
    One loop parses each segment name and ``np.load``s the file, so load
    and parse cannot be split inside the function without changing it. The
    replay times a parse-only pass (the parser call and the dictionary
    grouping), a load-only pass (``np.load`` of every file), and a
    reduce-only pass (``np.sum``, ``np.concatenate``, or the bin-edge
    equality check on the loaded arrays). ``residual`` is the end-to-end
    median minus the sum of the three replay medians: the cost of the
    function's own bookkeeping, or of interactions between the steps (such
    as cache effects) that the replay does not reproduce.
``ensemble``
    No load. The replay times grouping (one parser call per whole name)
    and reducing (the public ``ens_from_*`` function on each ensemble).
    Computing ``bin_center`` from the edges is in the end-to-end time only.
``ensemble_avg``
    No load and no grouping: the end-to-end time is the reduce time.

Each timed call runs ``--repeats`` times. The first run of the load step
reads files the operating system may not have cached yet, so the summary
gives it separately (``first``) besides the median, minimum, and maximum.

Memory
======
Peak memory is measured in a separate, untimed run of each stage, because
:mod:`tracemalloc` slows every allocation:

- ``tracemalloc_peak_bytes``: the peak of Python and NumPy allocations
  during the stage, above what was allocated when it started. NumPy
  reports its array buffers to :mod:`tracemalloc`; pandas' NumPy-backed
  columns are counted through NumPy.
- ``ru_maxrss_bytes``: the process's peak resident set size from
  :func:`resource.getrusage`, read after the stage. It never decreases, so
  it is the peak of the whole process up to that stage, imports included.
  Linux reports it in KiB and macOS in bytes; the script converts both to
  bytes.

Usage
=====
Run from a checkout or an environment where ``polyphys`` imports::

    python benchmarks/organizer_timing.py <probe_dir> \\
        --parser TransFociCub --group bug --prop gyrTMon \\
        --relation tseries --repeats 5

``<probe_dir>`` is the probe directory of one space, as
:func:`~polyphys.analyze.pipeline.reduce_space` takes it; ``--hierarchy``
is the glob under it (default ``*/*``, the research code's ``dir/file``
layout). Files are selected as ``reduce_space`` selects them: names ending
in ``.<prop>.npy`` (segments) or ``-<prop>.npy`` (wholes), in natural
order. Run ``--help`` for every option.

Complexity: the script holds one property's segments, wholes, and
ensembles at once, as the organizer does; time is ``repeats`` times the
stages' own time, plus one untimed run per stage for memory.
"""

from __future__ import annotations

import argparse
import gc
import glob
import json
import os
import platform
import re
import resource
import statistics
import sys
import time
import tracemalloc
from collections import Counter
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd

import polyphys
from polyphys.manage import organizer
from polyphys.manage import parser as parser_module

#: The ensemble_avg type each whole type maps to, as reduce_space maps it.
ENS_TYPES = {'vector': 'dataframe', 'matrix': 'ndarray',
             'bin_edge': 'bin_edge'}

#: The ens_from_* reducer each whole type uses inside ensemble().
ENS_REDUCERS: dict[str, Callable[..., Any]] = {
    'vector': organizer.ens_from_vec,
    'matrix': organizer.ens_from_mat_t,
    'bin_edge': organizer.ens_from_bin_edge,
}


def natural_key(path: str) -> list[str | int]:
    """Sort key that orders digit runs as numbers, so j2 precedes j10."""
    return [int(token) if token.isdigit() else token
            for token in re.split(r'(\d+)', os.path.basename(path))]


def select(observations: Sequence[str], prop: str) -> list[str]:
    """The ``.npy`` files of `prop`, in natural order, as reduce_space
    selects them."""
    suffixes = (f'-{prop}.npy', f'.{prop}.npy')
    return sorted((p for p in observations if p.endswith(suffixes)),
                  key=natural_key)


def summarize(samples: Sequence[float]) -> dict[str, float]:
    """Median, minimum, maximum, and first sample of `samples`, seconds."""
    return {
        'median_s': statistics.median(samples),
        'min_s': min(samples),
        'max_s': max(samples),
        'first_s': samples[0],
        'n': len(samples),
    }


def repeat(func: Callable[[], Any], repeats: int) -> tuple[list[float], Any]:
    """Run `func` `repeats` times; return the wall times and the last
    result."""
    times: list[float] = []
    result: Any = None
    for _ in range(repeats):
        result = None  # drop the previous result before the next run
        gc.collect()
        start = time.perf_counter()
        result = func()
        times.append(time.perf_counter() - start)
    return times, result


def maxrss_bytes() -> int:
    """The process's peak resident set size so far, in bytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak) if sys.platform == 'darwin' else int(peak) * 1024


def measure_memory(func: Callable[[], Any]) -> dict[str, int]:
    """Peak memory of one untimed run of `func`."""
    gc.collect()
    tracemalloc.start()
    base, _ = tracemalloc.get_traced_memory()
    result = func()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    del result
    return {'tracemalloc_peak_bytes': peak - base,
            'ru_maxrss_bytes': maxrss_bytes()}


def describe_arrays(arrays: Sequence[Any]) -> dict[str, Any]:
    """Counts of shapes and dtypes, and the total size, of `arrays`.

    DataFrames are described by their shape and their column dtypes.
    """
    shapes: Counter[str] = Counter()
    dtypes: Counter[str] = Counter()
    nbytes = 0
    for item in arrays:
        if isinstance(item, pd.DataFrame):
            shapes[str(item.shape)] += 1
            for dtype in item.dtypes:
                dtypes[str(dtype)] += 1
            nbytes += int(item.memory_usage(index=True, deep=False).sum())
        elif isinstance(item, dict):  # ens_avg_from_ndarray's result
            for value in item.values():
                shapes[str(np.shape(value))] += 1
                dtypes[str(np.asarray(value).dtype)] += 1
                nbytes += int(np.asarray(value).nbytes)
        else:
            array = np.asarray(item)
            shapes[str(array.shape)] += 1
            dtypes[str(array.dtype)] += 1
            nbytes += int(array.nbytes)
    return {'shapes': dict(shapes), 'dtypes': dict(dtypes),
            'total_bytes': nbytes}


def replay_whole(paths: Sequence[str], parser: type, group: str,
                 relation: str, repeats: int) -> dict[str, Any]:
    """Time the parse, load, and reduce steps of whole_from_segments
    separately, by repeating its calls outside it."""

    def parse() -> dict[str, list[str]]:
        grouped: dict[str, list[str]] = {}
        for path in paths:
            info = parser(path, 'segment', group, ispath=True)
            grouped.setdefault(getattr(info, 'whole'), []).append(path)
            getattr(info, 'segment_id')
        return grouped

    parse_times, grouped = repeat(parse, repeats)

    def load() -> list[np.ndarray]:
        return [np.load(path) for path in paths]

    load_times, _ = repeat(load, repeats)

    loaded = {name: [np.load(p) for p in members]
              for name, members in grouped.items()}

    def reduce() -> dict[str, np.ndarray]:
        if relation == 'histogram':
            return {n: np.sum(a, axis=0) for n, a in loaded.items()}
        if relation == 'tseries':
            return {n: np.concatenate(a, axis=0) for n, a in loaded.items()}
        out = {}
        for n, arrs in loaded.items():
            if not all(np.array_equal(arrs[0], a) for a in arrs[1:]):
                raise ValueError('Bin edges differ within a whole.')
            out[n] = arrs[0].copy()
        return out

    reduce_times, _ = repeat(reduce, repeats)
    segments_per_whole = Counter(len(m) for m in grouped.values())
    return {
        'parse_and_group': summarize(parse_times),
        'load': summarize(load_times),
        'reduce': summarize(reduce_times),
        'segments_per_whole': {str(k): v for k, v in
                               sorted(segments_per_whole.items())},
        'segment_arrays': describe_arrays(
            [a for arrs in loaded.values() for a in arrs]),
    }


def replay_ensemble(wholes: dict[str, Any], parser: type, group: str,
                    whole_type: str, repeats: int) -> dict[str, Any]:
    """Time the grouping and reducing steps of ensemble separately."""

    def group_wholes() -> dict[str, dict[str, Any]]:
        grouped: dict[str, dict[str, Any]] = {}
        for name, data in wholes.items():
            ens = getattr(parser(name, 'whole', group), 'ensemble_long')
            grouped.setdefault(ens, {})[name] = data
        return grouped

    group_times, grouped = repeat(group_wholes, repeats)
    reducer = ENS_REDUCERS[whole_type]
    reduce_times, _ = repeat(
        lambda: dict(map(reducer, grouped.items())), repeats)
    wholes_per_ens = Counter(len(m) for m in grouped.values())
    return {
        'parse_and_group': summarize(group_times),
        'reduce': summarize(reduce_times),
        'wholes_per_ensemble': {str(k): v for k, v in
                                sorted(wholes_per_ens.items())},
    }


def residual(end_to_end: dict[str, float],
             parts: Sequence[dict[str, float]]) -> float:
    """End-to-end median minus the sum of the parts' medians, seconds."""
    return end_to_end['median_s'] - sum(p['median_s'] for p in parts)


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Time the three stages and return the summary."""
    parser = getattr(parser_module, args.parser)
    pattern = os.path.join(args.probe_dir, args.hierarchy)
    observations = sorted(glob.glob(pattern))
    paths = select(observations, args.prop)
    if not paths:
        raise SystemExit(
            f"No '.{args.prop}.npy' or '-{args.prop}.npy' file under the "
            "given directory and hierarchy.")
    edge_paths = select(observations, args.edge_prop) if args.edge_prop \
        else []
    segment_args = [(p,) for p in paths]
    edge_args = [(p,) for p in edge_paths]
    summary: dict[str, Any] = {
        'issue': 'amirhs1/PolyPhys#69 step 1',
        'environment': {
            'python': platform.python_version(),
            'implementation': platform.python_implementation(),
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'polyphys': polyphys.__version__,
            'system': platform.system(),
            'machine': platform.machine(),
            'cpu_count': os.cpu_count(),
            'sched_affinity': (len(os.sched_getaffinity(0))
                               if hasattr(os, 'sched_getaffinity')
                               else None),
        },
        'arguments': {
            'parser': args.parser, 'group': args.group, 'prop': args.prop,
            'relation': args.relation, 'lineage': args.lineage,
            'edge_prop': args.edge_prop, 'repeats': args.repeats,
        },
        'counts': {'files_matched_by_glob': len(observations),
                   'property_files': len(paths),
                   'edge_files': len(edge_paths)},
        'stages': {},
    }
    stages: dict[str, dict[str, Any]] = summary['stages']
    stage: dict[str, Any]

    # Stage 1: segment -> whole (or whole files -> wholes)
    if args.lineage == 'segment':
        def stage1() -> dict[str, np.ndarray]:
            return organizer.whole_from_segments(
                args.prop, segment_args, parser, args.group, args.relation)
        e2e_times, wholes = repeat(stage1, args.repeats)
        stage = {'function': 'whole_from_segments',
                 'end_to_end': summarize(e2e_times)}
        stage.update(replay_whole(paths, parser, args.group, args.relation,
                                  args.repeats))
        stage['residual_median_s'] = residual(
            stage['end_to_end'],
            [stage['parse_and_group'], stage['load'], stage['reduce']])
    else:
        def stage1() -> dict[str, np.ndarray]:
            return organizer.whole_from_file(segment_args, parser,
                                             args.group)
        e2e_times, wholes = repeat(stage1, args.repeats)
        stage = {'function': 'whole_from_file',
                 'end_to_end': summarize(e2e_times),
                 'note': 'parse and load are one loop; not replayed'}
    stage.update(measure_memory(stage1))
    stage['wholes'] = len(wholes)
    stage['whole_arrays'] = describe_arrays(list(wholes.values()))
    stages['segment_to_whole'] = stage

    whole_edges = None
    if edge_paths:
        if args.lineage == 'segment':
            whole_edges = organizer.whole_from_segments(
                args.edge_prop, edge_args, parser, args.group, 'bin_edge')
        else:
            whole_edges = organizer.whole_from_file(edge_args, parser,
                                                    args.group)

    # Stage 2: whole -> ensemble
    whole_type = args.whole_type
    if whole_type == 'auto':
        if args.relation == 'bin_edge':
            whole_type = 'bin_edge'
        else:
            ndim = np.asarray(next(iter(wholes.values()))).ndim
            whole_type = 'vector' if ndim == 1 else 'matrix'
    summary['arguments']['whole_type'] = whole_type

    def stage2() -> dict[str, Any]:
        return organizer.ensemble(args.prop, wholes, parser, args.group,
                                  whole_type, whole_edges=whole_edges)

    e2e_times, ensembles = repeat(stage2, args.repeats)
    stage = {'function': 'ensemble', 'end_to_end': summarize(e2e_times)}
    stage.update(replay_ensemble(wholes, parser, args.group, whole_type,
                                 args.repeats))
    stage['residual_median_s'] = residual(
        stage['end_to_end'], [stage['parse_and_group'], stage['reduce']])
    stage.update(measure_memory(stage2))
    stage['ensembles'] = len(ensembles)
    stage['ensemble_arrays'] = describe_arrays(list(ensembles.values()))
    stages['whole_to_ensemble'] = stage

    # Stage 3: ensemble -> ensAvg
    ens_type = ENS_TYPES[whole_type]
    summary['arguments']['ens_type'] = ens_type

    def stage3() -> dict[str, Any]:
        return organizer.ensemble_avg(args.prop, ensembles, args.group,
                                      ens_type)  # type: ignore[arg-type]

    e2e_times, averages = repeat(stage3, args.repeats)
    stage = {'function': 'ensemble_avg', 'end_to_end': summarize(e2e_times),
             'note': 'no load or grouping step: end_to_end is the reduce'}
    stage.update(measure_memory(stage3))
    stage['averages'] = len(averages)
    stage['average_arrays'] = describe_arrays(list(averages.values()))
    stages['ensemble_to_ensAvg'] = stage

    total = sum(s['end_to_end']['median_s'] for s in stages.values())
    summary['total_median_s'] = total
    summary['share_of_total'] = {
        name: (s['end_to_end']['median_s'] / total if total else None)
        for name, s in stages.items()}
    return summary


def format_seconds(seconds: float) -> str:
    """Seconds in a readable unit."""
    if abs(seconds) < 1e-3:
        return f'{seconds * 1e6:.1f} us'
    if abs(seconds) < 1:
        return f'{seconds * 1e3:.2f} ms'
    return f'{seconds:.3f} s'


def format_bytes(n: int) -> str:
    """Bytes in a readable binary unit."""
    size = float(n)
    for unit in ('B', 'KiB', 'MiB', 'GiB'):
        if abs(size) < 1024 or unit == 'GiB':
            return f'{size:.1f} {unit}'
        size /= 1024
    return f'{size:.1f} GiB'


def to_text(summary: dict[str, Any]) -> str:
    """A plain-text rendering of `summary`."""
    lines = [f"PolyPhys organizer timing ({summary['issue']})"]
    env = summary['environment']
    lines.append(
        f"python {env['python']}, numpy {env['numpy']}, pandas "
        f"{env['pandas']}, polyphys {env['polyphys']}, {env['system']} "
        f"{env['machine']}, cpus {env['cpu_count']}"
        + (f" (affinity {env['sched_affinity']})"
           if env['sched_affinity'] is not None else ''))
    a = summary['arguments']
    lines.append(
        f"parser {a['parser']}, group {a['group']}, prop {a['prop']}, "
        f"relation {a['relation']}, lineage {a['lineage']}, whole_type "
        f"{a['whole_type']}, ens_type {a['ens_type']}, edge_prop "
        f"{a['edge_prop']}, repeats {a['repeats']}")
    c = summary['counts']
    lines.append(f"files: glob {c['files_matched_by_glob']}, property "
                 f"{c['property_files']}, edges {c['edge_files']}")
    lines.append('')
    header = f"{'step':34} {'median':>11} {'min':>11} {'max':>11} " \
             f"{'first':>11}"
    for name, stage in summary['stages'].items():
        lines.append(f"[{name}] {stage['function']}")
        lines.append(header)
        for step in ('end_to_end', 'parse_and_group', 'load', 'reduce'):
            if step in stage:
                t = stage[step]
                lines.append(
                    f"  {step:32} {format_seconds(t['median_s']):>11} "
                    f"{format_seconds(t['min_s']):>11} "
                    f"{format_seconds(t['max_s']):>11} "
                    f"{format_seconds(t['first_s']):>11}")
        if 'residual_median_s' in stage:
            lines.append(f"  {'residual (e2e - parts)':32} "
                         f"{format_seconds(stage['residual_median_s']):>11}")
        lines.append(
            f"  peak memory: tracemalloc "
            f"{format_bytes(stage['tracemalloc_peak_bytes'])}, process "
            f"maxrss {format_bytes(stage['ru_maxrss_bytes'])}")
        for key in ('segments_per_whole', 'wholes_per_ensemble'):
            if key in stage:
                lines.append(f"  {key}: {stage[key]}")
        for key in ('wholes', 'ensembles', 'averages'):
            if key in stage:
                lines.append(f"  {key}: {stage[key]}")
        for key in ('segment_arrays', 'whole_arrays', 'ensemble_arrays',
                    'average_arrays'):
            if key in stage:
                d = stage[key]
                lines.append(
                    f"  {key}: shapes {d['shapes']}, dtypes {d['dtypes']}, "
                    f"total {format_bytes(d['total_bytes'])}")
        if 'note' in stage:
            lines.append(f"  note: {stage['note']}")
        lines.append('')
    lines.append(f"total of stage medians: "
                 f"{format_seconds(summary['total_median_s'])}")
    shares = ', '.join(f"{k} {v:.0%}" for k, v in
                       summary['share_of_total'].items() if v is not None)
    lines.append(f"share of total: {shares}")
    return '\n'.join(lines)


def main(argv: Sequence[str] | None = None) -> None:
    """Parse the command line, run the timing, and print the summary."""
    cli = argparse.ArgumentParser(
        description='Time whole_from_segments -> ensemble -> ensemble_avg '
                    'on one property of one space (PolyPhys #69, step 1). '
                    'Prints timings and counts only; writes nothing.')
    cli.add_argument('probe_dir', help='the probe directory of one space')
    cli.add_argument('--parser', required=True,
                     help='parser class in polyphys.manage.parser, such as '
                          'TransFociCub')
    cli.add_argument('--group', required=True, help="particle group, 'bug'")
    cli.add_argument('--prop', required=True,
                     help="property in the file names, such as 'gyrTMon'")
    cli.add_argument('--relation', required=True,
                     choices=('tseries', 'histogram', 'bin_edge'),
                     help='how segments join into a whole')
    cli.add_argument('--lineage', default='segment',
                     choices=('segment', 'whole'),
                     help="lineage of the files (default 'segment')")
    cli.add_argument('--whole-type', default='auto',
                     choices=('auto', 'vector', 'matrix', 'bin_edge'),
                     help="ensemble() whole_type; 'auto' takes 'vector' for "
                          "1D wholes, 'matrix' for 2D, 'bin_edge' for that "
                          "relation")
    cli.add_argument('--edge-prop', default=None,
                     help="bin-edge property of a histogram, such as "
                          "'rEdgeMon'; passed to ensemble() as whole_edges")
    cli.add_argument('--hierarchy', default='*/*',
                     help="glob under probe_dir (default '*/*')")
    cli.add_argument('--repeats', type=int, default=5,
                     help='timed runs per step (default 5)')
    cli.add_argument('--format', default='text', choices=('text', 'json'),
                     help="output format (default 'text')")
    args = cli.parse_args(argv)
    if args.repeats < 1:
        cli.error('--repeats must be at least 1')
    summary = run(args)
    if args.format == 'json':
        print(json.dumps(summary, indent=2))
    else:
        print(to_text(summary))


if __name__ == '__main__':
    main()
