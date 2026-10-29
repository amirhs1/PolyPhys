"""\
==========================================================
:mod:`polyphys.manage.organizer`
==========================================================
The :mod:`polyphys.manage.organizer` module provides tools for organizing and
aggregating simulation data based on filename patterns, file structures, and
processing stages. It facilitates handling segmented data, grouping experiments
into ensembles, consolidating ensembles into spaces, and organizing spaces into
galaxies for high-level analysis.

The module is tailored for use with LAMMPS and similar simulation tools, where
it automates data consolidation across hierarchical levels. It is particularly
useful for organizing data generated in high-throughput simulations, enabling
efficient post-processing, ensemble averaging, and preparation of data for
statistical and comparative analysis.

Terminology
===========
The following terms are used throughout the module to structure and organize
simulation output data:

- **Experiment**: A complete simulation run, typically using tools like LAMMPS,
  representing a specific set of conditions.

- **Lineage**: A hierarchical structure for organizing data across different
  levels of aggregation:

    - **Segment**: A portion of data from a long simulation.
    - **Whole** (Instance): A complete dataset made by combining multiple
      segments.
    - **Ensemble**: A set of wholes with identical macroscopic properties but
      varied initial conditions (e.g., different random seeds).
    - **Space**: A group of ensembles with identical attributes except for one,
      which systematically varies across the space.
    - **Galaxy**: The highest level of aggregation, representing a collection
      of spaces with differing attributes.

- **Attributes**: Macroscopic characteristics such as box size or polymer
  topology that remain constant within a space.

- **Properties**: Physical features measured during the simulation, like
  density or radius of gyration, with optional associated measures (e.g.,
  autocorrelation).

- **Species**: The particle or molecular entity type, such as 'Mon' for
  monomers, each with specific behaviors or attributes as needed.

- **Phases**: Directories organized by processing phase, such as `logs`,
  `trjs` (trajectories), `probe` (probing data), `analysis` (post-processing),
  and `viz` (visualization-ready data).

- **Stages**: Aggregation levels within each phase, e.g., `wholeSim`, `ens`,
  or `ensAvg`.

Naming Conventions
==================
File and directory naming conventions in the `organizer` module provide
structure to data organization across phases and lineage levels. Each name
pattern encodes metadata about lineage, phase, group, property, and processing
stage. Example patterns include:

- `whole|segment.group.run.lammpstrj` for trajectory files.
- `lineage-phase-group-prop_stage.extension` for processed files.

Dependencies
============
- `numpy`: For numerical operations on arrays, such as summing or concatenating
  segments.
- `pandas`: For handling tabular data and DataFrames used in data aggregation.

Contents
========

Naming and persistence:
:func:`create_fullname`, :func:`save_artifact`, and :func:`make_database`.

Building wholes from their parts:
:func:`whole_from_segments`, :func:`whole_from_file`, and
:func:`whole_from_dist_mat_t`.

Grouping wholes into ensembles:
:func:`ens_from_bin_edge`, :func:`ens_from_vec`, :func:`ens_from_mat_t`,
:func:`ens_from_df`, and :func:`ensemble`.

Averaging over an ensemble:
:func:`ens_avg_from_bin_edge`, :func:`ens_avg_from_ndarray`,
:func:`ens_avg_from_df`, and :func:`ensemble_avg`.

Space-level collection and bookkeeping:
:func:`children_stamps`, :func:`parents_stamps`,
:func:`find_unique_properties`, :func:`space_tseries`, and
:func:`space_hists`.

Usage
=====
The :mod:`organizer` module streamlines organization and aggregation of
simulation data. Users can assemble segmented data, organize data into
ensembles, and generate ensemble-averaged data at various aggregation stages.
Its multi-level organization supports comparative analysis across varied
conditions within spaces and galaxies.

Examples
========
Two replicates (``ens1`` and ``ens2``) of one state point, each probed in
two segments (``j01`` and ``j02``), reduce to one ensemble-averaged time
series. Segment ``s`` of replicate ``e`` holds two frames of the value
``10 * e + s``:

>>> import tempfile
>>> from polyphys.manage.parser import SumRuleCyl
>>> tmp = tempfile.TemporaryDirectory()
>>> root = Path(tmp.name)
>>> state = 'N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000'
>>> paths = []
>>> for e in (1, 2):
...     for s in (1, 2):
...         path = root / f'{state}ens{e}.j0{s}.bug.gyrTMon.npy'
...         np.save(path, np.full(2, 10.0 * e + s))
...         paths.append(path)

Join each replicate's segments into a whole, collect the wholes of the
state point into an ensemble, and average over the replicates:

>>> wholes = whole_from_segments('gyrTMon', paths, SumRuleCyl, 'bug',
...                              'tseries')
>>> wholes[state + 'ens1']
array([11., 11., 12., 12.])
>>> ensembles = ensemble('gyrTMon', wholes, SumRuleCyl, 'bug', 'vector')
>>> averages = ensemble_avg('gyrTMon', ensembles, 'bug', 'dataframe')
>>> averages[state]
   gyrTMon-mean  gyrTMon-var  gyrTMon-sem
0          16.0         50.0          5.0
1          16.0         50.0          5.0
2          17.0         50.0          5.0
3          17.0         50.0          5.0
>>> tmp.cleanup()

Notes
=====
- Filenames must follow consistent patterns for the `organizer` module to
  accurately parse and organize data.
- This module supports only cylindrical, slit, or cubic geometries.
- Functions in this module rely on consistent naming conventions for groups,
  properties, species, and stages.
"""

import warnings
import logging
from typing import Any, Callable, Literal, Sequence, Iterator, TYPE_CHECKING
import os
import pathlib
from pathlib import Path
from glob import glob
import numpy as np
import pandas as pd
from .types import (
    Pathish,
    PathGroup,
    AnyPathIter,
    GroupT,
    PropertyT,
    ParserInstanceT,
    ParserClassT,
    StageT,
    PhaseT,
    WholeRelationT,
    PrimitiveLineageT,
)
from .utils import invalid_keyword, sort_filenames, round_to_multiple

if TYPE_CHECKING:
    # Unused here, but sphinx_autodoc_typehints runs this block to resolve
    # the "ParserBase" forward reference inside ParserClassT.
    from .parser import ParserBase  # noqa: F401

logger = logging.getLogger(__name__)


def create_fullname(name: str, group: GroupT, prop: PropertyT) -> str:
    """
    Create a structured filename based on a base name, particle group,
    and property.

    Parameters
    ----------
    name : str
        The base name.
    group : :py:data:`GroupT`
        Particle group (e.g., 'bug' or 'nucleoid').
    prop : :py:data:`PropertyT`
        Physical property name (e.g., 'density').

    Returns
    -------
    str
        A string that combines `name`, `group`, and `prop` in a
        hyphen-separated format.

    Examples
    --------
    >>> create_fullname('sample', 'bug', 'density')
    'sample-bug-density'
    """
    return "-".join([name, group, prop])


def save_artifact(
    filename: str,
    data: np.ndarray | pd.DataFrame | dict[str, np.ndarray],
    save_to: Pathish,
) -> None:
    """
    Save the `data` to a specified file format, allowing structured file
    naming.

    Parameters
    ----------
    filename : str
        The output file name.
    data : np.ndarray | pd.DataFrame | dict[str, np.ndarray]
        Data to be saved. Accepted types:

        - `np.ndarray`: Saved as ``<filename>.npy``.
        - `pd.DataFrame`: Saved as ``<filename>.csv``, without the index.
        - `dict` of `np.ndarray`: Each key must have the form
          ``<prop>-<measure>``, with exactly one hyphen, and each entry is
          saved as ``<filename>-<measure>.npy``.

    save_to : Pathish
        Path to an existing directory where the file is saved.

    Raises
    ------
    TypeError
        If `data` type does not match any of the supported formats.
    ValueError
        If `data` is a dict and one of its keys does not contain exactly one
        hyphen.

    Examples
    --------
    >>> import tempfile
    >>> tmp = tempfile.TemporaryDirectory()

    Save a numpy array to a ``.npy`` file:

    >>> save_artifact('w-bug-rhoMon', np.array([1.0, 2.0]), tmp.name)
    >>> np.load(Path(tmp.name) / 'w-bug-rhoMon.npy')
    array([1., 2.])

    Save a DataFrame to a ``.csv`` file:

    >>> save_artifact('w-bug-rhoMon', pd.DataFrame({'A': [1, 2]}), tmp.name)
    >>> pd.read_csv(Path(tmp.name) / 'w-bug-rhoMon.csv')
       A
    0  1
    1  2

    Save each entry of a dict to its own file, named by its measure:

    >>> save_artifact('w-bug-rhoMon', {'rhoMon-mean': np.array([1.5])},
    ...               tmp.name)
    >>> sorted(path.name for path in Path(tmp.name).iterdir())
    ['w-bug-rhoMon-mean.npy', 'w-bug-rhoMon.csv', 'w-bug-rhoMon.npy']
    >>> tmp.cleanup()
    """
    file_path = os.path.join(save_to, filename)

    # Save data based on type
    if isinstance(data, pd.DataFrame):
        data.to_csv(f"{file_path}.csv", index=False)
    elif isinstance(data, np.ndarray):
        np.save(f"{file_path}.npy", data)
    elif isinstance(data, dict):
        for prop_key, prop_data in data.items():
            _, prop_measure = prop_key.split("-")
            np.save(f"{file_path}-{prop_measure}.npy", prop_data)
    else:
        raise TypeError(
            f"Unsupported data type {type(data).__name__}."
            + "Expected pd.DataFrame, np.ndarray, or dict of np.ndarray."
        )


def make_database(
    old_database: Pathish,
    phase: PhaseT | str,
    stage: StageT | str,
    group: GroupT | None,
) -> str:
    """
    Create a new directory path based on the provided `old_database` path and
    specified parameters (`phase`, `group`, and `stage`). If the directory does
    not already exist, it will be created.

    The `old_database` is expected to follow a structured naming convention:
    `prefix-old_phase-old_group-old_stage`, where each part represents an
    aspect of the directory's content. This base structure helps generate the
    new path.

    The newly constructed directory name follows:
        `prefix-phase-group-stage`

    Parameters
    ----------
    old_database : Pathish
        Path to the reference directory whose structure will be used as a
        base. It need not exist.
    phase : {'simsAll', 'simsCont', 'logs', 'trjs', 'probes', 'analysis', \
'viz', 'galaxy', 'allInOne'}
        The new phase name for the directory, specifying its role in
        the workflow.
    stage : str
        The stage of processing or data type represented in the new
        directory, such as 'segment', 'wholeSim', 'ens', 'ensAvg', 'space',
        or 'galaxy'. It is not checked.
    group : GroupT or None
        Particle group related to the data in the directory; left out of the
        name when None.

    Returns
    -------
    str
        The path to the new directory, ending with ``/``. If the directory
        already exists, the same path is returned and an INFO record naming
        it is logged to the module's logger.

    Raises
    ------
    ValueError
        If `phase` is not one of the accepted values.

    Examples
    --------
    The new directory sits next to the old phase directory:

    - Input: ``root/parent1/.../parentN/old_phase/old_dir``
    - Result: ``root/parent1/.../parentN/phase/new_dir``

    >>> import tempfile
    >>> tmp = tempfile.TemporaryDirectory()
    >>> old = Path(tmp.name) / 'probes' / 'N20-bug-segment'
    >>> new = make_database(old, 'analysis', 'ens', 'bug')
    >>> Path(new).relative_to(tmp.name).as_posix()
    'analysis/N20-bug-ens'
    >>> Path(new).is_dir()
    True
    >>> tmp.cleanup()
    """
    invalid_keyword(
        phase,
        ["simsAll", "simsCont", "logs", "trjs", "probes", "analysis", "viz",
         "galaxy", "allInOne"],
    )

    old_path = pathlib.Path(old_database)

    # Common prefix derived from old_database to create new directory name
    prefix = old_path.parts[-1].split("*")[0].split("-")[0]
    new_directory = "-".join([part for part in [prefix, group, stage] if part])

    # Construct the new database path in the same parent directory level as
    # the old one
    new_database_parts = list(old_path.parts[:-2]) + [phase, new_directory]
    new_database = pathlib.Path(*new_database_parts)

    # Avoid double slashes at the start of the path
    if str(new_database).startswith("//"):
        new_database = pathlib.Path(str(new_database)[1:])

    # Create the new directory if it doesn't already exist
    try:
        new_database.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        logger.info("Reusing the existing directory %s.", new_database)

    return str(new_database) + "/"


def _as_path(item: Pathish | PathGroup) -> Path:
    """Normalize a single item into a pathlib.Path.

    Accepts either:
      - a bare pathish (str, Path, os.PathLike)
      - a sequence of pathish values (we take the first element)
    """
    # Bare pathish?
    if isinstance(item, (str, os.PathLike)):
        return Path(item)

    # Sequence of pathish (e.g., tuple[Path, ...] from sort_filenames)
    if isinstance(item, Sequence) and item:
        first = item[0]
        if isinstance(first, (str, os.PathLike)):
            return Path(first)

    raise TypeError(f"Unsupported path item: {type(item)!r} -> {item!r}")


def _normalize_paths(paths: AnyPathIter) -> Iterator[Path]:
    """Yield normalized Path objects from either a flat or grouped iterable."""
    for item in paths:
        # item is either Pathish or PathGroup; _as_path handles both
        yield _as_path(item)


def _assert_equal_bin_edges(arrays: Sequence[np.ndarray]) -> np.ndarray:
    """
    Ensure all input arrays of bin edges are identical.

    Parameters
    ----------
    arrays : Sequence[np.ndarray]
        Sequence of 1D arrays representing bin edges from different wholes
        or segments.

    Returns
    -------
    np.ndarray
        A copy of the first array if all arrays are identical.

    Raises
    ------
    ValueError
        If any array differs from the first one.
    """
    ref = arrays[0]
    if not all(np.array_equal(ref, a) for a in arrays[1:]):
        raise ValueError("Bin edges are not identical across all arrays.")
    return ref.copy()


def whole_from_segments(
    prop: PropertyT,
    segments: AnyPathIter,
    parser: ParserClassT,
    group: GroupT,
    relation: WholeRelationT,
    save_to: Pathish | None = None,
) -> dict[str, np.ndarray]:
    """
    Generate 'whole' arrays by combining data segments for a specified
    property of a particle `group`. The combination method is determined by
    the `relation` parameter.

    Parameters
    ----------
    prop: :py:data:`PropertyT`
        The physical property to process (e.g., 'density').
    segments : :py:data:`AnyPathIter`
        An iterable of paths to segment ``.npy`` files for `prop`, either
        bare paths or 1-tuples of a path, as
        :func:`~polyphys.manage.utils.sort_filenames` returns them. Each path
        points to a single segment file. Segments of different wholes may be
        interleaved, but the segments of each whole must come in strictly
        increasing ``segment_id`` order (the ``j<NN>`` keyword of their
        names).
    parser : :py:class:`ParserClassT`
        A parser class from :mod:`polyphys.manage.parser` that interprets
        file paths and extracts attributes such as the `whole` name and the
        ``segment_id``.
    group : :py:data:`GroupT`
        Type of the particle group.
    relation : :py:data:`WholeRelationT`
        Specifies how to combine segments into a whole. Accepted values:

        - 'histogram'
            'segment' is an N-dimensional histogram-like file, so the segments
            of a 'whole' are summed element-wise. In the two-dimensional
            space, a histogram has a (x_nbins, y_nbins) shape.

        - 'tseries'
            'segment' is a time-series-like file, so the segments of a 'whole'
            are concatenated along axis 0, in ``segment_id`` order.

        - 'bin_edge'
            'segment' is a bin-edge file. All the segments of a whole must
            hold identical edges, and one copy of them is the whole's
            bin-edge array.

    save_to : str, optional
        Directory path where the output files are saved, if specified. Each
        whole is saved as ``<whole>-<group>-<prop>.npy``.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary with whole names as keys and the combined whole data as
        numpy arrays, in the order in which the wholes first appear in
        `segments`.

    Raises
    ------
    ValueError
        If `relation` is not one of the accepted values; if the segments of
        a whole do not arrive in strictly increasing ``segment_id`` order,
        which includes a segment passed twice; or, for 'bin_edge', if the
        segments of a whole hold different edges.
    TypeError
        If an item of `segments` is neither a path nor a sequence whose first
        element is a path.

    Notes
    -----
    The order of a whole's segments is part of the data's meaning for
    'tseries', and a duplicated segment would be counted twice by
    'histogram'. This function only checks the order; sort the paths with
    :func:`~polyphys.manage.utils.sort_filenames` first, which orders the
    ``j<NN>`` keyword naturally.

    Let n_segments be the number of segment files and m the number of
    values in one segment. Time is O(n_segments * m). Every segment of
    every whole is loaded and held in memory before any whole is reduced,
    so peak memory is O(n_segments * m), not the size of one whole; the
    returned wholes add O(n_wholes * m) for 'histogram' and 'bin_edge' and
    O(n_segments * m) for 'tseries', where n_wholes is the number of
    wholes.

    Examples
    --------
    Two segments of one whole, written as the probe phase names them, join
    into one time series in ``segment_id`` order:

    >>> import tempfile
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> tmp = tempfile.TemporaryDirectory()
    >>> whole = ('N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000'
    ...          'adump5000ens1')
    >>> paths = []
    >>> for seg_id in (1, 2):
    ...     path = Path(tmp.name) / f'{whole}.j0{seg_id}.bug.gyrTMon.npy'
    ...     np.save(path, np.array([10.0, 11.0]) * seg_id)
    ...     paths.append(path)
    >>> wholes = whole_from_segments('gyrTMon', paths, SumRuleCyl, 'bug',
    ...                              'tseries')
    >>> wholes[whole]
    array([10., 11., 20., 22.])

    'histogram' sums the same segments element-wise:

    >>> whole_from_segments('gyrTMon', paths, SumRuleCyl, 'bug',
    ...                     'histogram')[whole]
    array([30., 33.])

    Segments out of order raise rather than giving a reversed series:

    >>> whole_from_segments('gyrTMon', paths[::-1], SumRuleCyl, 'bug',
    ...                     'tseries')  # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ValueError: The segments of whole '...ens1' are not in ... [2, 1]...
    >>> tmp.cleanup()
    """
    invalid_keyword(relation, ["histogram", "tseries", "bin_edge"])
    # Combine segments into wholes based on segment names and relation type,
    # checking that each whole's segment_id values strictly increase
    grouped_segments: dict[str, list[np.ndarray]] = {}
    segment_ids: dict[str, list[int]] = {}
    for seg_path in _normalize_paths(segments):
        seg_info = parser(seg_path, "segment", group, ispath=True)
        whole_name = getattr(seg_info, "whole")
        seg_id = getattr(seg_info, "segment_id")
        whole_ids = segment_ids.setdefault(whole_name, [])
        if whole_ids and seg_id <= whole_ids[-1]:
            raise ValueError(
                f"The segments of whole '{whole_name}' are not in strictly "
                f"increasing segment_id order: {whole_ids + [seg_id]}. Sort "
                "them with sort_filenames first."
            )
        whole_ids.append(seg_id)
        grouped_segments.setdefault(whole_name, []).append(np.load(seg_path))

    # Define reduction strategies for each relation type
    reducers: dict[str, Callable[[Sequence[np.ndarray]], np.ndarray]] = {
        "histogram": lambda arrs: np.sum(arrs, axis=0),
        "tseries": lambda arrs: np.concatenate(arrs, axis=0),
        "bin_edge": _assert_equal_bin_edges,
    }
    reducer = reducers[relation]

    # Apply the appropriate mapping function to merge each group of segments
    # into a whole.
    wholes: dict[str, np.ndarray] = \
        {name: reducer(vals) for name, vals in grouped_segments.items()}

    if save_to is not None:
        for name, data in wholes.items():
            save_artifact(create_fullname(name, group, prop), data, save_to)

    return wholes


def whole_from_file(
    whole_paths: AnyPathIter,
    parser: ParserClassT,
    group: GroupT,
) -> dict[str, np.ndarray]:
    """
    Load data from 'whole' ``.npy`` files of a particle `group`. Each file
    path in `whole_paths` points to a single whole file.

    Parameters
    ----------
    whole_paths : AnyPathIter
        Paths to whole files, either bare paths or 1-tuples of a path, as
        :func:`~polyphys.manage.utils.sort_filenames` returns them. Each file
        is loaded independently.
    parser : ParserClassT
        A parser class from :mod:`polyphys.manage.parser` that interprets
        file paths and extracts attributes such as the `whole` name.
    group : GroupT
        The particle group to which the data pertains.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary where keys are the `whole` names (derived from each file
        path) and values are numpy arrays containing the loaded data.

    Raises
    ------
    TypeError
        If an item of `whole_paths` is neither a path nor a sequence whose
        first element is a path.

    Notes
    -----
    This function assumes that each file in `whole_paths` contains data
    related to a complete simulation, referred to as a "whole." For further
    information on terminology, see the `organizer` module documentation.

    Let n_wholes be the number of files and m the number of values in one
    file. Time and memory are O(n_wholes * m): every whole is loaded and
    held in memory at once.

    Examples
    --------
    >>> import tempfile
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> tmp = tempfile.TemporaryDirectory()
    >>> whole = ('N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000'
    ...          'adump5000ens1')
    >>> path = Path(tmp.name) / f'{whole}-bug-gyrTMon.npy'
    >>> np.save(path, np.array([1.0, 2.0, 3.0]))
    >>> whole_from_file([path], SumRuleCyl, 'bug')[whole]
    array([1., 2., 3.])
    >>> tmp.cleanup()
    """
    wholes: dict[str, np.ndarray] = {}
    for whole_path in _normalize_paths(whole_paths):
        whole_info = parser(whole_path, "whole", group, ispath=True)
        wholes[getattr(whole_info, "whole")] = np.load(whole_path)

    return wholes


def whole_from_dist_mat_t(
    whole_paths: list[tuple[str]],
    parser: type[ParserInstanceT],
    group: GroupT,
    func: Callable[[str, ParserInstanceT],
                   tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]]
) -> tuple[
    dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, pd.DataFrame]
]:
    """
    Load pair distance matrix 'whole' data from given file paths, processing
    specific distance-related physical properties (i.e., frequencies, radial
    distribution functions, time series).

    Depending on the project, you must supply the `func` argument that knows
    how to process the distance matrix for your specific project.

    Parameters
    ----------
    whole_paths : list[tuple[str]]
        List of 1-tuples, as :func:`~polyphys.manage.utils.sort_filenames`
        returns them, each holding the path to a file representing a 'whole'
        dataset (typically a NumPy `.npy` distance matrix). Bare path
        strings are not supported: the first element of each item is taken
        as the path.
    parser : type[ParserInstanceT]
        A parser class used to create parser instances for each file. The
        parser provides attributes such as `whole` (name identifier) and
        simulation metadata.
    group : GroupT
        The particle group to which the data pertains (e.g., 'bug', 'all').
    func : Callable[[str, ParserInstanceT],
                    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]]
        Project-specific processing function. It must accept two arguments:

        - `whole_path` (str): path to the distance matrix file.
        - `whole_info` (ParserInstanceT): parser instance for the file.

        It must return a tuple of three pandas DataFrames corresponding to
        frequency/histogram data, radial distribution data, and time-series
        data.

    Returns
    -------
    tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame],
          dict[str, pd.DataFrame]]
        Three dictionaries with whole names as keys and DataFrames as values:

        - Frequency/histogram DataFrames
        - Radial distribution function (RDF) DataFrames
        - Time series DataFrames

    Notes
    -----
    Arguments are not validated here; errors raised by `parser` or `func`
    propagate.

    Let n_wholes be the number of files. `func` is called once per whole,
    and its three outputs for every whole are held in memory at once, so
    time and memory are n_wholes times those of one call of `func`.

    Examples
    --------
    A stand-in for a project-specific processor that summarizes each matrix
    by its mean:

    >>> import tempfile
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> def summarize(whole_path, whole_info):
    ...     frame = pd.DataFrame({'mean': [np.load(whole_path).mean()]})
    ...     return frame, frame, frame
    >>> tmp = tempfile.TemporaryDirectory()
    >>> whole = ('N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000'
    ...          'adump5000ens1')
    >>> path = Path(tmp.name) / f'{whole}-bug-distMatTFoci.npy'
    >>> np.save(path, np.array([[0.0, 2.0], [2.0, 0.0]]))
    >>> freqs, rdfs, tseries = whole_from_dist_mat_t(
    ...     [(str(path),)], SumRuleCyl, 'bug', summarize)
    >>> freqs[whole]
       mean
    0   1.0
    >>> tmp.cleanup()
    """
    wholes_freqs: dict[str, pd.DataFrame] = {}
    wholes_rdfs: dict[str, pd.DataFrame] = {}
    wholes_tseries: dict[str, pd.DataFrame] = {}
    for whole_path in whole_paths:
        whole_info = parser(whole_path[0], "whole", group, ispath=True)
        whole_freqs, whole_rdfs, whole_tseries = \
            func(whole_path[0], whole_info)
        whole_name = getattr(whole_info, "whole")
        wholes_freqs[whole_name] = whole_freqs
        wholes_rdfs[whole_name] = whole_rdfs
        wholes_tseries[whole_name] = whole_tseries
    return wholes_freqs, wholes_rdfs, wholes_tseries


def ens_from_bin_edge(
    ens_data: tuple[str, dict[str, np.ndarray]],
) -> tuple[str, np.ndarray]:
    """
    Reduce the bin edges of the 'whole' datasets in an 'ensemble' to the
    single array they share.

    Every 'whole' in an ensemble must be binned identically, so that the
    histograms paired with these edges stay aligned bin by bin. The edges are
    checked for equality with :func:`_assert_equal_bin_edges`, the check
    :func:`whole_from_segments` applies to segments, rather than merged.

    Parameters
    ----------
    ens_data : tuple[str, dict[str, np.ndarray]]
        A tuple where the first element is the ensemble name (str) and the
        second element is a non-empty dictionary of bin edge arrays from
        different 'whole' datasets within that ensemble. Within the
        dictionary, keys are whole names while values are 1D arrays of whole
        bin edges, each of shape ``(n_bins + 1,)``.

    Returns
    -------
    tuple[str, np.ndarray]
        A tuple where the first element is the ensemble name, and the second
        is a copy of the shared 1D array of bin edges, of shape
        ``(n_bins + 1,)``.

    Raises
    ------
    ValueError
        If the bin edges of any 'whole' differ from those of the first one,
        in length or in value.

    Notes
    -----
    Time and space complexity are O(n_wholes * n_edges) and O(n_edges).

    Examples
    --------
    >>> ens_from_bin_edge(
    ...     ('ens1', {'w1': np.array([0., 1., 2.]),
    ...               'w2': np.array([0., 1., 2.])})
    ... )
    ('ens1', array([0., 1., 2.]))

    Differing edges raise instead of being merged:

    >>> ens_from_bin_edge(
    ...     ('ens1', {'w1': np.array([0., 1., 2.]),
    ...               'w2': np.array([0., 1., 3.])})
    ... )
    Traceback (most recent call last):
        ...
    ValueError: Bin edges are not identical across all arrays.
    """
    ensemble_name, whole_edges = ens_data
    shared_bin_edges = _assert_equal_bin_edges(list(whole_edges.values()))
    return ensemble_name, shared_bin_edges


def ens_from_vec(
    ens_data: tuple[str, dict[str, np.ndarray]],
) -> tuple[str, pd.DataFrame]:
    """
    Convert a dictionary of 1D arrays (representing different 'whole'
    datasets) into a single ensemble-level DataFrame. This is useful when each
    'whole' dataset in the ensemble contains a vector (1D array) of
    measurements.

    Parameters
    ----------
    ens_data : tuple[str, dict[str, np.ndarray]]
        A tuple where the first element is the ensemble name, and the second
        is a dictionary mapping each 'whole' dataset name to a 1D numpy array
        of measurements.

    Returns
    -------
    tuple[str, pd.DataFrame]
        A tuple where the first element is the ensemble name, and the second
        is a pandas DataFrame. Each column represents a 'whole' dataset, in
        the order of the dictionary, and rows contain corresponding vector
        values from each dataset.

    Notes
    -----
    Let n_wholes be the number of wholes and m the length of one vector.
    Time and extra memory are O(n_wholes * m): the vectors are copied into
    the DataFrame.

    Examples
    --------
    >>> ens_name, ensemble_df = ens_from_vec(
    ...     ('ensemble1', {'whole1': np.array([0.1, 0.2, 0.3]),
    ...                    'whole2': np.array([0.4, 0.5, 0.6])})
    ... )
    >>> ens_name
    'ensemble1'
    >>> ensemble_df
       whole1  whole2
    0     0.1     0.4
    1     0.2     0.5
    2     0.3     0.6
    """
    ensemble_name, whole_vectors = ens_data
    ensemble_df = pd.DataFrame.from_dict(whole_vectors, orient="columns")
    return ensemble_name, ensemble_df


def ens_from_mat_t(
    ens_data: tuple[str, dict[str, np.ndarray]],
) -> tuple[str, np.ndarray]:
    """
    Combine multiple 2D arrays into a single 3D array for a specified
    'ensemble', where each 2D array corresponds to a 'whole' dataset in the
    ensemble.

    Parameters
    ----------
    ens_data : tuple[str, dict[str, np.ndarray]]
        A tuple where the first element is the ensemble name (str) and the
        second element is a dictionary with 'whole' names as keys and 2D numpy
        arrays as values.

    Returns
    -------
    tuple[str, np.ndarray]
        A tuple where the first element is the ensemble name and the second
        is a 3D numpy array of shape ``(n_wholes, n_rows, n_cols)``, formed
        by stacking the 2D arrays along the first axis in the order of the
        dictionary.

    Notes
    -----
    Let n_wholes be the number of wholes and m the number of elements of one
    matrix. Time and extra memory are O(n_wholes * m): the matrices are
    copied into the stacked array.

    Examples
    --------
    >>> ens_name, combined = ens_from_mat_t(
    ...     ('ensemble1', {'whole1': np.array([[1, 2], [3, 4]]),
    ...                    'whole2': np.array([[5, 6], [7, 8]])})
    ... )
    >>> combined.shape
    (2, 2, 2)
    >>> combined[1]
    array([[5, 6],
           [7, 8]])
    """
    ensemble_name, whole_matrices = ens_data
    combined_matrix = np.stack(list(whole_matrices.values()), axis=0)
    return ensemble_name, combined_matrix


def ens_from_df(
    ens_data: tuple[str, dict[str, pd.DataFrame]],
) -> tuple[str, pd.DataFrame]:
    """
    Average multiple DataFrames in an ensemble by grouping them into a single
    averaged DataFrame. Each DataFrame corresponds to a 'whole' dataset in the
    ensemble.

    In a 'whole' DataFrame, headers are 'elements' of a 'matrix' or 2D
    quantity and columns are the values of a given property.

    In each 'whole' DataFrame, headers represent 'elements' of a 2D quantity
    (matrix-like structure), and columns contain the values for a specified
    property. Each 'whole' DataFrame may be either a 'histogram' or a
    'time-series':

    - In a 'histogram' DataFrame, an additional column, `bin_center`, holds
      bin center values that are invariant under averaging and are retained
      in the output DataFrame. The length of histogram DataFrame is equal to
      the number of bins.

    - In a 'time-series' DataFrame, headers (columns) are elements of a matrix,
      and rows represent the values of each element over time. The length of
      timeseries DataFrame is equal to the number of time frames.


    Parameters
    ----------
    ens_data : tuple[str, dict[str, pd.DataFrame]]
        A tuple where the first element is the ensemble name and the second
        is a dictionary with 'whole' names as keys and DataFrames, sharing
        their columns and index labels, as values.

    Returns
    -------
    tuple[str, pd.DataFrame]
        A tuple with:

        - The ensemble name (str).
        - A DataFrame with the averaged values across all 'whole'
          DataFrames.

    Notes
    -----
    A clear distinction exists between ensembles of 'vector' or 'matrix' types
    and ensembles of 'dataframe' type. In vector or matrix types, each column
    in the resulting ensemble DataFrame is named after a 'whole' and represents
    that whole's values for a specific property. In contrast, for DataFrame
    type ensembles, the columns are elements of a 2D quantity (e.g., matrix or
    time-series values) and hold averaged values across all whole DataFrames.

    This function assumes that averaging non-element headers, such as
    'bin_center,' does not change their values, and these are carried forward
    without modification.

    Let n_wholes be the number of wholes, n_rows the rows and n_cols the
    columns of one DataFrame. Time and extra memory are
    O(n_wholes * n_rows * n_cols): the DataFrames are concatenated into one
    before they are grouped by index label and averaged.

    Examples
    --------
    >>> df1 = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    >>> df2 = pd.DataFrame({'A': [2, 3, 4], 'B': [5, 6, 7]})
    >>> ens_name, averaged_df = ens_from_df(
    ...     ('ensemble1', {'whole1': df1, 'whole2': df2})
    ... )
    >>> averaged_df
         A    B
    0  1.5  4.5
    1  2.5  5.5
    2  3.5  6.5
    """
    ensemble_name, whole_data = ens_data
    ensemble_df = pd.concat(list(whole_data.values())).groupby(level=0).mean()
    return ensemble_name, ensemble_df


def ensemble(
    prop: PropertyT,
    wholes: dict[str, np.ndarray] | dict[str, pd.DataFrame],
    parser: ParserClassT,
    group: GroupT,
    whole_type: Literal["vector", "matrix", "dataframe", "bin_edge"],
    whole_edges: dict[str, np.ndarray] | None = None,
    save_to: Pathish | None = None,
) -> dict[str, np.ndarray | pd.DataFrame]:
    """
    Generate an ensemble by  merging 'whole' data arrays or DataFrames
    representing a specified `prop_` of a particle `group` within a
    defined `geometry`. The `whole_type` determines the structure of each
    'whole' and the merging approach.

    Parameters
    ----------
    prop: PropertyT
        The physical property for which the ensemble is generated (e.g.,
        'density').
    wholes : dict[str, np.ndarray] | dict[str, pd.DataFrame]
        Dictionary containing 'whole' names as keys and their corresponding
        data as values, represented as numpy arrays or DataFrames.
    parser : ParserClassT
        An instance of the parser class to infer attributes like ensemble name
        from each 'whole'.
    group : GroupT
        Type of the particle group being processed.
    whole_type : {'vector', 'matrix', 'dataframe', 'bin_edge'}
        Specifies the type of each 'whole' and defines the merging method:

        - 'vector': A 1D numpy array, such as a time series or histogram.
        - 'matrix': A 2D numpy array, such as a gyration tensor.
        - 'dataframe': A DataFrame representing either a time series (with
          headers as matrix elements and rows as values over time) or a
          histogram (headers as elements, rows as counts or frequencies,
          and an additional 'bin_center' column for bin centers).
        - 'bin_edge': A 1D array representing bin edges. Assumes all bin edges
          are identical and selects one to represent the ensemble.

    whole_edges : dict[str, np.ndarray], optional
        Dictionary with 'whole' names as keys and 1D bin edge arrays, each of
        shape ``(n_bins + 1,)``, as values. Only valid when `whole_type` is
        'vector'. The wholes of an ensemble share their bin edges by design,
        so the edges of the first whole of each ensemble are used and the
        others are not compared. Each merged ensemble gets a 'bin_center'
        column of length ``n_bins`` holding the midpoints of those edges.
        Must contain every key of `wholes`.
    save_to : str, optional
        Directory path where output files are saved, if specified.

    Returns
    -------
    dict[str, np.ndarray | pd.DataFrame]
        Dictionary where keys are ensemble names and values are ensemble data
        as numpy arrays or DataFrames.

    Raises
    ------
    ValueError
        If any parameter has an invalid value; if `whole_edges` is given and
        `whole_type` is not 'vector'; or if the number of bins differs from
        the number of rows of the merged ensemble.

    Notes
    -----
    Different averaging methods apply based on `whole_type`. When `whole_type`
    is "DataFrame," the resulting DataFrame averages the values of each matrix
    element across wholes, retaining headers that are mean-invariant (e.g.,
    'bin_center'). For 'vector' and 'matrix' types, columns represent
    individual wholes in the ensemble DataFrame.

    The bin centre of bin ``i`` is the midpoint of its two edges,
    ``0.5 * (edges[i] + edges[i + 1])``. A centre is invariant under
    averaging over wholes, which is why :func:`ensemble_avg` excludes
    'bin_center' from its statistics by default.

    Let n_wholes be the number of wholes, m the number of values in one
    whole, n_ensembles the number of ensembles, and n_edges the length of
    one bin-edge array. Grouping parses each whole name once and merging
    copies every whole once, so time and extra memory are
    O(n_wholes * m): the caller's `wholes` and the merged ensembles are both
    held in memory. With `whole_edges`, computing the centres adds
    O(n_ensembles * n_edges) time and memory.

    With `save_to`, each ensemble is saved by :func:`save_artifact` under
    ``<ensemble_long>-<group>-<prop>``, as ``.csv`` for a DataFrame and
    ``.npy`` for an array.

    Examples
    --------
    Two replicates of one state point, with shared edges, merge into one
    ensemble that carries the edge midpoints as 'bin_center':

    >>> from polyphys.manage.parser import SumRuleCyl
    >>> base = ('N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000'
    ...         'adump10000ens')
    >>> wholes = {base + '1': np.array([1., 2., 3.]),
    ...           base + '2': np.array([3., 4., 5.])}
    >>> edges = {name: np.array([0., 1., 2., 4.]) for name in wholes}
    >>> ens = ensemble('rhoMon', wholes, SumRuleCyl, 'bug', 'vector',
    ...                whole_edges=edges)
    >>> ens_df = next(iter(ens.values()))
    >>> ens_df['bin_center'].tolist()
    [0.5, 1.5, 3.0]

    Without `whole_edges`, the ensemble has one column per whole only:

    >>> ens = ensemble('rhoMon', wholes, SumRuleCyl, 'bug', 'vector')
    >>> 'bin_center' in next(iter(ens.values())).columns
    False
    """
    invalid_keyword(whole_type, ["vector", "matrix", "dataframe", "bin_edge"])
    if whole_edges is not None and whole_type != "vector":
        raise ValueError(
            "'whole_edges' is only valid with whole_type 'vector', not "
            f"'{whole_type}'."
        )

    merging_func: dict[
        str, Callable[..., tuple[str, np.ndarray | pd.DataFrame]]
    ] = {
        "vector": ens_from_vec,
        "matrix": ens_from_mat_t,
        "dataframe": ens_from_df,
        "bin_edge": ens_from_bin_edge,
    }
    # Collecting wholes for each ensemble; the wholes of an ensemble share
    # their bin edges by design, so the centres come from its first whole
    ens_wholes: dict[str, dict[str, np.ndarray | pd.DataFrame]] = {}
    bin_centers: dict[str, np.ndarray] = {}
    for whole_name, whole_data in wholes.items():
        whole_info = parser(whole_name, "whole", group)
        ens_name = getattr(whole_info, "ensemble_long")
        ens_wholes.setdefault(ens_name, {})[whole_name] = whole_data
        if whole_edges is not None and ens_name not in bin_centers:
            edges = np.asarray(whole_edges[whole_name])
            bin_centers[ens_name] = 0.5 * (edges[:-1] + edges[1:])

    # Merging wholes into a single ensemble data object
    ensembles: dict[str, np.ndarray | pd.DataFrame] = dict(
        map(merging_func[whole_type], ens_wholes.items())
    )
    for ens_name, centers in bin_centers.items():
        ensembles[ens_name]["bin_center"] = centers

    # Save ensemble data object if save_to is specified
    if save_to:
        for ens_name, data in ensembles.items():
            ens_fullname = create_fullname(ens_name, group, prop)
            save_artifact(ens_fullname, data, save_to)

    return ensembles


def ens_avg_from_bin_edge(
    ens_data: np.ndarray | Sequence[np.ndarray],
) -> np.ndarray:
    """
    Return the bin edges of an 'ensemble' for its ensemble-averaged
    artifact.

    Bin edges are not averaged: every 'whole' in the ensemble must share
    them. A single 1D array, the form :func:`ens_from_bin_edge` produces and
    :func:`ensemble_avg` passes here, is returned as a copy. A sequence of
    arrays, one per 'whole', is checked for equality with
    :func:`_assert_equal_bin_edges`, and a copy of the shared array is
    returned.

    Parameters
    ----------
    ens_data : np.ndarray or Sequence[np.ndarray]
        Either the ensemble's 1D array of bin edges, of shape
        ``(n_bins + 1,)``, or a non-empty sequence of such arrays, one per
        'whole'. A 2D array of shape ``(n_wholes, n_bins + 1)`` counts as a
        sequence of its rows.

    Returns
    -------
    np.ndarray
        A copy of the shared 1D array of bin edges, of shape
        ``(n_bins + 1,)``.

    Raises
    ------
    ValueError
        If `ens_data` is a sequence and any array in it differs from the
        first one, in length or in value.

    Notes
    -----
    This function performs no averaging or statistical calculation. Time and
    space complexity are O(n_wholes * n_edges) and O(n_edges).

    Examples
    --------
    >>> ens_avg_from_bin_edge(np.array([0.0, 0.1, 0.2]))
    array([0. , 0.1, 0.2])
    >>> ens_avg_from_bin_edge(
    ...     [np.array([0.0, 0.1, 0.2]), np.array([0.0, 0.1, 0.2])]
    ... )
    array([0. , 0.1, 0.2])
    >>> ens_avg_from_bin_edge(
    ...     [np.array([0.0, 0.1, 0.2]), np.array([0.0, 0.1, 0.3])]
    ... )
    Traceback (most recent call last):
        ...
    ValueError: Bin edges are not identical across all arrays.
    """
    if isinstance(ens_data, np.ndarray) and ens_data.ndim == 1:
        return ens_data.copy()
    return _assert_equal_bin_edges([np.asarray(edges) for edges in ens_data])


def ens_avg_from_ndarray(
    ens_prop: PropertyT, ens_data: np.ndarray
) -> dict[str, np.ndarray]:
    """
    Compute ensemble statistics (mean, variance, and standard error of the
    mean)across an 'ensemble' 3D array where each slice along the first axis
    (axis=0) represents a 'whole' data matrix. The statistics are calculated
    along this axis.

    Parameters
    ----------
    ens_prop : PropertyT
        The name of the property to which this 'ensemble' data belongs. Used
        as a prefix for the computed statistics.
    ens_data : np.ndarray
        A 3D numpy array where each matrix slice along axis=0 represents a
        'whole' dataset.

    Returns
    -------
    dict[str, np.ndarray]
        A dictionary where keys are statistical measures (mean, variance, SEM)
        with `ens_prop` as a prefix while values are numpy arrays of the
        computed statistics across 'whole' matrices.

    Notes
    -----
    - The SEM is calculated with a degrees-of-freedom adjustment (`ddof=1`) for
      an unbiased estimate and is scaled by the square root of the number of
      wholes.
    - This function performs averaging along axis=0 for each element across
      the 'whole' data matrices.
    - Let n_wholes be the length of axis 0 and m the number of elements of
      one whole. Time is O(n_wholes * m). NumPy's ``var`` and ``std`` build
      temporaries of the size of `ens_data`, so peak extra memory is
      O(n_wholes * m); the three results take O(m).

    Examples
    --------
    >>> ens_data = np.array([
    ...     [[1.0, 2.0], [3.0, 4.0]],
    ...     [[2.0, 3.0], [4.0, 5.0]],
    ... ])
    >>> ens_avg = ens_avg_from_ndarray('gyrTMon', ens_data)
    >>> list(ens_avg)
    ['gyrTMon-mean', 'gyrTMon-var', 'gyrTMon-sem']
    >>> ens_avg['gyrTMon-mean']
    array([[1.5, 2.5],
           [3.5, 4.5]])
    >>> ens_avg['gyrTMon-sem']
    array([[0.5, 0.5],
           [0.5, 0.5]])
    """
    n_wholes_sqrt = ens_data.shape[0] ** 0.5
    ens_avg = {
        f"{ens_prop}-mean": np.mean(ens_data, axis=0),
        f"{ens_prop}-var": np.var(ens_data, axis=0, ddof=1),
        f"{ens_prop}-sem": np.std(ens_data, axis=0, ddof=1) / n_wholes_sqrt,
    }
    return ens_avg


def ens_avg_from_df(
    ens_prop: PropertyT, ens_data: pd.DataFrame, exclude: list[str]
) -> pd.DataFrame:
    """
    Calculate the mean, variance, and standard error of the mean (SEM) for
    columns representing 'whole' data within an 'ensemble' DataFrame. The
    resulting DataFrame retains excluded columns (e.g., bin_center), while
    replacing 'whole' columns with computed statistics.

    Parameters
    ----------
    ens_prop : PropertyT
        The name of the property to which this 'ensemble' data belongs, used
        as a prefix for the generated statistics columns.
    ens_data : pd.DataFrame
        A DataFrame containing ensemble data. Each 'whole' dataset is in its
        own column, and any columns specified in `exclude` remain unaffected.
    exclude : list of str
        A list of column names to exclude from averaging. For example,
        `exclude` might contain 'bin_center' for histogram-type ensembles.

    Returns
    -------
    pd.DataFrame
        A new DataFrame with the excluded columns of `ens_data`, followed by
        the columns ``<ens_prop>-mean``, ``<ens_prop>-var``, and
        ``<ens_prop>-sem``: the mean, variance, and SEM of each row across
        the 'whole' data columns. It has the same index as `ens_data`.

    Notes
    -----
    - `ens_data` is not modified; the statistics are written to a new
      DataFrame.
    - The function averages only the 'whole' columns, excluding any specified
      in `exclude`.
    - The standard error of the mean (SEM) is computed with a
      degrees-of-freedom adjustment (ddof=1) for an unbiased estimate.
    - Let n_rows be the number of rows of `ens_data` and n_wholes the number
      of 'whole' columns. Time is O(n_rows * n_wholes). Each statistic
      selects the 'whole' columns as a new DataFrame, so peak extra memory
      is O(n_rows * n_wholes), not O(n_rows); the result takes O(n_rows).

    Examples
    --------
    >>> ens_data = pd.DataFrame({
    ...     'bin_center': [1, 2, 3],
    ...     'whole1': [4.0, 5.0, 6.0],
    ...     'whole2': [5.0, 6.0, 7.0],
    ... })
    >>> ens_avg_from_df('density', ens_data, ['bin_center'])
       bin_center  density-mean  density-var  density-sem
    0           1           4.5          0.5          0.5
    1           2           5.5          0.5          0.5
    2           3           6.5          0.5          0.5
    >>> list(ens_data.columns)
    ['bin_center', 'whole1', 'whole2']
    """
    wholes = [col for col in ens_data.columns if col not in exclude]

    ens_avg = ens_data.drop(columns=wholes)
    ens_avg[ens_prop + "-mean"] = ens_data[wholes].mean(axis=1)
    ens_avg[ens_prop + "-var"] = ens_data[wholes].var(axis=1, ddof=1)
    ens_avg[ens_prop + "-sem"] = ens_data[wholes].sem(axis=1, ddof=1)
    return ens_avg


def ensemble_avg(
    prop: PropertyT,
    ensembles: dict[str, np.ndarray | pd.DataFrame],
    group: GroupT,
    ens_type: Literal["dataframe", "ndarray", "bin_edge"],
    exclude: list[str] | None = None,
    save_to: Pathish | None = None,
) -> dict[str, np.ndarray] | dict[str, pd.DataFrame]:
    """
    Generate ensemble-averaged data for a specified `prop_` by averaging
    'whole' data in each 'ensemble' DataFrame or array within `ensembles`.
    Columns listed in `exclude` are omitted from averaging.

    Parameters
    ----------
    prop: PropertyT
        The physical property for which ensemble averages are calculated.
    ensembles : dict[str, np.ndarray | pd.DataFrame]
        A dictionary where each key is an ensemble name and each value is
        either a DataFrame or ndarray representing 'whole' data.
    group : GroupT
        Type of the particle group.
    ens_type : {'dataframe', 'ndarray', 'bin_edge'}
        The data format for each ensemble:

        - 'dataframe': DataFrame format, with columns as 'whole' simulations
          or items in `exclude`.
        - 'ndarray': 3D ndarray format, with each slice along axis=0
          representing a 'whole' matrix.
        - 'bin_edge': 1D ndarray format for bin edges, with bin edges assumed
          to be identical across wholes.

    exclude : list[str], default None
        List of column names or items not included in ensemble averaging,
        such as 'bin_center' for histograms.
    save_to : str, optional
        Path to the directory where the averaged ensemble data will be saved.

    Returns
    -------
    dict[str, np.ndarray] | dict[str, pd.DataFrame]
        A dictionary where each key is an ensemble name, and each value is a
        DataFrame or ndarray of ensemble-averaged data.

    Raises
    ------
    ValueError
        If `ens_type` is not one of the accepted values.

    Notes
    -----
    With `save_to`, each result is saved by :func:`save_artifact` under
    ``<ensemble>-<group>-<prop>-ensAvg``: a DataFrame as ``.csv``, bin edges
    as ``.npy``, and the dictionary of an 'ndarray' ensemble as one
    ``...-ensAvg-<measure>.npy`` file per statistic.

    Let n_ensembles be the number of ensembles. Ensembles are averaged one
    after another, with the time and peak extra memory of
    :func:`ens_avg_from_df`, :func:`ens_avg_from_ndarray`, or
    :func:`ens_avg_from_bin_edge` for the largest one; the results for all
    n_ensembles ensembles are held in memory together.

    Examples
    --------
    >>> ensembles = {
    ...     'ensemble1': pd.DataFrame({'whole1': [1.0, 2.0],
    ...                                'whole2': [3.0, 4.0]}),
    ...     'ensemble2': pd.DataFrame({'whole1': [5.0, 6.0],
    ...                                'whole2': [5.0, 8.0]}),
    ... }
    >>> avg = ensemble_avg('rhoMon', ensembles, 'bug', 'dataframe')
    >>> avg['ensemble2']
       rhoMon-mean  rhoMon-var  rhoMon-sem
    0          5.0         0.0         0.0
    1          7.0         2.0         1.0
    """
    if exclude is None:
        exclude = ["bin_center"]

    # Compute ensemble averages
    ens_avgs: dict[str, np.ndarray] | dict[str, pd.DataFrame] = {}
    if ens_type == "bin_edge":
        for ens_name, ens_data in ensembles.items():
            ens_avg = ens_avg_from_bin_edge(ens_data)  # type: ignore
            ens_avgs[ens_name] = ens_avg  # type: ignore
    elif ens_type == "ndarray":
        for ens_name, ens_data in ensembles.items():
            ens_avg = ens_avg_from_ndarray(prop, ens_data)  # type: ignore
            ens_avgs[ens_name] = ens_avg  # type: ignore
    elif ens_type == "dataframe":
        for ens_name, ens_data in ensembles.items():
            ens_avg = ens_avg_from_df(prop, ens_data, exclude)  # type: ignore
            ens_avgs[ens_name] = ens_avg  # type: ignore
    else:
        invalid_keyword(ens_type, ["dataframe", "ndarray", "bin_edge"])

    # Save averaged ensemble data if save_to path is provided
    if save_to is not None:
        save_property = prop + "-ensAvg"
        for ens_name, data in ens_avgs.items():
            ens_fullname = create_fullname(ens_name, group, save_property)
            save_artifact(ens_fullname, data, save_to)

    return ens_avgs


def children_stamps(
    stamps: list[tuple[str]],
    lineage: PrimitiveLineageT,
    save_to: Pathish | None = None,
) -> pd.DataFrame:
    """
    Combine individual 'stamp' CSV files into a single DataFrame for a given
    `lineage` in a *space*. Optionally, saves the resulting
    DataFrame to a specified directory.

    Parameters
    ----------
    stamps : list[tuple[str]]
        List of 1-tuples, as :func:`~polyphys.manage.utils.sort_filenames`
        returns them, each holding the path to a 'stamp' CSV file. Rows are
        concatenated in this order.
    lineage : {'segment', 'whole'}
        Lineage type, either 'segment' or 'whole'.
    save_to : str, optional
        Directory path where the combined DataFrame is saved, if specified,
        as ``<space>-<group>-<lineage>-stamps.csv``, with ``space`` and
        ``group`` read from the first row.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing data from all the 'stamp' files, with
        a fresh ``0..n-1`` index.

    Raises
    ------
    ValueError
        If `lineage` is not one of the accepted values.

    Warns
    -----
    UserWarning
        If `lineage` is 'whole' and the 'segment' and 'segment_id' columns
        are dropped.

    Notes
    -----
    - It is assumed that all the stamp files have belong to the same *group* in
      a given *space*.
    - If `lineage` is 'whole', columns 'segment' and 'segment_id' are removed,
      as they are redundant for whole-lineage data, with a `UserWarning`.
      Either may be absent, as in stamps the probe phase writes for a whole;
      if neither is present, nothing is dropped and a DEBUG record is
      logged.
    - Let n_files be the number of stamp files, n_rows the total number of
      rows, and n_cols the number of columns. Time is O(n_rows * n_cols).
      Every file is read before the frames are concatenated, so peak memory
      is about twice the combined data, O(n_rows * n_cols).

    Examples
    --------
    >>> import tempfile
    >>> tmp = tempfile.TemporaryDirectory()
    >>> paths = []
    >>> for seg_id in (1, 2):
    ...     path = Path(tmp.name) / f'w1.j0{seg_id}-stamps.csv'
    ...     pd.DataFrame({'segment': [f'w1.j0{seg_id}'], 'whole': ['w1'],
    ...                   'segment_id': [seg_id], 'n_frames': [10]}
    ...                  ).to_csv(path, index=False)
    ...     paths.append((str(path),))
    >>> children_stamps(paths, 'segment')
      segment whole  segment_id  n_frames
    0  w1.j01    w1           1        10
    1  w1.j02    w1           2        10
    >>> tmp.cleanup()
    """
    invalid_keyword(lineage, ["segment", "whole"])

    # Load and concatenate all stamps
    space_stamps = pd.concat(
        [pd.read_csv(stamp[0]) for stamp in stamps], ignore_index=True
    )

    # Drop specific columns for 'whole' lineage
    cols_to_drop = ["segment", "segment_id"]
    if lineage == "whole":
        dropped = [col for col in cols_to_drop if col in space_stamps.columns]
        if dropped:
            space_stamps.drop(columns=dropped, inplace=True)
            warnings.warn(
                "'segment' and 'segment_id' columns are dropped when"
                " individual 'whole' stamps combined to create a single"
                " dataframe of 'whole' stamps by 'children_stamps'.",
                UserWarning,
            )
        else:
            logger.debug("No %s columns to drop from the whole stamps.",
                         cols_to_drop)

    # Save the resulting DataFrame to a CSV if specified
    if save_to is not None:
        space_name = str(space_stamps.loc[0, "space"])
        group = str(space_stamps.loc[0, "group"])
        filename = f"{space_name}-{group}-{lineage}-stamps.csv"
        space_stamps.to_csv(os.path.join(save_to, filename), index=False)

    return space_stamps


def parents_stamps(
    stamps: pd.DataFrame,
    lineage: PrimitiveLineageT,
    properties: dict[str, Callable[..., Any]] | None = None,
    save_to: Pathish | None = None,
) -> pd.DataFrame:
    """
    Aggregate data from child 'stamps' into a parent stamp by applying
    specified aggregation functions in a *space*. Optionally, saves the
    resulting DataFrame to a CSV file.

    Parameters
    ----------
    stamps : pd.DataFrame
        DataFrame containing child stamp data, the lineage of whose rows is
        `lineage`. It must have a 'whole' and a 'segment_id' column for
        'segment' stamps, and an 'ensemble_long' and an 'ensemble_id' column
        for 'whole' stamps, and an 'n_frames' column for both.
    lineage : {'segment', 'whole'}
        Lineage type of the child stamps, either 'segment' or 'whole'.
    properties : dict of str -> Callable, optional
        Dictionary specifying aggregation functions for specific columns,
        replacing the default 'last'.
    save_to : str, optional
        Directory path where the parent DataFrame is saved, if specified,
        as ``<space>-<group>-stamps-whole.csv`` for 'segment' stamps and
        ``<space>-<group>-stamps-ensAvg.csv`` for 'whole' stamps, with
        ``space`` and ``group`` read from the first row.

    Returns
    -------
    pd.DataFrame
        Aggregated parent DataFrame, one row per parent:

        - 'segment': one row per 'whole', with 'n_segments' (the number of
          segments) and 'n_frames' (their sum).
        - 'whole': one row per 'ensemble_long', with 'n_ensembles' (the
          number of wholes) and 'n_frames' (the last whole's); 'segment',
          'segment_id', and 'n_segments' are dropped if present.

        Every other column takes its last value in each group unless
        `properties` says otherwise; 'lineage_name' and the `lineage` column
        are dropped.

    Raises
    ------
    ValueError
        If `lineage` is not one of the accepted values.

    Notes
    -----
    - It is assumed that all the stamp files have belong to the same *group* in
      a given *space*.
    - If `lineage` is 'segment', stamps correspond to individual segments,
      each with a unique `segment_id`.
    - If `lineage` is 'whole', stamps correspond to whole simulations, each
      identified by an `ensemble_id`.
    - The 'lineage_name' column and the `lineage` column are not carried
      into the parent stamps. Either may be absent, as 'lineage_name' is
      from the stamps the probe phase writes; an absent one is skipped with
      a DEBUG record.
    - Let n_rows be the number of child stamps and n_cols the number of
      columns. The group-by is hash-based, so expected time is
      O(n_rows * n_cols); memory is one copy of the parent stamps.

    Examples
    --------
    >>> stamps = pd.DataFrame({
    ...     'lineage_name': ['w1.j01', 'w1.j02', 'w2.j01'],
    ...     'segment': ['w1.j01', 'w1.j02', 'w2.j01'],
    ...     'whole': ['w1', 'w1', 'w2'],
    ...     'segment_id': [1, 2, 1],
    ...     'n_frames': [10, 20, 30],
    ... })
    >>> parents_stamps(stamps, 'segment')
      whole  n_segments  n_frames
    0    w1           2        30
    1    w2           1        30
    """
    invalid_keyword(lineage, ["segment", "whole"])

    # Base aggregation functions for columns
    base_columns = list(stamps.columns)
    for column in ("lineage_name", lineage):
        if column in base_columns:
            base_columns.remove(column)
        else:
            logger.debug("Column %r is not among the stamps columns.", column)
    agg_funcs: dict[str, str | Callable[..., Any]] = {
        col: "last" for col in base_columns
    }
    if properties is not None:
        agg_funcs.update(properties)

    # Define grouping column and lineage-specific aggregations
    if lineage == "whole":
        parent_groupby = "ensemble_long"
        agg_funcs.update({"ensemble_id": "count", "n_frames": "last"})
    else:
        parent_groupby = "whole"
        agg_funcs.update({"segment_id": "count", "n_frames": "sum"})

    # Perform the groupby and aggregate
    agg_funcs.pop(parent_groupby)
    parent_stamps = stamps.groupby(parent_groupby).agg(agg_funcs).reset_index()

    # Rename and drop columns based on lineage
    if lineage == "whole":
        parent_stamps.rename(
            columns={"ensemble_id": "n_ensembles"}, inplace=True
        )
        # If the 'whole' stamps are generated directly in the 'probe' phase,
        # then 'segment' and 'segment_id' columns are "N/A" and are removed
        # from the list of stamps columns that are added to the parents
        # stamps.
        # There is no need to have the 'n_segment' column in the parent
        # stamps, so it is removed. The 'whole' stamps directly generated
        # in the 'probe' phase do not have such a column, but those generated
        # from 'segment' stamps have.
        # Dropping redundant columns silently:
        parent_stamps.drop(
            columns=["n_segments", "segment_id", "segment"],
            inplace=True,
            errors="ignore",
        )
    else:
        parent_stamps.rename(
            columns={"segment_id": "n_segments"}, inplace=True
        )

    if save_to is not None:
        space_name = str(parent_stamps.loc[0, "space"])
        group = str(parent_stamps.loc[0, "group"])
        file_suffix = "ensAvg" if lineage == "whole" else "whole"
        filename = f"{space_name}-{group}-stamps-{file_suffix}.csv"
        parent_stamps.to_csv(os.path.join(save_to, filename), index=False)

    return parent_stamps


def find_unique_properties(
    filepath: str,
    prop_idx: int,
    extensions: list[str],
    drop_properties: list[str] | None = None,
    sep: str = "-",
) -> tuple[list[str], list[str]]:
    """
    Extract unique physical properties and property-measures from filenames
    matched by a glob pattern. The function identifies unique segments in
    filenames based on specified extensions and index position, then sorts and
    returns them.

    Parameters
    ----------
    filepath : str
        The glob-friendly pattern used to locate filenames, e.g.,
        `path/to/files/*`.
    prop_idx : int
        The index position in the filename where the property or
        property-measure name starts after splitting by the separator.
    extensions : list[str]
        A list of suffixes that indicate the end of a property or
        property-measure name (e.g., `'-ensAvg'`, `'-ens'`, `'-whole'`). These
        are not file extensions like `".csv"` or `".npy"`.
    drop_properties : list[str] | None, default None
        A list of property names to ignore when determining unique properties
        and property-measures.
    sep : str, default '-'
        The separator used to split properties and measures within a filename.

    Returns
    -------
    tuple[list[str], list[str]]
        A tuple containing:

        - **uniq_props** (*list[str]*): A sorted list of unique physical
          properties.
        - **uniq_prop_measures** (*list[str]*): A sorted list of unique
          property-measures.

    Examples
    --------
    Three ensemble-averaged files of one ensemble; the property starts at
    index 2 after splitting the name by ``-``:

    >>> import tempfile
    >>> tmp = tempfile.TemporaryDirectory()
    >>> for name in ('N20-bug-gyrTMon-ensAvg.csv',
    ...              'N20-bug-gyrTMon-acf-ensAvg.csv',
    ...              'N20-bug-rhoMon-ensAvg.csv'):
    ...     (Path(tmp.name) / name).touch()
    >>> find_unique_properties(str(Path(tmp.name) / '*'), prop_idx=2,
    ...                        extensions=['-ensAvg'])
    (['gyrTMon', 'rhoMon'], ['gyrTMon-acf'])
    >>> tmp.cleanup()

    Notes
    -----
    - Ensure `prop_idx` accurately reflects the location in the filename where
      the property or measure name starts after splitting by `sep`.
    - This function assumes a consistent filename format where property-measure
      segments follow each other with defined separators and suffixes.
    - Every file is split by every entry of `extensions`, so each matched
      file should contain all of them; a file that lacks one contributes
      the rest of its basename, file extension included, under that entry.
    - Basenames are taken by splitting paths on ``/``.
    - Let n_files be the number of matched files, n_ext the number of
      `extensions`, and n_uniq the number of unique names. Time is
      O(n_files * n_ext + n_uniq log n_uniq); memory is O(n_files) for the
      list of matched paths.
    """
    props_measures = glob(filepath)
    uniq_prop_measures = set()
    for ext in extensions:
        for prop in props_measures:
            prop_name = sep.join(
                prop.split("/")[-1].split(ext)[0].split(sep)[prop_idx:]
            )
            uniq_prop_measures.add(prop_name)

    if drop_properties is not None:
        uniq_prop_measures.difference_update(drop_properties)

    # Extract unique properties (first segment before `sep`)
    # from property-measures
    uniq_props = {prop.split(sep)[0] for prop in uniq_prop_measures}

    # Remove any full properties that match unique property names from
    # prop-measures
    uniq_prop_measures.difference_update(uniq_props)

    # Sort and return the results as lists
    return sorted(uniq_props), sorted(uniq_prop_measures)


def space_tseries(
    input_database: str,
    prop: PropertyT,
    parser: ParserClassT,
    hierarchy: str,
    physical_attrs: list[str],
    group: GroupT,
    dump_key: str,
    divisor: float = 0.025,
    round_to: int = 3,
    is_save: bool | None = False,
    crowder_attr: str = "phi_bulk_c",
) -> pd.DataFrame:
    """
    Aggregate ensemble-averaged time-series data for a specified physical
    property across multiple files, adds specified physical attributes as
    columns, and concatenates into a single DataFrame.

    Parameters
    ----------
    input_database : str
        Path to the directory containing the ensemble-averaged time-series
        files, named ``<ensemble_long>-<group>-<prop>-ensAvg.csv``, ending
        with a path separator.
    prop : str
        Name of the physical property of interest.
    parser : ParserClassT
        Parser class to infer file-specific information from filenames or
        paths.
    hierarchy : str
        Pattern prefix for the time-series filenames (e.g., `"N*"`).
    physical_attrs : list[str]
        List of physical attributes to add as new columns in the output
        DataFrame.
    group : GroupT
        Particle group type.
    dump_key : str
        Dumping attribute/keyword used to extract dumping frequency from an
        instance of the `parser`.
    divisor : float, default 0.025
        Rounding step for the crowder volume fraction named by
        `crowder_attr`.
    round_to : int, default 3
        Number of decimal places for rounding the crowder volume fraction.
    is_save : bool, default False
        If True, saves the concatenated DataFrame to a CSV file named
        ``<space>-<group>-<prop>-space.csv`` in the ``analysis`` phase's
        ``space`` directory, which :func:`make_database` creates next to
        `input_database`.
    crowder_attr : str, default 'phi_bulk_c'
        Name of the `parser` attribute that holds the bulk volume fraction
        of crowders. Its value is added as a column of that name, whether or
        not it is listed in `physical_attrs`, and rounded into the
        ``phi_c_bulk_round`` column.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of all time-series with added physical
        attributes. Besides the columns of each file, it has ``t_index``
        (the frame index times the dump frequency named by `dump_key`),
        ``time`` (``t_index`` times the parser's ``dt``), one column per
        entry in `physical_attrs`, a `crowder_attr` column, and
        ``phi_c_bulk_round``.

    Examples
    --------
    One ensemble-averaged time series of two frames, dumped every 5000
    steps of 0.005 time units:

    >>> import tempfile
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> tmp = tempfile.TemporaryDirectory()
    >>> database = Path(tmp.name) / 'analysis' / 'N20-bug-ensAvg'
    >>> database.mkdir(parents=True)
    >>> state = 'N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000'
    >>> pd.DataFrame({'gyrTMon-mean': [1.0, 2.0]}).to_csv(
    ...     database / f'{state}-bug-gyrTMon-ensAvg.csv', index=False)
    >>> space = space_tseries(str(database) + '/', 'gyrTMon', SumRuleCyl,
    ...                       'N*', ['nmon'], 'bug', 'adump')
    >>> space[['t_index', 'time', 'gyrTMon-mean', 'nmon']]
       t_index  time  gyrTMon-mean  nmon
    0        0   0.0           1.0    20
    1     5000  25.0           2.0    20
    >>> tmp.cleanup()

    Notes
    -----
    - Files are selected by globbing ``input_database + hierarchy +
      '-<prop>-ensAvg.csv'`` and sorted alphanumerically by name.
    - Requires a parser class with methods to retrieve attribute information
      for each file.
    - Let n_files be the number of matched files, n_rows the total number
      of rows across them, and n_cols the number of columns of the result.
      Time is O(n_rows * n_cols + n_files log n_files). Every file is read
      and kept before the frames are concatenated, so peak memory is about
      twice the result, O(n_rows * n_cols).
    """
    prop_ext = f"-{prop}-ensAvg.csv"
    ens_avg_csvs = sort_filenames(
        glob(input_database + hierarchy + prop_ext), ["csv"]
    )
    prop_csvs = []

    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        prop_info = parser(
            ens_avg_csv[0], "ensemble_long", group, ispath=True
        )

        ens_avg.reset_index(inplace=True)
        ens_avg.rename(columns={"index": "t_index"}, inplace=True)

        # Calculate `t_index` and `time` columns based on `dumping_freq`
        ens_avg["t_index"] *= getattr(prop_info, dump_key)
        ens_avg["time"] = ens_avg["t_index"] * getattr(prop_info, "dt")

        # Add physical attributes
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(prop_info, attr_name)

        # Add and round the bulk volume fraction of crowders
        ens_avg[crowder_attr] = getattr(prop_info, crowder_attr)
        ens_avg["phi_c_bulk_round"] = ens_avg[crowder_attr].apply(
            lambda x: round_to_multiple(x, divisor, round_to)
        )

        prop_csvs.append(ens_avg)

    # Concatenate all time-series DataFrames
    prop_db = pd.concat(prop_csvs, axis=0)
    prop_db.reset_index(drop=True, inplace=True)

    # Optionally save to file
    if is_save:
        save_to_space = make_database(
            input_database, "analysis", stage="space", group=group
        )
        space = save_to_space.split("/")[-2].split("-")[0]
        filepath = (
            save_to_space + f"{create_fullname(space, group, prop)}-space.csv"
        )
        prop_db.to_csv(filepath, index=False)

    return prop_db


def space_hists(
    input_database: str,
    prop: PropertyT,
    parser: ParserClassT,
    hierarchy: str,
    physical_attrs: list[str],
    group: GroupT,
    bin_center: np.ndarray | None = None,
    normalize: bool | None = False,
    divisor: float = 0.025,
    round_to: int = 3,
    is_save: bool | None = False,
    crowder_attr: str = "phi_bulk_c",
) -> pd.DataFrame:
    """
    Aggregate ensemble-averaged histogram data for a specified physical
    property across multiple files, normalizes data if specified, and
    concatenates into a single DataFrame.

    Parameters
    ----------
    input_database : str
        Path to the directory containing the ensemble-averaged histogram
        files, named ``<ensemble_long>-<group>-<prop>-ensAvg.csv``, ending
        with a path separator.
    prop: str
        Name of the physical property of interest.
    parser : ParserClassT
        Parser class to infer file-specific information from filenames or
        paths.
    hierarchy : str
        Pattern prefix for the histogram filenames (e.g., `"N*"`).
    physical_attrs : list[str]
        List of physical attributes to add as new columns in the output
        DataFrame.
    group : GroupT
        Particle group type.
    bin_center : np.ndarray, optional
        Array of bin centers. If not provided, must be present in the
        DataFrames.
    normalize : bool, default False
        If True, normalizes the histogram data.
    divisor : float, default 0.025
        Rounding step for the crowder volume fraction named by
        `crowder_attr`.
    round_to : int, default 3
        Number of decimal places for rounding the crowder volume fraction.
    is_save : bool, default False
        If True, saves the concatenated DataFrame to a CSV file named
        ``<space>-<group>-<prop>-space.csv`` in the ``analysis`` phase's
        ``space`` directory, which :func:`make_database` creates next to
        `input_database`.
    crowder_attr : str, default 'phi_bulk_c'
        Name of the `parser` attribute that holds the bulk volume fraction
        of crowders. Its value is added as a column of that name, whether or
        not it is listed in `physical_attrs`, and rounded into the
        ``phi_c_bulk_round`` column.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of all histograms with added physical
        attributes: one column per entry in `physical_attrs`, a
        `crowder_attr` column, and ``phi_c_bulk_round``.

    Warns
    -----
    UserWarning
        If `normalize` is True and the ``<prop>-mean`` column of a file sums
        to zero; its ``<prop>-norm`` column is then set to zero.

    Examples
    --------
    One ensemble-averaged histogram of three bins, normalized to sum to one:

    >>> import tempfile
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> tmp = tempfile.TemporaryDirectory()
    >>> database = Path(tmp.name) / 'analysis' / 'N20-bug-ensAvg'
    >>> database.mkdir(parents=True)
    >>> state = 'N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000'
    >>> pd.DataFrame({'rhoMon-mean': [1.0, 2.0, 1.0]}).to_csv(
    ...     database / f'{state}-bug-rhoMon-ensAvg.csv', index=False)
    >>> space = space_hists(str(database) + '/', 'rhoMon', SumRuleCyl, 'N*',
    ...                     ['nmon'], 'bug',
    ...                     bin_center=np.array([0.5, 1.5, 2.5]),
    ...                     normalize=True)
    >>> space[['bin_center', 'rhoMon-mean', 'rhoMon-norm', 'nmon']]
       bin_center  rhoMon-mean  rhoMon-norm  nmon
    0         0.5          1.0         0.25    20
    1         1.5          2.0         0.50    20
    2         2.5          1.0         0.25    20
    >>> tmp.cleanup()

    Notes
    -----
    - Files are selected by globbing ``input_database + hierarchy +
      '-<prop>-ensAvg.csv'`` and sorted alphanumerically by name.
    - If `normalize` is True, histogram values will be scaled to sum to 1.
    - Let n_files be the number of matched files, n_rows the total number
      of rows across them, and n_cols the number of columns of the result.
      Time is O(n_rows * n_cols + n_files log n_files). Every file is read
      and kept before the frames are concatenated, so peak memory is about
      twice the result, O(n_rows * n_cols).
    - The `bin_center` should be provided if it is not available in the input
      data. When all the ensemble-averaged DataFrames have the same bin
      center, like 'clustersHistTFoci' or 'bondsHistTFoci' properties, the
      `bin_center` may not be provided.
    """
    prop_ext = f"-{prop}-ensAvg.csv"
    ens_avg_csvs_ungrouped = glob(input_database + hierarchy + prop_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs_ungrouped, formats=["csv"])
    prop_csvs = []

    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        prop_info = parser(
            ens_avg_csv[0], "ensemble_long", group, ispath=True
        )

        # Handle bin_center if provided
        if bin_center is not None:
            ens_avg["bin_center"] = bin_center.tolist()
            ens_avg["bin_center-norm"] = (
                ens_avg["bin_center"] / ens_avg["bin_center"].max()
            )

        # Normalize data if specified
        if normalize:
            normalizer = ens_avg[prop + "-mean"].sum()
            if normalizer != 0:
                ens_avg[prop + "-norm"] = ens_avg[prop + "-mean"] / normalizer
            else:
                warnings.warn(
                    "All values are zero; normalized values set to zero.",
                    UserWarning,
                )
                ens_avg[prop + "-norm"] = 0

        # Add physical attributes
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(prop_info, attr_name)

        # Add and round the bulk volume fraction of crowders
        ens_avg[crowder_attr] = getattr(prop_info, crowder_attr)
        ens_avg["phi_c_bulk_round"] = ens_avg[crowder_attr].apply(
            lambda x: round_to_multiple(x, divisor, round_to)
        )

        prop_csvs.append(ens_avg)

    # Concatenate all histograms
    prop_db = pd.concat(prop_csvs, axis=0)
    prop_db.reset_index(drop=True, inplace=True)

    # Optionally save to file
    if is_save:
        save_to_space = make_database(
            input_database, "analysis", stage="space", group=group
        )
        space = save_to_space.split("/")[-2].split("-")[0]
        filepath = (
            save_to_space + f"{create_fullname(space, group, prop)}-space.csv"
        )
        prop_db.to_csv(filepath, index=False)

    return prop_db
