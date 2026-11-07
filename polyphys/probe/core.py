"""
==========================================================
:mod:`polyphys.probe.core`
==========================================================

The :mod:`polyphys.probe.core` module holds what every probe shares: how a
study's LAMMPS data file and dump are opened as one MDAnalysis
:class:`~MDAnalysis.core.universe.Universe`, and the *stamps* file, the
one-row table of a probed artifact's lineage names, parser attributes, and
frame count that the organizer later combines across a space.

The research code opened every universe with the same four arguments and
wrote the stamps with the same columns in each of its probe functions;
:func:`make_universe` and :func:`write_stamps` are those two pieces, written
once.

Contents
========

- :func:`make_universe`, a universe of a LAMMPS data file and dump in the
  convention of the studies.
- :func:`write_stamps`, the stamps CSV of one probed artifact.
- :data:`UNIVERSE_KWARGS`, the fixed arguments of :func:`make_universe`.

Dependencies
============
- Python standard library: `copy`, `os`, `pathlib`, `types`, `typing`.
- `MDAnalysis <https://docs.mdanalysis.org/>`_ 2.8 or later, through the
  ``mda`` extra.
- `pandas <https://pandas.pydata.org/>`_.
- :mod:`polyphys.manage.parser`: :class:`~polyphys.manage.parser.ParserBase`.

Every probe analysis computes each part of a split run on its own shallow
copy of itself (``_compute_part``), so that the ``serial`` backend, which
hands every part the same analysis object, keeps each part's ``results``,
``frames``, and ``times`` apart, as the ``multiprocessing`` and ``dask``
backends do by pickling.

Examples
========
>>> import tempfile
>>> from polyphys.manage.parser import SumRuleCyl
>>> parser = SumRuleCyl('N4epsilon5r2.5lz10sig1nc0dt0.005bdump1adump1ens1',
...                     'whole', 'bug')
>>> with tempfile.TemporaryDirectory() as tmp:
...     path = write_stamps(parser, n_frames=5, save_to=tmp)
...     path.name
'N4epsilon5r2.5lz10sig1nc0dt0.005bdump1adump1ens1-bug-stamps.csv'

.. versionadded:: 0.7.0
"""
import copy
import os
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import Results

from ..manage.parser import ParserBase
from ..manage.types import Pathish

#: The fixed arguments of :func:`make_universe`: a LAMMPS data file as the
#: topology, a LAMMPS dump as the trajectory with its coordinates read as the
#: wrapped, unscaled ``x y z`` columns, and the ``Atoms`` lines of the data
#: file read as ``id resid type x y z``, the ``bond``, ``angle``, and
#: ``molecular`` atom styles whose molecule ID becomes the residue ID that
#: a selection such as ``'resid 1'`` picks.
#:
#: .. versionadded:: 0.7.0
UNIVERSE_KWARGS: Mapping[str, str] = MappingProxyType({
    'topology_format': 'DATA',
    'format': 'LAMMPSDUMP',
    'lammps_coordinate_convention': 'unscaled',
    'atom_style': 'id resid type x y z',
})


def make_universe(
    topology: Pathish,
    trajectory: Pathish,
    **kwargs: Any,
) -> mda.Universe:
    """
    Open a LAMMPS data file and dump as one MDAnalysis universe, as the
    studies read them.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    topology : str or pathlib.Path
        The LAMMPS data file, whose ``Atoms`` lines are ``atom-ID
        molecule-ID atom-type x y z`` and whose ``Masses`` section gives
        every atom type a mass. The masses weight the radius of gyration
        and the gyration tensor; without them MDAnalysis cannot form the
        tensor.
    trajectory : str or pathlib.Path
        The LAMMPS dump, in the ``atom`` style with unscaled, wrapped
        coordinates (``ITEM: ATOMS id type x y z``). MDAnalysis sorts the
        atoms of each frame by ``id``, so ``positions`` are in atom-ID
        order, the order along the chain.
    **kwargs
        Further arguments of :class:`~MDAnalysis.core.universe.Universe`,
        such as ``dt``, the time between frames, or ``in_memory``. The
        keys of :data:`UNIVERSE_KWARGS` are fixed and may not be given.

    Returns
    -------
    MDAnalysis.core.universe.Universe
        The universe, with the trajectory at its first frame.

    Raises
    ------
    ValueError
        If `kwargs` names a key of :data:`UNIVERSE_KWARGS`.

    Warns
    -----
    UserWarning
        From the dump reader, ``Reader has no dt information``, when `dt` is
        not given: the frame times are then the dump's time steps.

    Notes
    -----
    The dump's box is read per frame; the data file's box is not used after
    the topology is built. Opening a universe reads the whole data file and
    indexes the dump, :math:`O(N + F)` time for :math:`N` atoms and
    :math:`F` frames, and holds the topology, :math:`O(N)` space; frames
    are read one at a time unless `in_memory` is set.

    Examples
    --------
    A two-bead chain in a one-frame dump:

    >>> import tempfile
    >>> import numpy as np
    >>> from polyphys.simulate.core import Box
    >>> from polyphys.simulate.data import DataFile, write_data
    >>> pair = DataFile(
    ...     atom_style='bond', box=Box(lo=np.zeros(3), hi=np.full(3, 5.0)),
    ...     ids=np.array([1, 2]), types=np.array([1, 1]),
    ...     positions=np.array([[2.0, 2.5, 2.5], [3.0, 2.5, 2.5]]),
    ...     n_atom_types=1, masses=np.array([1.0]), mol_ids=np.array([1, 1]),
    ...     bonds=np.array([[1, 1, 2]]), n_bond_types=1)
    >>> dump = '''ITEM: TIMESTEP
    ... 0
    ... ITEM: NUMBER OF ATOMS
    ... 2
    ... ITEM: BOX BOUNDS pp pp pp
    ... 0 5
    ... 0 5
    ... 0 5
    ... ITEM: ATOMS id type x y z
    ... 2 1 3.0 2.5 2.5
    ... 1 1 2.0 2.5 2.5
    ... '''
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     data = write_data(pair, Path(tmp) / 'pair.data')
    ...     _ = (Path(tmp) / 'pair.lammpstrj').write_text(dump)
    ...     universe = make_universe(data, Path(tmp) / 'pair.lammpstrj',
    ...                              dt=1.0)
    ...     chain = universe.select_atoms('resid 1')
    ...     chain.n_atoms, universe.trajectory.n_frames, chain.ids.tolist()
    (2, 1, [1, 2])
    """
    fixed = sorted(set(kwargs) & set(UNIVERSE_KWARGS))
    if fixed:
        raise ValueError(
            f"{fixed} are fixed by make_universe and cannot be given.")
    return mda.Universe(os.fspath(topology), os.fspath(trajectory),
                        **UNIVERSE_KWARGS, **kwargs)


def write_stamps(
    parser: ParserBase,
    n_frames: int,
    save_to: Pathish,
) -> Path:
    """
    Write the stamps of one probed artifact: a one-row CSV of its lineage
    names, its parser attributes, and its number of frames.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    parser : :class:`~polyphys.manage.parser.ParserBase`
        The artifact's parser, at the ``segment`` or ``whole`` lineage,
        built from the trajectory's name or path.
    n_frames : int
        The number of frames probed, zero or more.
    save_to : str or pathlib.Path
        An existing directory; the file is
        ``<name>-<group>-stamps.csv`` in it, with the parser's `name` and
        `group`.

    Returns
    -------
    pathlib.Path
        The file written.

    Raises
    ------
    ValueError
        If `n_frames` is not a non-negative integer.

    Notes
    -----
    The columns are, in order, the lineage names of
    :attr:`~polyphys.manage.parser.ParserBase.lineage_genealogy`, the
    ``group``, every name in
    :attr:`~polyphys.manage.parser.ParserBase.attributes`, and
    ``n_frames``; the row holds the parser's values. The research code's
    ``stamps_report`` wrote the same columns without ``group``, which
    :func:`~polyphys.manage.organizer.children_stamps` reads from the first
    row when it saves. Time and space are :math:`O(a)` for :math:`a`
    attributes.

    Examples
    --------
    >>> import tempfile
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> parser = SumRuleCyl(
    ...     'N4epsilon5r2.5lz10sig1nc0dt0.005bdump1adump1ens1.j01', 'segment',
    ...     'bug')
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     stamps = pd.read_csv(write_stamps(parser, 3, tmp))
    >>> stamps.columns[:3].tolist(), stamps.columns[-2:].tolist()
    (['segment', 'whole', 'ensemble_long'], ['rho_bulk_c', 'n_frames'])
    >>> int(stamps.loc[0, 'nmon']), int(stamps.loc[0, 'n_frames'])
    (4, 3)
    """
    if isinstance(n_frames, bool) or not isinstance(n_frames, int) \
            or n_frames < 0:
        raise ValueError(
            f"'n_frames' must be a non-negative integer, not {n_frames!r}.")
    columns = [*parser.lineage_genealogy, 'group', *parser.attributes,
               'n_frames']
    values = [
        *(getattr(parser, lineage) for lineage in parser.lineage_genealogy),
        parser.group,
        *(getattr(parser, attr) for attr in parser.attributes),
        n_frames,
    ]
    stamps = pd.DataFrame([values], columns=columns)
    path = Path(save_to) / f"{parser.name}-{parser.group}-stamps.csv"
    stamps.to_csv(path, index=False)
    return path


def _compute_part(
    analysis: AnalysisBase,
    indexed_frames: np.ndarray,
    verbose: bool | None = None,
    *,
    progressbar_kwargs: dict[str, Any] | None = None,
) -> AnalysisBase:
    """
    Compute one part of a run on a shallow copy of `analysis` with empty
    results, and return the copy.

    MDAnalysis 2.8 and later split a run's frames into parts and call
    ``_compute`` once per part, then merge the returned objects' ``results``,
    ``frames``, and ``times`` with the analysis's ``_get_aggregator``. The
    ``multiprocessing`` and ``dask`` backends pickle the analysis for each
    part; the ``serial`` backend passes the analysis itself to every part,
    so without a copy each part overwrites the one before it and the merge
    repeats the last part ``n_parts`` times.

    Parameters
    ----------
    analysis : MDAnalysis.analysis.base.AnalysisBase
        The analysis whose ``_compute`` calls this, with its frames set up
        by ``run``.
    indexed_frames : numpy.ndarray
        The part's frames, shape ``(n, 2)``: the run's iteration index and
        the trajectory frame of each.
    verbose : bool or None, default None
        Whether to show a progress bar; ``None`` keeps the analysis's own
        setting.
    progressbar_kwargs : dict or None, default None
        Arguments of :class:`MDAnalysis.lib.log.ProgressBar`.

    Returns
    -------
    MDAnalysis.analysis.base.AnalysisBase
        The copy, holding this part's ``results``, ``frames``, and
        ``times``; `analysis` is unchanged.

    Notes
    -----
    The copy shares the analysis's atom groups, trajectory, and settings;
    only ``results`` is new, and ``_compute`` rebinds the part's frames,
    times, and sliced trajectory on the copy. The copy costs :math:`O(a)`
    time and space for :math:`a` attributes, independent of the frames and
    atoms.

    Examples
    --------
    Two parts of a four-frame run, computed serially, stay apart:

    >>> import MDAnalysis as mda
    >>> from polyphys.probe.chain_size import ChainSize
    >>> frames = np.zeros((4, 2, 3))
    >>> frames[:, 1, 2] = [1.0, 2.0, 3.0, 4.0]
    >>> universe = mda.Universe.empty(2, trajectory=True)
    >>> universe.add_TopologyAttr('masses', np.ones(2))
    >>> universe.load_new(frames, format='MEMORY')
    <Universe with 2 atoms>
    >>> analysis = ChainSize(universe.atoms, ('rfloryT',))
    >>> analysis.run(backend='serial', n_parts=2).results['rfloryT']
    array([1., 2., 3., 4.])
    """
    part = copy.copy(analysis)
    part.results = Results()
    return AnalysisBase._compute(part, indexed_frames, verbose,
                                 progressbar_kwargs=progressbar_kwargs)
