"""Tests for :mod:`polyphys.probe.core`.

The inputs are tier A: a data file and a dump of one linear chain written
by the ``make_chain_files`` fixture, and parsers built from hand-written
names. Expected values are the positions the fixture wrote and the
attributes the parser sets.

The real-frame tests open the tier B files of the studies' wholes
(``polyphys/tests/probe/conftest.py``, "Real frames"). Their expected
values are the dump headers and coordinates as the package's own LAMMPS
reader, :func:`polyphys.simulate.dump.iter_dump`, reads them. The tests
that run a study's probe spec on these frames are in the polyphys-lab
repository, with the studies' parsers.
"""
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import polyphys
from polyphys.manage.organizer import children_stamps
from polyphys.tests.synthetic import ChainCyl
from polyphys.simulate.dump import iter_dump

pytest.importorskip("MDAnalysis")

from polyphys.probe import make_universe, write_stamps  # noqa: E402
from polyphys.probe.binders import BinderStatistics  # noqa: E402
from polyphys.probe.chain_size import ChainSize  # noqa: E402
from polyphys.probe.core import (  # noqa: E402
    _ClosingMultiprocessing, _ProbeAnalysis)
from polyphys.probe.histograms import SpatialHistograms  # noqa: E402
from polyphys.probe.pairs import PairDistances  # noqa: E402

# The dump reader warns that the frames carry no time step; the tests that
# give no dt measure per frame.
pytestmark = pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")


@pytest.fixture
def chain_frames() -> np.ndarray:
    """Five frames of a 10-bead chain along x, moved 0.1 along z per frame."""
    first = np.zeros((10, 3))
    first[:, 0] = np.arange(10) + 5.0
    first[:, 1:] = 10.0
    return np.stack([first + [0.0, 0.0, 0.1 * k] for k in range(5)])


@pytest.fixture
def whole_parser() -> ChainCyl:
    """A ChainCyl parser of a whole, from a hand-written name."""
    return ChainCyl(
        'N10r2.5lz10ac1nc0dt0.005bdump1adump1ens1', 'whole', 'bug')


# make_universe


def test_make_universe_reads_the_chain_and_its_frames(
        make_chain_files, chain_frames):
    """The universe has the data file's atoms, all in residue 1, and the
    dump's frames."""
    data_path, dump_path = make_chain_files(chain_frames)

    universe = make_universe(data_path, dump_path)

    assert (universe.atoms.n_atoms, universe.trajectory.n_frames) == (10, 5)
    assert_array_equal(universe.select_atoms('resid 1').ids, np.arange(1, 11))


def test_make_universe_positions_follow_the_dump_frame_by_frame(
        make_chain_files, chain_frames):
    """Frame k of the trajectory holds the dump's frame k, in atom-ID
    order."""
    data_path, dump_path = make_chain_files(chain_frames)
    universe = make_universe(data_path, dump_path)

    read = np.stack([universe.atoms.positions.copy()
                     for _ in universe.trajectory])

    # 1e-5: the dump holds six decimals and MDAnalysis reads float32.
    assert_allclose(read, chain_frames, atol=1e-5)


def test_make_universe_passes_extra_keywords_through(make_chain_files,
                                                     chain_frames):
    """A keyword such as dt reaches the MDAnalysis universe."""
    data_path, dump_path = make_chain_files(chain_frames)

    universe = make_universe(data_path, dump_path, dt=2.0)

    assert universe.trajectory.dt == 2.0


def test_make_universe_rejects_a_fixed_keyword(make_chain_files,
                                               chain_frames):
    """A keyword that make_universe fixes cannot be overridden."""
    data_path, dump_path = make_chain_files(chain_frames)

    with pytest.raises(ValueError, match=r"\['format'\] are fixed"):
        make_universe(data_path, dump_path, format='DCD')


# write_stamps


def test_write_stamps_columns_are_genealogy_group_attributes_and_frames(
        whole_parser, tmp_path):
    """The stamps have the lineage names, group, exactly the parser's
    attributes in order, then n_frames, with the parser's values."""
    path = write_stamps(whole_parser, 7, tmp_path)

    stamps = pd.read_csv(path)

    assert path.name == f'{whole_parser.name}-bug-stamps.csv'
    assert stamps.columns.tolist() == (
        ['whole', 'ensemble_long', 'ensemble', 'space', 'group']
        + whole_parser.attributes + ['n_frames'])
    assert stamps['whole'].tolist() == [whole_parser.name]
    assert stamps['group'].tolist() == ['bug']
    assert stamps['n_frames'].tolist() == [7]
    for attr in whole_parser.attributes:
        assert stamps.loc[0, attr] == getattr(whole_parser, attr), attr


def test_write_stamps_reads_back_through_children_stamps(whole_parser,
                                                         tmp_path):
    """The organizer combines two written stamps into one row per whole."""
    second = ChainCyl(
        'N10r2.5lz10ac1nc0dt0.005bdump1adump1ens2', 'whole', 'bug')
    paths = [(str(write_stamps(parser, n, tmp_path)),)
             for parser, n in ((whole_parser, 3), (second, 4))]

    combined = children_stamps(paths, 'whole')

    assert combined['ensemble_id'].tolist() == [1, 2]
    assert combined['n_frames'].tolist() == [3, 4]


@pytest.mark.parametrize('n_frames', [
    pytest.param(-1, id='negative'),
    pytest.param(2.5, id='float'),
    pytest.param(True, id='bool'),
])
def test_write_stamps_rejects_a_bad_frame_count(whole_parser, tmp_path,
                                                n_frames):
    """n_frames must be a non-negative integer."""
    with pytest.raises(ValueError, match="non-negative integer"):
        write_stamps(whole_parser, n_frames, tmp_path)


# The multiprocessing backend


def test_multiprocessing_run_leaves_no_leaked_semaphore():
    """A run on the multiprocessing backend lets its workers exit, so the
    resource tracker finds no leaked semaphore at shutdown."""
    script = textwrap.dedent("""\
        import numpy as np
        import MDAnalysis as mda
        from polyphys.probe.chain_size import ChainSize

        if __name__ == '__main__':
            frames = np.zeros((4, 2, 3))
            frames[:, 1, 2] = [1.0, 2.0, 3.0, 4.0]
            universe = mda.Universe.empty(2, trajectory=True)
            universe.add_TopologyAttr('masses', np.ones(2))
            universe.load_new(frames, format='MEMORY')
            ChainSize(universe.atoms, ('rfloryT',)).run(
                backend='multiprocessing', n_workers=2, n_parts=2)
        """)
    # The child and its workers import this checkout's package.
    env = dict(os.environ,
               PYTHONPATH=str(Path(polyphys.__file__).parents[1]))

    completed = subprocess.run(
        [sys.executable, '-W', 'default', '-c', script], env=env,
        capture_output=True, text=True, timeout=120, check=False)

    assert completed.returncode == 0, completed.stderr
    assert 'leaked semaphore' not in completed.stderr


def test_closing_multiprocessing_apply_passes_a_worker_error_on():
    """An error raised in a worker reaches the caller, after the pool is
    terminated and joined."""
    backend = _ClosingMultiprocessing(n_workers=1)

    # int of a non-numeric string raises ValueError in the worker.
    with pytest.raises(ValueError, match="invalid literal for int"):
        backend.apply(int, ['1', 'x'])


@pytest.mark.parametrize('analysis_class', [
    pytest.param(BinderStatistics, id='BinderStatistics'),
    pytest.param(ChainSize, id='ChainSize'),
    pytest.param(PairDistances, id='PairDistances'),
    pytest.param(SpatialHistograms, id='SpatialHistograms'),
])
def test_every_probe_analysis_runs_on_the_closing_pool(analysis_class):
    """Every probe analysis takes its run, and so its multiprocessing
    backend, from the shared probe base."""
    assert analysis_class.run is _ProbeAnalysis.run


# Real frames


def test_make_universe_opens_every_real_polymer_dump(
        real_polymer_whole, real_files, real_universe):
    """A polymer dump opens with its topology, with the dump's atoms and
    all 50 of its frames."""
    files = real_files(real_polymer_whole)
    n_atoms = next(iter_dump(files.polymer)).n_atoms

    universe = real_universe(files.topology, files.polymer)

    assert (universe.atoms.n_atoms, universe.trajectory.n_frames) == (
        n_atoms, 50)


def test_make_universe_opens_every_real_window(real_window_whole, real_files,
                                               real_universe):
    """An ``all`` window opens with its whole's data file, with the dump's
    atoms and every frame of the window."""
    files = real_files(real_window_whole)
    n_atoms = next(iter_dump(files.window)).n_atoms

    universe = real_universe(files.data, files.window)

    assert (universe.atoms.n_atoms, universe.trajectory.n_frames) == (
        n_atoms, real_window_whole.n_window_frames)


def test_make_universe_positions_are_the_dump_coordinates_from_the_origin(
        real_polymer_whole, real_files, real_universe, trajectory_arrays):
    """Every frame's positions are the dump's coordinates, sorted by atom
    ID, less the box's lower corner."""
    files = real_files(real_polymer_whole)
    expected = np.stack([frame.positions() - frame.box.lo
                         for frame in iter_dump(files.polymer)])

    positions = trajectory_arrays(
        real_universe(files.topology, files.polymer))['positions']

    # float32 positions of at most about 165 carry an absolute round-off
    # of about 1e-5.
    assert_allclose(positions, expected, rtol=0, atol=2e-5)
