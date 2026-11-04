"""Tests for :mod:`polyphys.probe.core`.

The inputs are tier A: a data file and a dump of one linear chain written
by the ``make_chain_files`` fixture, and parsers built from hand-written
names. Expected values are the positions the fixture wrote and the
attributes the parser sets.
"""
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.manage.organizer import children_stamps
from polyphys.manage.parser import SumRuleCyl

pytest.importorskip("MDAnalysis")

from polyphys.probe import make_universe, write_stamps  # noqa: E402

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
def whole_parser() -> SumRuleCyl:
    """A SumRuleCyl parser of a whole, from a hand-written name."""
    return SumRuleCyl(
        'N10epsilon5r2.5lz10sig1nc0dt0.005bdump1adump1ens1', 'whole', 'bug')


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
    second = SumRuleCyl(
        'N10epsilon5r2.5lz10sig1nc0dt0.005bdump1adump1ens2', 'whole', 'bug')
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
