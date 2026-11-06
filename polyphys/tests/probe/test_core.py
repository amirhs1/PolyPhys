"""Tests for :mod:`polyphys.probe.core`.

The inputs are tier A: a data file and a dump of one linear chain written
by the ``make_chain_files`` fixture, and parsers built from hand-written
names. Expected values are the positions the fixture wrote and the
attributes the parser sets.

The real-frame tests read the tier B and tier C files of the studies'
wholes (``polyphys/tests/probe/conftest.py``, "Real frames"). Their
expected values are the dump headers and coordinates as the package's own
LAMMPS readers, :func:`polyphys.simulate.dump.iter_dump` and
:func:`polyphys.simulate.data.read_data`, read them, and the particle
counts in the wholes' names.
"""
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.manage.organizer import children_stamps
from polyphys.manage.parser import SumRuleCyl
from polyphys.simulate.data import read_data
from polyphys.simulate.dump import iter_dump

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


def test_spec_selections_pick_the_species_of_a_real_polymer_dump(
        real_polymer_whole, real_files, real_universe, study_spec):
    """Each selection of the study's probe spec picks, in the polymer dump,
    the atoms the whole's name gives its species."""
    probe = study_spec(real_polymer_whole.study).PROBE
    files = real_files(real_polymer_whole)
    universe = real_universe(files.topology, files.polymer)

    counts = {species: universe.select_atoms(selection).n_atoms
              for species, selection in probe.selections.items()}

    assert counts == real_polymer_whole.polymer_counts()


def test_spec_selections_pick_the_species_of_a_real_window(
        real_window_whole, real_files, real_universe, study_spec):
    """Each selection of the study's probe spec picks, in the ``all``
    window, the atoms the whole's name gives its species, crowders
    included."""
    probe = study_spec(real_window_whole.study).PROBE
    files = real_files(real_window_whole)
    universe = real_universe(files.data, files.window)

    counts = {species: universe.select_atoms(selection).n_atoms
              for species, selection in probe.selections.items()}

    assert counts == real_window_whole.window_counts()


@pytest.mark.parametrize('study, whole, expected', [
    pytest.param('hns_cub', 'N200kbmm2nh30ac2l25epshc1nc10743ens1',
                 {'nmon': 200, 'nhns': 30, 'dcrowd': 2.0, 'lcube': 50.0,
                  'eps_hc': 1.0, 'ncrowd': 10743, 'ensemble_id': 1},
                 id='HnsCub'),
    pytest.param('hns_cyl', 'N200kbmm2r4.5nh30ac2lz75epshc1nc386ens2',
                 {'nmon': 200, 'nhns': 30, 'dcrowd': 2.0, 'dcyl': 8.0,
                  'lcyl': 150.0, 'ncrowd': 386, 'ensemble_id': 2},
                 id='HnsCyl'),
    pytest.param('trans_foci_cub',
                 'al6nl5ml216ns400ac1nc112776l27dt0.005bdump2000adump5000'
                 'ens1',
                 {'dmon_large': 6.0, 'nmon_large': 5, 'mmon_large': 216.0,
                  'nmon_small': 400, 'ncrowd': 112776, 'lcube': 54.0},
                 id='TransFociCub'),
    pytest.param('trans_foci_cyl',
                 'epss5epsl5r10.5al6nl5ml216ns400ac1nc0lz82.5dt0.005'
                 'bdump2000adump5000ens2',
                 {'dmon_large': 6.0, 'nmon_small': 400, 'ncrowd': 0,
                  'dcyl': 20.0, 'lcyl': 165.0, 'ensemble_id': 2},
                 id='TransFociCyl'),
])
def test_write_stamps_of_a_real_whole_carry_its_name(
        study_spec, tmp_path, study, whole, expected):
    """The stamps of a real whole hold the values its name encodes and the
    50 frames of its polymer dump."""
    spec = study_spec(study).SPEC
    parser = spec.parser(f"{whole}.ring", 'whole', spec.group)

    stamps = pd.read_csv(write_stamps(parser, 50, tmp_path))

    assert {key: stamps.loc[0, key] for key in expected} == expected
    # A whole's lineage name never carries the topology suffix ``.ring``
    # (ParserBase._set_parents).
    assert (stamps.loc[0, 'whole'], int(stamps.loc[0, 'n_frames'])) == (
        whole, 50)


@pytest.mark.private_data
def test_borrowed_topology_is_the_crowded_whole_own_polymer_topology(
        real_crowded_whole, real_files, private_window_files):
    """The polymer atoms of a crowded whole's own data file have the IDs,
    types, molecule IDs, masses, and bonds of the data file its polymer
    dump borrows."""
    _, own_path = private_window_files(real_crowded_whole)
    borrowed = read_data(real_files(real_crowded_whole).topology)
    own = read_data(own_path)
    n_polymer = borrowed.n_atoms
    in_polymer = np.all(own.bonds[:, 1:] <= n_polymer, axis=1)

    def as_set(bonds):
        # Each run's write_data numbers its bonds in its own order, so the
        # bonds are compared as a set: sorted by type and atoms.
        return bonds[np.lexsort((bonds[:, 2], bonds[:, 1], bonds[:, 0]))]

    assert_array_equal(own.ids[:n_polymer], borrowed.ids)
    assert_array_equal(own.types[:n_polymer], borrowed.types)
    assert_array_equal(own.mol_ids[:n_polymer], borrowed.mol_ids)
    assert_array_equal(own.masses[borrowed.types - 1],
                       borrowed.masses[borrowed.types - 1])
    assert_array_equal(as_set(own.bonds[in_polymer]),
                       as_set(borrowed.bonds))
