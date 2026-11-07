"""Tests for :mod:`polyphys.probe.chain_size`.

The inputs are tier A: chains of ten beads at hand-picked positions, written
as a data file and a dump by the ``make_chain_universe`` fixture. Expected
values are analytic: the end-to-end distance and the span of a straight
chain by hand, the radius of gyration of equally spaced beads from the
variance of a discrete uniform distribution, and the asphericity and shape
parameter of a cube's vertices (0, an isotropic gyration tensor) and of
collinear beads (1 and 2, eigenvalues (l, 0, 0)) from the formulas of Dima
and Thirumalai, J. Phys. Chem. B 108, 6564 (2004), as MDAnalysis
implements them.

The real-frame tests read the tier B and tier C files of the studies'
wholes (``polyphys/tests/probe/conftest.py``, "Real frames"). Their expected
values are properties or other files of the same run: the chain measured in
the ``all`` window at a time step the polymer dump shares, which holds the
same configuration; the FENE bound on a bond, :math:`R_0 = 1.5 \\sigma`
(``docs/source/design.md``), far below half the box, so a wrapped chain
whose bonds are all shorter than half the box is whole and unwrapping it
changes nothing; and the first level of the blocking transform of
Flyvbjerg and Petersen, J. Chem. Phys. 91, 461 (1989), whose standard
error is the naive one of the uncorrelated estimate.
"""
import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_equal,
    assert_equal,
)

from polyphys.analyze.stats import block_average
from polyphys.manage.parser import SumRuleCyl
from polyphys.spec import (
    CHAIN_SIZE_OBSERVABLES,
    ProbeSpec,
    ProjectSpec,
)

pytest.importorskip("MDAnalysis")

from polyphys.probe import ChainSize  # noqa: E402

# The dump reader warns that the frames carry no time step; the tests
# measure per frame and give no dt.
pytestmark = pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")

# Tolerances: MDAnalysis holds positions as float32, so a value carries a
# relative round-off of about 1e-7; 1e-6 leaves room for a few operations.
RTOL = 1e-6
ATOL = 1e-6


@pytest.fixture
def straight_chain() -> np.ndarray:
    """Five frames of ten beads one unit apart along x, from x = 5, moved
    0.1 along z per frame."""
    first = np.zeros((10, 3))
    first[:, 0] = np.arange(10) + 5.0
    first[:, 1:] = 10.0
    return np.stack([first + [0.0, 0.0, 0.1 * k] for k in range(5)])


@pytest.fixture
def cube_vertices() -> np.ndarray:
    """One frame of the eight vertices of a cube of side 2 about (10, 10,
    10), plus two beads at its center so the chain has ten beads."""
    corners = np.array([[sx, sy, sz] for sx in (-1.0, 1.0)
                        for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)])
    return (np.vstack([corners, np.zeros((2, 3))]) + 10.0)[np.newaxis]


# Hand values


def test_end_to_end_of_a_straight_chain(make_chain_universe, straight_chain):
    """Ten beads one unit apart along x span 9 from first to last, in every
    frame."""
    universe = make_chain_universe(straight_chain)

    analysis = ChainSize(universe.atoms, ('rfloryT',)).run()

    assert_allclose(analysis.results['rfloryT'], np.full(5, 9.0), rtol=RTOL)


def test_fsd_is_the_span_along_the_axis(make_chain_universe):
    """Beads at z = 10, 10.5, ..., 14.5 span 4.5 along z, unchanged when
    the chain moves."""
    first = np.full((10, 3), 10.0)
    first[:, 2] = 10.0 + 0.5 * np.arange(10)
    frames = np.stack([first, first + [1.0, -1.0, 2.0]])
    universe = make_chain_universe(frames)

    analysis = ChainSize(universe.atoms, ('fsdT',), axis=2).run()

    assert_allclose(analysis.results['fsdT'], [4.5, 4.5], rtol=RTOL)


def test_transverse_size_of_beads_on_a_circle(make_chain_universe):
    """Ten beads on a circle of radius 2 in the xy plane have transverse
    size 4 across z: twice the mean distance from their center."""
    angles = 2 * np.pi * np.arange(10) / 10
    first = np.full((10, 3), 10.0)
    first[:, 0] += 2.0 * np.cos(angles)
    first[:, 1] += 2.0 * np.sin(angles)
    universe = make_chain_universe(first[np.newaxis])

    analysis = ChainSize(universe.atoms, ('transSizeT',), axis=2).run()

    assert_allclose(analysis.results['transSizeT'], [4.0], rtol=RTOL)


def test_radius_of_gyration_of_equally_spaced_beads(make_chain_universe,
                                                    straight_chain):
    """N equal beads one unit apart on a line have Rg^2 = (N^2 - 1) / 12."""
    universe = make_chain_universe(straight_chain)

    analysis = ChainSize(universe.atoms, ('gyrT',)).run()

    assert_allclose(analysis.results['gyrT'], np.full(5, np.sqrt(99 / 12)),
                    rtol=RTOL)


# Property tests of the gyration-tensor shapes


def test_asphericity_is_zero_for_a_cube_and_one_for_a_line(
        make_chain_universe, cube_vertices, straight_chain):
    """A cube's vertices have an isotropic gyration tensor, asphericity 0;
    collinear beads have eigenvalues (l, 0, 0), asphericity 1."""
    cube = make_chain_universe(cube_vertices, name='cube')
    line = make_chain_universe(straight_chain, name='line')

    on_cube = ChainSize(cube.atoms, ('asphericityT',)).run()
    on_line = ChainSize(line.atoms, ('asphericityT',)).run()

    assert_allclose(on_cube.results['asphericityT'], [0.0], atol=ATOL)
    assert_allclose(on_line.results['asphericityT'], np.ones(5), rtol=RTOL)


def test_shape_parameter_is_zero_for_a_cube_and_two_for_a_line(
        make_chain_universe, cube_vertices, straight_chain):
    """The shape parameter 27 prod(l_k - mean) / (sum l_k)^3 is 0 for an
    isotropic tensor and 2 for eigenvalues (l, 0, 0)."""
    cube = make_chain_universe(cube_vertices, name='cube')
    line = make_chain_universe(straight_chain, name='line')

    on_cube = ChainSize(cube.atoms, ('shapeT',)).run()
    on_line = ChainSize(line.atoms, ('shapeT',)).run()

    assert_allclose(on_cube.results['shapeT'], [0.0], atol=ATOL)
    assert_allclose(on_line.results['shapeT'], np.full(5, 2.0), rtol=RTOL)


def test_principal_axes_are_orthonormal_with_the_line_last(
        make_chain_universe, straight_chain):
    """A straight chain along x has a (3, 3) orthonormal frame per frame
    whose last axis, the smallest moment of inertia, is along x."""
    universe = make_chain_universe(straight_chain)

    axes = ChainSize(universe.atoms, ('principalT',)).run().results[
        'principalT']

    assert axes.shape == (5, 3, 3)
    assert_allclose(axes @ axes.transpose(0, 2, 1),
                    np.broadcast_to(np.eye(3), (5, 3, 3)), atol=ATOL)
    assert_allclose(np.abs(axes[:, 2, 0]), np.ones(5), rtol=RTOL)


def test_unwrap_reaches_the_asphericity(make_chain_universe):
    """A straight diagonal chain wrapped across a periodic face is straight
    again when unwrapped, asphericity 1, and the whole-chain check leaves
    an unwrapped asphericity alone."""
    first = np.full((10, 3), 10.0)
    first[:, 0] = np.arange(16.0, 26.0) % 20.0
    first[:, 1] = np.arange(5.0, 15.0)
    universe = make_chain_universe(first[np.newaxis], boundary='pp pp pp')

    unwrapped = ChainSize(universe.atoms, ('asphericityT',),
                          unwrap=True).run()

    assert_allclose(unwrapped.results['asphericityT'], [1.0], rtol=RTOL)


@pytest.fixture
def split_chain() -> np.ndarray:
    """Two frames of a four-bead chain with unit bonds along x in a box of
    10: whole at x = 1.5 to 4.5, then wrapped at x = 8.5, 9.5, 0.5, 1.5,
    its third bond split by the boundary."""
    frames = np.full((2, 4, 3), 5.0)
    frames[0, :, 0] = [1.5, 2.5, 3.5, 4.5]
    frames[1, :, 0] = [8.5, 9.5, 0.5, 1.5]
    return frames


@pytest.mark.parametrize('observables, unwrap', [
    pytest.param(('gyrT',), True, id='gyration-with-unwrap'),
    pytest.param(('asphericityT',), False, id='asphericity-without-unwrap'),
    pytest.param(('asphericityT', 'shapeT'), True,
                 id='shape-beside-unwrapped-asphericity'),
])
def test_chain_split_by_the_boundary_raises(make_chain_universe,
                                            split_chain, observables,
                                            unwrap):
    """A frame whose chain bond spans half the box or more stops the run,
    naming the frame, for every observable that reads the wrapped
    positions."""
    universe = make_chain_universe(split_chain, box=((0.0, 10.0),) * 3,
                                   boundary='pp pp pp')
    analysis = ChainSize(universe.atoms, observables, unwrap=unwrap)

    with pytest.raises(ValueError,
                       match="frame 1: 1 of 3 chain bonds span half the "
                             "box or more along x"):
        analysis.run()


def test_whole_chain_past_the_box_edge_gives_its_radius_of_gyration(
        make_chain_universe):
    """A whole chain that reaches past the box edge passes the check and
    gives the analytic radius of gyration of four unit-spaced beads,
    sqrt(5) / 2."""
    frames = np.full((1, 4, 3), 5.0)
    frames[0, :, 0] = [8.5, 9.5, 10.5, 11.5]
    universe = make_chain_universe(frames, box=((0.0, 10.0),) * 3,
                                   boundary='pp pp pp')

    analysis = ChainSize(universe.atoms, ('gyrT',)).run()

    assert_allclose(analysis.results['gyrT'], [np.sqrt(5.0) / 2.0],
                    rtol=RTOL)


def test_chain_without_bonds_is_not_checked():
    """Without bonds in the topology the check cannot run, and the wrapped
    end-to-end distance is measured as it stands, 7 for beads at
    x = 8.5 and 1.5."""
    import MDAnalysis as mda

    universe = mda.Universe.empty(4, trajectory=True)
    universe.add_TopologyAttr('masses', np.ones(4))
    universe.dimensions = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    universe.atoms.positions = np.array(
        [[8.5, 5.0, 5.0], [9.5, 5.0, 5.0], [0.5, 5.0, 5.0], [1.5, 5.0, 5.0]])

    analysis = ChainSize(universe.atoms, ('rfloryT',)).run()

    assert_allclose(analysis.results['rfloryT'], [7.0], rtol=RTOL)


# The analysis contract


def test_every_declared_observable_is_measured(make_chain_universe,
                                               straight_chain):
    """Every name in CHAIN_SIZE_OBSERVABLES gives one array per frame."""
    universe = make_chain_universe(straight_chain)

    analysis = ChainSize(universe.atoms, CHAIN_SIZE_OBSERVABLES).run()

    assert tuple(analysis.results) == CHAIN_SIZE_OBSERVABLES
    assert {analysis.results[obs].shape[0]
            for obs in CHAIN_SIZE_OBSERVABLES} == {5}


def test_run_slices_the_frames(make_chain_universe, straight_chain):
    """run(start, stop, step) measures only those frames, in order."""
    universe = make_chain_universe(straight_chain)

    analysis = ChainSize(universe.atoms, ('fsdT',)).run(
        start=1, stop=5, step=2)

    assert_array_equal(analysis.frames, [1, 3])
    # The chain's z span is 0: its ten beads share one z in every frame.
    assert_allclose(analysis.results['fsdT'], [0.0, 0.0], atol=ATOL)
    assert analysis.n_frames == 2


def test_parallel_run_matches_the_serial_run(make_chain_universe,
                                             straight_chain):
    """Two workers over two parts give the serial arrays, scalar and
    matrix series alike."""
    universe = make_chain_universe(straight_chain)
    observables = ('rfloryT', 'gyrT', 'principalT')
    serial = ChainSize(universe.atoms, observables).run()

    parallel = ChainSize(universe.atoms, observables).run(
        backend='multiprocessing', n_workers=2, n_parts=2)

    for obs in observables:
        assert_array_equal(parallel.results[obs], serial.results[obs])


@pytest.mark.parametrize('n_parts', [
    pytest.param(2, id='two-parts'),
    pytest.param(3, id='three-parts'),
])
def test_serial_parts_match_the_multiprocessing_run(
        make_chain_universe, straight_chain, n_parts):
    """The serial backend split into parts gives the frames, times, and
    arrays of the multiprocessing backend over the same parts: the result
    does not depend on the backend."""
    universe = make_chain_universe(straight_chain)
    observables = ('rfloryT', 'gyrT', 'principalT')
    parallel = ChainSize(universe.atoms, observables).run(
        backend='multiprocessing', n_workers=2, n_parts=n_parts)

    serial = ChainSize(universe.atoms, observables).run(
        backend='serial', n_parts=n_parts)

    assert_array_equal(serial.frames, parallel.frames)
    assert_array_equal(serial.times, parallel.times)
    assert_equal(dict(serial.results), dict(parallel.results))


def test_from_spec_selects_the_chain_and_carries_the_spec(
        make_chain_universe, straight_chain):
    """from_spec picks the spec's chain selection and its observables,
    chain label, and unwrap."""
    project = ProjectSpec(
        parser=SumRuleCyl, group='bug', divisor=0.025,
        crowder_attr='phi_bulk_c',
        physical_attrs=('space', 'nmon'),
        species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
        tseries_properties=('fsdTMon', 'gyrTMon'))
    spec = ProbeSpec(project, {'Mon': 'resid 1', 'Crd': 'resid 0'}, 'Mon',
                     ('fsdT', 'gyrT'), unwrap=True)
    universe = make_chain_universe(straight_chain)

    analysis = ChainSize.from_spec(universe, spec)

    assert (analysis.atomgroup.n_atoms, analysis.observables,
            analysis.chain, analysis.unwrap) == (
        10, ('fsdT', 'gyrT'), 'Mon', True)


@pytest.mark.parametrize('kwargs, match', [
    pytest.param({'observables': ()}, "must not be empty",
                 id='no-observables'),
    pytest.param({'observables': ('gyrT', 'radius')},
                 r"names \['radius'\], which are not in",
                 id='unknown-observable'),
    pytest.param({'observables': ('gyrT', 'gyrT')}, "repeats a name",
                 id='repeated-observable'),
    pytest.param({'chain': ''}, "'chain' must be a non-empty",
                 id='empty-chain'),
    pytest.param({'axis': 3}, "'axis' must be 0, 1, or 2", id='bad-axis'),
])
def test_rejects_invalid_arguments(make_chain_universe, straight_chain,
                                   kwargs, match):
    """An empty, unknown, or repeated observable, an empty chain label, or
    an axis outside 0 to 2 raises at construction."""
    universe = make_chain_universe(straight_chain)
    fields = {'observables': ('gyrT',), **kwargs}

    with pytest.raises(ValueError, match=match):
        ChainSize(universe.atoms, **fields)


def test_rejects_a_chain_of_one_bead(make_chain_universe, straight_chain):
    """A chain needs at least two beads."""
    universe = make_chain_universe(straight_chain)

    with pytest.raises(ValueError, match="at least two beads"):
        ChainSize(universe.select_atoms('id 1'), ('gyrT',))


# save


def test_save_writes_one_file_per_observable(make_chain_universe,
                                             straight_chain, tmp_path):
    """save writes <name>-<observable><chain>.npy per observable, holding
    the results."""
    universe = make_chain_universe(straight_chain)
    analysis = ChainSize(universe.atoms, ('rfloryT', 'principalT'),
                         chain='Mon').run()

    paths = analysis.save(tmp_path, 'N10ens1-bug')

    assert [path.name for path in paths] == [
        'N10ens1-bug-rfloryTMon.npy', 'N10ens1-bug-principalTMon.npy']
    assert_array_equal(np.load(paths[0]), analysis.results['rfloryT'])
    assert_array_equal(np.load(paths[1]), analysis.results['principalT'])


def test_save_before_run_raises(make_chain_universe, straight_chain,
                                tmp_path):
    """save needs the results of a run."""
    universe = make_chain_universe(straight_chain)
    analysis = ChainSize(universe.atoms, ('rfloryT',))

    with pytest.raises(RuntimeError, match=r"Call run\(\) before save"):
        analysis.save(tmp_path, 'N10ens1-bug')


# Real frames


def _at_shared_steps(polymer, window, probe, trajectory_arrays):
    """The spec's observables of the chain in the polymer dump and in the
    window, at the time steps both hold, each flattened and joined in the
    spec's order; and the number of those steps."""
    polymer_steps = trajectory_arrays(polymer)['steps']
    window_steps = trajectory_arrays(window)['steps']
    shared = np.intersect1d(polymer_steps, window_steps)
    values = [
        ChainSize.from_spec(universe, probe).run(
            frames=np.flatnonzero(np.isin(steps, shared))).results
        for universe, steps in ((polymer, polymer_steps),
                                (window, window_steps))]
    return ([np.concatenate([result[obs].ravel()
                             for obs in probe.observables])
             for result in values], shared.size)


def test_real_chain_size_agrees_between_the_polymer_dump_and_the_window(
        real_window_whole, real_files, real_universe, study_spec,
        trajectory_arrays):
    """At every time step the polymer dump and the ``all`` window share,
    each chain-size observable of the spec is the same in both."""
    probe = study_spec(real_window_whole.study).PROBE
    files = real_files(real_window_whole)
    polymer = real_universe(files.topology, files.polymer)
    window = real_universe(files.data, files.window)

    (from_polymer, from_window), n_shared = _at_shared_steps(
        polymer, window, probe, trajectory_arrays)

    assert n_shared == real_window_whole.n_shared_frames
    assert_allclose(from_polymer, from_window, rtol=RTOL, atol=ATOL)


@pytest.mark.private_data
def test_real_chain_size_agrees_with_the_crowded_window(
        real_crowded_whole, real_files, private_window_files, real_universe,
        study_spec, trajectory_arrays):
    """In a crowded whole, the polymer dump opened with the borrowed
    topology and the ``all`` window opened with the whole's own data file
    give the same chain-size observables at every shared time step."""
    probe = study_spec(real_crowded_whole.study).PROBE
    files = real_files(real_crowded_whole)
    window_path, data_path = private_window_files(real_crowded_whole)
    polymer = real_universe(files.topology, files.polymer)
    window = real_universe(data_path, window_path)

    (from_polymer, from_window), n_shared = _at_shared_steps(
        polymer, window, probe, trajectory_arrays)

    assert n_shared == real_crowded_whole.n_shared_frames
    assert_allclose(from_polymer, from_window, rtol=RTOL, atol=ATOL)


def test_real_chain_bonds_are_shorter_than_half_the_periodic_box(
        real_polymer_whole, real_files, real_universe, study_spec,
        trajectory_arrays):
    """In every frame of a polymer dump, every bond of the chain spans less
    than half the box along each periodic direction: the wrapping splits no
    chain. (Free H-NS proteins do cross the boundary of the periodic cube;
    their core-patch bonds are not checked.)"""
    probe = study_spec(real_polymer_whole.study).PROBE
    files = real_files(real_polymer_whole)
    universe = real_universe(files.topology, files.polymer)
    arrays = trajectory_arrays(universe)
    first, second = universe.select_atoms(
        probe.selections[probe.chain]).bonds.indices.T
    periodic = {'cubic': [0, 1, 2], 'cylindrical': [2]}[probe.geometry]

    spans = np.abs(arrays['positions'][:, first] -
                   arrays['positions'][:, second])[..., periodic]
    halves = 0.5 * arrays['boxes'][:, None, periodic]

    assert np.count_nonzero(spans >= halves) == 0


def test_real_asphericity_does_not_depend_on_unwrap(
        real_polymer_whole, real_files, real_universe, study_spec):
    """On real chains that no boundary splits, the asphericity is the same
    whether or not the chain is unwrapped first, in every study."""
    probe = study_spec(real_polymer_whole.study).PROBE
    files = real_files(real_polymer_whole)
    chain = real_universe(files.topology, files.polymer).select_atoms(
        probe.selections[probe.chain])

    wrapped, unwrapped = (
        ChainSize(chain, ('asphericityT',), unwrap=unwrap).run().results[
            'asphericityT'] for unwrap in (False, True))

    assert_allclose(unwrapped, wrapped, rtol=RTOL, atol=ATOL)


def test_block_average_of_a_real_gyration_series_starts_at_the_naive_error(
        real_polymer_whole, real_files, real_universe, study_spec):
    """Block averaging runs on the 50-frame radius-of-gyration series of a
    real chain, with finite errors at every level and the naive standard
    error at the first."""
    probe = study_spec(real_polymer_whole.study).PROBE
    files = real_files(real_polymer_whole)
    chain = real_universe(files.topology, files.polymer).select_atoms(
        probe.selections[probe.chain])
    series = ChainSize(chain, ('gyrT',)).run().results['gyrT']

    result = block_average(series)

    assert np.all(np.isfinite(result.sem))
    # Both are float64 reductions of the same 50 values.
    assert result.sem[0] == pytest.approx(
        series.std(ddof=1) / np.sqrt(series.size), rel=1e-12)
