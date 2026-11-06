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
"""
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.manage.parser import SumRuleCyl
from polyphys.manage.projects import (
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
    again only when unwrapped: asphericity 1 with unwrap, below 1
    without, since the wrapped beads lie on two parallel segments."""
    first = np.full((10, 3), 10.0)
    first[:, 0] = np.arange(16.0, 26.0) % 20.0
    first[:, 1] = np.arange(5.0, 15.0)
    universe = make_chain_universe(first[np.newaxis], boundary='pp pp pp')

    unwrapped = ChainSize(universe.atoms, ('asphericityT',),
                          unwrap=True).run()
    wrapped = ChainSize(universe.atoms, ('asphericityT',),
                        unwrap=False).run()

    assert_allclose(unwrapped.results['asphericityT'], [1.0], rtol=RTOL)
    assert wrapped.results['asphericityT'][0] < 0.99


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


def test_from_spec_selects_the_chain_and_carries_the_spec(
        make_chain_universe, straight_chain):
    """from_spec picks the spec's chain selection and its observables,
    chain label, and unwrap."""
    project = ProjectSpec(
        parser=SumRuleCyl, group='bug', divisor=0.025,
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
