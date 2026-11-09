"""Tests for :mod:`polyphys.probe.histograms`.

The inputs are tier A: two species of hand-placed beads in a box of side 20
about (10, 10, 10), written as a data file and a dump by the
``make_chain_universe`` fixture, binned with hand edges of unit or
half-unit width. Expected values are analytic: the bin of each bead by
hand from its recentered coordinates; the conservation of each species'
particle number per frame; the per-block means from the block sums
divided by the frame counts; and the edges of ``histogram_edges`` from the
research code's bin sizes over the ranges the maintainer set on
2026-11-04 (the research probe's ``trans_foci_cyl_all`` and
``trans_foci_cub_all``, extended to the wall surface and the cube's
corner).

The real-frame tests bin the tier B and tier C frames of the studies'
wholes (``polyphys/tests/probe/conftest.py``, "Real frames") with each
study's own probe spec and parser. Their expected values are conservation:
every position of a real frame falls inside the edges, so each frame's
counts sum to the species' particle number in the whole's name, and the run's
counts to the number of frames times it.
"""
import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_equal,
    assert_equal,
)

from polyphys.tests.synthetic import ChainCyl, RingCub
from polyphys.spec import ProbeSpec, ProjectSpec

pytest.importorskip("MDAnalysis")

from polyphys.probe import SpatialHistograms, histogram_edges  # noqa: E402

# The dump reader warns that the frames carry no time step; the tests
# measure per frame and give no dt.
pytestmark = pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")

# Tolerance for the per-block means and the edges, both a few float64
# operations on exact values.
RTOL = 1e-12
ATOL = 1e-12


@pytest.fixture
def center() -> float:
    """The box center along each axis, for a box from 0 to 20."""
    return 10.0


@pytest.fixture
def two_species(center) -> np.ndarray:
    """Five frames of six beads: beads 1-3 (type 1, ``Mon``) and 4-6 (type 2,
    ``Crd``), the first of each species moved 0.5 along z per frame.

    About the center ``c``: bead 1 at ``(c + 2.5, c, c)``, bead 2 at
    ``(c, c + 1, c + 3.2)``, bead 3 at ``(c - 1.5, c - 1.5, c - 2)``, bead 4
    at ``(c + 0.5, c + 0.5, c + 0.5)``, bead 5 at ``(c - 3, c, c - 7)``, and
    bead 6 at ``(c, c - 2.5, c + 6)``.
    """
    first = np.array([
        [2.5, 0.0, 0.0], [0.0, 1.0, 3.2], [-1.5, -1.5, -2.0],
        [0.5, 0.5, 0.5], [-3.0, 0.0, -7.0], [0.0, -2.5, 6.0]]) + center
    frames = np.stack([first.copy() for _ in range(5)])
    frames[:, 0, 2] += 0.5 * np.arange(5)
    frames[:, 3, 2] += 0.5 * np.arange(5)
    return frames


@pytest.fixture
def types() -> list[int]:
    """The atom type of each of the six beads."""
    return [1, 1, 1, 2, 2, 2]


@pytest.fixture
def cylinder_edges() -> dict[str, np.ndarray]:
    """Hand edges of a cylinder of radius 4 and length 20: unit bins along
    r, x, and y, half-unit bins along z, 72 bins along theta."""
    return {
        'r': np.arange(0.0, 4.5, 1.0),
        'z': np.linspace(-10.0, 10.0, 41),
        'theta': np.linspace(-np.pi, np.pi, 73),
        'x': np.arange(-4.0, 4.5, 1.0),
        'y': np.arange(-4.0, 4.5, 1.0),
    }


@pytest.fixture
def cube_edges() -> dict[str, np.ndarray]:
    """Hand edges of a cube of side 20: unit bins along r to the corner and
    along x, y, and z."""
    return {
        'r': np.arange(0.0, 18.5, 1.0),
        'x': np.arange(-10.0, 10.5, 1.0),
        'y': np.arange(-10.0, 10.5, 1.0),
        'z': np.arange(-10.0, 10.5, 1.0),
    }


@pytest.fixture
def make_analysis(make_chain_universe, two_species, types, cylinder_edges):
    """Factory: the cylinder analysis of the two species, with any argument
    overridden; `positions` replaces the frames."""

    def _make(positions=None, **overrides):
        universe = make_chain_universe(
            two_species if positions is None else positions,
            atom_types=types, masses=(1.0, 1.0))
        groups = {'Mon': universe.select_atoms('type 1'),
                  'Crd': universe.select_atoms('type 2')}
        kwargs = {'geometry': 'cylindrical', 'block_frames': 2,
                  'recenter': 'half'}
        kwargs.update(overrides)
        return SpatialHistograms(groups, cylinder_edges, **kwargs)

    return _make


# Hand bins


@pytest.mark.parametrize('key, filled', [
    pytest.param('r', {1: 1, 2: 2}, id='r-1-2.5-and-2.12'),
    pytest.param('theta', {36: 1, 54: 1, 9: 1},
                 id='theta-0-pi-over-2-and-minus-3pi-over-4'),
    pytest.param('z', {20: 1, 26: 1, 16: 1}, id='z-0-3.2-and-minus-2'),
])
def test_histograms_bin_each_bead_by_its_recentered_coordinate(
        make_analysis, key, filled):
    """In frame 0 the three Mon beads, at (2.5, 0, 0), (0, 1, 3.2), and
    (-1.5, -1.5, -2) from the center, fill the r, theta, and z bins
    computed by hand from the edges and nothing else."""
    analysis = make_analysis().run()

    row = analysis.results[f'{key}HistTMon'][0]

    expected = np.zeros(row.size, dtype=np.int64)
    for index, count in filled.items():
        expected[index] = count
    assert_array_equal(row, expected)


def test_histogram_t_tracks_a_moving_bead(make_analysis):
    """A bead moved 0.5 along z per frame moves one z bin per frame, so the
    per-frame rows differ as the bead does."""
    analysis = make_analysis().run()

    rows = analysis.results['zHistTMon']

    assert_array_equal(rows[:, 20:25].diagonal(), np.ones(5, dtype=int))


@pytest.mark.parametrize('plane, cell, mean', [
    pytest.param('xy', (6, 4), 1.0, id='xy-of-(2.5,0)'),
    pytest.param('xz', (6, 20), 0.5, id='xz-of-(2.5,0)'),
    pytest.param('yz', (4, 20), 0.5, id='yz-of-(2.5,0)'),
])
def test_maps_put_a_bead_in_its_cell_with_rows_along_the_first_axis(
        make_analysis, plane, cell, mean):
    """Bead 1 at (2.5, 0, 0) falls in the cell indexed by its first-axis
    coordinate then its second, in the first block's mean map: in both
    frames of the block in the xy plane, in frame 0 only in the planes with
    z, where it moves from the half-unit bin [0, 0.5) to the next."""
    analysis = make_analysis().run()

    block = analysis.results[f'{plane}HistBlockMon'][0]

    assert block[cell] == pytest.approx(mean, abs=ATOL)


# Conservation and consistency


@pytest.mark.parametrize('key', ['r', 'z', 'theta'])
def test_each_frame_conserves_the_particle_number(make_analysis, key):
    """Every frame's counts along a direction sum to the species' three
    beads, for both species."""
    analysis = make_analysis().run()

    assert_array_equal(analysis.results[f'{key}HistTMon'].sum(axis=1),
                       np.full(5, 3))
    assert_array_equal(analysis.results[f'{key}HistTCrd'].sum(axis=1),
                       np.full(5, 3))


@pytest.mark.parametrize('key', ['r', 'z', 'theta'])
def test_summed_histogram_is_the_column_sum_of_the_frames(make_analysis, key):
    """The frame-summed histogram equals the per-frame rows summed down the
    columns, int64."""
    analysis = make_analysis().run()

    summed = analysis.results[f'{key}HistCrd']

    assert summed.dtype == np.int64
    assert_array_equal(summed, analysis.results[f'{key}HistTCrd'].sum(axis=0))


@pytest.mark.parametrize('plane', ['xy', 'xz', 'yz'])
def test_block_means_times_block_frames_give_the_summed_map(make_analysis,
                                                            plane):
    """Each block's mean map times its frame count, summed over the blocks,
    is the frame-summed map, whose total is the beads times the frames."""
    analysis = make_analysis().run()

    means = analysis.results[f'{plane}HistBlockMon']
    counts = analysis.results['blockFrames']

    assert_array_equal(counts, [2, 2, 1])
    assert_allclose((means * counts[:, None, None]).sum(axis=0),
                    analysis.results[f'{plane}HistMon'], rtol=RTOL, atol=ATOL)
    assert analysis.results[f'{plane}HistMon'].sum() == pytest.approx(
        15.0, abs=ATOL)


def test_blocks_follow_the_frames_of_the_run(make_analysis):
    """A run over four of the five frames with blocks of two gives two full
    blocks, numbered over the frames visited."""
    analysis = make_analysis().run(start=1, stop=5)

    assert_array_equal(analysis.results['blockFrames'], [2, 2])
    assert analysis.results['xyHistBlockMon'].shape == (2, 8, 8)


def test_parallel_run_matches_serial(make_analysis):
    """Two workers over two parts, the split falling inside a block, give
    the same per-frame rows, block means, and frame counts as one worker."""
    serial = make_analysis().run()
    parallel = make_analysis().run(backend='multiprocessing', n_workers=2,
                                   n_parts=2)

    for key in serial.results:
        assert_allclose(parallel.results[key], serial.results[key],
                        rtol=RTOL, atol=ATOL)


@pytest.mark.parametrize('n_parts', [
    pytest.param(2, id='two-parts'),
    pytest.param(3, id='three-parts'),
])
def test_serial_parts_match_the_multiprocessing_run(make_analysis, n_parts):
    """The serial backend split into parts gives the frames, times,
    per-frame rows, block means, and frame counts of the multiprocessing
    backend over the same parts: the result does not depend on the
    backend."""
    parallel = make_analysis().run(backend='multiprocessing', n_workers=2,
                                   n_parts=n_parts)

    serial = make_analysis().run(backend='serial', n_parts=n_parts)

    assert_array_equal(serial.frames, parallel.frames)
    assert_array_equal(serial.times, parallel.times)
    assert_equal(dict(serial.results), dict(parallel.results))


def test_empty_group_gives_zero_histograms(make_chain_universe, two_species,
                                           types, cylinder_edges):
    """A species a whole lacks yields all-zero rows of the right shape."""
    universe = make_chain_universe(two_species, atom_types=types,
                                   masses=(1.0, 1.0))
    groups = {'Mon': universe.select_atoms('type 1'),
              'Hns': universe.select_atoms('type 3')}

    analysis = SpatialHistograms(groups, cylinder_edges,
                                 geometry='cylindrical', block_frames=2).run()

    assert_array_equal(analysis.results['rHistTHns'], np.zeros((5, 4), int))
    assert_array_equal(analysis.results['rHistTMon'].sum(axis=1),
                       np.full(5, 3))


# Recentering


@pytest.mark.parametrize('recenter', [
    pytest.param('half', id='half-the-box'),
    pytest.param((10.0, 10.0, 10.0), id='the-same-triple'),
])
def test_recenter_half_and_its_triple_agree(make_analysis, recenter):
    """Half of a box of side 20 and the triple (10, 10, 10) give the same
    bins: the Mon r row [0, 1, 2, 0]."""
    analysis = make_analysis(recenter=recenter).run()

    assert_array_equal(analysis.results['rHistTMon'][0], [0, 1, 2, 0])


def test_recenter_half_undoes_the_readers_translation(
        make_chain_universe, cylinder_edges):
    """For a box from -10 to 10, which the dump reader moves to the origin,
    'half' restores the written coordinates: a bead written at x = 2.5 is
    in r bin 2 and x bin 6."""
    positions = np.array([[[2.5, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    universe = make_chain_universe(positions, box=((-10.0, 10.0),) * 3)

    analysis = SpatialHistograms(
        {'Mon': universe.atoms}, cylinder_edges, geometry='cylindrical',
        block_frames=1, recenter='half').run()

    assert analysis.results['rHistTMon'][0, 2] == 1
    assert analysis.results['xyHistBlockMon'][0, 6, 4] == pytest.approx(
        1.0, abs=ATOL)


def test_recenter_none_keeps_the_readers_frame(make_chain_universe,
                                               cube_edges):
    """With 'none', positions stay in the reader's frame, from 0 to 20, so
    a bead written at (12.5, 10, 10) lies beyond cube edges centered on the
    origin and the run raises instead of binning it."""
    positions = np.array([[[12.5, 10.0, 10.0], [10.0, 10.0, 10.0]]])
    universe = make_chain_universe(positions)

    analysis = SpatialHistograms(
        {'Mon': universe.atoms}, cube_edges, geometry='cubic',
        block_frames=1, recenter='none')

    with pytest.raises(ValueError, match="positions fall outside"):
        analysis.run()


# Out-of-range positions


def test_position_outside_a_direction_raises(make_analysis, two_species):
    """A bead beyond the last r edge aborts the run, naming the count."""
    positions = two_species.copy()
    positions[0, 0, 0] = 10.0 + 4.5  # r = 4.5 > 4

    with pytest.raises(ValueError, match="1 of 3 positions fall outside"):
        make_analysis(positions=positions).run()


def test_position_outside_the_cartesian_edges_raises(
        make_chain_universe, two_species, types, cylinder_edges):
    """A bead inside every one-dimensional range but beyond a Cartesian
    edge of the maps aborts the run, naming the axis and the frame."""
    edges = dict(cylinder_edges)
    edges['x'] = np.arange(-2.0, 2.5, 1.0)  # bead 1 at x = 2.5 is beyond
    universe = make_chain_universe(two_species, atom_types=types,
                                   masses=(1.0, 1.0))

    analysis = SpatialHistograms(
        {'Mon': universe.select_atoms('type 1')}, edges,
        geometry='cylindrical', block_frames=2)

    with pytest.raises(ValueError,
                       match="1 of 3 Mon positions fall outside the x edges"):
        analysis.run()


@pytest.mark.parametrize('position, axis, expected_bin', [
    pytest.param([20.03, 10.0, 10.0], 'x', 0, id='beyond-hi-x'),
    pytest.param([10.0, -0.02, 10.0], 'y', 19, id='below-lo-y'),
    pytest.param([10.0, 10.0, 20.04], 'z', 0, id='beyond-hi-z'),
])
def test_cube_bins_a_position_past_the_box_as_its_periodic_image(
        make_chain_universe, cube_edges, position, axis, expected_bin):
    """A position a little past a periodic face, as a LAMMPS dump writes
    between reneighborings, is binned as its periodic image."""
    positions = np.array([[position, [10.0, 10.0, 10.0]]])
    universe = make_chain_universe(positions, boundary='pp pp pp')

    analysis = SpatialHistograms(
        {'Mon': universe.atoms}, cube_edges, geometry='cubic',
        block_frames=1).run()

    # One bead at the center and one at the periodic image: x = 0.03 and
    # z = 0.04 fall in the first unit bin, y = 19.98 in the last.
    plane = {'x': 'xyHistMon', 'y': 'xyHistMon', 'z': 'xzHistMon'}[axis]
    marginal = analysis.results[plane].sum(axis={'x': 1, 'y': 0,
                                                 'z': 0}[axis])
    assert marginal[expected_bin] == 1
    assert int(analysis.results['rHistTMon'].sum()) == 2


def test_frame_without_a_box_raises(cube_edges):
    """A frame that carries no box cannot be wrapped, and the run says
    so."""
    import MDAnalysis as mda

    universe = mda.Universe.empty(1, trajectory=True)
    universe.atoms.positions = np.array([[1.0, 2.0, 3.0]])
    analysis = SpatialHistograms({'Mon': universe.atoms}, cube_edges,
                                 geometry='cubic', block_frames=1,
                                 recenter='none')

    with pytest.raises(ValueError, match="carries no box"):
        analysis.run()


def test_cylinder_bins_a_position_past_the_axis_end_as_its_periodic_image(
        make_analysis, two_species):
    """A position a little past the end of the periodic cylinder axis is
    binned as its periodic image, while the walls still bound x and y."""
    positions = two_species.copy()
    positions[0, 1, 2] = 20.03  # z past the box [0, 20]: image at 0.03

    analysis = make_analysis(positions=positions).run()

    # Recentered by half the box, the image is at z = -9.97, in the first
    # half-unit bin of [-10, 10].
    assert analysis.results['zHistTMon'][0, 0] == 1
    assert int(analysis.results['zHistTMon'][0].sum()) == 3


# Cube geometry


def test_cube_bins_r_from_the_box_center(make_chain_universe, cube_edges):
    """In a cube, r is the distance from the box center: a bead at
    (c + 3, c + 4, c) is in r bin 5, and only r is binned."""
    positions = np.array([[[13.0, 14.0, 10.0], [10.0, 10.0, 10.0]]])
    universe = make_chain_universe(positions)

    analysis = SpatialHistograms(
        {'Mon': universe.atoms}, cube_edges, geometry='cubic',
        block_frames=1).run()

    assert_array_equal(analysis.results['rHistTMon'][0, [0, 5]], [1, 1])
    assert 'zHistTMon' not in analysis.results


# Validation


@pytest.mark.parametrize('overrides, match', [
    pytest.param({'geometry': 'spherical'}, "'geometry' must be one of",
                 id='unknown-geometry'),
    pytest.param({'block_frames': 0}, "'block_frames' must be a positive",
                 id='zero-block-frames'),
    pytest.param({'block_frames': 2.0}, "'block_frames' must be a positive",
                 id='float-block-frames'),
    pytest.param({'recenter': 'center'}, "'recenter' must be",
                 id='unknown-recenter'),
    pytest.param({'recenter': (1.0, 2.0)}, "'recenter' must be",
                 id='short-recenter'),
])
def test_invalid_arguments_raise(make_analysis, overrides, match):
    """A bad geometry, block length, or recentering raises at construction."""
    with pytest.raises(ValueError, match=match):
        make_analysis(**overrides)


@pytest.mark.parametrize('edit, match', [
    pytest.param(lambda edges: edges.pop('theta'), r"lacks \['theta'\]",
                 id='missing-direction'),
    pytest.param(lambda edges: edges.pop('x'), r"lacks \['x'\]",
                 id='missing-cartesian-axis'),
    pytest.param(lambda edges: edges.update(r=np.array([0.0, 2.0, 1.0])),
                 "edges along 'r' must be", id='decreasing-edges'),
    pytest.param(lambda edges: edges.update(r=np.array([0.0])),
                 "edges along 'r' must be", id='single-edge'),
])
def test_invalid_edges_raise(make_chain_universe, two_species, cylinder_edges,
                             edit, match):
    """Edges that lack a key the geometry needs, or are not increasing with
    at least two values, raise at construction."""
    edges = dict(cylinder_edges)
    edit(edges)
    universe = make_chain_universe(two_species)

    with pytest.raises(ValueError, match=match):
        SpatialHistograms({'Mon': universe.atoms}, edges,
                          geometry='cylindrical', block_frames=2)


def test_empty_groups_raise(make_chain_universe, two_species, cylinder_edges):
    """No species at all raises."""
    make_chain_universe(two_species)

    with pytest.raises(ValueError, match="'groups' must map"):
        SpatialHistograms({}, cylinder_edges, geometry='cylindrical',
                          block_frames=2)


def test_groups_from_two_universes_raise(make_chain_universe, two_species,
                                         cylinder_edges):
    """Species from different universes cannot share one trajectory."""
    first = make_chain_universe(two_species, name='first')
    second = make_chain_universe(two_species, name='second')

    with pytest.raises(ValueError, match="one universe"):
        SpatialHistograms({'Mon': first.atoms, 'Crd': second.atoms},
                          cylinder_edges, geometry='cylindrical',
                          block_frames=2)


# histogram_edges and from_spec


@pytest.fixture
def cylinder_spec():
    """A ChainCyl probe spec with the research code's bin sizes: 0.1 of
    the smaller diameter along r, x, and y, half of it along z, pi/36 along
    theta."""
    project = ProjectSpec(
        parser=ChainCyl, group='bug', divisor=0.025,
        crowder_attr='phi_bulk_c',
        physical_attrs=('space', 'nmon', 'dcyl', 'dcrowd'),
        species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
        tseries_properties=('gyrTMon',))

    def tenth(info):
        return 0.1 * min(info.dmon, info.dcrowd)

    def half(info):
        return 0.5 * min(info.dmon, info.dcrowd)

    return ProbeSpec(
        project, {'Mon': 'type 1', 'Crd': 'type 2'}, 'Mon', ('gyrT',),
        unwrap=False, histogram_species=('Mon', 'Crd'), block_frames=2,
        bin_sizes={'r': tenth, 'z': half, 'theta': np.pi / 36, 'x': tenth,
                   'y': tenth})


def test_histogram_edges_of_a_cylinder(cylinder_spec):
    """For dcyl = 4 and lcyl = 20 with unit beads: r has 25 bins of 0.1 to
    the wall surface 2.5, z 40 bins of 0.5 over [-10, 10], theta 72 bins
    over [-pi, pi], and x and y 50 bins of 0.1 over [-2.5, 2.5]."""
    whole = ChainCyl(
        'N10r2.5lz10ac1nc0dt0.005bdump1adump1ens1', 'whole', 'all')

    edges = histogram_edges(cylinder_spec, whole)

    assert list(edges) == ['r', 'z', 'theta', 'x', 'y']
    assert_allclose(edges['r'], np.linspace(0.0, 2.5, 26), rtol=RTOL,
                    atol=ATOL)
    assert_allclose(edges['z'], np.linspace(-10.0, 10.0, 41), rtol=RTOL,
                    atol=ATOL)
    assert_allclose(edges['theta'], np.linspace(-np.pi, np.pi, 73),
                    rtol=RTOL, atol=ATOL)
    assert_allclose(edges['x'], np.linspace(-2.5, 2.5, 51), rtol=RTOL,
                    atol=ATOL)
    assert_array_equal(edges['y'], edges['x'])


def test_histogram_edges_scale_with_the_whole(cylinder_spec):
    """A crowder of diameter 0.5 halves the r bin size, so the same cylinder
    gets 50 r bins."""
    whole = ChainCyl(
        'N10r2.5lz10ac0.5nc4dt0.005bdump1adump1ens1', 'whole',
        'all')

    edges = histogram_edges(cylinder_spec, whole)

    assert edges['r'].size == 51


def test_histogram_edges_of_a_cube():
    """For lcube = 90: r runs to the corner, sqrt(3) * 45, in bins of 0.3,
    extended upward to a whole bin; x, y, and z have 60 bins of 1.5 over
    [-45, 45]."""
    project = ProjectSpec(
        parser=RingCub, group='bug', divisor=0.025,
        crowder_attr='phi_bulk_c',
        physical_attrs=('space', 'nmon_small', 'lcube'),
        species_sizes={'Mon': 'dmon_small', 'Foci': 'dmon_large',
                       'Crd': 'dcrowd'},
        tseries_properties=('gyrTMon',))
    spec = ProbeSpec(
        project, {'Mon': 'resid 1', 'Crd': 'resid 0'}, 'Mon', ('gyrT',),
        unwrap=True, histogram_species=('Mon',), block_frames=2,
        bin_sizes={'r': 0.3, 'x': 1.5, 'y': 1.5, 'z': 1.5})
    whole = RingCub(
        'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4', 'whole',
        'all')

    edges = histogram_edges(spec, whole)

    corner = np.sqrt(3.0) * 45.0
    assert list(edges) == ['r', 'x', 'y', 'z']
    assert edges['r'][0] == 0.0
    assert edges['r'].size == int(np.ceil(corner / 0.3)) + 1
    assert edges['r'][-1] == pytest.approx(0.3 * np.ceil(corner / 0.3),
                                           abs=1e-9)
    assert_allclose(edges['z'], np.linspace(-45.0, 45.0, 61), rtol=RTOL,
                    atol=ATOL)


def test_histogram_edges_need_every_key(cylinder_spec):
    """A spec lacking a bin size the geometry needs raises, naming it."""
    spec = ProbeSpec(
        cylinder_spec.project, dict(cylinder_spec.selections), 'Mon',
        ('gyrT',), unwrap=False, histogram_species=('Mon',), block_frames=2,
        bin_sizes={'r': 0.1, 'z': 0.5, 'theta': np.pi / 36})
    whole = ChainCyl(
        'N10r2.5lz10ac1nc0dt0.005bdump1adump1ens1', 'whole', 'all')

    with pytest.raises(ValueError, match=r"lacks \['x', 'y'\]"):
        histogram_edges(spec, whole)


def test_from_spec_selects_the_species_and_carries_the_spec(
        make_chain_universe, two_species, types, cylinder_spec):
    """from_spec selects each histogram species, builds the whole's edges,
    and takes the geometry, block length, and recentering from the spec."""
    universe = make_chain_universe(two_species, atom_types=types,
                                   masses=(1.0, 1.0))
    whole = ChainCyl(
        'N3r2.5lz10ac1nc3dt0.005bdump1adump1ens1', 'whole', 'all')

    analysis = SpatialHistograms.from_spec(universe, cylinder_spec, whole)

    assert [group.n_atoms for group in analysis.groups.values()] == [3, 3]
    assert list(analysis.groups) == ['Mon', 'Crd']
    assert (analysis.geometry, analysis.block_frames, analysis.recenter) == (
        'cylindrical', 2, 'half')
    assert analysis.edges['r'].size == 26


def test_from_spec_needs_histogram_species(make_chain_universe, two_species,
                                           cylinder_spec):
    """A spec without histogram species cannot build the analysis."""
    spec = ProbeSpec(cylinder_spec.project, dict(cylinder_spec.selections),
                     'Mon', ('gyrT',), unwrap=False)
    universe = make_chain_universe(two_species)
    whole = ChainCyl(
        'N3r2.5lz10ac1nc3dt0.005bdump1adump1ens1', 'whole', 'all')

    with pytest.raises(ValueError, match="names no histogram species"):
        SpatialHistograms.from_spec(universe, spec, whole)


# Saving


def test_save_writes_the_research_names(make_analysis, tmp_path):
    """save writes the edges, summed and per-frame histograms per species
    and direction, the Cartesian edges without a species, the summed and
    per-block maps, and the block frame counts, under the research names."""
    analysis = make_analysis().run()
    folder = tmp_path / 'probes'
    folder.mkdir()

    paths = analysis.save(folder, 'N6ens1-all')

    names = [path.name for path in paths]
    assert names[:9] == [
        'N6ens1-all-rEdgeMon.npy', 'N6ens1-all-rHistMon.npy',
        'N6ens1-all-rHistTMon.npy', 'N6ens1-all-zEdgeMon.npy',
        'N6ens1-all-zHistMon.npy', 'N6ens1-all-zHistTMon.npy',
        'N6ens1-all-thetaEdgeMon.npy', 'N6ens1-all-thetaHistMon.npy',
        'N6ens1-all-thetaHistTMon.npy']
    assert names[18:21] == ['N6ens1-all-xEdge.npy', 'N6ens1-all-yEdge.npy',
                            'N6ens1-all-zEdge.npy']
    assert names[21:23] == ['N6ens1-all-xyHistMon.npy',
                            'N6ens1-all-xyHistBlockMon.npy']
    assert names[-1] == 'N6ens1-all-blockFrames.npy'
    assert len(names) == 34
    assert_array_equal(np.load(paths[2]), analysis.results['rHistTMon'])
    assert_array_equal(np.load(paths[-1]), [2, 2, 1])


def test_save_before_run_raises(make_analysis, tmp_path):
    """save needs the results of run."""
    analysis = make_analysis()

    with pytest.raises(RuntimeError, match="Call run"):
        analysis.save(tmp_path, 'N6ens1-all')


# Real frames


@pytest.fixture
def real_histograms(study_spec, real_universe):
    """Factory: :meth:`SpatialHistograms.from_spec` run on a real dump with
    its topology, with the edges of the whole's parser in `group`."""

    def _run(whole, topology, trajectory, group):
        module = study_spec(whole.study)
        parser = module.SPEC.parser(f"{whole.whole}.ring", 'whole', group)
        universe = real_universe(topology, trajectory)
        return SpatialHistograms.from_spec(universe, module.PROBE,
                                           parser).run()

    return _run


@pytest.fixture
def directions():
    """The one-dimensional histograms of each geometry."""
    return {'cubic': ('r',), 'cylindrical': ('r', 'z', 'theta')}


def _counted(analysis, species, directions):
    """Per species and direction: the distinct per-frame sums, and the sum
    over the run."""
    return {
        f"{d}{s}": (np.unique(analysis.results[f"{d}HistT{s}"].sum(axis=1))
                    .tolist(), int(analysis.results[f"{d}Hist{s}"].sum()))
        for s in species for d in directions}


def test_real_window_histograms_conserve_every_species(
        real_window_whole, real_files, real_histograms, study_spec,
        directions):
    """Every position of a real ``all`` window falls inside the edges: each
    frame counts every particle of each species once, and the run counts it
    once per frame."""
    probe = study_spec(real_window_whole.study).PROBE
    files = real_files(real_window_whole)
    counts = real_window_whole.window_counts()

    analysis = real_histograms(real_window_whole, files.data, files.window,
                               'all')

    assert _counted(analysis, probe.histogram_species,
                    directions[probe.geometry]) == {
        f"{d}{s}": ([counts[s]], counts[s] * real_window_whole.n_window_frames)
        for s in probe.histogram_species
        for d in directions[probe.geometry]}


def test_real_polymer_histograms_conserve_every_species(
        real_polymer_whole, real_files, real_histograms, study_spec,
        directions):
    """Every position of a real polymer dump falls inside the edges of its
    whole: each frame counts every particle of each species once, and the
    crowders, absent from the dump, not at all."""
    probe = study_spec(real_polymer_whole.study).PROBE
    files = real_files(real_polymer_whole)
    counts = real_polymer_whole.polymer_counts()

    analysis = real_histograms(real_polymer_whole, files.topology,
                               files.polymer, real_polymer_whole.group)

    assert _counted(analysis, probe.histogram_species,
                    directions[probe.geometry]) == {
        f"{d}{s}": ([counts[s]], counts[s] * 50)
        for s in probe.histogram_species
        for d in directions[probe.geometry]}


def test_real_window_histograms_have_one_int64_row_per_frame_and_bin(
        real_window_whole, real_files, real_histograms, study_spec,
        directions):
    """Each per-frame histogram of a real window is int64 with one row per
    frame and one column per bin of the whole's edges."""
    module = study_spec(real_window_whole.study)
    probe = module.PROBE
    files = real_files(real_window_whole)
    edges = histogram_edges(probe, module.SPEC.parser(
        f"{real_window_whole.whole}.ring", 'whole', 'all'))

    analysis = real_histograms(real_window_whole, files.data, files.window,
                               'all')

    assert {f"{d}{s}": (analysis.results[f"{d}HistT{s}"].dtype,
                        analysis.results[f"{d}HistT{s}"].shape)
            for s in probe.histogram_species
            for d in directions[probe.geometry]} == {
        f"{d}{s}": (np.dtype(np.int64),
                    (real_window_whole.n_window_frames, edges[d].size - 1))
        for s in probe.histogram_species
        for d in directions[probe.geometry]}


@pytest.mark.private_data  # pragma: no cover - needs data
def test_crowded_window_histograms_conserve_every_species(
        real_crowded_whole, private_window_files, real_histograms, study_spec,
        directions):
    """Every position of a crowded ``all`` window, crowders at the wall and
    in the corners included, falls inside the edges: each frame counts every
    particle of each species once."""
    probe = study_spec(real_crowded_whole.study).PROBE
    window, data = private_window_files(real_crowded_whole)
    counts = real_crowded_whole.window_counts()

    analysis = real_histograms(real_crowded_whole, data, window, 'all')

    assert _counted(analysis, probe.histogram_species,
                    directions[probe.geometry]) == {
        f"{d}{s}": ([counts[s]],
                    counts[s] * real_crowded_whole.n_window_frames)
        for s in probe.histogram_species
        for d in directions[probe.geometry]}
