"""Tests for :mod:`polyphys.probe.pairs`.

The inputs are tier A: a chain of ten beads of which three, the foci, sit
at hand-picked positions, written as a data file and a dump by the
``make_chain_universe`` fixture. Expected values are analytic: the sides
of a 3-4-5 triangle, the contacts, bonds, and clusters they give at a
cutoff of 3.5 by hand, the minimum-image distance across a periodic
boundary by hand, and the backbone indices of the foci from their places
in the chain.

The real-frame tests run the TransFoci studies' probe specs on the tier B
polymer dumps of their wholes (``polyphys/tests/probe/conftest.py``, "Real
frames"). Their expected values are the foci's places along the ring in the
data file, read with :func:`polyphys.simulate.data.read_data`; the
minimum-image distances of :func:`MDAnalysis.lib.distances.distance_array`;
and the conservation of the foci among the bonds and the clusters.
"""
import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_equal,
    assert_equal,
)

from polyphys.manage.parser import TransFociCyl
from polyphys.manage.projects import ProbeSpec, ProjectSpec
from polyphys.simulate.data import read_data

pytest.importorskip("MDAnalysis")

from polyphys.probe import PairDistances  # noqa: E402

# The dump reader warns that the frames carry no time step; the tests
# measure per frame and give no dt.
pytestmark = pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")

# Tolerances: MDAnalysis holds positions as float32, so a distance carries a
# relative round-off of about 1e-7; 1e-6 leaves room for a few operations.
RTOL = 1e-6
ATOL = 1e-6


@pytest.fixture
def types() -> list[int]:
    """The atom type of each of ten beads: the foci, type 2, are beads 3,
    6, and 9 (LAMMPS IDs), at backbone indices 2, 5, and 8."""
    return [1, 1, 2, 1, 1, 2, 1, 1, 2, 1]


@pytest.fixture
def triangle(types) -> np.ndarray:
    """Two frames of ten beads about (10, 10, 10): the foci at the corners
    of a 3-4-5 triangle, (0, 0, 0), (3, 0, 0), and (3, 4, 0) from the
    center, the other beads at the center; in frame 1 the third focus moves
    to (3, 3, 0), 3 from the second and sqrt(18) from the first."""
    first = np.full((10, 3), 10.0)
    first[2] += [0.0, 0.0, 0.0]
    first[5] += [3.0, 0.0, 0.0]
    first[8] += [3.0, 4.0, 0.0]
    second = first.copy()
    second[8] = [13.0, 13.0, 10.0]
    return np.stack([first, second])


@pytest.fixture
def make_analysis(make_chain_universe, triangle, types):
    """Factory: the analysis of the three foci of the triangle frames, with
    any argument overridden; `positions` replaces the frames."""

    def _make(positions=None, **overrides):
        universe = make_chain_universe(
            triangle if positions is None else positions,
            atom_types=types, masses=(1.0, 1.0))
        kwargs = {'cutoff': 3.5, 'genomic_positions': [2, 5, 8]}
        kwargs.update(overrides)
        return PairDistances(universe.select_atoms('type 2'), **kwargs)

    return _make


# Hand values


def test_distance_matrix_holds_the_upper_triangle_and_the_genomic_diagonal(
        make_analysis):
    """Frame 0 gives the triangle's sides 3, 5, and 4 above the diagonal,
    the backbone indices 2, 5, 8 on it, and zeros below."""
    analysis = make_analysis().run()

    matrix = analysis.results['distMatTFoci']

    assert matrix.shape == (2, 3, 3)
    assert matrix.dtype == np.float64
    assert_allclose(matrix[0], [[2.0, 3.0, 5.0],
                                [0.0, 5.0, 4.0],
                                [0.0, 0.0, 8.0]], rtol=RTOL, atol=ATOL)
    assert_allclose(matrix[1], [[2.0, 3.0, np.sqrt(18.0)],
                                [0.0, 5.0, 3.0],
                                [0.0, 0.0, 8.0]], rtol=RTOL, atol=ATOL)


def test_contacts_bonds_and_clusters_at_the_cutoff(make_analysis):
    """At a cutoff of 3.5, foci 1 and 2 (3 apart) touch and focus 3 is
    alone in frame 0: one focus with no contact and two with one each,
    bonds [1, 2, 0], and one cluster of size 1 and one of size 2."""
    analysis = make_analysis().run()

    contacts = analysis.results['directContactsMatTFoci'][0]

    assert contacts.dtype == np.bool_
    assert_array_equal(contacts, [[True, True, False],
                                  [True, True, False],
                                  [False, False, True]])
    assert_array_equal(analysis.results['bondsHistTFoci'][0], [1, 2, 0])
    assert_array_equal(analysis.results['clustersHistTFoci'][0],
                       [0, 1, 1, 0])


def test_frame_1_forms_a_chain_of_three(make_analysis):
    """In frame 1, focus 3 touches focus 2 (3 apart) but not focus 1
    (sqrt(18) apart): the bonds are [0, 2, 1], one cluster of size 3."""
    analysis = make_analysis().run()

    assert_array_equal(analysis.results['bondsHistTFoci'][1], [0, 2, 1])
    assert_array_equal(analysis.results['clustersHistTFoci'][1],
                       [0, 0, 0, 1])


@pytest.mark.parametrize('inclusive, bonds', [
    pytest.param(True, [1, 2, 0], id='at-the-cutoff-touches'),
    pytest.param(False, [3, 0, 0], id='at-the-cutoff-does-not'),
])
def test_inclusive_decides_a_pair_at_the_cutoff(make_analysis, inclusive,
                                                bonds):
    """A pair exactly at the cutoff, 3 apart with cutoff 3, is in contact
    only when the cutoff is inclusive."""
    analysis = make_analysis(cutoff=3.0, inclusive=inclusive).run()

    assert_array_equal(analysis.results['bondsHistTFoci'][0], bonds)


# Minimum image


@pytest.mark.parametrize('periodic, distance, contact', [
    pytest.param((False, False, True), 1.0, True, id='z-periodic'),
    pytest.param(None, 19.0, False, id='no-minimum-image'),
])
def test_minimum_image_along_the_periodic_axis(make_chain_universe, types,
                                               periodic, distance, contact):
    """Two foci at z = 0.5 and 19.5 in a box of length 20 are 1 apart with
    the minimum image along z and 19 without."""
    positions = np.full((1, 10, 3), 10.0)
    positions[0, 2] = [10.0, 10.0, 0.5]
    positions[0, 5] = [10.0, 10.0, 19.5]
    positions[0, 8] = [15.0, 10.0, 10.0]
    universe = make_chain_universe(positions, atom_types=types,
                                   masses=(1.0, 1.0))

    analysis = PairDistances(universe.select_atoms('type 2'), 1.5,
                             [2, 5, 8], periodic=periodic).run()

    assert analysis.results['distMatTFoci'][0, 0, 1] == pytest.approx(
        distance, rel=RTOL)
    assert analysis.results['directContactsMatTFoci'][0, 0, 1] is np.bool_(
        contact)


def test_numpy_boolean_periodic_flags_act_as_booleans(make_analysis):
    """Periodic flags given as numpy booleans give the run of the same
    flags given as Python booleans."""
    expected = make_analysis(periodic=(False, False, True)).run()

    analysis = make_analysis(
        periodic=(np.bool_(False), np.bool_(False), np.bool_(True))).run()

    assert analysis.periodic == (False, False, True)
    assert_equal(dict(analysis.results), dict(expected.results))


def test_periodic_direction_needs_a_box():
    """A frame without box dimensions cannot take the minimum image, and
    the run raises ValueError, as the other probes do."""
    import MDAnalysis as mda

    universe = mda.Universe.empty(2, trajectory=True)
    universe.atoms.positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    analysis = PairDistances(universe.atoms, 1.5, [0, 1],
                             periodic=(False, False, True))

    with pytest.raises(ValueError, match="carries no box"):
        analysis.run()


# Parallel


def test_parallel_run_matches_serial(make_analysis):
    """Two workers over two parts give the same four arrays as one."""
    serial = make_analysis().run()
    parallel = make_analysis().run(backend='multiprocessing', n_workers=2,
                                   n_parts=2)

    for key in serial.results:
        assert_array_equal(parallel.results[key], serial.results[key])


def test_serial_parts_match_the_multiprocessing_run(make_analysis):
    """The serial backend split into two parts, one frame each, gives the
    frames, times, and arrays of the multiprocessing backend over the same
    parts: the result does not depend on the backend."""
    parallel = make_analysis().run(backend='multiprocessing', n_workers=2,
                                   n_parts=2)

    serial = make_analysis().run(backend='serial', n_parts=2)

    assert_array_equal(serial.frames, parallel.frames)
    assert_array_equal(serial.times, parallel.times)
    assert_equal(dict(serial.results), dict(parallel.results))


# from_spec


@pytest.fixture
def foci_spec():
    """A TransFociCyl probe spec whose pair species is the foci, with the
    research code's cutoff, the large-monomer and crowder diameters summed."""
    project = ProjectSpec(
        parser=TransFociCyl, group='bug', divisor=0.025,
        physical_attrs=('space', 'nmon_small', 'nmon_large', 'dcyl'),
        species_sizes={'Mon': 'dmon_small', 'Foci': 'dmon_large',
                       'Crd': 'dcrowd'},
        tseries_properties=('gyrTMon',))

    def cutoff(info):
        return info.dmon_large + info.dcrowd

    return ProbeSpec(
        project, {'Mon': 'resid 1', 'Foci': 'type 2'}, 'Mon', ('gyrT',),
        unwrap=False, pair_species='Foci', contact_cutoff=cutoff)


@pytest.fixture
def foci_whole():
    """A TransFociCyl whole with dmon_large 2 and dcrowd 1.5, so the cutoff
    is 3.5."""
    return TransFociCyl(
        'epss5epsl5r5al2nl3ml8ns7ac1.5nc0lz10dt0.005bdump1adump1ens1',
        'whole', 'bug')


def test_from_spec_takes_backbone_indices_cutoff_and_periodicity(
        make_chain_universe, triangle, types, foci_spec, foci_whole):
    """from_spec selects the foci, places them at their 0-based positions
    among the chain's atoms, resolves the cutoff for the whole, and uses
    the minimum image along the cylinder axis."""
    universe = make_chain_universe(triangle, atom_types=types,
                                   masses=(1.0, 1.0))

    analysis = PairDistances.from_spec(universe, foci_spec, foci_whole)

    assert analysis.atomgroup.n_atoms == 3
    assert_array_equal(analysis.genomic_positions, [2, 5, 8])
    assert (analysis.cutoff, analysis.species, analysis.periodic,
            analysis.inclusive) == (3.5, 'Foci', (False, False, True), True)


def test_from_spec_without_minimum_image(make_chain_universe, triangle, types,
                                         foci_spec, foci_whole):
    """A spec with minimum_image False gives no periodic direction."""
    spec = ProbeSpec(
        foci_spec.project, dict(foci_spec.selections), 'Mon', ('gyrT',),
        unwrap=False, pair_species='Foci',
        contact_cutoff=foci_spec.contact_cutoff, minimum_image=False)
    universe = make_chain_universe(triangle, atom_types=types,
                                   masses=(1.0, 1.0))

    analysis = PairDistances.from_spec(universe, spec, foci_whole)

    assert analysis.periodic is None


def test_from_spec_needs_a_pair_species(make_chain_universe, triangle, types,
                                        foci_spec, foci_whole):
    """A spec without a pair species cannot build the analysis."""
    spec = ProbeSpec(foci_spec.project, dict(foci_spec.selections), 'Mon',
                     ('gyrT',), unwrap=False)
    universe = make_chain_universe(triangle, atom_types=types,
                                   masses=(1.0, 1.0))

    with pytest.raises(ValueError, match="names no pair species"):
        PairDistances.from_spec(universe, spec, foci_whole)


def test_from_spec_rejects_foci_outside_the_chain(make_chain_universe,
                                                  triangle, types, foci_spec,
                                                  foci_whole):
    """A pair species whose atoms are not among the chain's raises, naming
    their IDs."""
    spec = ProbeSpec(
        foci_spec.project, {'Mon': 'type 1', 'Foci': 'type 2'}, 'Mon',
        ('gyrT',), unwrap=False, pair_species='Foci',
        contact_cutoff=foci_spec.contact_cutoff)
    universe = make_chain_universe(triangle, atom_types=types,
                                   masses=(1.0, 1.0))

    with pytest.raises(ValueError, match=r"IDs \[3, 6, 9\] are not among"):
        PairDistances.from_spec(universe, spec, foci_whole)


# Validation


@pytest.mark.parametrize('overrides, match', [
    pytest.param({'cutoff': 0.0}, "'cutoff' must be a finite positive",
                 id='zero-cutoff'),
    pytest.param({'cutoff': float('inf')}, "'cutoff' must be a finite",
                 id='infinite-cutoff'),
    pytest.param({'genomic_positions': [2, 5]},
                 "one non-negative integer per particle, 3 in all",
                 id='too-few-positions'),
    pytest.param({'genomic_positions': [2, -5, 8]},
                 "one non-negative integer per particle",
                 id='negative-position'),
    pytest.param({'genomic_positions': [2.0, 5.0, 8.0]},
                 "one non-negative integer per particle",
                 id='float-positions'),
    pytest.param({'genomic_positions': [2, 2, 8]}, "repeats an index",
                 id='repeated-position'),
    pytest.param({'species': ''}, "'species' must be a non-empty",
                 id='empty-species'),
    pytest.param({'periodic': (True, False)}, "'periodic' must be None",
                 id='short-periodic'),
    pytest.param({'periodic': (1, 0, 1)}, "'periodic' must be None",
                 id='non-boolean-periodic'),
])
def test_invalid_arguments_raise(make_analysis, overrides, match):
    """A bad cutoff, backbone index list, species, or periodicity raises at
    construction."""
    with pytest.raises(ValueError, match=match):
        make_analysis(**overrides)


def test_fewer_than_two_particles_raise(make_chain_universe, triangle):
    """One particle has no pair."""
    universe = make_chain_universe(triangle)

    with pytest.raises(ValueError, match="at least two particles"):
        PairDistances(universe.atoms[:1], 3.5, [0])


# Saving


def test_save_writes_the_research_names(make_analysis, tmp_path):
    """save writes the four arrays as the research code named them, with
    the species label, and each file holds its array."""
    analysis = make_analysis().run()
    folder = tmp_path / 'probes'
    folder.mkdir()

    paths = analysis.save(folder, 'ns7nl3ens1-bug')

    assert [path.name for path in paths] == [
        'ns7nl3ens1-bug-distMatTFoci.npy',
        'ns7nl3ens1-bug-directContactsMatTFoci.npy',
        'ns7nl3ens1-bug-bondsHistTFoci.npy',
        'ns7nl3ens1-bug-clustersHistTFoci.npy']
    assert_array_equal(np.load(paths[0]), analysis.results['distMatTFoci'])
    assert_array_equal(np.load(paths[3]),
                       analysis.results['clustersHistTFoci'])


def test_save_before_run_raises(make_analysis, tmp_path):
    """save needs the results of run."""
    analysis = make_analysis()

    with pytest.raises(RuntimeError, match="Call run"):
        analysis.save(tmp_path, 'ns7nl3ens1-bug')


# Real frames


@pytest.fixture
def real_pairs(real_files, real_universe, study_spec):
    """Factory: the study's :meth:`PairDistances.from_spec` run on a real
    TransFoci polymer dump, with the universe it read."""

    def _run(whole):
        module = study_spec(whole.study)
        files = real_files(whole)
        universe = real_universe(files.topology, files.polymer)
        parser = module.SPEC.parser(f"{whole.whole}.ring", 'whole',
                                    whole.group)
        return (PairDistances.from_spec(universe, module.PROBE,
                                        parser).run(), universe)

    return _run


def test_real_foci_diagonal_holds_their_places_along_the_ring(
        real_foci_whole, real_files, real_pairs):
    """Every frame's diagonal holds the 0-based places of the large
    monomers along the ring, as the data file orders its atoms."""
    data = read_data(real_files(real_foci_whole).topology)
    expected = np.flatnonzero(data.types[data.mol_ids == 1] == 2)

    analysis, _ = real_pairs(real_foci_whole)

    assert_array_equal(
        np.diagonal(analysis.results['distMatTFoci'], axis1=1, axis2=2),
        np.tile(expected, (50, 1)))


def test_real_foci_distances_are_the_minimum_image_distances(
        real_foci_whole, real_pairs, study_spec, trajectory_arrays,
        reference_distances):
    """Above the diagonal, every frame holds the foci's minimum-image
    distances along the geometry's periodic directions."""
    probe = study_spec(real_foci_whole.study).PROBE

    analysis, universe = real_pairs(real_foci_whole)

    arrays = trajectory_arrays(universe)
    foci = arrays['positions'][:, analysis.atomgroup.indices]
    expected = reference_distances(foci, foci, arrays['boxes'],
                                   probe.geometry)
    upper = np.triu_indices(analysis.atomgroup.n_atoms, k=1)
    # float32 positions of at most about 165 carry a round-off of about
    # 1e-5 in a separation.
    assert_allclose(analysis.results['distMatTFoci'][:, upper[0], upper[1]],
                    expected[:, upper[0], upper[1]], rtol=0, atol=5e-5)


def test_real_foci_bonds_and_clusters_conserve_the_foci(real_foci_whole,
                                                        real_pairs):
    """In every frame, the bonds histogram counts each large monomer once
    and the clusters histogram, weighted by cluster size, does too."""
    analysis, _ = real_pairs(real_foci_whole)
    clusters = analysis.results['clustersHistTFoci']
    sizes = np.arange(clusters.shape[1])

    assert_array_equal(analysis.results['bondsHistTFoci'].sum(axis=1),
                       np.full(50, real_foci_whole.n_foci))
    assert_array_equal(clusters @ sizes, np.full(50, real_foci_whole.n_foci))
