"""Tests for :mod:`polyphys.probe.binders`.

The inputs are tier A: a ring of ten monomers one unit apart along x and
four two-patch binders at hand-picked positions over five frames, written
as a data file and a dump by the ``make_binder_universe`` fixture. Expected
values are worked by hand from the counting rules of
:func:`polyphys.analyze.clusters.binding_statistics`, as the clusters guide
states them: a patch 0.3 off a monomer touches it at a cutoff of 0.5, a
patch at the midpoint of two neighbors is exactly 0.5 from both, the
single-patch rule keeps the lower index, cis is a genomic distance of at
most ``cis_threshold`` bonds, and the loop-length row is the bincount of
the bridged pairs' distances. The minimum-image distances are checked
against ``MDAnalysis.lib.distances.distance_array``.

The real-frame tests run the Hns studies' probe specs on the tier B polymer
dumps of their wholes (``polyphys/tests/probe/conftest.py``, "Real
frames"). Their expected values are the identities the counts obey by
definition, with the binder numbers of the wholes' names and two patches per
H-NS core, and the minimum-image distances of ``distance_array``.
"""
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_equal,
    assert_equal,
)

from polyphys.analyze.clusters import threshold_sweep
from polyphys.manage.parser import HnsCub, HnsCyl, ParserBase
from polyphys.manage.projects import (
    BINDER_SERIES,
    BinderSpec,
    ProbeSpec,
    ProjectSpec,
)
from polyphys.simulate.core import Box
from polyphys.simulate.data import DataFile, write_data

pytest.importorskip("MDAnalysis")

import MDAnalysis as mda  # noqa: E402
from MDAnalysis.lib.distances import distance_array  # noqa: E402

from polyphys.probe import BinderStatistics  # noqa: E402

# The dump reader warns that the frames carry no time step; the tests
# count per frame and give no dt.
pytestmark = pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")

# Tolerances for the distances: MDAnalysis holds positions as float32, so a
# distance carries a relative round-off of about 1e-7, and distance_array
# may sum the squares in float32; 1e-5 leaves room for both.
RTOL = 1e-5
ATOL = 1e-6


@pytest.fixture
def make_binder_universe(
    tmp_path: Path, make_lammps_dump: Callable[..., Path]
) -> Callable[..., Any]:
    """Factory: the universe of a chain and its binders' patches, written
    as a ``bond``-style data file (type 1 and molecule 1 for the chain,
    type 2 and molecule ``2 + k`` for the patches of binder ``k``) and a
    dump of the frames.

    Takes `chain` of shape ``(n_frames, n_monomers, 3)`` and `patches` of
    shape ``(n_frames, n_patches, 3)``, `box` as three ``(lo, hi)`` pairs,
    and `boundary` as the dump's three flags.
    """

    def _make(
        chain: Any,
        patches: Any,
        *,
        box: Sequence[Sequence[float]] = ((0.0, 20.0),) * 3,
        boundary: str = "pp pp pp",
        name: str = "binders",
    ) -> Any:
        from polyphys.probe import make_universe

        chain = np.asarray(chain, dtype=float)
        patches = np.asarray(patches, dtype=float)
        n_monomers, n_patches = chain.shape[1], patches.shape[1]
        frames = np.concatenate([chain, patches], axis=1)
        types = np.array([1] * n_monomers + [2] * n_patches)
        mol_ids = np.array([1] * n_monomers
                           + [2 + k // 2 for k in range(n_patches)])
        data = DataFile(
            atom_style="bond",
            box=Box(lo=np.array([lo for lo, _ in box], dtype=float),
                    hi=np.array([hi for _, hi in box], dtype=float)),
            ids=np.arange(1, n_monomers + n_patches + 1), types=types,
            positions=frames[0], n_atom_types=2,
            masses=np.array([1.0, 1.0]), mol_ids=mol_ids,
            bonds=np.array([[1, i, i + 1] for i in range(1, n_monomers)]),
            n_bond_types=1)
        data_path = write_data(data, tmp_path / f"{name}.data")
        dump_path = make_lammps_dump(
            name=f"{name}.lammpstrj", positions=frames, box=box,
            boundary=boundary, atom_types=types.tolist())
        return make_universe(data_path, dump_path)

    return _make


@pytest.fixture
def ring() -> np.ndarray:
    """One frame of ten monomers at x = 1, ..., 10 and y = z = 10."""
    chain = np.full((10, 3), 10.0)
    chain[:, 0] = np.arange(10) + 1.0
    return chain


@pytest.fixture
def far() -> np.ndarray:
    """Eight patch positions far from the ring, all distinct."""
    patches = np.full((8, 3), 15.0)
    patches[:, 1] = 12.0 + 0.5 * np.arange(8)
    return patches


@pytest.fixture
def binding_frames(ring: np.ndarray, far: np.ndarray) -> tuple[np.ndarray,
                                                               np.ndarray]:
    """Five frames of the ring and four binders A (patches 0, 1), B (2, 3),
    C (4, 5), and D (6, 7), each patch 0.3 off its monomer:

    0. A free; B dangling from monomer 2; C bridging 4 and 6 (cis, two
       bonds); D bridging 0 and 5 (trans, five bonds).
    1. As frame 0, but D's first patch at the midpoint of monomers 0 and
       1, exactly 0.5 from both.
    2. All free.
    3. C and D both bridging 4 and 6.
    4. C with both patches on monomer 4; the rest free.
    """
    chain = np.stack([ring] * 5)
    patches = np.stack([far] * 5)
    patches[0, 2] = [3.0, 10.3, 10.0]
    patches[0, 4] = [5.0, 10.3, 10.0]
    patches[0, 5] = [7.0, 10.3, 10.0]
    patches[0, 6] = [1.0, 10.3, 10.0]
    patches[0, 7] = [6.0, 9.7, 10.0]
    patches[1] = patches[0]
    patches[1, 6] = [1.5, 10.0, 10.0]
    patches[3, 4] = [5.0, 10.3, 10.0]
    patches[3, 5] = [7.0, 10.3, 10.0]
    patches[3, 6] = [5.0, 9.7, 10.0]
    patches[3, 7] = [7.0, 9.7, 10.0]
    patches[4, 4] = [5.0, 10.3, 10.0]
    patches[4, 5] = [5.0, 9.7, 10.0]
    return chain, patches


@pytest.fixture
def hns_project() -> Callable[[type[ParserBase]], ProjectSpec]:
    """Factory: a project spec of an H-NS parser class that lists the eight
    binding series."""

    def _make(parser: type[ParserBase]) -> ProjectSpec:
        return ProjectSpec(
            parser=parser, group='nucleoid', divisor=0.01,
            physical_attrs=('space', 'nmon', 'nhns'),
            species_sizes={'Mon': 'dmon', 'Hns': 'dhns'},
            tseries_properties=(
                'gyrTMon', 'nBoundTHnsPatch', 'nFreeTHnsPatch',
                'nEngagedTHnsPatch', 'nFreeTHnsCore', 'nBridgeTHnsCore',
                'nDangleTHnsCore', 'nCisTHnsCore', 'nTransTHnsCore'))

    return _make


# Hand values


@pytest.mark.parametrize('series, expected', [
    pytest.param('n_patch_monomer_contacts', [5, 6, 0, 4, 2],
                 id='contacts'),
    pytest.param('n_patches_free', [3, 3, 8, 4, 6], id='free-patches'),
    pytest.param('n_patches_engaged', [5, 5, 0, 4, 2],
                 id='engaged-patches'),
    pytest.param('n_binders_free', [1, 1, 4, 2, 3], id='free-binders'),
    pytest.param('n_binders_bridging', [2, 2, 0, 2, 1], id='bridging'),
    pytest.param('n_binders_dangling', [1, 1, 0, 0, 0], id='dangling'),
    pytest.param('n_pairs_cis', [1, 1, 0, 1, 0], id='cis'),
    pytest.param('n_pairs_trans', [1, 1, 0, 0, 0], id='trans'),
])
def test_counts_follow_the_binding_rules(make_binder_universe,
                                         binding_frames, series, expected):
    """Each series holds one int64 count per frame: the midpoint patch of
    frame 1 adds one contact but no engaged patch, two binders on one pair
    bridge one cis pair, and two patches on one monomer bridge no pair."""
    universe = make_binder_universe(*binding_frames)

    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', cis_threshold=4).run()

    assert analysis.results[series].dtype == np.int64
    assert_array_equal(analysis.results[series], expected)


def test_loop_length_rows_count_bridged_pairs_by_distance(
        make_binder_universe, binding_frames):
    """A ring of ten has six loop lengths; each row counts the frame's
    distinct bridged pairs at their genomic distance, with column 0 empty."""
    universe = make_binder_universe(*binding_frames)

    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', cis_threshold=4).run()

    histogram = analysis.results['loop_length_histogram']
    assert histogram.dtype == np.int64
    assert_array_equal(histogram, [[0, 0, 1, 0, 0, 1],
                                   [0, 0, 1, 0, 0, 1],
                                   [0, 0, 0, 0, 0, 0],
                                   [0, 0, 1, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0]])


@pytest.mark.parametrize('inclusive, expected', [
    # Inclusive: the midpoint patch touches monomers 0 and 1, D bridges 0
    # and 5 through the lower index. Strict: it touches neither, D dangles.
    pytest.param(True, {'n_patch_monomer_contacts': 6,
                        'n_binders_dangling': 1, 'n_binders_bridging': 2,
                        'n_pairs_trans': 1}, id='inclusive'),
    pytest.param(False, {'n_patch_monomer_contacts': 4,
                         'n_binders_dangling': 2, 'n_binders_bridging': 1,
                         'n_pairs_trans': 0}, id='strict'),
])
def test_inclusive_decides_a_pair_at_the_cutoff(make_binder_universe,
                                                binding_frames, inclusive,
                                                expected):
    """A patch exactly at the cutoff from two neighbors is in contact with
    both when inclusive and with neither when strict."""
    universe = make_binder_universe(*binding_frames)

    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', inclusive=inclusive).run(start=1,
                                                              stop=2)

    assert {key: int(analysis.results[key][0]) for key in expected} \
        == expected


@pytest.mark.parametrize('topology, cis, trans, row', [
    # Monomers 0 and 8 are two bonds apart the short way round a ring of
    # ten and eight bonds apart on a linear chain.
    pytest.param('ring', 1, 0, [0, 0, 1, 0, 0, 0], id='ring'),
    pytest.param('linear', 0, 1, [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                 id='linear'),
])
def test_topology_sets_the_genomic_distance(make_binder_universe, ring, far,
                                            topology, cis, trans, row):
    """A bridge between monomers 0 and 8 is cis at two bonds on a ring and
    trans at eight bonds on a linear chain, whose histogram is longer."""
    patches = far.copy()
    patches[6] = [1.0, 10.3, 10.0]
    patches[7] = [9.0, 10.3, 10.0]
    universe = make_binder_universe(ring[np.newaxis], patches[np.newaxis])

    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology=topology, cis_threshold=4).run()

    assert (int(analysis.results['n_pairs_cis'][0]),
            int(analysis.results['n_pairs_trans'][0])) == (cis, trans)
    assert_array_equal(analysis.results['loop_length_histogram'], [row])


@pytest.mark.parametrize('periodic, contacts', [
    pytest.param((True, True, True), 1, id='all-periodic'),
    pytest.param((False, False, True), 0, id='axis-only'),
    pytest.param((False, False, False), 0, id='none'),
])
def test_minimum_image_follows_the_periodic_directions(
        make_binder_universe, ring, far, periodic, contacts):
    """A patch 0.4 across the x face from monomer 0 touches it only when x
    is periodic; a cylinder's z-only image leaves them 19.6 apart."""
    chain = ring.copy()
    chain[:, 0] -= 0.7
    patches = far.copy()
    patches[0] = [19.9, 10.0, 10.0]
    universe = make_binder_universe(chain[np.newaxis], patches[np.newaxis])

    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', periodic=periodic).run()

    assert int(analysis.results['n_patch_monomer_contacts'][0]) == contacts


@pytest.mark.parametrize('periodic, box', [
    pytest.param((False, False, False), False, id='no-box'),
    pytest.param((True, True, True), True, id='full-box'),
])
def test_distances_match_mdanalysis(make_binder_universe, ring, far,
                                    periodic, box):
    """The stored distance matrix equals distance_array without a box when
    nothing is periodic and with the frame's box when everything is."""
    rng = np.random.default_rng(7)
    chain = ring + rng.uniform(-0.5, 0.5, size=ring.shape)
    patches = rng.uniform(0.0, 20.0, size=far.shape)
    universe = make_binder_universe(chain[np.newaxis], patches[np.newaxis])
    monomers = universe.select_atoms('resid 1')
    patch_atoms = universe.select_atoms('type 2')

    analysis = BinderStatistics(
        monomers, patch_atoms, cutoff=0.5, topology='ring',
        periodic=periodic, store_distances=True).run()

    expected = distance_array(
        monomers.positions, patch_atoms.positions,
        box=universe.trajectory.ts.dimensions if box else None)
    assert analysis.results['distances'].shape == (1, 10, 8)
    assert analysis.results['distances'].dtype == np.float64
    assert_allclose(analysis.results['distances'][0], expected, rtol=RTOL,
                    atol=ATOL)


# The analysis contract


def test_run_slices_the_frames(make_binder_universe, binding_frames):
    """run(start, stop, step) counts only those frames, in order."""
    universe = make_binder_universe(*binding_frames)

    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring').run(start=1, stop=5, step=2)

    assert_array_equal(analysis.frames, [1, 3])
    assert_array_equal(analysis.results['n_binders_bridging'], [2, 2])
    assert_array_equal(analysis.results['n_binders_dangling'], [1, 0])
    assert analysis.n_frames == 2


def test_parallel_run_matches_the_serial_run(make_binder_universe,
                                             binding_frames):
    """Two workers over two parts give the serial arrays: the series, the
    histogram rows, and the distance matrices."""
    universe = make_binder_universe(*binding_frames)
    serial = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', store_distances=True).run()

    parallel = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', store_distances=True).run(
            backend='multiprocessing', n_workers=2, n_parts=2)

    assert parallel.keys == serial.keys
    for key in serial.keys:
        assert_array_equal(parallel.results[key], serial.results[key])


@pytest.mark.parametrize('n_parts', [
    pytest.param(2, id='two-parts'),
    pytest.param(3, id='three-parts'),
])
def test_serial_parts_match_the_multiprocessing_run(
        make_binder_universe, binding_frames, n_parts):
    """The serial backend split into parts gives the frames, times, and
    arrays of the multiprocessing backend over the same parts: the result
    does not depend on the backend."""
    universe = make_binder_universe(*binding_frames)
    chain = universe.select_atoms('resid 1')
    patches = universe.select_atoms('type 2')
    parallel = BinderStatistics(
        chain, patches, cutoff=0.5, topology='ring',
        store_distances=True).run(
            backend='multiprocessing', n_workers=2, n_parts=n_parts)

    serial = BinderStatistics(
        chain, patches, cutoff=0.5, topology='ring',
        store_distances=True).run(backend='serial', n_parts=n_parts)

    assert_array_equal(serial.frames, parallel.frames)
    assert_array_equal(serial.times, parallel.times)
    assert_equal(dict(serial.results), dict(parallel.results))


def test_stored_distances_swept_at_the_spec_reproduce_the_series(
        make_binder_universe, binding_frames):
    """The threshold sweep of the stored matrices at the run's cutoff, cis
    threshold, and comparison gives the frame means of the saved series and
    histogram rows, so the sweep's grid point at the spec reproduces the
    probe's counts."""
    universe = make_binder_universe(*binding_frames)
    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', inclusive=True, cis_threshold=4,
        store_distances=True).run()

    counts, histogram = threshold_sweep(
        analysis.results['distances'], [0.5], [4], 'ring', inclusive=True)

    # Both sides are an integer sum over five frames divided by five.
    for key in BINDER_SERIES:
        assert counts.loc[0, key] == pytest.approx(
            analysis.results[key].mean(), abs=1e-12)
    assert_allclose(histogram[0],
                    analysis.results['loop_length_histogram'].mean(axis=0),
                    rtol=1e-12, atol=1e-12)


def test_keys_list_the_series_then_the_histogram_then_the_distances(
        make_binder_universe, binding_frames):
    """keys is the eight series in BINDER_SERIES order, the histogram, and,
    with store_distances, the distances."""
    universe = make_binder_universe(*binding_frames)
    chain, patches = (universe.select_atoms('resid 1'),
                      universe.select_atoms('type 2'))

    without = BinderStatistics(chain, patches, cutoff=0.5, topology='ring')
    with_distances = BinderStatistics(chain, patches, cutoff=0.5,
                                      topology='ring', store_distances=True)

    assert without.keys == (*BINDER_SERIES, 'loop_length_histogram')
    assert with_distances.keys == (*BINDER_SERIES, 'loop_length_histogram',
                                   'distances')


@pytest.mark.parametrize('parser, periodic', [
    pytest.param(HnsCub, (True, True, True), id='HnsCub'),
    pytest.param(HnsCyl, (False, False, True), id='HnsCyl'),
])
def test_from_spec_selects_the_groups_and_carries_the_spec(
        make_binder_universe, binding_frames, hns_project, parser,
        periodic):
    """from_spec picks the spec's chain and patch selections, its binder
    spec, the project's topology, and the geometry's periodic directions."""
    spec = ProbeSpec(
        hns_project(parser), {'Mon': 'resid 1', 'HnsPatch': 'type 2'},
        'Mon', ('gyrT',), unwrap=True,
        binders=BinderSpec('HnsPatch', 'HnsCore', cutoff=0.685,
                           inclusive=True, cis_threshold=4,
                           store_distances=True))
    universe = make_binder_universe(*binding_frames)

    analysis = BinderStatistics.from_spec(universe, spec)

    assert (analysis.chain.n_atoms, analysis.patches.n_atoms) == (10, 8)
    assert (analysis.cutoff, analysis.inclusive, analysis.cis_threshold,
            analysis.patches_per_binder, analysis.store_distances) == (
        0.685, True, 4, 2, True)
    assert (analysis.topology, analysis.periodic) == ('ring', periodic)
    assert (analysis.chain_label, analysis.patch_label,
            analysis.binder_label) == ('Mon', 'HnsPatch', 'HnsCore')


def test_from_spec_needs_a_binder_spec(make_binder_universe, binding_frames,
                                       hns_project):
    """A probe spec without binders cannot build the analysis."""
    spec = ProbeSpec(hns_project(HnsCub),
                     {'Mon': 'resid 1', 'HnsPatch': 'type 2'}, 'Mon',
                     ('gyrT',), unwrap=True)
    universe = make_binder_universe(*binding_frames)

    with pytest.raises(ValueError, match="has no 'binders'"):
        BinderStatistics.from_spec(universe, spec)


@pytest.mark.parametrize('kwargs, match', [
    pytest.param({'patches_per_binder': 3}, "must be 2, not 3",
                 id='three-patches'),
    pytest.param({'cutoff': 0.0}, "'cutoff' must be a finite positive",
                 id='zero-cutoff'),
    pytest.param({'cutoff': float('inf')},
                 "'cutoff' must be a finite positive", id='infinite-cutoff'),
    pytest.param({'topology': 'star'}, "'topology' must be one of",
                 id='bad-topology'),
    pytest.param({'cis_threshold': -1}, "'cis_threshold' must be",
                 id='negative-cis'),
    pytest.param({'periodic': (True, True)}, "'periodic' must be three",
                 id='two-flags'),
    pytest.param({'periodic': (1, 0, 0)}, "'periodic' must be three",
                 id='int-flags'),
    pytest.param({'chain_label': ''}, "'chain_label' must be a non-empty",
                 id='empty-label'),
])
def test_rejects_invalid_arguments(make_binder_universe, binding_frames,
                                   kwargs, match):
    """A valence other than two, a non-positive or infinite cutoff, an
    unknown topology, a negative cis threshold, malformed periodic flags,
    or an empty label raises at construction."""
    universe = make_binder_universe(*binding_frames)
    fields = {'cutoff': 0.5, 'topology': 'ring', **kwargs}

    with pytest.raises(ValueError, match=match):
        BinderStatistics(universe.select_atoms('resid 1'),
                         universe.select_atoms('type 2'), **fields)


@pytest.mark.parametrize('chain, patches, match', [
    pytest.param('id 1', 'type 2', "at least two monomers",
                 id='one-monomer'),
    pytest.param('resid 1', 'type 3', "positive multiple of 2",
                 id='no-patches'),
    pytest.param('resid 1', 'id 11 12 13', "positive multiple of 2",
                 id='odd-patches'),
])
def test_rejects_malformed_groups(make_binder_universe, binding_frames,
                                  chain, patches, match):
    """A chain of one monomer, no patches, or an odd number of patches
    raises at construction."""
    universe = make_binder_universe(*binding_frames)

    with pytest.raises(ValueError, match=match):
        BinderStatistics(universe.select_atoms(chain),
                         universe.select_atoms(patches), cutoff=0.5,
                         topology='ring')


def test_rejects_patches_from_another_universe(make_binder_universe,
                                               binding_frames):
    """The patches must come from the chain's universe."""
    universe = make_binder_universe(*binding_frames)
    other = make_binder_universe(*binding_frames, name='other')

    with pytest.raises(ValueError, match="universe of 'chain'"):
        BinderStatistics(universe.select_atoms('resid 1'),
                         other.select_atoms('type 2'), cutoff=0.5,
                         topology='ring')


def test_periodic_direction_needs_a_box():
    """A frame without box dimensions cannot take the minimum image."""
    universe = mda.Universe.empty(4, trajectory=True)
    universe.atoms.positions = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.3, 0.0], [5.0, 5.0, 5.0]])
    analysis = BinderStatistics(universe.atoms[:2], universe.atoms[2:],
                                cutoff=0.5, topology='linear',
                                periodic=(True, True, True))

    with pytest.raises(ValueError, match="carries no box"):
        analysis.run()


# save


@pytest.mark.parametrize('store_distances, names', [
    pytest.param(False, [
        'nBoundTHnsPatch', 'nFreeTHnsPatch', 'nEngagedTHnsPatch',
        'nFreeTHnsCore', 'nBridgeTHnsCore', 'nDangleTHnsCore',
        'nCisTHnsCore', 'nTransTHnsCore', 'loopLengthHistTMon'],
        id='without-distances'),
    pytest.param(True, [
        'nBoundTHnsPatch', 'nFreeTHnsPatch', 'nEngagedTHnsPatch',
        'nFreeTHnsCore', 'nBridgeTHnsCore', 'nDangleTHnsCore',
        'nCisTHnsCore', 'nTransTHnsCore', 'loopLengthHistTMon',
        'distMatTMonHnsPatch'], id='with-distances'),
])
def test_save_writes_the_research_names(make_binder_universe,
                                        binding_frames, tmp_path,
                                        store_distances, names):
    """save writes <name>-<series>.npy per key under the research code's
    names, holding the results."""
    universe = make_binder_universe(*binding_frames)
    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', store_distances=store_distances).run()

    paths = analysis.save(tmp_path, 'N10ens1-nucleoid')

    assert [path.name for path in paths] == [
        f'N10ens1-nucleoid-{name}.npy' for name in names]
    for path, key in zip(paths, analysis.keys):
        assert_array_equal(np.load(path), analysis.results[key])


def test_save_uses_the_labels(make_binder_universe, binding_frames,
                              tmp_path):
    """The chain, patch, and binder labels name the saved files."""
    universe = make_binder_universe(*binding_frames)
    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring', store_distances=True, chain_label='Dna',
        patch_label='Site', binder_label='Prot').run()

    paths = analysis.save(tmp_path, 'N10ens1-nucleoid')

    assert [path.name for path in paths][3:] == [
        'N10ens1-nucleoid-nFreeTProt.npy', 'N10ens1-nucleoid-nBridgeTProt.npy',
        'N10ens1-nucleoid-nDangleTProt.npy', 'N10ens1-nucleoid-nCisTProt.npy',
        'N10ens1-nucleoid-nTransTProt.npy',
        'N10ens1-nucleoid-loopLengthHistTDna.npy',
        'N10ens1-nucleoid-distMatTDnaSite.npy']
    assert paths[0].name == 'N10ens1-nucleoid-nBoundTSite.npy'


def test_save_before_run_raises(make_binder_universe, binding_frames,
                                tmp_path):
    """save needs the results of a run."""
    universe = make_binder_universe(*binding_frames)
    analysis = BinderStatistics(
        universe.select_atoms('resid 1'), universe.select_atoms('type 2'),
        cutoff=0.5, topology='ring')

    with pytest.raises(RuntimeError, match=r"Call run\(\) before save"):
        analysis.save(tmp_path, 'N10ens1-nucleoid')


# Real frames


@pytest.fixture
def real_binders(real_files, real_universe, study_spec):
    """Factory: the study's :meth:`BinderStatistics.from_spec` on a real Hns
    polymer dump, with the universe it reads; not yet run."""

    def _build(whole):
        files = real_files(whole)
        universe = real_universe(files.topology, files.polymer)
        return (BinderStatistics.from_spec(
            universe, study_spec(whole.study).PROBE), universe)

    return _build


def test_real_binding_counts_obey_their_identities(real_hns_whole,
                                                   real_binders):
    """In every real frame, each patch is free or engaged, each binder is
    free, dangling, or bridging, the loop histogram counts every distinct
    cis or trans pair, no more pairs than bridges form, and engaged patches
    make at least one contact each."""
    analysis, _ = real_binders(real_hns_whole)
    n_hns = real_hns_whole.n_hns

    results = analysis.run().results

    pairs = results['n_pairs_cis'] + results['n_pairs_trans']
    assert_array_equal(results['n_patches_free']
                       + results['n_patches_engaged'], np.full(50, 2 * n_hns))
    assert_array_equal(results['n_binders_free']
                       + results['n_binders_dangling']
                       + results['n_binders_bridging'], np.full(50, n_hns))
    assert_array_equal(results['loop_length_histogram'].sum(axis=1), pairs)
    assert np.all(pairs <= results['n_binders_bridging'])
    assert np.all(results['n_patch_monomer_contacts']
                  >= results['n_patches_engaged'])


def test_real_binder_distances_are_the_minimum_image_distances(
        real_hns_whole, real_binders, study_spec, trajectory_arrays,
        reference_distances):
    """The stored monomer-to-patch distances of every real frame are the
    minimum-image distances along the geometry's periodic directions."""
    probe = study_spec(real_hns_whole.study).PROBE
    analysis, universe = real_binders(real_hns_whole)

    distances = analysis.run().results['distances']

    arrays = trajectory_arrays(universe)
    expected = reference_distances(
        arrays['positions'][:, analysis.chain.indices],
        arrays['positions'][:, analysis.patches.indices], arrays['boxes'],
        probe.geometry)
    # float32 positions of at most about 150 carry a round-off of about
    # 1e-5 in a separation.
    assert_allclose(distances, expected, rtol=0, atol=5e-5)


def test_real_whole_without_hns_has_no_binders_to_count(real_nh0_whole,
                                                        real_binders):
    """A real Hns whole without H-NS proteins has no patches, which the
    binder analysis rejects."""
    with pytest.raises(ValueError, match="positive multiple of 2"):
        real_binders(real_nh0_whole)
