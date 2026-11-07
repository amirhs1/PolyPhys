"""
Tests for :mod:`polyphys.analyze.clusters`.

Every expected value is independent of the code under test. It comes from
one of these sources, named in the test that uses it:

H.  A hand count on a graph or contact matrix small enough to draw.

R1. The Boolean transitive closure of an adjacency matrix, computed in
    :func:`_warshall_closure` by Warshall's algorithm (S. Warshall, "A theorem
    on Boolean matrices", J. ACM 9, 11 (1962)): ``C[i, j]`` is true when a
    path joins ``i`` to ``j``. Two particles are in one cluster exactly when
    ``C[i, j]`` holds, so it is the reference for the cluster labels.

R2. Brute force: a double loop over the definition, in
    :func:`_brute_force_distances` (minimum image, one pair at a time),
    :func:`_walk_genomic_distance` (steps along the backbone, each way round
    a ring), :func:`_binder_links_by_definition` (the same-monomer and
    neighboring-monomer rules, one pair of binders at a time), and
    :func:`_sweep_frame_by_frame` (the per-frame kernels called at one grid
    point and one frame at a time, then averaged).

R3. :func:`numpy.histogram`, an independent implementation of binning.

P.  A property: conservation (sizes sum to the particle count, the handshake
    lemma, unit normalization), symmetry, idempotence, or invariance under
    relabeling.

Random inputs are drawn by Hypothesis from a seed, then built with
``numpy.random.default_rng(seed)``, so every example is reproducible from its
seed. Float comparisons use ``FLOAT64_RTOL`` where two computations differ
only in the order of float64 operations; the one binning tolerance is derived
where it is used.
"""

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, strategies as st
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.analyze.clusters import (
    BindingStats,
    FrameBinding,
    binder_contacts,
    binding_statistics,
    cluster_labels,
    contact_matrix,
    count_bonds,
    count_clusters,
    enforce_single_patch_contact,
    find_direct_contacts,
    foci_pair_info,
    genomic_distance,
    monomer_binder_contacts,
    pair_distance_histograms,
    pair_distances,
    pair_separation_density,
    split_binders,
    threshold_sweep,
)

# Relative and absolute tolerance for results that differ only by float64
# round-off over a handful of operations: far above the 2.2e-16 epsilon.
FLOAT64_RTOL = 1e-12
FLOAT64_ATOL = 1e-12


# ---------------------------------------------------------------------------
# Independent references (R1, R2)
# ---------------------------------------------------------------------------


def _warshall_closure(adjacency: np.ndarray) -> np.ndarray:
    """Boolean transitive closure of `adjacency` by Warshall's algorithm."""
    closure = adjacency.astype(bool) | np.eye(len(adjacency), dtype=bool)
    for k in range(len(adjacency)):
        closure |= closure[:, k, None] & closure[None, k, :]
    return closure


def _brute_force_distances(positions: np.ndarray, box: list) -> np.ndarray:
    """Pair distances with the minimum image, one pair and axis at a time."""
    n_particles = len(positions)
    dist = np.zeros((n_particles, n_particles))
    for i in range(n_particles):
        for j in range(n_particles):
            delta = positions[i] - positions[j]
            for axis, length in enumerate(box):
                if 0 < length < np.inf:
                    delta[axis] -= length * np.round(delta[axis] / length)
            dist[i, j] = np.sqrt(np.sum(delta**2))
    return dist


def _walk_genomic_distance(i: int, j: int, topology: str, n: int) -> int:
    """Bonds between monomers `i` and `j`, walking the backbone stepwise."""
    forward = 0
    k = i
    while k != j:
        k = (k + 1) % n
        forward += 1
    if topology == 'linear':
        return forward if i <= j else _walk_genomic_distance(j, i, 'linear', n)
    backward = 0
    k = i
    while k != j:
        k = (k - 1) % n
        backward += 1
    return min(forward, backward)


def _binder_links_by_definition(contacts: np.ndarray, topology: str):
    """Binder links from the definition, one pair of binders at a time."""
    n_monomers, n_binders = contacts.shape
    linked = np.eye(n_binders, dtype=bool)
    for a in range(n_binders):
        for b in range(n_binders):
            if np.any(contacts[:, a] & contacts[:, b]):
                linked[a, b] = True
    for k in range(n_monomers - 1):
        on_pair = np.flatnonzero(contacts[k] | contacts[k + 1])
        for a in on_pair:
            linked[a, on_pair] = True
    if topology == 'ring':
        on_pair = np.flatnonzero(contacts[0] | contacts[-1])
        for a in on_pair:
            linked[a, on_pair] = True
    return linked


def _sweep_frame_by_frame(
    distances: np.ndarray, cutoff: float, cis_threshold: int,
    topology: str, inclusive: bool,
) -> tuple[dict, np.ndarray]:
    """The frame means at one grid point, one frame at a time, from the
    per-frame kernels: the eight counts and the loop-length histogram."""
    n_monomers = distances.shape[1]
    stats = BindingStats(n_monomers, topology, cis_threshold)
    for matrix in distances:
        stats.update(find_direct_contacts(matrix, cutoff, inclusive))
    counts = {name: series.mean()
              for name, series in stats.as_dict().items()}
    return counts, stats.loop_length_histogram / len(distances)


def _graph_from_edges(n_particles: int, edges: list) -> np.ndarray:
    """Symmetric adjacency matrix with a true diagonal from an edge list."""
    adjacency = np.eye(n_particles, dtype=bool)
    for i, j in edges:
        adjacency[i, j] = adjacency[j, i] = True
    return adjacency


# ---------------------------------------------------------------------------
# Strategies and fixtures
# ---------------------------------------------------------------------------


@st.composite
def contact_graphs(draw, max_particles: int = 40) -> np.ndarray:
    """A random symmetric contact matrix with a true diagonal."""
    n_particles = draw(st.integers(min_value=1, max_value=max_particles))
    probability = draw(st.floats(min_value=0.0, max_value=1.0))
    seed = draw(st.integers(min_value=0, max_value=2**32 - 1))
    rng = np.random.default_rng(seed)
    adjacency = rng.random((n_particles, n_particles)) < probability
    adjacency |= adjacency.T
    np.fill_diagonal(adjacency, True)
    return adjacency


@st.composite
def monomer_binder_matrices(draw) -> np.ndarray:
    """A random monomer-binder contact matrix."""
    n_monomers = draw(st.integers(min_value=1, max_value=12))
    n_binders = draw(st.integers(min_value=1, max_value=8))
    probability = draw(st.floats(min_value=0.0, max_value=0.6))
    seed = draw(st.integers(min_value=0, max_value=2**32 - 1))
    rng = np.random.default_rng(seed)
    return rng.random((n_monomers, n_binders)) < probability


@st.composite
def distance_stacks(draw) -> np.ndarray:
    """Random monomer-to-patch distance matrices over a few frames, with
    the distances spread across the sweep grid of 0.3 to 1.0."""
    n_frames = draw(st.integers(min_value=1, max_value=4))
    n_monomers = draw(st.integers(min_value=1, max_value=8))
    n_binders = draw(st.integers(min_value=1, max_value=4))
    seed = draw(st.integers(min_value=0, max_value=2**32 - 1))
    rng = np.random.default_rng(seed)
    return rng.uniform(0.3, 1.0, size=(n_frames, n_monomers, 2 * n_binders))


@pytest.fixture
def sweep_distances() -> np.ndarray:
    """Two frames of four monomers and two binders on a linear chain, every
    distance 2.0 but these; hand-counted in the sweep tests (H).

    Frame 0: patch 0 is 0.5 from monomer 0, patch 1 is 0.6 from monomer 2,
    patch 2 is exactly 0.685 from monomer 1, patch 3 is 0.7 from monomer 3.
    Frame 1: patches 0 and 1 are 0.5 and 0.6 from monomer 0; binder 1 is far.
    So at 0.55 binder 0 dangles in both frames; from 0.6 it bridges monomers
    0 and 2 (two bonds) in frame 0 and holds monomer 0 with both patches in
    frame 1, a bridge with no pair; binder 1 engages one patch only at an
    inclusive 0.685.
    """
    distances = np.full((2, 4, 4), 2.0)
    distances[0, 0, 0] = 0.5
    distances[0, 2, 1] = 0.6
    distances[0, 1, 2] = 0.685
    distances[0, 3, 3] = 0.7
    distances[1, 0, 0] = 0.5
    distances[1, 0, 1] = 0.6
    return distances


@pytest.fixture
def regression_graph() -> np.ndarray:
    """The 7-node graph on which the single-pass merge stopped short (H)."""
    return _graph_from_edges(7, [(0, 6), (1, 6), (2, 3), (3, 4), (4, 6)])


@pytest.fixture
def hns_toy() -> np.ndarray:
    """Six monomers, four two-patch binders; hand-counted below (H).

    Binder 0: patch 0 on monomer 0, patch 1 on monomer 2 (bridge, 2 bonds on
    a ring of six, 2 on a line). Binder 1: patch 2 on monomers 3 and 4 (an
    artifact), patch 3 free (dangling). Binder 2: free. Binder 3: patch 6 on
    monomer 5, patch 7 on monomer 0 (bridge, 1 bond on the ring, 5 on a
    line).
    """
    contacts = np.zeros((6, 8), dtype=int)
    contacts[0, 0] = contacts[2, 1] = 1
    contacts[3, 2] = contacts[4, 2] = 1
    contacts[5, 6] = contacts[0, 7] = 1
    return contacts


# ---------------------------------------------------------------------------
# Clusters
# ---------------------------------------------------------------------------


def test_count_clusters_reaches_the_transitive_closure(regression_graph):
    """A merge order that defeated the single pass gives one 6 and one 1."""
    hist = count_clusters(regression_graph)
    assert_array_equal(hist, [0, 1, 0, 0, 0, 0, 1, 0])  # H


@pytest.mark.parametrize(
    ('n_particles', 'edges', 'expected'),
    [
        pytest.param(5, [], [0, 5, 0, 0, 0, 0], id='no-contacts'),
        pytest.param(
            4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            [0, 0, 0, 0, 1], id='complete-k4',
        ),
        pytest.param(5, [(0, 1), (1, 2), (2, 3)], [0, 1, 0, 0, 1, 0],
                     id='path-plus-isolated'),
        pytest.param(
            6, [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)],
            [0, 0, 0, 2, 0, 0, 0], id='two-triangles',
        ),
        pytest.param(4, [(0, 1), (0, 2), (0, 3)], [0, 0, 0, 0, 1],
                     id='star'),
        pytest.param(5, [(4, 3), (3, 2), (2, 1), (1, 0)],
                     [0, 0, 0, 0, 0, 1], id='chain-listed-backwards'),
    ],
)
def test_count_clusters_known_sizes(n_particles, edges, expected):
    """Hand-drawn graphs give their hand-counted size histograms (H)."""
    hist = count_clusters(_graph_from_edges(n_particles, edges))
    assert_array_equal(hist, expected)


@given(adjacency=contact_graphs())
def test_contact_matrix_matches_boolean_closure(adjacency):
    """The contact matrix equals the Warshall transitive closure (R1)."""
    assert_array_equal(contact_matrix(adjacency), _warshall_closure(adjacency))


@given(adjacency=contact_graphs())
def test_cluster_labels_partition_matches_boolean_closure(adjacency):
    """Two particles share a label exactly when the closure joins them (R1)."""
    n_clusters, labels = cluster_labels(adjacency)
    closure = _warshall_closure(adjacency)
    assert_array_equal(labels[:, None] == labels[None, :], closure)
    assert n_clusters == len(np.unique(labels))


@given(adjacency=contact_graphs())
def test_count_clusters_sizes_sum_to_n_particles(adjacency):
    """Cluster sizes weighted by their counts add up to the particles (P)."""
    hist = count_clusters(adjacency)
    sizes = np.arange(len(hist))
    assert hist.shape == (len(adjacency) + 1,)
    assert hist[0] == 0
    assert int(np.sum(sizes * hist)) == len(adjacency)


@given(adjacency=contact_graphs())
def test_contact_matrix_is_a_closure(adjacency):
    """The contact matrix is symmetric, reflexive, and its own closure (P)."""
    closure = contact_matrix(adjacency)
    assert_array_equal(closure, closure.T)
    assert np.all(np.diagonal(closure))
    assert_array_equal(closure.astype(int) @ closure.astype(int) > 0, closure)


@given(adjacency=contact_graphs(), seed=st.integers(0, 2**32 - 1))
def test_count_clusters_is_invariant_under_relabeling(adjacency, seed):
    """Permuting the particles leaves the size histogram unchanged (P)."""
    order = np.random.default_rng(seed).permutation(len(adjacency))
    permuted = adjacency[np.ix_(order, order)]
    assert_array_equal(count_clusters(permuted), count_clusters(adjacency))


def test_cluster_labels_of_no_particles():
    """An empty contact matrix has no clusters and no labels."""
    n_clusters, labels = cluster_labels(np.zeros((0, 0), dtype=bool))
    assert n_clusters == 0
    assert labels.shape == (0,)


@pytest.mark.parametrize(
    'bad',
    [
        pytest.param(np.ones((2, 3), dtype=bool), id='not-square'),
        pytest.param(np.array([[1, 1], [0, 1]]), id='asymmetric'),
        pytest.param(np.ones(3, dtype=bool), id='one-dimensional'),
    ],
)
def test_cluster_labels_rejects_a_bad_matrix(bad):
    """A non-square, asymmetric, or non-matrix input raises."""
    with pytest.raises(ValueError, match='square|symmetric|two-dimensional'):
        cluster_labels(bad)


# ---------------------------------------------------------------------------
# Bonds and direct contacts
# ---------------------------------------------------------------------------


def test_count_bonds_known_degrees():
    """A star of four has one particle with three bonds and three with one."""
    star = _graph_from_edges(4, [(0, 1), (0, 2), (0, 3)])
    assert_array_equal(count_bonds(star), [0, 3, 0, 1])  # H


@given(adjacency=contact_graphs())
def test_count_bonds_satisfies_the_handshake_lemma(adjacency):
    """Degrees sum to twice the contacts and counts to the particles (P)."""
    hist = count_bonds(adjacency)
    n_particles = len(adjacency)
    n_contacts = (int(adjacency.sum()) - n_particles) // 2
    assert hist.shape == (n_particles,)
    assert int(hist.sum()) == n_particles
    assert int(np.sum(np.arange(n_particles) * hist)) == 2 * n_contacts


def test_count_bonds_rejects_an_asymmetric_matrix():
    """A directed contact matrix has no well-defined degree and raises."""
    with pytest.raises(ValueError, match='symmetric'):
        count_bonds(np.array([[1, 1], [0, 1]]))


@pytest.mark.parametrize(
    ('inclusive', 'expected'),
    [
        pytest.param(True, [[1, 1], [1, 1]], id='inclusive'),
        pytest.param(False, [[1, 0], [0, 1]], id='strict'),
    ],
)
def test_find_direct_contacts_at_the_cutoff(inclusive, expected):
    """A pair exactly at the cutoff is a contact only when inclusive."""
    dist = np.array([[0.0, 1.5], [1.5, 0.0]])
    contacts = find_direct_contacts(dist, 1.5, inclusive=inclusive)
    assert contacts.dtype == np.bool_
    assert_array_equal(contacts, np.array(expected, dtype=bool))


@pytest.mark.parametrize('cutoff', [0.0, -1.0])
def test_find_direct_contacts_rejects_a_nonpositive_cutoff(cutoff):
    """A zero or negative cutoff raises."""
    with pytest.raises(ValueError, match='cutoff must be positive'):
        find_direct_contacts(np.zeros((2, 2)), cutoff)


# ---------------------------------------------------------------------------
# Pair distances
# ---------------------------------------------------------------------------


def test_pair_distances_known_values():
    """A 3-4-5 right triangle gives its side lengths (H)."""
    dist = pair_distances([[0.0, 0.0], [3.0, 0.0], [3.0, 4.0]])
    assert_allclose(dist, [[0, 3, 5], [3, 0, 4], [5, 4, 0]],
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)


@pytest.mark.parametrize(
    ('box', 'expected'),
    [
        pytest.param([10.0, 0.0], 0.2, id='periodic'),
        pytest.param([0.0, 0.0], 9.8, id='length-zero-is-not-periodic'),
        pytest.param([np.inf, 0.0], 9.8, id='length-inf-is-not-periodic'),
        pytest.param(None, 9.8, id='no-box'),
    ],
)
def test_pair_distances_applies_the_minimum_image(box, expected):
    """Points at 0.1 and 9.9 are 0.2 apart only through a periodic edge (H)."""
    dist = pair_distances([[0.1, 0.0], [9.9, 0.0]], box=box)
    assert dist[0, 1] == pytest.approx(expected, abs=FLOAT64_ATOL)


@pytest.mark.parametrize(
    'box',
    [
        pytest.param([10.0, 10.0, 10.0], id='cubic'),
        pytest.param([0.0, 0.0, 10.0], id='z-only'),
    ],
)
@pytest.mark.parametrize('seed', [1, 2, 3])
def test_pair_distances_matches_brute_force(box, seed):
    """The vectorized distances equal a pair-by-pair minimum image (R2)."""
    positions = np.random.default_rng(seed).uniform(-5.0, 15.0, (25, 3))
    dist = pair_distances(positions, box=box)
    assert_allclose(dist, _brute_force_distances(positions, box),
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)
    assert_allclose(dist, dist.T, rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)


@pytest.mark.parametrize(
    ('positions', 'box', 'message'),
    [
        pytest.param([1.0, 2.0], None, 'n_particles, n_dims',
                     id='one-dimensional-positions'),
        pytest.param([[0.0, 0.0]], [10.0], 'one length per dimension',
                     id='wrong-box-length'),
        pytest.param([[0.0, 0.0]], [-1.0, 0.0], 'non-negative',
                     id='negative-box'),
    ],
)
def test_pair_distances_rejects_bad_shapes(positions, box, message):
    """Positions without two axes or a box of the wrong shape raise."""
    with pytest.raises(ValueError, match=message):
        pair_distances(positions, box=box)


# ---------------------------------------------------------------------------
# Genomic distance and foci pairs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('topology', ['linear', 'ring'])
@pytest.mark.parametrize(
    ('n_backbone', 'i', 'j'),
    [
        pytest.param(10, 0, 9, id='ends'),
        pytest.param(10, 3, 3, id='same-monomer'),
        pytest.param(10, 2, 7, id='half-way'),
        pytest.param(9, 1, 6, id='odd-ring'),
        pytest.param(10, 8, 1, id='reversed-order'),
    ],
)
def test_genomic_distance_matches_walking_the_chain(
    topology, n_backbone, i, j
):
    """The genomic distance is the shortest walk along the backbone (R2)."""
    expected = _walk_genomic_distance(i, j, topology, n_backbone)
    assert_array_equal(genomic_distance([i], [j], topology, n_backbone),
                       [expected])


def test_genomic_distance_on_a_linear_chain_runs():
    """A linear chain gives |i - j|; the research code raised here (H)."""
    assert_array_equal(genomic_distance([0, 1], [3, 9], 'linear', 10), [3, 8])


def test_genomic_distance_rejects_an_unknown_topology():
    """A topology other than linear or ring raises."""
    with pytest.raises(ValueError, match="topology must be 'linear' or"):
        genomic_distance([0], [1], 'star', 5)


def test_genomic_distance_rejects_an_index_beyond_the_chain():
    """An index at or past n_backbone raises."""
    with pytest.raises(ValueError, match='indices in 0, ..., 4'):
        genomic_distance([0], [5], 'ring', 5)


@pytest.mark.parametrize(
    ('topology', 'expected'),
    [
        pytest.param('linear', [[0, 8, 8]], id='linear'),
        pytest.param('ring', [[0, 8, 2]], id='ring'),
    ],
)
def test_foci_pair_info_respects_the_topology(topology, expected):
    """Foci at 0 and 8 of 10 are 8 bonds apart on a line, 2 on a ring (H)."""
    assert_array_equal(foci_pair_info([0, 8], 10, topology), expected)


def test_foci_pair_info_lists_pairs_in_upper_triangle_order():
    """Three foci give three rows in numpy.triu_indices order (H)."""
    info = foci_pair_info([0, 4, 8], 10, 'ring')
    assert_array_equal(info, [[0, 4, 4], [0, 8, 2], [4, 8, 4]])


# ---------------------------------------------------------------------------
# Pair separation histograms and density
# ---------------------------------------------------------------------------


def test_pair_distance_histograms_known_counts():
    """Hand-placed separations land in their bins; the last edge counts (H)."""
    counts, edges = pair_distance_histograms(
        [[0.5, 1.5, 1.5, 3.0], [2.0, 2.0, 2.0, 2.0]], bin_size=1.0, r_max=3.0
    )
    assert_array_equal(counts, [[1, 2, 1], [0, 0, 4]])
    assert_array_equal(edges, [0.0, 1.0, 2.0, 3.0])


@pytest.mark.parametrize('seed', [11, 12, 13])
def test_pair_distance_histograms_matches_numpy_histogram(seed):
    """Counts agree with numpy.histogram on the same edges (R3)."""
    series = np.random.default_rng(seed).uniform(0.0, 7.5, (4, 300))
    counts, edges = pair_distance_histograms(series, bin_size=0.5, r_max=7.5)
    reference = [np.histogram(row, bins=edges)[0] for row in series]
    assert_array_equal(counts, reference)
    assert_array_equal(counts.sum(axis=1), [300, 300, 300, 300])


def test_pair_distance_histograms_rounds_the_grid_up_to_whole_bins():
    """A range that is not a multiple of the bin size gets one more bin."""
    counts, edges = pair_distance_histograms([0.0], bin_size=1.0, r_max=2.5)
    assert_array_equal(edges, [0.0, 1.0, 2.0, 3.0])
    assert counts.shape == (3,)


@pytest.mark.parametrize(
    ('bin_size', 'r_max', 'n_bins'),
    [
        pytest.param(0.3, 2.1, 7, id='ratio-above-whole'),
        pytest.param(0.03, 15.9, 530, id='edge-below-r-max'),
    ],
)
def test_pair_distance_histograms_counts_whole_bins_without_round_off(
    bin_size, r_max, n_bins
):
    """A range that is a whole number of bins in decimal arithmetic gets
    exactly that many, ending at `r_max`, although its float64 ratio or
    product misses the whole number by an ulp."""
    counts, edges = pair_distance_histograms([r_max], bin_size, r_max)
    assert counts.shape == (n_bins,)
    assert edges[-1] == r_max
    assert counts[-1] == 1


def test_pair_distance_histograms_rejects_a_separation_beyond_the_grid():
    """A separation above the last edge raises instead of being dropped."""
    with pytest.raises(ValueError, match='must lie in'):
        pair_distance_histograms([0.5, 3.5], bin_size=1.0, r_max=3.0)


@pytest.mark.parametrize(
    ('bin_size', 'r_max', 'message'),
    [
        pytest.param(0.0, 3.0, 'bin_size must be positive', id='bin-size'),
        pytest.param(1.0, 0.0, 'r_max must be positive', id='r-max'),
    ],
)
def test_pair_distance_histograms_rejects_a_nonpositive_grid(
    bin_size, r_max, message
):
    """A zero bin size or range raises."""
    with pytest.raises(ValueError, match=message):
        pair_distance_histograms([0.5], bin_size=bin_size, r_max=r_max)


def test_pair_separation_density_known_value():
    """All frames in one shell give one over that shell's volume (H)."""
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    density = pair_separation_density([0, 0, 4], edges)
    shell = 4.0 * np.pi / 3.0 * (27.0 - 8.0)
    assert_allclose(density, [0.0, 0.0, 1.0 / shell],
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)


@pytest.mark.parametrize('seed', [21, 22])
def test_pair_separation_density_integrates_to_one_and_keeps_the_moment(seed):
    """The density integrates to 1 and reproduces the mean square (P)."""
    bin_size, r_max = 0.25, 8.0
    series = np.random.default_rng(seed).uniform(0.0, r_max, (3, 400))
    counts, edges = pair_distance_histograms(series, bin_size, r_max)
    density = pair_separation_density(counts, edges)
    shells = 4.0 * np.pi / 3.0 * np.diff(edges**3)
    centers = 0.5 * (edges[1:] + edges[:-1])
    assert_allclose(np.sum(density * shells, axis=1), [1.0, 1.0, 1.0],
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)
    # Replacing each separation r by its bin center changes r**2 by at most
    # (bin_size / 2) * 2 r_max, so the moment agrees to bin_size * r_max.
    assert_allclose(np.sum(density * shells * centers**2, axis=1),
                    np.mean(series**2, axis=1), atol=bin_size * r_max)


def test_pair_separation_density_rejects_an_empty_row():
    """A pair with no counts has no density and raises."""
    with pytest.raises(ValueError, match='positive total'):
        pair_separation_density([[1, 0], [0, 0]], [0.0, 1.0, 2.0])


@pytest.mark.parametrize(
    'edges',
    [
        pytest.param([0.0, 1.0], id='wrong-length'),
        pytest.param([0.0, 2.0, 1.0], id='not-increasing'),
        pytest.param([-1.0, 0.0, 1.0], id='negative'),
    ],
)
def test_pair_separation_density_rejects_bad_edges(edges):
    """Edges of the wrong length, order, or sign raise."""
    with pytest.raises(ValueError, match='edges must'):
        pair_separation_density([1, 1], edges)


# ---------------------------------------------------------------------------
# Binding: contact-matrix reshaping
# ---------------------------------------------------------------------------


def test_enforce_single_patch_contact_keeps_the_lowest_monomer():
    """A patch on monomers 1 and 2 keeps monomer 1; a triple keeps one (H)."""
    contacts = np.array([[0, 1, 0],
                         [1, 0, 1],
                         [1, 0, 1],
                         [0, 0, 1]])
    assert_array_equal(
        enforce_single_patch_contact(contacts).astype(int),
        [[0, 1, 0], [1, 0, 1], [0, 0, 0], [0, 0, 0]],
    )


@given(contacts=monomer_binder_matrices())
def test_enforce_single_patch_contact_is_idempotent_with_one_per_column(
    contacts
):
    """Each column keeps at most one contact; a second pass is a no-op (P)."""
    once = enforce_single_patch_contact(contacts)
    assert np.all(once.sum(axis=0) <= 1)
    assert_array_equal(once.sum(axis=0) > 0, contacts.sum(axis=0) > 0)
    assert_array_equal(enforce_single_patch_contact(once), once)


def test_monomer_binder_contacts_collapses_each_binder_s_patches():
    """A binder touches a monomer when any of its patches does (H)."""
    contacts = np.array([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 0, 1]])
    assert_array_equal(monomer_binder_contacts(contacts, 2).astype(int),
                       [[1, 0], [1, 0], [0, 1]])


def test_monomer_binder_contacts_rejects_a_wrong_patch_count():
    """A patch count that does not divide the columns raises."""
    with pytest.raises(ValueError, match='must be a multiple of'):
        monomer_binder_contacts(np.zeros((2, 3)), 2)


def test_split_binders_by_the_number_of_monomers_bound():
    """One contact dangles, two bridge, zero or three go to neither (H)."""
    contacts = np.array([[0, 1, 1, 1],
                         [0, 0, 1, 1],
                         [0, 0, 0, 1]])
    dangling, bridging = split_binders(contacts)
    assert_array_equal(dangling.astype(int),
                       [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    assert_array_equal(bridging.astype(int),
                       [[0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 0]])


@pytest.mark.parametrize('topology', ['linear', 'ring'])
@given(contacts=monomer_binder_matrices())
def test_binder_contacts_matches_the_definition(topology, contacts):
    """Binders are linked exactly as the pairwise definition says (R2)."""
    assert_array_equal(binder_contacts(contacts, topology),
                       _binder_links_by_definition(contacts, topology))


@pytest.mark.parametrize(
    ('topology', 'expected'),
    [
        pytest.param('ring', [[1, 1], [1, 1]], id='ring-wraps'),
        pytest.param('linear', [[1, 0], [0, 1]], id='linear-does-not'),
    ],
)
def test_binder_contacts_wraps_only_on_a_ring(topology, expected):
    """Binders on the first and last monomer are linked only on a ring (H)."""
    contacts = np.array([[1, 0], [0, 0], [0, 0], [0, 1]])
    assert_array_equal(binder_contacts(contacts, topology).astype(int),
                       expected)


# ---------------------------------------------------------------------------
# Binding statistics
# ---------------------------------------------------------------------------


def test_binding_statistics_counts_the_toy_on_a_ring(hns_toy):
    """The toy gives the hand-counted patch, binder, and pair counts (H)."""
    frame = binding_statistics(hns_toy, 'ring', cis_threshold=1)
    assert isinstance(frame, FrameBinding)
    assert frame.n_patch_monomer_contacts == 6
    assert frame.n_patches_engaged == 5
    assert frame.n_patches_free == 3
    assert frame.n_binders_free == 1
    assert frame.n_binders_dangling == 1
    assert frame.n_binders_bridging == 2
    assert frame.n_pairs_cis == 1
    assert frame.n_pairs_trans == 1
    assert_array_equal(frame.bridged_pairs, [[0, 2, 2], [0, 5, 1]])


def test_binding_statistics_counts_the_toy_on_a_linear_chain(hns_toy):
    """On a line the 5-0 bridge is 5 bonds long, so both bridges are trans."""
    frame = binding_statistics(hns_toy, 'linear', cis_threshold=1)
    assert (frame.n_binders_free, frame.n_binders_dangling,
            frame.n_binders_bridging) == (1, 1, 2)
    assert (frame.n_pairs_cis, frame.n_pairs_trans) == (0, 2)
    assert_array_equal(frame.bridged_pairs, [[0, 2, 2], [0, 5, 5]])  # H


def test_binding_statistics_counts_a_shared_pair_once():
    """Two binders bridging the same monomers give one bridged pair (H)."""
    contacts = np.array([[1, 0, 1, 0],
                         [0, 1, 0, 1]])
    frame = binding_statistics(contacts, 'linear', cis_threshold=1)
    assert frame.n_binders_bridging == 2
    assert frame.n_pairs_cis == 1
    assert_array_equal(frame.bridged_pairs, [[0, 1, 1]])


def test_binding_statistics_bridge_on_one_monomer_has_no_pair():
    """A binder with both patches on one monomer bridges but joins no pair."""
    frame = binding_statistics([[1, 1]], 'linear', cis_threshold=1)
    assert frame.n_binders_bridging == 1
    assert (frame.n_pairs_cis, frame.n_pairs_trans) == (0, 0)
    assert frame.bridged_pairs.shape == (0, 3)


def test_binding_statistics_rejects_an_odd_number_of_patches():
    """A matrix with an odd number of columns is not two patches per binder."""
    with pytest.raises(ValueError, match='two columns per binder'):
        binding_statistics(np.zeros((2, 3)), 'ring', cis_threshold=1)


def test_binding_statistics_rejects_a_negative_threshold():
    """A negative cis threshold raises."""
    with pytest.raises(ValueError, match='cis_threshold must be non-negative'):
        binding_statistics(np.zeros((2, 2)), 'ring', cis_threshold=-1)


def test_binding_stats_accumulates_frames_and_loop_lengths(hns_toy):
    """Two toy frames and an empty one give three rows and twice the loops."""
    stats = BindingStats(n_monomers=6, topology='ring', cis_threshold=1)
    stats.update(hns_toy)
    stats.update(hns_toy)
    stats.update(np.zeros((6, 8), dtype=int))
    series = stats.as_dict()
    assert stats.n_frames == 3
    assert_array_equal(series['n_binders_bridging'], [2, 2, 0])
    assert_array_equal(series['n_binders_free'], [1, 1, 4])
    assert_array_equal(series['n_pairs_cis'], [1, 1, 0])
    assert_array_equal(series['n_pairs_trans'], [1, 1, 0])
    assert_array_equal(stats.loop_length_histogram, [0, 2, 2, 0])  # H


@pytest.mark.parametrize(
    ('topology', 'expected_length'),
    [
        pytest.param('linear', 6, id='linear'),
        pytest.param('ring', 4, id='ring'),
    ],
)
def test_binding_stats_histogram_spans_the_longest_loop(
    topology, expected_length
):
    """The histogram has one bin per possible genomic distance."""
    stats = BindingStats(n_monomers=6, topology=topology, cis_threshold=1)
    assert stats.loop_length_histogram.shape == (expected_length,)
    assert stats.n_frames == 0


def test_binding_stats_rejects_a_frame_from_a_longer_chain(hns_toy):
    """A frame whose loop is longer than the chain allows raises."""
    stats = BindingStats(n_monomers=4, topology='linear', cis_threshold=1)
    frame = binding_statistics(hns_toy, 'linear', cis_threshold=1)
    with pytest.raises(ValueError, match='beyond the chain'):
        stats.append(frame)


def test_binding_stats_rejects_a_matrix_with_the_wrong_monomer_count():
    """A contact matrix with another number of rows raises on update."""
    stats = BindingStats(n_monomers=4, topology='linear', cis_threshold=1)
    with pytest.raises(ValueError, match='n_monomers = 4'):
        stats.update(np.zeros((6, 2), dtype=int))


@pytest.mark.parametrize(
    ('n_monomers', 'topology', 'cis_threshold', 'message'),
    [
        pytest.param(0, 'ring', 1, 'n_monomers must be positive',
                     id='no-monomers'),
        pytest.param(4, 'star', 1, "topology must be 'linear'",
                     id='topology'),
        pytest.param(4, 'ring', -1, 'cis_threshold must be non-negative',
                     id='threshold'),
    ],
)
def test_binding_stats_rejects_bad_parameters(
    n_monomers, topology, cis_threshold, message
):
    """An empty chain, an unknown topology, or a negative threshold raises."""
    with pytest.raises(ValueError, match=message):
        BindingStats(n_monomers=n_monomers, topology=topology,
                     cis_threshold=cis_threshold)


# ---------------------------------------------------------------------------
# Shape and range checks on the remaining inputs
# ---------------------------------------------------------------------------


def test_find_direct_contacts_rejects_a_non_matrix():
    """A one-dimensional distance input raises."""
    with pytest.raises(ValueError, match='two-dimensional'):
        find_direct_contacts(np.zeros(3), 1.0)


def test_count_bonds_rejects_a_non_square_matrix():
    """A rectangular matrix has no per-particle degree and raises."""
    with pytest.raises(ValueError, match='square'):
        count_bonds(np.ones((2, 3), dtype=bool))


def test_genomic_distance_rejects_an_empty_chain():
    """A chain of no monomers raises."""
    with pytest.raises(ValueError, match='n_backbone must be positive'):
        genomic_distance([0], [0], 'ring', 0)


def test_foci_pair_info_rejects_a_matrix_of_positions():
    """Foci positions must be one-dimensional."""
    with pytest.raises(ValueError, match='one-dimensional'):
        foci_pair_info([[0, 1], [2, 3]], 10, 'ring')


def test_pair_distance_histograms_rejects_three_dimensions():
    """A separation array with more than two axes raises."""
    with pytest.raises(ValueError, match='n_pairs, n_frames'):
        pair_distance_histograms(np.zeros((2, 2, 2)), 1.0, 3.0)


def test_monomer_binder_contacts_rejects_a_nonpositive_patch_count():
    """Zero patches per binder raises."""
    with pytest.raises(ValueError, match='patches_per_binder must be'):
        monomer_binder_contacts(np.zeros((2, 2)), 0)


# ---------------------------------------------------------------------------
# Threshold sweep
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('cutoff', 'cis_threshold', 'expected'),
    [
        # (contacts, engaged, free patches, free, dangling, bridging,
        #  cis, trans), each the mean over the two frames.
        pytest.param(0.55, 1, (1.0, 1.0, 3.0, 1.0, 1.0, 0.0, 0.0, 0.0),
                     id='dangling-below-0.6'),
        pytest.param(0.65, 1, (2.0, 2.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.5),
                     id='trans-at-cis-1'),
        pytest.param(0.65, 2, (2.0, 2.0, 2.0, 1.0, 0.0, 1.0, 0.5, 0.0),
                     id='cis-at-cis-2'),
        pytest.param(0.685, 2, (2.5, 2.5, 1.5, 0.5, 0.5, 1.0, 0.5, 0.0),
                     id='binder-1-engages-at-the-cutoff'),
    ],
)
def test_threshold_sweep_counts_the_hand_frames(
    sweep_distances, cutoff, cis_threshold, expected
):
    """Each row holds the frame means of the eight hand counts (H)."""
    counts, _ = threshold_sweep(
        sweep_distances, [0.55, 0.65, 0.685], [1, 2], 'linear',
        inclusive=True)
    row = counts[(counts['cutoff'] == cutoff)
                 & (counts['cis_threshold'] == cis_threshold)]
    assert len(row) == 1
    # The means of small integers over two frames are exact in float64.
    assert_allclose(row.iloc[0, 2:].to_numpy(dtype=float), expected,
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)


@pytest.mark.parametrize(
    ('inclusive', 'expected'),
    [
        pytest.param(True, 2.5, id='inclusive'),
        pytest.param(False, 2.0, id='strict'),
    ],
)
def test_threshold_sweep_inclusive_decides_a_pair_at_the_cutoff(
    sweep_distances, inclusive, expected
):
    """A patch exactly at the cutoff is a contact only when inclusive, so
    the mean contact count at 0.685 is 2.5 inclusive and 2.0 strict (H)."""
    counts, _ = threshold_sweep(
        sweep_distances, [0.685], [1], 'linear', inclusive=inclusive)
    assert counts.loc[0, 'n_patch_monomer_contacts'] == pytest.approx(
        expected, abs=FLOAT64_ATOL)


@given(distances=distance_stacks(), inclusive=st.booleans(),
       topology=st.sampled_from(['linear', 'ring']))
def test_threshold_sweep_matches_the_per_frame_kernels(
    distances, inclusive, topology
):
    """Every row and histogram row equals the per-frame kernels run at that
    grid point and averaged over the frames (R2)."""
    cutoffs = [0.5, 0.7, 0.9]
    cis_thresholds = [1, 3]
    counts, histogram = threshold_sweep(
        distances, cutoffs, cis_thresholds, topology, inclusive=inclusive)
    for row, (_, record) in enumerate(counts.iterrows()):
        expected, loops = _sweep_frame_by_frame(
            distances, record['cutoff'], int(record['cis_threshold']),
            topology, inclusive)
        for name, value in expected.items():
            # Both sides are an integer sum divided by n_frames.
            assert record[name] == pytest.approx(value, abs=FLOAT64_ATOL)
        assert_allclose(histogram[row], loops,
                        rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)


@given(distances=distance_stacks(), inclusive=st.booleans())
def test_threshold_sweep_contacts_grow_with_the_cutoff(distances, inclusive):
    """Along an increasing grid the contacts and engaged patches never fall
    and the free binders never rise (P: monotonicity of ``dist <= cutoff``)."""
    counts, _ = threshold_sweep(
        distances, [0.4, 0.6, 0.8, 1.0], [2], 'ring', inclusive=inclusive)
    assert np.all(np.diff(counts['n_patch_monomer_contacts']) >= 0)
    assert np.all(np.diff(counts['n_patches_engaged']) >= 0)
    assert np.all(np.diff(counts['n_binders_free']) <= 0)


@given(distances=distance_stacks(),
       topology=st.sampled_from(['linear', 'ring']))
def test_threshold_sweep_cis_threshold_only_splits_the_pairs(
    distances, topology
):
    """Across the cis thresholds of one cutoff, the six other counts, the
    sum of cis and trans pairs, and the histogram row are the same (P)."""
    counts, histogram = threshold_sweep(
        distances, [0.6], [0, 1, 2, 5], topology, inclusive=True)
    fixed = ['n_patch_monomer_contacts', 'n_patches_engaged',
             'n_patches_free', 'n_binders_free', 'n_binders_dangling',
             'n_binders_bridging']
    assert counts[fixed].nunique().eq(1).all()
    pairs = counts['n_pairs_cis'] + counts['n_pairs_trans']
    # The two means split one integer sum, so their sum differs across rows
    # only by the round-off of two float64 divisions.
    assert_allclose(pairs, pairs.iloc[0], rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)
    assert_array_equal(histogram, np.broadcast_to(histogram[0],
                                                  histogram.shape))


@pytest.mark.parametrize(
    ('topology', 'n_lengths'),
    [
        pytest.param('linear', 4, id='linear'),
        pytest.param('ring', 3, id='ring'),
    ],
)
def test_threshold_sweep_histogram_counts_the_bridged_pairs(
    sweep_distances, topology, n_lengths
):
    """The histogram has one column per genomic distance of the topology,
    an empty column 0, and rows that sum to the cis and trans pairs (P)."""
    counts, histogram = threshold_sweep(
        sweep_distances, [0.55, 0.65], [1, 2], topology, inclusive=True)
    assert histogram.dtype == np.float64
    assert histogram.shape == (4, n_lengths)
    assert_array_equal(histogram[:, 0], 0.0)
    assert_allclose(histogram.sum(axis=1),
                    counts['n_pairs_cis'] + counts['n_pairs_trans'],
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)


def test_threshold_sweep_keeps_the_grid_order_without_repeats(
    sweep_distances
):
    """Rows follow the cutoffs as given, cis thresholds inside each, once
    each; the research script's cumulative tables repeated earlier rows."""
    counts, histogram = threshold_sweep(
        sweep_distances, [0.65, 0.55], [4, 2], 'linear', inclusive=True)
    assert list(counts.columns[:2]) == ['cutoff', 'cis_threshold']
    assert_allclose(counts['cutoff'], [0.65, 0.65, 0.55, 0.55],
                    rtol=FLOAT64_RTOL)
    assert_array_equal(counts['cis_threshold'], [4, 2, 4, 2])
    assert counts.shape == (4, 10)
    assert histogram.shape[0] == 4


def test_threshold_sweep_accepts_whole_float_cis_thresholds():
    """A cis threshold given as a whole float counts as that integer."""
    distances = np.full((2, 3, 2), 2.0)
    distances[0, 0, 0] = 0.5
    distances[0, 2, 1] = 0.6
    as_ints, _ = threshold_sweep(distances, [0.65], [1, 2], 'linear',
                                 inclusive=True)
    as_floats, _ = threshold_sweep(distances, [0.65], [1.0, 2.0], 'linear',
                                   inclusive=True)
    pd.testing.assert_frame_equal(as_floats, as_ints)


@pytest.mark.parametrize(
    ('distances', 'cutoffs', 'cis_thresholds', 'topology', 'message'),
    [
        pytest.param(np.ones((3, 4)), [0.5], [1], 'linear',
                     'n_frames, n_monomers, n_patches', id='two-dimensional'),
        pytest.param(np.ones((0, 3, 4)), [0.5], [1], 'linear',
                     'at least one frame', id='no-frames'),
        pytest.param(np.ones((1, 3, 4)), [], [1], 'linear',
                     'cutoffs must not be empty', id='no-cutoffs'),
        pytest.param(np.ones((1, 3, 4)), [0.5, 0.0], [1], 'linear',
                     'cutoff must be positive', id='zero-cutoff'),
        pytest.param(np.ones((1, 3, 4)), [0.5], [], 'linear',
                     'cis_thresholds must not be empty', id='no-cis'),
        pytest.param(np.ones((1, 3, 4)), [0.5], [1, -1], 'linear',
                     'non-negative', id='negative-cis'),
        pytest.param(np.ones((1, 3, 3)), [0.5], [1], 'linear',
                     'two patches per binder', id='odd-patches'),
        pytest.param(np.ones((1, 3, 4)), [0.5], [2.7, 0.5], 'linear',
                     'whole number of bonds', id='fractional-cis'),
        pytest.param(np.ones((1, 3, 4)), [0.5], [1], 'star',
                     'topology', id='bad-topology'),
    ],
)
def test_threshold_sweep_rejects_bad_input(
    distances, cutoffs, cis_thresholds, topology, message
):
    """A malformed array, grid, or topology raises before any counting."""
    with pytest.raises(ValueError, match=message):
        threshold_sweep(distances, cutoffs, cis_thresholds, topology,
                        inclusive=True)
