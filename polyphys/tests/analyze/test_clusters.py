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
    a ring), and :func:`_binder_links_by_definition` (the same-monomer and
    neighbouring-monomer rules, one pair of binders at a time).

R3. :func:`numpy.histogram`, an independent implementation of binning.

P.  A property: conservation (sizes sum to the particle count, the handshake
    lemma, unit normalization), symmetry, idempotence, or invariance under
    relabelling.

Random inputs are drawn by Hypothesis from a seed, then built with
``numpy.random.default_rng(seed)``, so every example is reproducible from its
seed. Float comparisons use ``FLOAT64_RTOL`` where two computations differ
only in the order of float64 operations; the one binning tolerance is derived
where it is used.
"""

import numpy as np
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


@pytest.fixture
def regression_graph() -> np.ndarray:
    """The 7-node graph on which the single-pass merge stopped short (H)."""
    return _graph_from_edges(7, [(0, 6), (1, 6), (2, 3), (3, 4), (4, 6)])


@pytest.fixture
def hns_toy() -> np.ndarray:
    """Six monomers, four two-patch binders; hand-counted below (H).

    Binder 0: patch 0 on monomer 0, patch 1 on monomer 2 (bridge, 2 bonds on
    a ring of six, 2 on a line). Binder 1: patch 2 on monomers 3 and 4 (an
    artefact), patch 3 free (dangling). Binder 2: free. Binder 3: patch 6 on
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
def test_count_clusters_is_invariant_under_relabelling(adjacency, seed):
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
    centres = 0.5 * (edges[1:] + edges[:-1])
    assert_allclose(np.sum(density * shells, axis=1), [1.0, 1.0, 1.0],
                    rtol=FLOAT64_RTOL, atol=FLOAT64_ATOL)
    # Replacing each separation r by its bin centre changes r**2 by at most
    # (bin_size / 2) * 2 r_max, so the moment agrees to bin_size * r_max.
    assert_allclose(np.sum(density * shells * centres**2, axis=1),
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
