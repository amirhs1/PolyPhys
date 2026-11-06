"""
Tests for :mod:`polyphys.analyze.density`.

Every expected value is independent of the code under test. It comes from
one of these sources, named in the test that uses it:

A1. Closed-form volumes: a spherical shell ``4 pi (r2**3 - r1**3) / 3``, a
    cylindrical shell ``pi L (r2**2 - r1**2)``, a slab ``pi D**2 dz / 4``,
    and a bead ``pi a**3 / 6``.

A2. The spherical segment and the sphere-sphere intersection volume
    (Weisstein, MathWorld, "Spherical Segment" eq. 14 and "Sphere-Sphere
    Intersection" eq. 16), through :mod:`polyphys.analyze.measurer`, whose
    own tests hold their analytic values.

A3. The volume of a ball of radius ``R_b`` inside a coaxial cylinder of
    radius ``R < R_b``: ``int_0^R 2 sqrt(R_b**2 - s**2) 2 pi s ds =
    4 pi [R_b**3 - (R_b**2 - R**2)**1.5] / 3``.

A4. The averaged slab shares of a bead of diameter 1 in slabs of width 1:
    the cap beyond the next plane has height ``t = x - 1/2`` for a center at
    ``x > 1/2``, so the share into a neighbor is ``int_0^{1/2} pi t**2 (3/2
    - t) / 3 dt = pi / 64`` and the share into the bead's own slab is
    ``pi / 6 - 2 pi / 64 = 13 pi / 96``.

A5. An ideal gas: ``N`` points uniform in a volume ``V`` give a bin count
    that is Poisson with mean ``N V_i / V``, so the density ``N / V`` holds
    within ``5 sqrt(N V_i / V) / V_i`` in every bin.

R1. The research code's volume shares, restated in
    :func:`_sciops_concentric_shares`: for a bead at each bin center, the
    difference of the sphere-sphere volumes at consecutive edges.

R2. A Monte Carlo estimate of a ball-in-cylinder volume from
    ``numpy.random.default_rng(seed)`` with the seed shown, within five
    binomial standard deviations.

R3. The research code's axial slab volume. Its ``SpatialDistribution``
    integrates ``pi dcyl**2 / 4`` along ``z`` but passes the cylinder length
    ``lcyl`` as the ``dcyl`` argument, so each of its slabs has volume
    ``pi lcyl**2 dz / 4``. Its axial ``rho`` and ``phi`` therefore equal the
    closed-form ones (A1) times ``(dcyl / lcyl)**2``; its radial ones carry
    no such factor.

P.  A property: conservation of the total bead volume, the symmetry
    ``V_i S_ij = V_j S_ji``, translation invariance of the periodic slab
    matrix, a uniform fluid returning ``rho v``, or row-wise independence.

Float comparisons use ``FLOAT64_ATOL`` where two computations differ only
in the order of float64 operations; every other tolerance is derived where
it is used.
"""

import json
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.analyze.density import (
    Profiles,
    bin_volumes,
    density_profiles,
    normalize_profile,
    number_density,
    share_matrix,
    sphere_cylinder_intersection,
    volume_fraction,
)
from polyphys.analyze.measurer import (
    sphere_sphere_intersection,
    spherical_segment,
)

# Absolute tolerance for results that differ only by float64 round-off over
# sums of a few hundred terms: far above the 2.2e-16 epsilon, far below any
# value the tests compare.
FLOAT64_ATOL = 1e-12


def _bead_volume(diameter: float) -> float:
    """Return ``pi a**3 / 6`` (A1)."""
    return np.pi * diameter**3 / 6.0


# ---------------------------------------------------------------------------
# Reference computations
# ---------------------------------------------------------------------------


def _sciops_concentric_shares(edges: np.ndarray, r_bead: float) -> np.ndarray:
    """Return the research code's concentric shares, bead at each center (R1).

    ``distributions.py`` built, for the bead centered at ``centers[c]``, the
    share of bin ``k`` as the sphere-sphere volume at edge ``k + 1`` minus
    the one at edge ``k``, with the sphere-sphere formula standing in for
    the sphere-cylinder one along ``r`` in a cylinder.
    """
    centers = 0.5 * (edges[:-1] + edges[1:])
    shares = np.zeros((centers.size, centers.size))
    for c, center in enumerate(centers):
        previous = 0.0
        for k in range(1, edges.size):
            current = sphere_sphere_intersection(r_bead, edges[k], center)
            shares[c, k - 1] = current - previous
            previous = current
    return shares


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def make_uniform_counts():
    """Return a factory that bins ideal-gas points in one geometry (A5).

    ``make_uniform_counts(geometry, direction, edges, n_points, seed)``
    draws points uniform in the ball, the cylinder, or the periodic slab the
    edges span and returns their histogram along the direction.
    """

    def factory(geometry, direction, edges, n_points, seed):
        rng = np.random.default_rng(seed)
        outer = edges[-1]
        if (geometry, direction) == ('cubic', 'r'):
            coordinate = outer * np.cbrt(rng.uniform(size=n_points))
        elif (geometry, direction) == ('cylindrical', 'r'):
            coordinate = outer * np.sqrt(rng.uniform(size=n_points))
        else:
            coordinate = rng.uniform(edges[0], outer, size=n_points)
        counts, _ = np.histogram(coordinate, bins=edges)
        return counts

    return factory


# ---------------------------------------------------------------------------
# bin_volumes (T1, T2, D6)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'geometry, direction, kwargs, expected',
    [
        pytest.param(
            'cubic', 'r', {}, [4 * np.pi / 3, 28 * np.pi / 3],
            id='spherical-shells',
        ),
        pytest.param(
            'cylindrical', 'r', {'lcyl': 3.0}, [3 * np.pi, 9 * np.pi],
            id='cylindrical-shells',
        ),
        pytest.param(
            'cylindrical', 'z', {'dcyl': 4.0}, [4 * np.pi, 4 * np.pi],
            id='slabs',
        ),
    ],
)
def test_bin_volumes_match_closed_forms(geometry, direction, kwargs, expected):
    """Each scheme returns the closed-form shell or slab volume (A1)."""
    volumes = bin_volumes([0.0, 1.0, 2.0], geometry, direction, **kwargs)
    assert_allclose(volumes, expected, rtol=1e-15)  # float64 round-off


@pytest.mark.parametrize(
    'geometry, direction, kwargs, total',
    [
        pytest.param('cubic', 'r', {}, 4 * np.pi * 5.0**3 / 3,
                     id='ball'),
        pytest.param('cylindrical', 'r', {'lcyl': 7.0},
                     np.pi * 7.0 * 5.0**2, id='cylinder'),
        pytest.param('cylindrical', 'z', {'dcyl': 3.0},
                     np.pi * 3.0**2 * 5.0 / 4, id='slab-stack'),
    ],
)
def test_bin_volumes_sum_to_the_domain(geometry, direction, kwargs, total):
    """Unequal bins' volumes telescope to the whole domain's volume (A1)."""
    edges = np.array([0.0, 0.3, 1.1, 2.0, 3.7, 5.0])
    volumes = bin_volumes(edges, geometry, direction, **kwargs)
    assert volumes.sum() == pytest.approx(total, rel=1e-14)  # round-off


def test_bin_volumes_warns_beyond_half_the_cube():
    """Shells past half the cube side plus one bin raise a warning."""
    edges = np.arange(0.0, 8.5, 0.5)  # outermost 8 > 10/2 + 0.5
    with pytest.warns(UserWarning, match='beyond half the box side'):
        volumes = bin_volumes(edges, 'cubic', 'r', lcube=10.0)
    assert volumes.shape == (16,)


def test_bin_volumes_is_silent_within_half_the_cube():
    """Shells that end within half the side plus one bin raise no warning."""
    edges = np.arange(0.0, 5.51, 0.5)  # outermost 5.5 == 10/2 + 0.5
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        bin_volumes(edges, 'cubic', 'r', lcube=10.0)


# ---------------------------------------------------------------------------
# sphere_cylinder_intersection (T9, T10)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'r_cylinder, d, expected',
    [
        pytest.param(2.0, 0.5, _bead_volume(2.0), id='ball-inside'),
        pytest.param(2.0, 3.0, 0.0, id='ball-outside'),
        pytest.param(0.0, 0.0, 0.0, id='empty-cylinder'),
        pytest.param(
            0.6, 0.0, 4 * np.pi * (1 - (1 - 0.36) ** 1.5) / 3,
            id='on-axis-narrow-cylinder',
        ),
    ],
)
def test_sphere_cylinder_intersection_has_the_analytic_limits(
    r_cylinder, d, expected
):
    """The ball-in-cylinder volume matches its closed-form cases (A1, A3)."""
    volume = sphere_cylinder_intersection(1.0, r_cylinder, d)
    # The quadrature agrees with an adaptive reference to about 1e-8 of the
    # ball volume; the limits are reached exactly.
    assert volume == pytest.approx(expected, abs=1e-8)


def test_sphere_cylinder_intersection_tends_to_the_half_space():
    """A cylinder wide against the bead cuts it like a plane: half its volume.

    The planar limit is approached at first order in ``a / R`` (about
    ``0.09 a / R`` of the bead volume), so at ``R = 2000 a`` the two agree
    to better than ``1e-4`` relative.
    """
    volume = sphere_cylinder_intersection(0.5, 2000.0, 2000.0)
    assert volume == pytest.approx(_bead_volume(1.0) / 2, rel=1e-4)


def test_sphere_cylinder_intersection_agrees_with_monte_carlo():
    """The quadrature agrees with a seeded Monte Carlo estimate (R2)."""
    r_ball, r_cylinder, d = 0.5, 0.7, 0.4
    n_points = 400_000
    rng = np.random.default_rng(20261103)
    cube = rng.uniform(-r_ball, r_ball, size=(n_points, 3))
    in_ball = np.linalg.norm(cube, axis=1) < r_ball
    points = cube[in_ball]
    inside = np.hypot(points[:, 0] + d, points[:, 1]) < r_cylinder
    fraction = inside.mean()
    estimate = fraction * _bead_volume(2 * r_ball)
    # Five binomial standard deviations of the inside fraction.
    sigma = 5 * np.sqrt(fraction * (1 - fraction) / points.shape[0])
    volume = sphere_cylinder_intersection(r_ball, r_cylinder, d)
    assert volume == pytest.approx(estimate, abs=sigma * _bead_volume(1.0))


@pytest.mark.parametrize(
    'r_sphere, d, message',
    [
        pytest.param(0.0, 0.5, 'r_sphere must be positive', id='zero-ball'),
        pytest.param(1.0, -0.1, 'd must be non-negative', id='negative-d'),
    ],
)
def test_sphere_cylinder_intersection_rejects_bad_input(r_sphere, d, message):
    """A non-positive ball radius or a negative distance raises ValueError."""
    with pytest.raises(ValueError, match=message):
        sphere_cylinder_intersection(r_sphere, 1.0, d)


# ---------------------------------------------------------------------------
# share_matrix: analytic entries (T7, T8, T11)
# ---------------------------------------------------------------------------


def test_share_matrix_center_slab_entries_are_spherical_segments():
    """Centered slab shares are the segments between the slab planes (A2)."""
    edges = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    r_bead = 0.75
    share = share_matrix(edges, 2 * r_bead, 'cylindrical', 'z',
                         shares='center').toarray()
    center = 1.5
    expected = [
        spherical_segment(r_bead, edges[k] - center, edges[k + 1] - center)
        for k in range(3)
    ]
    assert_allclose(share[1, :3], expected, atol=FLOAT64_ATOL)
    assert share[1, 3] == 0.0


def test_share_matrix_average_slab_entries_match_the_cap_integral():
    """Averaged slab shares: ``13 pi / 96`` home, ``pi / 64`` beside (A4)."""
    share = share_matrix([0.0, 1.0, 2.0, 3.0, 4.0], 1.0, 'cylindrical',
                         'z').toarray()
    assert_allclose(
        share[1], [np.pi / 64, 13 * np.pi / 96, np.pi / 64, 0.0],
        atol=FLOAT64_ATOL,
    )


@pytest.mark.parametrize(
    'geometry, kwargs',
    [
        pytest.param('cubic', {}, id='spherical'),
        pytest.param('cylindrical', {'intersection': 'exact'},
                     id='cylindrical-exact'),
        pytest.param('cylindrical', {'intersection': 'sphere'},
                     id='cylindrical-sphere'),
    ],
)
def test_share_matrix_keeps_a_small_bead_in_its_first_shell(geometry, kwargs):
    """A centered bead reaching no edge of its shell lies wholly in it (A1).

    With edges 0, 1, 2 and a bead of radius 0.25, the bead at center 0.5
    spans 0.25 to 0.75 and the one at center 1.5 spans 1.25 to 1.75.
    """
    share = share_matrix([0.0, 1.0, 2.0], 0.5, geometry, 'r',
                         shares='center', **kwargs).toarray()
    v = _bead_volume(0.5)
    assert_allclose(share, [[v, 0.0], [0.0, v]],
                    atol=1e-8)  # the exact kernel's quadrature level


@pytest.mark.parametrize(
    'geometry',
    [pytest.param('cubic', id='cubic-r'),
     pytest.param('cylindrical', id='cylindrical-r')],
)
def test_share_matrix_center_sphere_reproduces_the_research_code(geometry):
    """The center mode with the sphere kernel equals the SciOps loop (R1)."""
    edges = np.arange(0.0, 6.01, 0.25)
    share = share_matrix(edges, 1.0, geometry, 'r', shares='center',
                         intersection='sphere').toarray()
    assert_allclose(share, _sciops_concentric_shares(edges, 0.5),
                    atol=FLOAT64_ATOL)


# ---------------------------------------------------------------------------
# share_matrix: properties (T4, T5, T6)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'geometry, direction, shares',
    [
        pytest.param('cubic', 'r', 'average', id='spherical-average'),
        pytest.param('cylindrical', 'r', 'average',
                     id='cylindrical-average'),
        pytest.param('cylindrical', 'z', 'average', id='slab-average'),
        pytest.param('cubic', 'r', 'center', id='spherical-center'),
        pytest.param('cylindrical', 'r', 'center', id='cylindrical-center'),
        pytest.param('cylindrical', 'z', 'center', id='slab-center'),
    ],
)
def test_share_matrix_rows_sum_to_the_bead_volume(geometry, direction, shares):
    """Every bin whose bead stays in the domain shares out exactly ``v`` (P).

    Along ``z`` every bin qualifies through the periodic images; along
    ``r`` the bins whose outer edge plus the bead radius stays within the
    last edge do.
    """
    edges = np.arange(0.0, 4.01, 0.25) if direction == 'r' else \
        np.arange(-2.0, 2.01, 0.25)
    share = share_matrix(edges, 1.0, geometry, direction, shares=shares)
    row_sums = np.asarray(share.sum(axis=1)).ravel()
    inside = edges[1:] + 0.5 <= edges[-1] if direction == 'r' else \
        np.ones(edges.size - 1, dtype=bool)
    assert_allclose(row_sums[inside], _bead_volume(1.0), rtol=1e-8)


@pytest.mark.parametrize(
    'geometry, direction, kwargs',
    [
        pytest.param('cubic', 'r', {}, id='spherical'),
        pytest.param('cylindrical', 'r', {'lcyl': 3.0}, id='cylindrical'),
        pytest.param('cylindrical', 'z', {'dcyl': 2.0}, id='slab'),
    ],
)
def test_share_matrix_average_is_volume_symmetric(geometry, direction, kwargs):
    """Averaged shares satisfy ``V_i S_ij = V_j S_ji`` in every scheme (P)."""
    edges = np.arange(0.0, 3.01, 0.25) if direction == 'r' else \
        np.arange(-1.5, 1.51, 0.25)
    volumes = bin_volumes(edges, geometry, direction, **kwargs)
    share = share_matrix(edges, 1.0, geometry, direction).toarray()
    weighted = volumes[:, None] * share
    # The exact cylinder kernel is a quadrature: 1e-8 of a bead volume.
    assert_allclose(weighted, weighted.T, atol=1e-8 * _bead_volume(1.0))


@pytest.mark.parametrize(
    'geometry, direction, kwargs',
    [
        pytest.param('cubic', 'r', {}, id='spherical'),
        pytest.param('cylindrical', 'r', {'lcyl': 3.0}, id='cylindrical'),
        pytest.param('cylindrical', 'z', {'dcyl': 2.0}, id='slab'),
    ],
)
def test_volume_fraction_conserves_the_bead_volume(geometry, direction,
                                                   kwargs):
    """``sum_i V_i phi_i`` equals the number of beads times ``v`` (P, A1).

    One bead sits in each bin whose bead stays inside the domain; the
    averaged shares return their total volume to round-off.
    """
    edges = np.arange(0.0, 4.01, 0.25) if direction == 'r' else \
        np.arange(-2.0, 2.01, 0.25)
    volumes = bin_volumes(edges, geometry, direction, **kwargs)
    inside = edges[1:] + 0.5 <= edges[-1] if direction == 'r' else \
        np.ones(edges.size - 1, dtype=bool)
    counts = inside.astype(float)
    share = share_matrix(edges, 1.0, geometry, direction)
    phi = volume_fraction(number_density(counts, volumes), share)
    total = (volumes * phi).sum()
    assert total == pytest.approx(counts.sum() * _bead_volume(1.0),
                                  rel=1e-8)  # quadrature level


def test_share_matrix_wraps_a_bead_across_the_periodic_boundary():
    """A bead in the last slab reaches into the first as much as into its
    inner neighbor, and the matrix is invariant under a cyclic shift (P)."""
    edges = np.arange(0.0, 4.01, 0.5)
    share = share_matrix(edges, 1.0, 'cylindrical', 'z').toarray()
    assert share[7, 0] > 0.0
    assert share[7, 0] == pytest.approx(share[7, 6], abs=FLOAT64_ATOL)
    shifted = np.roll(np.roll(share, 3, axis=0), 3, axis=1)
    assert_allclose(shifted, share, atol=FLOAT64_ATOL)


@pytest.mark.parametrize(
    'geometry, direction, kwargs, shares',
    [
        pytest.param('cubic', 'r', {}, 'average', id='spherical-average'),
        pytest.param('cylindrical', 'r', {'lcyl': 3.0}, 'average',
                     id='cylindrical-average'),
        pytest.param('cylindrical', 'z', {'dcyl': 2.0}, 'average',
                     id='slab-average'),
        pytest.param('cubic', 'r', {}, 'center', id='spherical-center'),
        pytest.param('cylindrical', 'z', {'dcyl': 2.0}, 'center',
                     id='slab-center'),
    ],
)
def test_volume_fraction_of_a_uniform_fluid_is_rho_v(geometry, direction,
                                                     kwargs, shares):
    """A constant density gives ``phi = rho v`` wherever the bead fits (P)."""
    edges = np.arange(0.0, 4.01, 0.25) if direction == 'r' else \
        np.arange(-2.0, 2.01, 0.25)
    rho = np.full(edges.size - 1, 0.3)
    share = share_matrix(edges, 1.0, geometry, direction, shares=shares)
    phi = volume_fraction(rho, share)
    inside = edges[1:] + 0.5 <= edges[-1] if direction == 'r' else \
        np.ones(edges.size - 1, dtype=bool)
    assert_allclose(phi[inside], 0.3 * _bead_volume(1.0), rtol=1e-8)


# ---------------------------------------------------------------------------
# number_density and the ideal gas (T3, T13, T15)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'geometry, direction, kwargs',
    [
        pytest.param('cubic', 'r', {}, id='ball'),
        pytest.param('cylindrical', 'r', {'lcyl': 2.0}, id='cylinder'),
        pytest.param('cylindrical', 'z', {'dcyl': 2.0}, id='slab-stack'),
    ],
)
def test_number_density_of_an_ideal_gas_is_flat(make_uniform_counts, geometry,
                                                direction, kwargs):
    """Uniform points give ``rho = N / V`` flat to counting noise (A5)."""
    edges = np.arange(0.0, 5.01, 0.5) if direction == 'r' else \
        np.arange(-2.5, 2.51, 0.5)
    n_points = 200_000
    counts = make_uniform_counts(geometry, direction, edges, n_points,
                                 seed=20261103)
    volumes = bin_volumes(edges, geometry, direction, **kwargs)
    rho = number_density(counts, volumes, n_particles=n_points)
    bulk = n_points / volumes.sum()
    # Five Poisson standard deviations of each bin's count, per volume.
    tolerance = 5 * np.sqrt(bulk * volumes) / volumes
    assert np.all(np.abs(rho - bulk) <= tolerance)


def test_number_density_divides_by_frames_and_volume():
    """``rho`` is the count per frame per bin volume (A1)."""
    rho = number_density([6, 12], [2.0, 3.0], n_frames=3)
    assert_allclose(rho, [1.0, 4.0 / 3.0], rtol=1e-15)  # round-off


def test_number_density_treats_rows_as_frames():
    """A ``(n_frames, n_bins)`` input gives each row's own density (P)."""
    counts = np.array([[2, 4], [6, 0]])
    rho = number_density(counts, [1.0, 2.0])
    assert_allclose(rho, [[2.0, 2.0], [6.0, 0.0]], rtol=1e-15)
    assert_allclose(rho[1], number_density(counts[1], [1.0, 2.0]),
                    rtol=1e-15)


def test_number_density_raises_when_a_bead_is_missing():
    """Counts not summing to ``n_frames * n_particles`` raise ValueError."""
    with pytest.raises(ValueError, match='a bead was lost'):
        number_density([2, 3], [1.0, 1.0], n_frames=1, n_particles=6)


@pytest.mark.parametrize(
    'counts, volumes, kwargs, message',
    [
        pytest.param([[[1]]], [1.0], {}, 'one- or two-dimensional',
                     id='3-d-counts'),
        pytest.param([1, 2], [1.0], {}, 'volumes must have shape',
                     id='shape-mismatch'),
        pytest.param([1, -2], [1.0, 1.0], {}, 'non-negative',
                     id='negative-count'),
        pytest.param([1, 2], [1.0, 0.0], {}, 'volumes must be positive',
                     id='zero-volume'),
        pytest.param([1, 2], [1.0, 1.0], {'n_frames': 0},
                     'n_frames must be positive', id='zero-frames'),
    ],
)
def test_number_density_rejects_bad_input(counts, volumes, kwargs, message):
    """Mismatched shapes, negative counts, or bad volumes raise ValueError."""
    with pytest.raises(ValueError, match=message):
        number_density(counts, volumes, **kwargs)


def test_volume_fraction_treats_rows_as_frames():
    """A ``(n_frames, n_bins)`` density gives each row its own ``phi`` (P)."""
    share = share_matrix([0.0, 1.0, 2.0, 3.0, 4.0], 1.0, 'cylindrical', 'z')
    rho = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 2.0, 0.0]])
    phi = volume_fraction(rho, share)
    assert phi.shape == (2, 4)
    assert_allclose(phi[1], volume_fraction(rho[1], share),
                    atol=FLOAT64_ATOL)


@pytest.mark.parametrize(
    'rho, message',
    [
        pytest.param(np.ones((2, 2, 4)), 'one- or two-dimensional',
                     id='3-d-rho'),
        pytest.param(np.ones(3), 'has 3 bins but share has 4',
                     id='bin-mismatch'),
    ],
)
def test_volume_fraction_rejects_bad_input(rho, message):
    """A density of the wrong rank or length raises ValueError."""
    share = share_matrix([0.0, 1.0, 2.0, 3.0, 4.0], 1.0, 'cylindrical', 'z')
    with pytest.raises(ValueError, match=message):
        volume_fraction(rho, share)


# ---------------------------------------------------------------------------
# normalize_profile (T12)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'mode, kwargs, expected',
    [
        pytest.param('none', {}, [1.0, 3.0], id='none'),
        pytest.param('sum', {}, [0.25, 0.75], id='sum'),
        pytest.param('integral', {'volumes': [2.0, 2.0]}, [0.125, 0.375],
                     id='integral'),
    ],
)
def test_normalize_profile_modes(mode, kwargs, expected):
    """``'sum'`` sums to 1, ``'integral'`` integrates to 1 against the volumes,
    and ``'none'`` returns the input (A1)."""
    out = normalize_profile([1.0, 3.0], mode, **kwargs)
    assert_allclose(out, expected, rtol=1e-15)  # round-off


def test_normalize_profile_normalizes_each_row():
    """A two-dimensional profile is normalized row by row (P)."""
    out = normalize_profile([[1.0, 3.0], [2.0, 2.0]], 'sum')
    assert_allclose(out, [[0.25, 0.75], [0.5, 0.5]], rtol=1e-15)


@pytest.mark.parametrize(
    'profile, mode, kwargs, message',
    [
        pytest.param([1.0], 'max', {}, 'mode must be one of', id='bad-mode'),
        pytest.param([1.0], 'integral', {}, 'volumes are needed',
                     id='missing-volumes'),
        pytest.param([1.0, 2.0], 'integral', {'volumes': [1.0]},
                     'volumes must have shape', id='volume-shape'),
        pytest.param([0.0, 0.0], 'sum', {}, 'zero normalizer',
                     id='zero-sum'),
        pytest.param([[[1.0]]], 'sum', {}, 'one- or two-dimensional',
                     id='3-d-profile'),
    ],
)
def test_normalize_profile_rejects_bad_input(profile, mode, kwargs, message):
    """A bad mode, missing volumes, or a zero normalizer raises ValueError."""
    with pytest.raises(ValueError, match=message):
        normalize_profile(profile, mode, **kwargs)


# ---------------------------------------------------------------------------
# density_profiles (D4, D5, T12)
# ---------------------------------------------------------------------------


def test_density_profiles_uses_the_effective_slab_diameter():
    """Along ``z`` the slab area is ``pi (D - a)**2 / 4``, so three beads of
    diameter 1 per unit slab in ``D = 3`` give ``rho = 3 / pi`` and
    ``phi = 1/2`` (A1)."""
    out = density_profiles([3, 3, 3, 3], [0.0, 1.0, 2.0, 3.0, 4.0], 1.0,
                           'cylindrical', 'z', dcyl=3.0, n_particles=12)
    assert isinstance(out, Profiles)
    assert_allclose(out.volumes, np.pi, rtol=1e-15)
    assert_allclose(out.rho, 3 / np.pi, rtol=1e-15)
    assert_allclose(out.phi, 0.5, rtol=1e-14)


def test_density_profiles_divides_by_frames():
    """Counts summed over ``n_frames`` give the per-frame density (A1)."""
    out = density_profiles([6, 6, 6, 6], [0.0, 1.0, 2.0, 3.0, 4.0], 1.0,
                           'cylindrical', 'z', dcyl=3.0, n_frames=2)
    assert_allclose(out.rho, 3 / np.pi, rtol=1e-15)


def test_density_profiles_integral_normalization_integrates_to_one():
    """``normalize='integral'`` makes both profiles integrate to 1 (A1)."""
    out = density_profiles([1, 4, 2, 5], [0.0, 1.0, 2.0, 3.0, 4.0], 1.0,
                           'cylindrical', 'z', dcyl=3.0, normalize='integral')
    assert (out.rho * out.volumes).sum() == pytest.approx(1.0, rel=1e-14)
    assert (out.phi * out.volumes).sum() == pytest.approx(1.0, rel=1e-14)


def test_density_profiles_rejects_a_bead_wider_than_the_cylinder():
    """``dcyl`` no larger than the bead diameter raises ValueError."""
    with pytest.raises(ValueError, match='must exceed the bead diameter'):
        density_profiles([1, 1], [0.0, 1.0, 2.0], 1.0, 'cylindrical', 'z',
                         dcyl=1.0)


# ---------------------------------------------------------------------------
# Input validation shared by the kernels (T14)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'edges, geometry, direction, kwargs, message',
    [
        pytest.param([1.0], 'cubic', 'r', {}, 'at least two values',
                     id='one-edge'),
        pytest.param([[0.0, 1.0]], 'cubic', 'r', {}, 'one-dimensional',
                     id='2-d-edges'),
        pytest.param([0.0, np.inf], 'cubic', 'r', {}, 'finite',
                     id='infinite-edge'),
        pytest.param([0.0, 1.0, 1.0], 'cubic', 'r', {}, 'strictly increasing',
                     id='repeated-edge'),
        pytest.param([-1.0, 1.0], 'cubic', 'r', {}, 'cannot be negative',
                     id='negative-radius'),
        pytest.param([0.0, 1.0], 'cubic', 'z', {}, 'is not supported',
                     id='cubic-z'),
        pytest.param([0.0, 1.0], 'slit', 'r', {}, 'is not supported',
                     id='slit'),
        pytest.param([0.0, 1.0], 'cylindrical', 'r', {}, 'lcyl is needed',
                     id='missing-lcyl'),
        pytest.param([0.0, 1.0], 'cylindrical', 'z', {'dcyl': -2.0},
                     'dcyl must be positive', id='negative-dcyl'),
    ],
)
def test_bin_volumes_rejects_bad_input(edges, geometry, direction, kwargs,
                                       message):
    """Bad edges, an unsupported pair, or no box length raise ValueError."""
    with pytest.raises(ValueError, match=message):
        bin_volumes(edges, geometry, direction, **kwargs)


@pytest.mark.parametrize(
    'bead_diameter, direction, kwargs, message',
    [
        pytest.param(0.0, 'r', {}, 'bead_diameter must be positive',
                     id='zero-bead'),
        pytest.param(5.0, 'z', {}, 'does not fit the periodic length',
                     id='bead-wider-than-period'),
        pytest.param(1.0, 'r', {'shares': 'mean'}, 'shares must be one of',
                     id='bad-shares'),
        pytest.param(1.0, 'r', {'intersection': 'plane'},
                     'intersection must be one of', id='bad-intersection'),
    ],
)
def test_share_matrix_rejects_bad_input(bead_diameter, direction, kwargs,
                                        message):
    """A bad bead size or an unknown mode raises ValueError."""
    with pytest.raises(ValueError, match=message):
        share_matrix([0.0, 1.0, 2.0, 3.0, 4.0], bead_diameter, 'cylindrical',
                     direction, **kwargs)


def test_share_matrix_returns_a_sparse_band():
    """The share matrix is sparse, with about ``a / h + 2`` entries per row."""
    edges = np.arange(-50.0, 50.01, 0.1)
    share = share_matrix(edges, 1.0, 'cylindrical', 'z')
    assert share.shape == (1000, 1000)
    assert share.nnz <= 1000 * 14
    assert_array_equal(np.diff(share.indptr) >= 11, True)


# ---------------------------------------------------------------------------
# Tier C: the research code's output on one real whole (written blind)
# ---------------------------------------------------------------------------


@pytest.mark.private_data
@pytest.mark.parametrize(
    ('direction', 'research_scale'),
    [
        pytest.param('r', lambda meta: 1.0, id='r'),
        # R3: the research code's z slabs have lcyl, not dcyl, as diameter.
        pytest.param('z', lambda meta: (meta['dcyl'] / meta['lcyl']) ** 2,
                     id='z'),
    ],
)
def test_center_sphere_profiles_match_the_research_code(
    private_data_dir, direction, research_scale
):
    """The research code's ``rho`` and ``phi`` of one SciOps whole equal the
    center-sphere ones times its known axial slab factor (R1, R3).

    Expected layout under ``POLYPHYS_DATA_DIR/density/<direction>/``,
    prepared by the maintainer from one whole: ``edges.npy`` and
    ``counts.npy`` (the whole's bin edges and summed histogram),
    ``rho_sciops.npy`` and ``phi_sciops.npy`` (its ``distributions_generator``
    output), and ``meta.json`` with ``geometry``, ``direction``,
    ``bead_diameter``, ``lcyl``, and ``dcyl``. The research code divided by
    the quadrature of a polynomial, never by the number of frames, and used
    the sphere-sphere kernel, which is how the kernels are called here.
    """
    folder = private_data_dir / 'density' / direction
    meta = json.loads((folder / 'meta.json').read_text())
    edges = np.load(folder / 'edges.npy')
    counts = np.load(folder / 'counts.npy')
    volumes = bin_volumes(
        edges, meta['geometry'], meta['direction'],
        lcyl=meta['lcyl'], dcyl=meta['dcyl'],
    )
    share = share_matrix(edges, meta['bead_diameter'], meta['geometry'],
                         meta['direction'], shares='center',
                         intersection='sphere')
    scale = research_scale(meta)
    rho = number_density(counts, volumes)
    phi = volume_fraction(rho, share)
    # scipy's quadrature of a polynomial agrees with the closed form to
    # about 1e-10 relative.
    assert_allclose(scale * rho, np.load(folder / 'rho_sciops.npy'),
                    rtol=1e-8)
    assert_allclose(scale * phi, np.load(folder / 'phi_sciops.npy'),
                    rtol=1e-8)
