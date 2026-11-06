"""
Tests for :mod:`polyphys.analyze.correlations`.

Every expected value is independent of the code under test. It comes from
one of these sources, named in the test that uses it:

A1. AR(1), ``x_t = phi x_{t-1} + e_t`` with ``|phi| < 1``, independent
    innovations, and a stationary start: multiplying by ``x_{t-k}`` and
    taking expectations gives ``gamma_k = phi gamma_{k-1}``, so the
    autocorrelation is ``rho_k = phi**k``.

A2. The unit bond vectors of a regular planar polygon of ``N`` sides turn by
    ``2 pi / N`` from one bond to the next, so two bonds ``k`` apart have
    cosine ``cos(2 pi k / N)``.

A3. Every bond of a straight rod is parallel to every other, so the bond
    cosine is 1 at every lag.

A4. A pure exponential ``A exp(-x/tau) + B`` sampled without noise is in the
    model family, so least squares recovers its parameters exactly (to
    round-off).

R1. Brute force: the defining sums, one lag at a time, in
    :func:`_acf_by_definition` (the ``1/N`` autocovariance of the
    NIST/SEMATECH e-Handbook, section 1.3.3.1) and
    :func:`_bond_cosines_by_definition` (a double loop over bond pairs).

R2. Bartlett's variance as the e-Handbook states it and ``statsmodels``
    computes it, ``(1 + 2 sum_{j=1}^{k-1} r_j^2) / N`` for ``k >= 2``,
    ``1/N`` at ``k = 1``, and 0 at ``k = 0``, written out in
    :func:`_bartlett_half_width`.

P.  A property: shape, normalization, invariance under affine rescaling of
    a series, rigid motion of a chain, relabeling of a ring, or
    column-wise independence.

Random inputs come from ``numpy.random.default_rng(seed)`` with the seed
shown, so every example is reproducible. Float comparisons use
``FLOAT64_ATOL`` where two computations differ only in the order of float64
operations; the one statistical tolerance is derived where it is used.
"""

import logging

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.optimize import OptimizeWarning
from scipy.stats import norm

from polyphys.analyze.correlations import (
    DECAY_MODELS,
    AcfResult,
    DecayFit,
    acf,
    bartlett_half_width,
    bond_cosine_correlation,
    bond_vectors,
    fit_decay,
    initial_guess,
)

# Absolute tolerance for results that differ only by float64 round-off over
# an FFT of a few thousand points: far above the 2.2e-16 epsilon, far below
# any value the tests compare.
FLOAT64_ATOL = 1e-10


# ---------------------------------------------------------------------------
# Reference computations
# ---------------------------------------------------------------------------


def _acf_by_definition(x: np.ndarray, nlags: int) -> np.ndarray:
    """Return the ``1/N`` sample autocorrelation by its defining sum (R1)."""
    n = len(x)
    centered = x - x.mean()
    autocov = np.array(
        [np.sum(centered[: n - k] * centered[k:]) / n
         for k in range(nlags + 1)]
    )
    return autocov / autocov[0]


def _bartlett_half_width(
    values: np.ndarray, n: int, alpha: float
) -> np.ndarray:
    """Return Bartlett's half-widths for the given autocorrelations (R2)."""
    z = norm.ppf(1.0 - alpha / 2.0)
    out = np.empty(len(values))
    for k in range(len(values)):
        if k == 0:
            var = 0.0
        else:
            var = (1.0 + 2.0 * np.sum(values[1:k] ** 2)) / n
        out[k] = z * np.sqrt(var)
    return out


def _bond_cosines_by_definition(
    units: np.ndarray, topology: str
) -> np.ndarray:
    """Return the mean bond cosine at every lag by a double loop (R1)."""
    n = len(units)
    out = np.empty(n)
    for k in range(n):
        if topology == 'ring':
            pairs = [units[i] @ units[(i + k) % n] for i in range(n)]
        else:
            pairs = [units[i] @ units[i + k] for i in range(n - k)]
        out[k] = np.mean(pairs)
    return out


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def make_ar1():
    """Return a factory for a stationary AR(1) series with unit innovations."""

    def _make(phi: float, n_frames: int, seed: int) -> np.ndarray:
        rng = np.random.default_rng(seed)
        noise = rng.standard_normal(n_frames)
        x = np.empty(n_frames, dtype=np.float64)
        x[0] = noise[0] / np.sqrt(1.0 - phi ** 2)  # stationary start
        for t in range(1, n_frames):
            x[t] = phi * x[t - 1] + noise[t]
        return x

    return _make


@pytest.fixture
def make_noise():
    """Return a factory for i.i.d. standard normal series, ``(n,)`` or
    ``(n, m)``."""

    def _make(shape, seed: int) -> np.ndarray:
        return np.random.default_rng(seed).standard_normal(shape)

    return _make


@pytest.fixture
def make_polygon():
    """Return a factory for the vertices of a regular planar polygon."""

    def _make(n_sides: int, radius: float = 1.0) -> np.ndarray:
        angles = 2.0 * np.pi * np.arange(n_sides) / n_sides
        return radius * np.column_stack(
            [np.cos(angles), np.sin(angles), np.zeros(n_sides)]
        )

    return _make


@pytest.fixture
def make_random_walk():
    """Return a factory for a random walk of unit steps in three dimensions,
    with any leading frames axis."""

    def _make(n_atoms: int, seed: int, n_frames: int | None = None):
        rng = np.random.default_rng(seed)
        shape = (n_atoms, 3) if n_frames is None else (n_frames, n_atoms, 3)
        steps = rng.standard_normal(shape)
        steps /= np.linalg.norm(steps, axis=-1, keepdims=True)
        return np.cumsum(steps, axis=-2)

    return _make


# ---------------------------------------------------------------------------
# acf
# ---------------------------------------------------------------------------


def test_acf_lag_zero_is_one_with_zero_band(make_noise):
    """The autocorrelation at lag 0 is exactly 1 with a zero half-width,
    and the result has ``nlags + 1`` rows (P)."""
    x = make_noise(100, seed=1)
    result = acf(x, nlags=5)
    assert isinstance(result, AcfResult)
    assert result.values.shape == (6,)
    assert result.half_width is not None
    assert result.half_width.shape == (6,)
    assert result.values[0] == 1.0
    assert result.half_width[0] == 0.0


def test_acf_matches_defining_sum(make_ar1):
    """The FFT autocorrelation equals the ``1/N`` defining sum at every lag
    on an odd, non-power-of-two length (R1)."""
    x = make_ar1(0.5, n_frames=257, seed=2)
    result = acf(x, nlags=40, alpha=None)
    assert_allclose(
        result.values, _acf_by_definition(x, 40), atol=FLOAT64_ATOL
    )


def test_acf_of_ar1_is_phi_to_k_within_bartlett_band(make_ar1):
    """The AR(1) autocorrelation ``phi**k`` lies inside Bartlett's band at
    lags 1 to 10 (A1)."""
    phi, n_frames, n_lags = 0.6, 2 ** 15, 10
    # Bonferroni over the 10 lags: a 1 % family-wise level, so that the
    # check holds for almost every seed rather than one chosen seed.
    alpha = 0.01 / n_lags
    x = make_ar1(phi, n_frames, seed=3)
    result = acf(x, nlags=n_lags, alpha=alpha)
    assert result.half_width is not None
    lags = np.arange(1, n_lags + 1)
    deviation = np.abs(result.values[1:] - phi ** lags)
    assert np.all(deviation <= result.half_width[1:])


def test_acf_of_white_noise_is_within_root_n_band(make_noise):
    """White noise has autocorrelations inside ``z / sqrt(N)`` at every lag
    checked, the band the e-Handbook gives for randomness (P, R2)."""
    n_frames, n_lags = 2 ** 14, 10
    alpha = 0.01 / n_lags  # Bonferroni over the lags, as above
    x = make_noise(n_frames, seed=4)
    result = acf(x, nlags=n_lags, alpha=alpha)
    band = norm.ppf(1.0 - alpha / 2.0) / np.sqrt(n_frames)
    assert np.all(np.abs(result.values[1:]) <= band)


def test_acf_band_is_bartlett_formula():
    """The half-width at lag ``k`` is ``z sqrt((1 + 2 sum_{j<k} r_j^2)/N)``
    with ``1/N`` at lag 1 and 0 at lag 0, on a six-value series (R2)."""
    x = np.array([2.0, 4.0, 1.0, 5.0, 3.0, 6.0])
    result = acf(x, nlags=4, alpha=0.1)
    assert result.half_width is not None
    expected = _bartlett_half_width(result.values, n=6, alpha=0.1)
    assert_allclose(result.half_width, expected, atol=FLOAT64_ATOL)


def test_bartlett_half_width_rebuilds_the_band_of_acf(make_ar1):
    """The band rebuilt from a stored autocorrelation, the frame count, and
    the level equals the one ``acf`` computed with them."""
    x = make_ar1(0.5, 256, 11)
    result = acf(x, nlags=10, alpha=0.1)

    rebuilt = bartlett_half_width(result.values, 256, alpha=0.1)

    assert_allclose(rebuilt, result.half_width, atol=FLOAT64_ATOL)
    assert rebuilt[0] == 0.0


def test_bartlett_half_width_treats_columns_independently(make_noise):
    """Each column's band is the band of that column alone (R2)."""
    x = make_noise((128, 2), seed=5)
    values = acf(x, nlags=6, alpha=None).values

    band = bartlett_half_width(values, 128, alpha=0.05)

    for col in range(2):
        assert_allclose(band[:, col],
                        _bartlett_half_width(values[:, col], 128, 0.05),
                        atol=FLOAT64_ATOL)


@pytest.mark.parametrize("values, n_frames, alpha, match", [
    pytest.param(np.ones((2, 2, 2)), 8, 0.05, "one- or two-dimensional",
                 id="three-dimensional"),
    pytest.param([1.0, 0.5], 1, 0.05, "at least 2", id="one-frame"),
    pytest.param([1.0, 0.5], 8, 1.0, r"alpha must be in \(0, 1\)",
                 id="alpha-one"),
])
def test_bartlett_half_width_rejects_bad_input(values, n_frames, alpha,
                                               match):
    """Bad shapes, frame counts, and levels raise ``ValueError``."""
    with pytest.raises(ValueError, match=match):
        bartlett_half_width(values, n_frames, alpha)


def test_acf_treats_columns_as_independent_series(make_noise):
    """A ``(n, m)`` input gives, column by column, the result of ``m``
    one-dimensional calls (P)."""
    x = make_noise((200, 3), seed=5)
    result = acf(x, nlags=7)
    assert result.half_width is not None
    assert result.values.shape == (8, 3)
    for j in range(3):
        single = acf(x[:, j], nlags=7)
        assert single.half_width is not None
        assert_allclose(result.values[:, j], single.values, atol=FLOAT64_ATOL)
        assert_allclose(
            result.half_width[:, j], single.half_width, atol=FLOAT64_ATOL
        )


def test_acf_is_invariant_under_positive_affine_rescaling(make_noise):
    """Scaling a series by a positive factor and shifting it leaves its
    autocorrelation unchanged (P)."""
    x = make_noise(300, seed=6)
    reference = acf(x, nlags=10, alpha=None).values
    rescaled = acf(3.5 * x - 20.0, nlags=10, alpha=None).values
    assert_allclose(rescaled, reference, atol=FLOAT64_ATOL)


def test_acf_skips_band_when_alpha_is_none(make_noise):
    """``alpha=None`` returns no half-width and records no level (P)."""
    result = acf(make_noise(50, seed=7), nlags=3, alpha=None)
    assert result.half_width is None
    assert result.alpha is None


@pytest.mark.parametrize(
    'x, nlags, alpha, match',
    [
        pytest.param(np.ones(10), 2, 0.05, 'constant', id='constant-series'),
        pytest.param(
            np.arange(10.0), 10, 0.05, 'nlags must be', id='nlags-too-large'
        ),
        pytest.param(
            np.arange(10.0), -1, 0.05, 'nlags must be', id='nlags-negative'
        ),
        pytest.param(
            np.array([1.0, np.nan, 2.0]), 1, 0.05, 'finite', id='nan-value'
        ),
        pytest.param(
            np.arange(10.0), 2, 1.0, 'alpha must be', id='alpha-one'
        ),
        pytest.param(
            np.arange(10.0), 2, 0.0, 'alpha must be', id='alpha-zero'
        ),
        pytest.param(np.ones((2, 2, 2)), 1, 0.05, 'dimensional', id='3-d'),
        pytest.param(np.array([1.0]), 0, 0.05, 'two frames', id='one-frame'),
    ],
)
def test_acf_rejects_bad_input(x, nlags, alpha, match):
    """A constant series, a lag outside ``[0, N-1]``, a non-finite value, a
    level outside ``(0, 1)``, or a bad shape raises ``ValueError``."""
    with pytest.raises(ValueError, match=match):
        acf(x, nlags=nlags, alpha=alpha)


# ---------------------------------------------------------------------------
# fit_decay and initial_guess
# ---------------------------------------------------------------------------


def test_fit_decay_recovers_pure_exponential():
    """A noise-free ``exp(-x/tau)`` is recovered to 1e-10 from a start 50 %
    off (A4)."""
    x = np.arange(0.0, 40.0)
    fit = fit_decay(x, np.exp(-x / 3.7), 'exp', p0=[5.55])
    assert fit.converged
    assert fit.params[0] == pytest.approx(3.7, rel=1e-10)


def test_fit_decay_recovers_stretched_amp_offset():
    """A noise-free ``A exp(-(x/tau)^beta) + B`` is recovered to 1e-8 in
    every parameter (A4)."""
    x = np.arange(0.0, 80.0)
    truth = np.array([5.0, 0.7, 2.0, 0.1])
    y = DECAY_MODELS['stretched_amp_offset'].function(x, *truth)
    fit = fit_decay(x, y, 'stretched_amp_offset', p0=[4.0, 1.0, 1.5, 0.0])
    assert fit.converged
    assert fit.names == ('tau', 'beta', 'amp', 'offset')
    assert_allclose(fit.params, truth, rtol=1e-8)
    assert_allclose(fit.evaluate(x), y, rtol=1e-8)


@pytest.mark.parametrize('model', list(DECAY_MODELS))
def test_fit_decay_returns_one_value_per_parameter(model):
    """Every model's fit returns parameters, errors, a starting point, and a
    covariance sized by its parameter count, in its parameter order (P)."""
    x = np.arange(0.0, 30.0)
    y = np.exp(-x / 6.0)
    fit = fit_decay(x, y, model)
    n_params = len(DECAY_MODELS[model].names)
    assert isinstance(fit, DecayFit)
    assert fit.model == model
    assert fit.names == DECAY_MODELS[model].names
    assert fit.params.shape == (n_params,)
    assert fit.stderr.shape == (n_params,)
    assert fit.p0.shape == (n_params,)
    assert fit.pcov.shape == (n_params, n_params)


def test_fit_decay_reports_failure_without_raising(caplog):
    """A fit the optimizer gives up on returns ``converged=False`` with NaN
    parameters of the model's length and a logged warning (P)."""
    x = np.arange(0.0, 80.0)
    y = 2.0 * np.exp(-((x / 5.0) ** 0.7)) + 0.1
    with caplog.at_level(logging.WARNING, logger='polyphys'):
        fit = fit_decay(x, y, 'stretched_amp_offset', maxfev=1)
    assert fit.converged is False
    assert fit.params.shape == (4,)
    assert np.all(np.isnan(fit.params))
    assert np.all(np.isnan(fit.stderr))
    assert np.all(np.isnan(fit.pcov))
    assert np.all(np.isnan(fit.evaluate(x)))
    assert 'did not converge' in caplog.text


def test_fit_decay_stderr_is_root_of_covariance_diagonal():
    """The standard error of each parameter is the square root of the
    covariance diagonal, not the variance itself (P)."""
    x = np.arange(0.0, 60.0)
    rng = np.random.default_rng(8)
    y = np.exp(-x / 7.0) + 0.02 * rng.standard_normal(x.size)
    fit = fit_decay(x, y, 'stretched_amp', p0=[7.0, 1.0, 1.0])
    assert fit.converged
    assert_allclose(fit.stderr ** 2, np.diag(fit.pcov), rtol=1e-12)
    assert np.all(fit.stderr > 0.0)


def test_fit_decay_records_the_start_it_used():
    """With ``p0=None`` the fit starts from :func:`initial_guess` and
    records that start; an explicit ``p0`` is recorded as given (P)."""
    x = np.arange(0.0, 40.0)
    y = 2.0 * np.exp(-x / 4.0) + 0.5
    default = fit_decay(x, y, 'stretched_amp_offset')
    explicit = fit_decay(x, y, 'stretched_amp_offset', p0=[3.0, 1.0, 1.0, 0.0])
    assert_array_equal(
        default.p0, initial_guess(x, y, 'stretched_amp_offset')
    )
    assert_array_equal(explicit.p0, [3.0, 1.0, 1.0, 0.0])


def test_initial_guess_is_within_a_factor_of_two_and_fit_converges_from_it():
    """On a noise-free curve every guess lies within a factor of 2 of the
    truth, and the fit from it recovers the truth (A4)."""
    x = np.arange(0.0, 100.0)
    truth = np.array([5.0, 1.0, 2.0, 0.3])
    y = DECAY_MODELS['stretched_amp_offset'].function(x, *truth)
    guess = initial_guess(x, y, 'stretched_amp_offset')
    ratio = guess / truth
    assert np.all((ratio >= 0.5) & (ratio <= 2.0))
    fit = fit_decay(x, y, 'stretched_amp_offset', p0=guess)
    assert fit.converged
    assert_allclose(fit.params, truth, rtol=1e-8)


@pytest.mark.parametrize(
    'model, expected',
    [
        pytest.param('exp', [4.0], id='exp'),
        pytest.param('stretched', [4.0, 1.0], id='stretched'),
        pytest.param('stretched_amp', [4.0, 1.0, 1.0], id='stretched-amp'),
    ],
)
def test_initial_guess_reads_tau_at_one_over_e(model, expected):
    """The decay-time guess is the first ``x`` at which a unit exponential
    is at ``1/e``, which is ``tau`` itself when it is sampled, with ``beta``
    1 and amplitude 1 for a model without an offset (A4)."""
    x = np.arange(0.0, 20.0)
    guess = initial_guess(x, np.exp(-x / 4.0), model)
    assert_allclose(guess, expected, atol=FLOAT64_ATOL)


def test_initial_guess_uses_unit_amplitude_for_flat_curve():
    """A curve whose first value equals its tail gives an amplitude guess
    of 1, not 0, so the fit can start (P)."""
    x = np.arange(0.0, 20.0)
    guess = initial_guess(x, np.full(20, 0.3), 'stretched_amp_offset')
    assert guess[2] == 1.0


def test_initial_guess_rejects_empty_curve():
    """An empty curve raises ``ValueError`` instead of indexing nothing."""
    with pytest.raises(ValueError, match='empty'):
        initial_guess([], [], 'exp')


def test_fit_decay_marks_unestimable_covariance_as_nan(caplog):
    """A fit with as many points as parameters converges but gets NaN
    standard errors and a logged warning, since SciPy cannot estimate the
    covariance (P)."""
    x = np.array([2.0])
    y = np.exp(-x / 3.0)
    with caplog.at_level(logging.WARNING, logger='polyphys'):
        with pytest.warns(OptimizeWarning, match='Covariance'):
            fit = fit_decay(x, y, 'exp', p0=[2.5])
    assert fit.converged
    assert fit.params[0] == pytest.approx(3.0, rel=1e-8)
    assert np.all(np.isnan(fit.stderr))
    assert 'covariance' in caplog.text


def test_initial_guess_never_returns_zero_tau():
    """A curve already at ``1/e`` at ``x = 0`` gets the first positive ``x``
    as its decay time, never 0 (P)."""
    x = np.array([0.0, 0.5, 1.0, 1.5])
    guess = initial_guess(x, np.full(4, 0.1), 'exp')
    assert guess[0] == 0.5


@pytest.mark.parametrize(
    'kwargs, match',
    [
        pytest.param(
            {'model': 'mono_exp'}, 'model must be one of', id='unknown-model'
        ),
        pytest.param(
            {'model': 'exp', 'p0': [1.0, 2.0]}, 'p0 must have', id='p0-length'
        ),
        pytest.param(
            {'model': 'exp', 'y': np.ones(5)}, 'same length', id='length'
        ),
        pytest.param(
            {'model': 'exp', 'y': np.full(10, np.nan)}, 'finite', id='nan'
        ),
        pytest.param(
            {'model': 'exp', 'x': np.ones((2, 5)), 'y': np.ones((2, 5))},
            'one-dimensional',
            id='2-d',
        ),
        pytest.param(
            {'model': 'stretched_amp_offset', 'x': np.arange(3.0),
             'y': np.ones(3)},
            'only 3 point',
            id='too-few-points',
        ),
    ],
)
def test_fit_decay_rejects_bad_input(kwargs, match):
    """An unknown model, a wrong-length ``p0``, mismatched or non-finite
    curves, a 2-D curve, or fewer points than parameters raises
    ``ValueError``."""
    args = {'x': np.arange(10.0), 'y': np.exp(-np.arange(10.0))}
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        fit_decay(**args)


# ---------------------------------------------------------------------------
# bond_vectors and bond_cosine_correlation
# ---------------------------------------------------------------------------


def test_bond_cosine_of_straight_rod_is_one():
    """A straight rod has bond cosine 1 at every lag (A3)."""
    rod = np.outer(np.arange(10.0), [0.3, -0.4, 1.2])
    assert_allclose(
        bond_cosine_correlation(rod, 'linear'), np.ones(9), atol=FLOAT64_ATOL
    )


@pytest.mark.parametrize('n_sides', [pytest.param(12, id='12-gon'),
                                     pytest.param(7, id='7-gon')])
def test_bond_cosine_of_regular_ring_is_cosine_of_turning_angle(
    make_polygon, n_sides
):
    """A regular ring of ``N`` bonds has ``C(k) = cos(2 pi k / N)`` (A2)."""
    lags = np.arange(n_sides)
    result = bond_cosine_correlation(make_polygon(n_sides), 'ring')
    assert_allclose(
        result, np.cos(2.0 * np.pi * lags / n_sides), atol=FLOAT64_ATOL
    )


@pytest.mark.parametrize('topology', ['linear', 'ring'])
def test_bond_cosine_matches_double_loop(make_random_walk, topology):
    """The FFT bond cosine equals the mean over pairs by a double loop at
    every lag, for both topologies (R1)."""
    positions = make_random_walk(n_atoms=33, seed=9)
    bonds = bond_vectors(positions, topology)
    units = bonds / np.linalg.norm(bonds, axis=-1, keepdims=True)
    assert_allclose(
        bond_cosine_correlation(positions, topology),
        _bond_cosines_by_definition(units, topology),
        atol=FLOAT64_ATOL,
    )


def test_bond_cosine_is_invariant_under_rigid_motion(make_random_walk):
    """Rotating and translating a chain leaves its bond cosines unchanged
    (P)."""
    positions = make_random_walk(n_atoms=20, seed=10)
    rng = np.random.default_rng(11)
    rotation, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    moved = positions @ rotation.T + np.array([5.0, -2.0, 7.5])
    assert_allclose(
        bond_cosine_correlation(moved, 'linear'),
        bond_cosine_correlation(positions, 'linear'),
        atol=FLOAT64_ATOL,
    )


def test_bond_cosine_keeps_leading_frames_axis(make_random_walk):
    """A ``(n_frames, n_atoms, 3)`` input gives one row per frame, equal to
    the single-frame results (P)."""
    frames = make_random_walk(n_atoms=15, seed=12, n_frames=4)
    result = bond_cosine_correlation(frames, 'ring')
    assert result.shape == (4, 15)
    for f in range(4):
        assert_allclose(
            result[f], bond_cosine_correlation(frames[f], 'ring'),
            atol=FLOAT64_ATOL,
        )


def test_bond_cosine_of_ring_is_symmetric_and_relabeling_invariant(
    make_random_walk,
):
    """On a ring ``C(k) = C(N - k)``, and starting the labels at another
    monomer gives the same curve (P)."""
    positions = make_random_walk(n_atoms=21, seed=13)
    result = bond_cosine_correlation(positions, 'ring')
    assert_allclose(result[1:], result[1:][::-1], atol=FLOAT64_ATOL)
    relabeled = bond_cosine_correlation(np.roll(positions, 5, axis=0), 'ring')
    assert_allclose(relabeled, result, atol=FLOAT64_ATOL)


def test_bond_cosine_runs_on_linear_chain():
    """A linear chain of four beads gives three lags without error, which
    the research ``bond_info`` could not (P)."""
    positions = np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1, 0], [2.0, 1, 0]])
    result = bond_cosine_correlation(positions, 'linear')
    assert result.shape == (3,)
    assert_allclose(result, [1.0, 0.0, 1.0], atol=FLOAT64_ATOL)


@pytest.mark.parametrize(
    'topology, expected',
    [
        pytest.param(
            'linear', [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]], id='linear'
        ),
        pytest.param(
            'ring',
            [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]],
            id='ring',
        ),
    ],
)
def test_bond_vectors_differences_along_chain(topology, expected):
    """Bond ``i`` is ``r[i+1] - r[i]``, and a ring adds ``r[0] - r[-1]``
    as its last bond (hand count on a unit square)."""
    square = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    assert_array_equal(bond_vectors(square, topology), expected)


@pytest.mark.parametrize(
    'positions, topology, match',
    [
        pytest.param(
            np.array([[0.0, 0, 0], [0.0, 0, 0], [1.0, 0, 0]]),
            'linear',
            'positive length',
            id='coincident-beads',
        ),
        pytest.param(np.eye(3), 'helix', "must be 'linear' or 'ring'",
                     id='bad-topology'),
        pytest.param(np.arange(6.0), 'linear', 'dimension', id='1-d'),
        pytest.param(np.eye(2), 'ring', 'at least 3 atoms', id='ring-of-two'),
        pytest.param(np.ones((1, 3)), 'linear', 'at least 2 atoms',
                     id='one-atom'),
    ],
)
def test_bond_cosine_rejects_bad_input(positions, topology, match):
    """Coincident bonded beads, an unknown topology, a 1-D input, or a chain
    too short for its topology raises ``ValueError``."""
    with pytest.raises(ValueError, match=match):
        bond_cosine_correlation(positions, topology)
