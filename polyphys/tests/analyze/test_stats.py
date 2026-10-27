"""
Tests for :mod:`polyphys.analyze.stats`.

Expected values come from three analytic results, derived here so that no
expected value depends on the code under test.

R1. For a stationary series with autocorrelation function ``rho(t)``, the
    statistical inefficiency is ``s = sum_{t=-inf}^{inf} rho(t)``, twice the
    integrated autocorrelation time (Wolff, Comput. Phys. Commun. 156, 143
    (2004), Eq. (23) and (25)). For independent frames ``s = 1`` at every
    block size.

R2. AR(1), ``x_t = phi x_{t-1} + e_t`` with unit-variance innovations and a
    stationary start, has ``rho(t) = phi**|t|`` and ``var(x) = 1/(1-phi^2)``.
    The variance of the mean of ``B`` consecutive values is
    ``(var(x)/B) [1 + 2 sum_{t=1}^{B-1} (1 - t/B) phi**t]``, so the
    inefficiency at block size ``B`` is the bracket, and its ``B -> inf``
    limit is ``(1+phi)/(1-phi)`` by the geometric series in R1. Finite-``n``
    corrections are of order ``s/n``, below 3e-4 for the sizes used here.

R3. A block-constant series, ``x_t = c_{t // L}`` with i.i.d. ``c_i`` and
    ``L`` dividing ``n``: for ``B <= L`` every block mean is one of the
    ``c_i`` and both ``ddof=1`` variances are sums over the same repeated
    values, so exactly ``s(B) = B (n-1) / (n-B)``. For ``B >= L`` the block
    means are means of ``B/L`` independent ``c_i``, so ``E[var_B] = var(c)
    L/B`` and ``E[var_0] = var(c) (n-L)/(n-1)``, giving ``L (n-1)/(n-L)`` as
    the ratio of expectations; the two forms agree at ``B = L``.

Statistical tolerances are five times the leading-order Gaussian error of a
variance estimated from ``m`` block means, ``sqrt(2/(m-1))`` relative (Lee et
al., Phys. Rev. E 83, 066706 (2011), Eq. (15)), checked only at levels with at
least 32 blocks, where that Gaussian approximation holds. With about one
hundred such comparisons in this file, the false-failure probability is below
1e-4, and every input is seeded, so a run is deterministic.
"""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.analyze.measurer import sem
from polyphys.analyze.stats import (
    BlockingResult,
    block_average,
    select_plateau,
    sem_blocked,
    statistical_inefficiency,
)

# Relative tolerance for results that differ only by float64 round-off over
# at most ~2 n additions: far above the 2.2e-16 epsilon, far below any
# statistical effect.
FLOAT64_RTOL = 1e-12
# Gaussian-error multiple and the smallest block count at which it is used;
# see the module docstring.
SIGMA = 5.0
MIN_BLOCKS_FOR_GAUSSIAN = 32


def ar1_inefficiency(phi: float, block_size: np.ndarray) -> np.ndarray:
    """Return the AR(1) inefficiency at each block size, from R2."""
    out = np.empty(len(block_size), dtype=np.float64)
    for i, size in enumerate(block_size):
        t = np.arange(1, size)
        out[i] = 1.0 + 2.0 * np.sum((1.0 - t / size) * phi ** t)
    return out


def block_constant_inefficiency(
    period: int, n_frames: int, block_size: np.ndarray
) -> np.ndarray:
    """Return the block-constant inefficiency at each block size, from R3."""
    effective = np.minimum(block_size, period)
    return effective * (n_frames - 1) / (n_frames - effective)


def relative_gaussian_error(n_blocks: np.ndarray) -> np.ndarray:
    """Return the relative error of a variance from ``n_blocks`` samples."""
    return np.sqrt(2.0 / (n_blocks - 1.0))


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
def make_block_constant():
    """Return a factory for a series constant over blocks of ``period``."""

    def _make(period: int, n_frames: int, seed: int) -> np.ndarray:
        rng = np.random.default_rng(seed)
        return np.repeat(rng.standard_normal(n_frames // period), period)

    return _make


@pytest.fixture
def make_table():
    """Return a factory for a hand-built :class:`BlockingResult`.

    Only ``sem`` (and optionally ``sem_err``) carry the test's values; the
    other fields are filled consistently with the given ``n_frames``.
    """

    def _make(
        sem_values, n_frames: int, sem_err_values=None
    ) -> BlockingResult:
        sem_arr = np.asarray(sem_values, dtype=np.float64)
        n_levels = len(sem_arr)
        block_size = 2 ** np.arange(n_levels, dtype=np.int64)
        n_blocks = n_frames // block_size
        var = sem_arr ** 2 * n_blocks
        if sem_err_values is None:
            sem_err = sem_arr / np.sqrt(2.0 * (n_blocks - 1))
        else:
            sem_err = np.asarray(sem_err_values, dtype=np.float64)
        inefficiency = block_size * var / var[0]
        return BlockingResult(
            n_frames=n_frames,
            block_size=block_size,
            n_blocks=n_blocks,
            mean=np.zeros(n_levels),
            var=var,
            sem=sem_arr,
            sem_err=sem_err,
            statistical_inefficiency=inefficiency,
            inefficiency_err=inefficiency * np.sqrt(2.0 / (n_blocks - 1)),
        )

    return _make


# --- the table ---------------------------------------------------------------


@pytest.mark.parametrize(
    "n_frames,n_levels",
    [
        pytest.param(64, 5, id="power-of-two-64"),
        pytest.param(100, 5, id="odd-levels-100"),
        pytest.param(127, 5, id="just-below-128"),
        pytest.param(128, 6, id="power-of-two-128"),
        pytest.param(1000, 8, id="1000"),
    ],
)
def test_block_average_level_layout(n_frames, n_levels):
    """Levels have block sizes 1, 2, 4, ... and stop at min_blocks blocks."""
    # n_levels is floor(log2(n_frames // 4)) + 1, worked by hand per case.
    x = np.random.default_rng(0).standard_normal(n_frames)

    table = block_average(x, min_blocks=4)

    assert table.n_levels == n_levels
    assert_array_equal(table.block_size, 2 ** np.arange(n_levels))
    assert_array_equal(table.n_blocks, n_frames // table.block_size)
    assert table.n_blocks[-1] >= 4
    assert table.n_blocks[-1] // 2 < 4


def test_block_average_level_zero_is_naive_estimate():
    """Level 0 reproduces the naive SEM and an inefficiency of exactly 1."""
    x = np.random.default_rng(1).standard_normal(1000)

    table = block_average(x)

    assert table.sem[0] == pytest.approx(sem(x), rel=FLOAT64_RTOL)
    # var / var is exactly 1 in IEEE arithmetic.
    assert table.statistical_inefficiency[0] == pytest.approx(1.0, abs=0.0)


def test_block_average_preserves_mean_without_dropped_frames():
    """With a power-of-two frame count the mean is the same at every level."""
    x = np.random.default_rng(2).standard_normal(2 ** 10)

    table = block_average(x)

    assert_allclose(
        table.mean, np.full(table.n_levels, x.mean()), rtol=FLOAT64_RTOL
    )


@pytest.mark.parametrize(
    "scale,offset",
    [
        pytest.param(-3.0, 7.0, id="negative-scale"),
        pytest.param(0.5, 7.0, id="shrink"),
    ],
)
def test_block_average_is_affine_equivariant(scale, offset):
    """Scaling a series scales its SEMs; the inefficiency does not change."""
    x = np.random.default_rng(3).standard_normal(2 ** 10)

    base = block_average(x)
    moved = block_average(scale * x + offset)

    assert_allclose(moved.sem, abs(scale) * base.sem, rtol=FLOAT64_RTOL)
    assert_allclose(
        moved.statistical_inefficiency,
        base.statistical_inefficiency,
        rtol=FLOAT64_RTOL,
    )


# --- known answers -----------------------------------------------------------


def test_block_average_iid_inefficiency_is_one_at_every_level():
    """Independent frames give an inefficiency of 1 at every level (R1)."""
    n_frames = 2 ** 16
    x = np.random.default_rng(4).standard_normal(n_frames)

    table = block_average(x)

    checked = table.n_blocks >= MIN_BLOCKS_FOR_GAUSSIAN
    sigma = relative_gaussian_error(table.n_blocks[checked])
    pulls = (table.statistical_inefficiency[checked] - 1.0) / sigma
    assert_allclose(pulls, np.zeros_like(pulls), rtol=0.0, atol=SIGMA)


@pytest.mark.parametrize(
    "phi,seed",
    [
        pytest.param(-0.5, 5, id="anticorrelated"),
        pytest.param(0.5, 6, id="phi-0.5"),
        pytest.param(0.9, 7, id="phi-0.9"),
    ],
)
def test_block_average_ar1_inefficiency_matches_analytic(
    make_ar1, phi, seed
):
    """An AR(1) series gives the analytic inefficiency at each level (R2)."""
    n_frames = 2 ** 16
    x = make_ar1(phi, n_frames, seed)

    table = block_average(x)

    checked = table.n_blocks >= MIN_BLOCKS_FOR_GAUSSIAN
    expected = ar1_inefficiency(phi, table.block_size[checked])
    sigma = expected * relative_gaussian_error(table.n_blocks[checked])
    pulls = (table.statistical_inefficiency[checked] - expected) / sigma
    assert_allclose(pulls, np.zeros_like(pulls), rtol=0.0, atol=SIGMA)


@pytest.mark.parametrize(
    "phi,seed,method",
    [
        pytest.param(0.5, 8, "lee", id="phi-0.5-lee"),
        pytest.param(0.9, 9, "lee", id="phi-0.9-lee"),
        pytest.param(0.5, 8, "sem_err", id="phi-0.5-sem_err"),
        pytest.param(0.9, 9, "sem_err", id="phi-0.9-sem_err"),
    ],
)
def test_sem_blocked_ar1_matches_analytic_error(make_ar1, phi, seed, method):
    """The blocked SEM of an AR(1) series is the analytic error at the level
    each rule selects, and a long series raises no warning (R2)."""
    n_frames = 2 ** 16
    x = make_ar1(phi, n_frames, seed)
    variance = 1.0 / (1.0 - phi ** 2)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        table = block_average(x)
        level = select_plateau(table, method=method)
        estimate = sem_blocked(x, method=method)

    assert level >= 0
    block_size = table.block_size[level]
    expected = np.sqrt(
        ar1_inefficiency(phi, np.array([block_size]))[0]
        * variance / n_frames
    )
    tolerance = SIGMA * expected * relative_gaussian_error(
        np.array([table.n_blocks[level]])
    )[0] / 2.0  # the SEM is a square root, so half the relative error
    assert estimate == pytest.approx(expected, abs=tolerance)


@pytest.mark.parametrize(
    "period,seed",
    [
        pytest.param(8, 10, id="period-8"),
        pytest.param(32, 11, id="period-32"),
    ],
)
def test_block_average_block_constant_below_period_is_exact(
    make_block_constant, period, seed
):
    """Below the period a block-constant series has an exact inefficiency
    B (n-1)/(n-B) (R3)."""
    n_frames = 2 ** 14
    x = make_block_constant(period, n_frames, seed)

    table = block_average(x)

    below = table.block_size <= period
    expected = block_constant_inefficiency(
        period, n_frames, table.block_size[below]
    )
    assert_allclose(
        table.statistical_inefficiency[below], expected, rtol=FLOAT64_RTOL
    )


@pytest.mark.parametrize(
    "period,seed",
    [
        pytest.param(8, 10, id="period-8"),
        pytest.param(32, 11, id="period-32"),
    ],
)
def test_block_average_block_constant_plateaus_at_period(
    make_block_constant, period, seed
):
    """At and above the period the inefficiency plateaus at L (n-1)/(n-L)
    within statistical error (R3)."""
    n_frames = 2 ** 14
    x = make_block_constant(period, n_frames, seed)

    table = block_average(x)

    checked = (table.block_size >= period) & (
        table.n_blocks >= MIN_BLOCKS_FOR_GAUSSIAN
    )
    expected = block_constant_inefficiency(
        period, n_frames, table.block_size[checked]
    )
    sigma = expected * relative_gaussian_error(table.n_blocks[checked])
    pulls = (table.statistical_inefficiency[checked] - expected) / sigma
    assert_allclose(pulls, np.zeros_like(pulls), rtol=0.0, atol=SIGMA)


# --- contracts: shapes, dtypes, independence of observables ------------------


def test_block_average_treats_trailing_axes_independently(make_ar1):
    """Each observable of an n-D input gets the table of its own 1-D series."""
    first = make_ar1(0.3, 2 ** 10, 12)
    second = make_ar1(0.8, 2 ** 10, 13)
    stacked = np.column_stack([first, second]).reshape(2 ** 10, 2, 1)

    table = block_average(stacked)

    assert_allclose(
        table.sem[:, 0, 0], block_average(first).sem, rtol=FLOAT64_RTOL
    )
    assert_allclose(
        table.statistical_inefficiency[:, 1, 0],
        block_average(second).statistical_inefficiency,
        rtol=FLOAT64_RTOL,
    )


@pytest.mark.parametrize(
    "trailing_shape",
    [
        pytest.param((), id="scalar-series"),
        pytest.param((3,), id="vector-series"),
        pytest.param((2, 3), id="matrix-series"),
    ],
)
def test_block_average_shapes_follow_the_input(trailing_shape):
    """Per-level arrays are (n_levels, *trailing); counts are (n_levels,)."""
    n_frames = 256
    x = np.random.default_rng(14).standard_normal((n_frames, *trailing_shape))

    table = block_average(x)

    assert table.observable_shape == trailing_shape
    assert table.block_size.shape == (table.n_levels,)
    assert table.n_blocks.shape == (table.n_levels,)
    assert table.mean.shape == (table.n_levels, *trailing_shape)
    assert table.sem.shape == (table.n_levels, *trailing_shape)
    assert table.sem_err.shape == (table.n_levels, *trailing_shape)
    assert table.statistical_inefficiency.shape == (
        table.n_levels, *trailing_shape
    )
    assert table.inefficiency_err.shape == (table.n_levels, *trailing_shape)


def test_select_plateau_returns_int_for_scalar_series():
    """A scalar series gives a plain int level and float wrapped estimates."""
    x = np.random.default_rng(15).standard_normal(2 ** 12)

    level = select_plateau(block_average(x), method="lee")

    assert isinstance(level, int)
    assert isinstance(sem_blocked(x, method="lee"), float)
    assert isinstance(statistical_inefficiency(x, method="lee"), float)


def test_select_plateau_returns_array_per_observable():
    """An n-D series gives one level and one estimate per observable."""
    x = np.random.default_rng(16).standard_normal((2 ** 12, 2, 3))

    level = select_plateau(block_average(x), method="lee")

    assert level.shape == (2, 3)
    assert level.dtype == np.int64
    assert sem_blocked(x, method="lee").shape == (2, 3)
    assert statistical_inefficiency(x, method="lee").shape == (2, 3)


@pytest.mark.parametrize(
    "dtype",
    [
        pytest.param(np.float32, id="float32"),
        pytest.param(np.int64, id="int64"),
        pytest.param(np.float64, id="float64"),
    ],
)
def test_block_average_output_dtypes(dtype):
    """Statistics are float64 and counts int64 whatever the input dtype."""
    x = np.random.default_rng(17).integers(0, 100, 512).astype(dtype)

    table = block_average(x)

    assert table.block_size.dtype == np.int64
    assert table.n_blocks.dtype == np.int64
    assert table.mean.dtype == np.float64
    assert table.var.dtype == np.float64
    assert table.sem.dtype == np.float64
    assert table.sem_err.dtype == np.float64
    assert table.statistical_inefficiency.dtype == np.float64
    assert table.inefficiency_err.dtype == np.float64


def test_block_average_leaves_input_unchanged():
    """The input array is not modified by the in-place blocking buffer."""
    x = np.random.default_rng(18).standard_normal((512, 2))
    before = x.copy()

    block_average(x)

    assert_array_equal(x, before)


# --- the selection rules on hand-built tables --------------------------------


@pytest.mark.parametrize(
    "sem_ratios,n_frames,expected_level",
    [
        # n = 1024: a level qualifies when B**3 > 2048 r**4; with r = 1
        # everywhere that is B >= 16, so the run from the top ends at k = 4.
        pytest.param([1, 1, 1, 1, 1, 1, 1, 1], 1024, 4, id="flat-table"),
        # n = 512, r = 0.5 for k <= 3: at k = 2, B**3 = 64 equals
        # 2 * 512 * 0.5**4 = 64, so strict '>' stops the run at k = 3; a
        # non-strict rule would continue to k = 2.
        pytest.param([1, 0.5, 0.5, 0.5, 1, 1, 1], 512, 3, id="strict"),
        # n = 4096: a level qualifies when B**3 > 8192 r**4, so with r = 1
        # that is B >= 32. r = 10 at k = 5 breaks the run: k = 7, 6 qualify,
        # k = 5 does not, so the rule stops at k = 6 although k = 5 would
        # qualify with r = 1 and a smallest-anywhere rule would pick it.
        pytest.param([1, 1, 1, 1, 1, 10, 1, 1], 4096, 6, id="broken-run"),
    ],
)
def test_select_plateau_lee_scans_down_from_largest_block(
    make_table, sem_ratios, n_frames, expected_level
):
    """The Lee rule takes the smallest block of the qualifying run that
    includes the largest block, with a strict inequality."""
    table = make_table(sem_ratios, n_frames)

    level = select_plateau(table, method="lee")

    assert level == expected_level


def test_select_plateau_lee_reports_no_plateau_when_top_level_fails(
    make_table
):
    """A top level that fails the Lee inequality gives -1 with a warning."""
    # n = 1024, k = 7: B**3 = 2**21 is not above 2048 * 10**4.
    table = make_table([1, 1, 1, 1, 1, 1, 1, 10], 1024)

    with pytest.warns(RuntimeWarning, match="No plateau found"):
        level = select_plateau(table, method="lee")

    assert level == -1


def test_select_plateau_lee_warns_when_fewer_than_fifty_blocks_remain(
    make_table
):
    """Selecting a block size at or above n/50 raises Lee's reliability
    warning while still returning the level."""
    # n = 1024, r = 10 for k <= 4: the run from the top stops at k = 5,
    # B = 32, and 50 * 32 >= 1024.
    table = make_table([1, 10, 10, 10, 10, 1, 1, 1], 1024)

    with pytest.warns(RuntimeWarning, match="not below n_frames / 50"):
        level = select_plateau(table, method="lee")

    assert level == 5


@pytest.mark.parametrize(
    "sem_values,sem_err_values,expected_level",
    [
        # Level 3 is the first whose SEM equals the next two within 0.05.
        pytest.param(
            [1.0, 1.5, 1.9, 2.0, 2.0, 2.0, 2.0, 2.0],
            [0.05] * 8,
            3,
            id="first-flat-level",
        ),
        # Level 2 has a large error bar of its own, but the rule compares
        # against the error bars of levels 3 and 4, which it fails.
        pytest.param(
            [1.0, 1.5, 1.9, 2.0, 2.0, 2.0, 2.0, 2.0],
            [0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05],
            3,
            id="next-levels-error-bars",
        ),
    ],
)
def test_select_plateau_sem_err_takes_first_level_within_next_two_errors(
    make_table, sem_values, sem_err_values, expected_level
):
    """The sem_err rule picks the first level whose SEM lies within the
    error bar of each of the next two levels."""
    table = make_table(sem_values, 1024, sem_err_values)

    level = select_plateau(table, method="sem_err")

    assert level == expected_level


@pytest.mark.parametrize(
    "sem_values",
    [
        pytest.param([1.0, 2.0, 4.0, 8.0, 16.0], id="still-rising"),
        pytest.param([1.0, 1.0], id="two-levels-only"),
    ],
)
def test_select_plateau_sem_err_reports_no_plateau(make_table, sem_values):
    """A table that never flattens, or has fewer than three levels, gives
    -1 with a warning."""
    table = make_table(sem_values, 1024, [0.01] * len(sem_values))

    with pytest.warns(RuntimeWarning, match="No plateau found"):
        level = select_plateau(table, method="sem_err")

    assert level == -1


@pytest.mark.parametrize("method", ["lee", "sem_err"])
def test_wrappers_return_nan_for_a_drifting_series(method):
    """A series with a linear drift has no plateau: the wrappers give NaN."""
    x = np.arange(1024, dtype=np.float64)

    with pytest.warns(RuntimeWarning, match="No plateau found"):
        estimate = sem_blocked(x, method=method)

    assert np.isnan(estimate)


@pytest.mark.parametrize("method", ["lee", "sem_err"])
def test_wrappers_return_the_selected_row_of_the_table(make_ar1, method):
    """sem_blocked and statistical_inefficiency equal the table entries at
    the level select_plateau chooses."""
    x = make_ar1(0.7, 2 ** 14, 19)
    table = block_average(x)
    level = select_plateau(table, method=method)

    blocked_sem = sem_blocked(x, method=method)
    inefficiency = statistical_inefficiency(x, method=method)

    assert blocked_sem == pytest.approx(table.sem[level], rel=FLOAT64_RTOL)
    assert inefficiency == pytest.approx(
        table.statistical_inefficiency[level], rel=FLOAT64_RTOL
    )


# --- invalid input -----------------------------------------------------------


@pytest.mark.parametrize(
    "x,min_blocks,match",
    [
        pytest.param(3.0, 4, "got a scalar", id="scalar"),
        pytest.param([1.0, np.nan, 2.0, 3.0], 4, "finite", id="nan"),
        pytest.param([1.0, np.inf, 2.0, 3.0], 4, "finite", id="inf"),
        pytest.param([1.0, 2.0, 3.0], 4, "at least min_blocks", id="short"),
        pytest.param([1.0, 2.0, 3.0, 4.0], 1, "at least 2", id="min-blocks"),
        pytest.param(
            np.column_stack([np.arange(8.0), np.ones(8)]),
            4,
            "constant series",
            id="constant-observable",
        ),
    ],
)
def test_block_average_rejects_invalid_input(x, min_blocks, match):
    """Scalars, non-finite values, too few frames, a min_blocks below 2, and
    a constant observable raise ValueError."""
    with pytest.raises(ValueError, match=match):
        block_average(x, min_blocks=min_blocks)


def test_select_plateau_rejects_unknown_method(make_table):
    """An unknown method name raises ValueError naming the choices."""
    table = make_table([1, 1, 1, 1, 1], 256)

    with pytest.raises(ValueError, match="choose 'lee' or 'sem_err'"):
        select_plateau(table, method="bogus")


def test_sem_blocked_rejects_unknown_method():
    """The wrappers pass an unknown method name through as ValueError."""
    x = np.random.default_rng(20).standard_normal(256)

    with pytest.raises(ValueError, match="Unknown plateau method"):
        sem_blocked(x, method="bogus")
