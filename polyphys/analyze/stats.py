"""
==========================================================
:mod:`polyphys.analyze.stats`
==========================================================


The :mod:`polyphys.analyze.stats` module estimates the uncertainty of a
time average over correlated molecular-dynamics frames with the blocking
method of Flyvbjerg and Petersen.

Consecutive frames of a trajectory are not independent. The textbook
standard error of the mean, :math:`\\sigma / \\sqrt{n}`, where :math:`n` is
the number of frames and :math:`\\sigma` the standard deviation of the
observable over them, counts every frame as a fresh observation and so
understates the true error, often by a large factor. Blocking averages
neighbouring frames in pairs, repeatedly, and recomputes the standard error at
each level. Once the blocks are longer than the correlation time, the lag over
which the observable's autocorrelation decays, the block means are effectively
independent and the estimate stops growing: that plateau is the error bar to
report. The ratio of the plateau variance of the mean to the naive one, the
*statistical inefficiency* :math:`s`, says how many frames one independent
sample is worth; the true standard error is :math:`\\sigma \\sqrt{s / n}`.

The module distinguishes two kinds of average that look alike in a table:

- a *time average* over the frames of one run (one ``whole``), whose frames
  are correlated and whose error needs this module;
- an *ensemble average* over independent replicates (one ``ensemble``), whose
  members are independent runs and whose error is the naive standard error of
  the replicate means, as :func:`polyphys.analyze.measurer.sem` computes it.

Notation
========

:math:`x_t`
    The value of one observable at frame :math:`t = 1, \\dots, n`.
:math:`n`
    The number of frames, ``n_frames``.
:math:`\\sigma^2`
    The variance of :math:`x_t` over the frames, with ``ddof=1``.
:math:`k`
    The blocking level, :math:`k = 0, 1, 2, \\dots`.
:math:`B`
    The block size in frames at level :math:`k`, :math:`B = 2^k`.
:math:`m`
    The number of blocks at that level, :math:`\\lfloor n / B \\rfloor`.
:math:`\\mathrm{SE}(B)`
    The standard error of the mean estimated from the :math:`m` block means;
    :math:`\\mathrm{SE}(1) = \\sigma / \\sqrt{n}` is the naive estimate.
:math:`s(B)`
    The statistical inefficiency estimated at block size :math:`B`: the
    variance of the block means times :math:`B`, divided by :math:`\\sigma^2`.
    When :math:`n` is a power of two, :math:`s(B) = (\\mathrm{SE}(B) /
    \\mathrm{SE}(1))^2`.
:math:`s`
    The statistical inefficiency, the plateau of :math:`s(B)` at block sizes
    beyond the correlation time.
:math:`\\rho(t)`
    The normalized autocorrelation function of :math:`x_t` at a lag of
    :math:`t` frames, with :math:`\\rho(0) = 1`; for a stationary series
    :math:`s = \\sum_{t=-\\infty}^{\\infty} \\rho(t)`.
:math:`\\tau_\\mathrm{int}`
    The integrated autocorrelation time in Wolff's (2004) convention,
    :math:`\\tau_\\mathrm{int} = s / 2`, in frames.

Contents
========

The full blocking table:
:func:`block_average` returns a :class:`BlockingResult`, one row per blocking
level, so the caller can read the plateau directly, as Flyvbjerg and Petersen
do.

Choosing the plateau:
:func:`select_plateau` applies one of two rules to a table, chosen with the
``method`` keyword; there is no default, because the rule is an estimator
choice.

``'lee'``
    The block size satisfies :math:`B^3 > 2 n (\\mathrm{SE}(B) /
    \\mathrm{SE}(1))^4`, scanning from the largest block downward, as Lee et
    al. (2011), Eq. (14), and the ``pyblock`` package implement it. Wolff
    (2004), Eq. (47), derives the same balance between the systematic and
    the statistical error of the error estimate, with a constant four times
    smaller in :math:`B^3`, because his correlation time is
    :math:`\\tau_\\mathrm{int} = s / 2` where Lee et al. use :math:`s`.
``'sem_err'``
    The first level whose standard error lies within one error bar of each of
    the next two levels. A convenience heuristic without a published source;
    it tends to pick the start of the plateau and so to err on the low side.

Under either rule, :func:`select_plateau` warns when the series fails the
acceptance criterion of Lee et al. (2011): the block size the ``'lee'`` rule
selects must leave more than 50 blocks, which amounts to :math:`n / s > 500`,
more than about 500 independent samples.

Convenience wrappers:
:func:`sem_blocked` and :func:`statistical_inefficiency` run both steps and
return the selected row's value.

Conventions
===========

- Frames run along axis 0 of the input, shape ``(n_frames, ...)``. Every
  trailing axis is an independent observable: a per-frame histogram
  ``(n_frames, n_bins)`` gives one error bar per bin.
- Block sizes are powers of two, ``1, 2, 4, ...``. Level ``k`` has block size
  ``2**k``; when a level has an odd number of blocks, the last one is dropped
  before the next pairing, as Allen and Tildesley's example code does.
- Variances use ``ddof=1``. The error of a standard error,
  :math:`\\mathrm{SE} / \\sqrt{2 (m - 1)}` for :math:`m` blocks, and the error
  of the inefficiency, :math:`s \\sqrt{2 / (m - 1)}`, are the leading-order
  Gaussian estimates (Lee et al. 2011, Eq. (15)).
- Results are in the units of the input; the statistical inefficiency is
  dimensionless, in frames per independent sample.
- The input must be stationary: discard the equilibration transient before
  calling. A drift or a rare switch between metastable states shows up as a
  standard error that keeps rising with block size, and the selectors then
  report that no plateau was found rather than pick a level.

References
==========

- H. Flyvbjerg and H. G. Petersen, "Error estimates on averages of correlated
  data", J. Chem. Phys. 91, 461-466 (1989).
  https://doi.org/10.1063/1.457480
- U. Wolff (ALPHA Collaboration), "Monte Carlo errors with less errors",
  Comput. Phys. Commun. 156, 143-153 (2004); erratum 176, 383 (2007).
  https://doi.org/10.1016/S0010-4655(03)00467-3
- R. M. Lee, G. J. Conduit, N. Nemec, P. López Ríos, and N. D. Drummond,
  "Strategies for improving the efficiency of quantum Monte Carlo
  calculations", Phys. Rev. E 83, 066706 (2011); erratum 85, 029908 (2012).
  https://doi.org/10.1103/PhysRevE.83.066706
- M. P. Allen and D. J. Tildesley, Computer Simulation of Liquids, 2nd ed.
  (Oxford University Press, 2017); example code ``error_calc.py`` at
  https://github.com/Allen-Tildesley/examples/blob/master/python_examples/error_calc.py
- J. S. Spencer, ``pyblock``: ``pyblock.blocking.find_optimal_block``.
  https://github.com/jsspencer/pyblock

.. versionadded:: 0.5.0
"""

import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt

__all__ = [
    'BlockingResult',
    'PlateauMethod',
    'block_average',
    'select_plateau',
    'sem_blocked',
    'statistical_inefficiency',
]

PlateauMethod = Literal['lee', 'sem_err']
"""The plateau-selection rules :func:`select_plateau` accepts."""

#: Lee et al. (2011) accept a blocking analysis only when their optimal block
#: size is below ``n_frames / LEE_MIN_BLOCK_RATIO``, that is, when more than
#: this many blocks remain at that level. :func:`select_plateau` applies the
#: check to the level the ``'lee'`` rule selects, whichever rule is asked for.
LEE_MIN_BLOCK_RATIO = 50


@dataclass(frozen=True)
class BlockingResult:
    """
    The blocking table of one time series, one row per blocking level.

    Level ``k`` averages the frames in contiguous blocks of ``2**k`` frames
    and treats the block means as the sample. Every per-level array has the
    level along axis 0; the trailing axes are those of one frame of the input,
    one entry per observable.

    Attributes
    ----------
    n_frames : int
        Number of frames in the input, the sample size at level 0.
    block_size : numpy.ndarray
        Frames per block at each level, ``2**k``; shape ``(n_levels,)``,
        ``int64``.
    n_blocks : numpy.ndarray
        Blocks at each level, ``n_frames // block_size``; shape
        ``(n_levels,)``, ``int64``.
    mean : numpy.ndarray
        Mean of the block means, in the units of the input; shape
        ``(n_levels, *observable_shape)``, ``float64``. It changes across
        levels only when a trailing frame is dropped at an odd level.
    var : numpy.ndarray
        Variance of the block means with ``ddof=1``, in squared input units.
    sem : numpy.ndarray
        Standard error of the mean at each level, ``sqrt(var / n_blocks)``,
        in the units of the input. Level 0 is the naive estimate.
    sem_err : numpy.ndarray
        Leading-order error of ``sem``, ``sem / sqrt(2 (n_blocks - 1))``.
    statistical_inefficiency : numpy.ndarray
        ``block_size * var / var[0]``, dimensionless: frames per independent
        sample, as estimated at that level. Level 0 is exactly 1.
    inefficiency_err : numpy.ndarray
        Leading-order error of the inefficiency,
        ``statistical_inefficiency * sqrt(2 / (n_blocks - 1))``. It neglects
        the error of ``var[0]``.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.stats import block_average
    >>> table = block_average(np.repeat([1.0, 3.0, 2.0, 4.0], 4))
    >>> table.n_levels, table.observable_shape
    (3, ())
    >>> table.block_size
    array([1, 2, 4])

    .. versionadded:: 0.5.0
    """

    n_frames: int
    block_size: npt.NDArray[np.int64]
    n_blocks: npt.NDArray[np.int64]
    mean: npt.NDArray[np.float64]
    var: npt.NDArray[np.float64]
    sem: npt.NDArray[np.float64]
    sem_err: npt.NDArray[np.float64]
    statistical_inefficiency: npt.NDArray[np.float64]
    inefficiency_err: npt.NDArray[np.float64]

    @property
    def n_levels(self) -> int:
        """Number of blocking levels, the length of axis 0."""
        return int(self.block_size.shape[0])

    @property
    def observable_shape(self) -> tuple[int, ...]:
        """Shape of one frame of the input: ``()`` for a scalar series."""
        return tuple(self.mean.shape[1:])


def _as_series(x: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """
    Copy `x` to a finite ``float64`` array with at least one axis.

    Raises
    ------
    ValueError
        If `x` is a scalar or holds a NaN or an infinity.
    """
    data = np.array(x, dtype=np.float64)
    if data.ndim == 0:
        raise ValueError(
            'Input must have frames along axis 0; got a scalar.'
        )
    if not np.all(np.isfinite(data)):
        raise ValueError('Input must be finite; it holds NaN or inf.')
    return data


def block_average(x: npt.ArrayLike, min_blocks: int = 4) -> BlockingResult:
    """
    Compute the Flyvbjerg-Petersen blocking table of a time series.

    Starting from the frames themselves (block size 1), the series is
    replaced by the means of consecutive pairs, and the mean, variance,
    standard error, and statistical inefficiency are recorded at each level,
    until the next halving would leave fewer than `min_blocks` blocks. The
    full table is option (a) of the plateau choice: the caller reads where
    the standard error stops growing.

    Parameters
    ----------
    x : array_like
        Time series with frames along axis 0, shape ``(n_frames, ...)``. Every
        trailing axis is an independent observable. Any numeric dtype; the
        computation is in ``float64`` and `x` is not modified.
    min_blocks : int, default 4
        Smallest number of blocks kept; at least 2. The last level has at
        least this many blocks. With 4 blocks the error of the standard error
        is already 41% of its value, so levels this deep are mostly a guard
        against reading noise as a plateau.

    Returns
    -------
    BlockingResult
        The table; see its attributes for shapes and units.

    Raises
    ------
    ValueError
        If `x` is a scalar, holds NaN or inf, has fewer than `min_blocks`
        frames, or has an observable that is constant across frames (its
        inefficiency would be 0/0); or if `min_blocks` is below 2.

    Notes
    -----
    For a series :math:`x_t` of :math:`n` frames, level :math:`k = 0, 1,
    \\dots` has block size :math:`B_k = 2^k` frames, :math:`m_k = \\lfloor n
    / B_k \\rfloor` blocks, and block means :math:`\\bar{x}_j^{(k)}`, the
    average of :math:`x_t` over block :math:`j = 1, \\dots, m_k`. With
    :math:`\\langle \\bar{x}^{(k)} \\rangle` the mean of the block means, the
    table records

    .. math::

        \\mathrm{var}_k = \\frac{1}{m_k - 1} \\sum_{j=1}^{m_k}
            \\left(\\bar{x}_j^{(k)}
            - \\langle \\bar{x}^{(k)} \\rangle\\right)^2,
        \\quad
        \\mathrm{SE}_k = \\sqrt{\\mathrm{var}_k / m_k},
        \\quad
        s_k = B_k \\, \\mathrm{var}_k / \\mathrm{var}_0 ,

    the variance of the block means, the standard error of the mean, and the
    statistical inefficiency at level :math:`k`; :math:`\\mathrm{var}_0` is
    the variance of the frames themselves. As :math:`B_k` grows past the
    correlation time, :math:`\\mathrm{SE}_k` reaches a plateau at the true
    standard error of the time average and :math:`s_k` at the statistical
    inefficiency :math:`s = \\sum_{t=-\\infty}^{\\infty} \\rho(t)`, where
    :math:`\\rho(t)` is the normalized autocorrelation function of
    :math:`x_t` at a lag of :math:`t` frames. This is twice the integrated
    autocorrelation time :math:`\\tau_\\mathrm{int}` of Wolff (2004),
    Eq. (23) and (25).

    When a level has an odd number of blocks, the trailing block is left out
    of the next pairing, so the frames used shrink by at most one block per
    level and ``mean`` can drift slightly between levels.

    With :math:`n_\\mathrm{obs}` the number of observables, the product of
    the trailing dimensions of `x` (1 for a 1-D series), time complexity is
    :math:`O(n \\, n_\\mathrm{obs})`, since the levels touch :math:`n + n/2
    + n/4 + \\dots < 2n` frames in all. Space is one ``float64`` copy of the
    input, halved in place, plus the table of :math:`O(\\log_2 n \\,
    n_\\mathrm{obs})` entries.

    References
    ----------
    H. Flyvbjerg and H. G. Petersen, J. Chem. Phys. 91, 461 (1989),
    https://doi.org/10.1063/1.457480. Implementation after M. P. Allen and
    D. J. Tildesley, Computer Simulation of Liquids, 2nd ed. (2017), example
    code ``error_calc.py``.

    Examples
    --------
    A signal of :math:`n = 16` frames that holds each value for :math:`L = 4`
    frames has an exact table: for every block size :math:`B \\le L` the block
    means are the held values themselves, so the inefficiency is
    :math:`s(B) = B (n - 1) / (n - B)`, here 1, 30/14 and 5.

    >>> import numpy as np
    >>> from polyphys.analyze.stats import block_average
    >>> x = np.repeat([1.0, 3.0, 2.0, 4.0], 4)
    >>> table = block_average(x, min_blocks=4)
    >>> table.block_size, table.n_blocks
    (array([1, 2, 4]), array([16,  8,  4]))
    >>> table.statistical_inefficiency.round(4)
    array([1.    , 2.1429, 5.    ])
    >>> table.sem.round(4)
    array([0.2887, 0.4226, 0.6455])

    A two-column input gives one table column per observable:

    >>> table = block_average(np.column_stack([x, 2.0 * x]))
    >>> table.observable_shape, table.sem.shape
    ((2,), (3, 2))
    >>> np.allclose(table.sem[:, 1], 2.0 * table.sem[:, 0])
    True

    .. versionadded:: 0.5.0
    """
    if min_blocks < 2:
        raise ValueError(
            f'min_blocks must be at least 2; got {min_blocks}.'
        )
    data = _as_series(x)
    n_frames = int(data.shape[0])
    if n_frames < min_blocks:
        raise ValueError(
            f'Input has {n_frames} frames but min_blocks is {min_blocks}; '
            'at least min_blocks frames are needed.'
        )
    var_naive = data.var(axis=0, ddof=1)
    if np.any(var_naive == 0.0):
        raise ValueError(
            'Every observable must vary across frames; a constant series has '
            'an undefined statistical inefficiency.'
        )

    block_sizes: list[int] = []
    block_counts: list[int] = []
    means: list[npt.NDArray[np.float64]] = []
    variances: list[npt.NDArray[np.float64]] = []
    n_blocks = n_frames
    block_size = 1
    while True:
        blocks = data[:n_blocks]
        block_sizes.append(block_size)
        block_counts.append(n_blocks)
        means.append(blocks.mean(axis=0))
        variances.append(blocks.var(axis=0, ddof=1))
        if n_blocks // 2 < min_blocks:
            break
        n_blocks //= 2
        # The blocking transformation: pair consecutive blocks, dropping a
        # trailing unpaired one, and overwrite the front of the buffer.
        data[:n_blocks] = 0.5 * (
            data[0:2 * n_blocks:2] + data[1:2 * n_blocks:2]
        )
        block_size *= 2

    size_arr = np.asarray(block_sizes, dtype=np.int64)
    count_arr = np.asarray(block_counts, dtype=np.int64)
    mean = np.stack(means)
    var = np.stack(variances)
    # Broadcast the per-level counts against the trailing observable axes.
    level_shape = (-1,) + (1,) * (var.ndim - 1)
    m = count_arr.reshape(level_shape).astype(np.float64)
    b = size_arr.reshape(level_shape).astype(np.float64)
    sem = np.sqrt(var / m)
    sem_err = sem / np.sqrt(2.0 * (m - 1.0))
    inefficiency = b * var / var_naive
    inefficiency_err = inefficiency * np.sqrt(2.0 / (m - 1.0))
    return BlockingResult(
        n_frames=n_frames,
        block_size=size_arr,
        n_blocks=count_arr,
        mean=mean,
        var=var,
        sem=sem,
        sem_err=sem_err,
        statistical_inefficiency=inefficiency,
        inefficiency_err=inefficiency_err,
    )


def _select_lee(result: BlockingResult) -> npt.NDArray[np.int64]:
    """
    Apply the rule of Lee et al. (2011), Eq. (14), to a blocking table.

    Level ``k`` qualifies when ``B_k**3 > 2 n (SE_k / SE_0)**4``. Scanning
    from the largest block downward, the selected level is the last
    qualifying one before the first failure; ``-1`` marks an observable whose
    largest level does not qualify. :func:`select_plateau` checks the result
    against the acceptance criterion of the same paper.
    """
    level_shape = (-1,) + (1,) * (result.sem.ndim - 1)
    b = result.block_size.reshape(level_shape).astype(np.float64)
    ratio = result.sem / result.sem[0]
    qualifies = b ** 3 > 2.0 * result.n_frames * ratio ** 4
    # Length of the run of qualifying levels that starts at the top level.
    run_from_top = np.cumprod(qualifies[::-1], axis=0)
    run_length = run_from_top.sum(axis=0)
    return np.where(
        run_length > 0, result.n_levels - run_length, -1
    ).astype(np.int64)


def _select_sem_err(result: BlockingResult) -> npt.NDArray[np.int64]:
    """
    Select the first level whose SEM agrees with the next two within error.

    Level ``k`` qualifies when ``|SE_{k+j} - SE_k| <= sem_err_{k+j}`` for
    ``j = 1, 2``; ``-1`` marks an observable with no such level.
    """
    sem = result.sem
    err = result.sem_err
    if result.n_levels < 3:
        return np.full(result.observable_shape, -1, dtype=np.int64)
    near_next = np.abs(sem[1:-1] - sem[:-2]) <= err[1:-1]
    near_second = np.abs(sem[2:] - sem[:-2]) <= err[2:]
    qualifies = near_next & near_second
    found = qualifies.any(axis=0)
    first = qualifies.argmax(axis=0)
    return np.where(found, first, -1).astype(np.int64)


def select_plateau(
    result: BlockingResult,
    *,
    method: PlateauMethod,
) -> int | npt.NDArray[np.int64]:
    """
    Pick the blocking level that estimates the plateau, by a named rule.

    Parameters
    ----------
    result : BlockingResult
        A table from :func:`block_average`.
    method : {'lee', 'sem_err'}
        The selection rule; required, because it is an estimator choice.

        ``'lee'``
            The largest-block-downward scan of Lee et al. (2011): a level
            qualifies when its block size :math:`B`, in frames, satisfies
            :math:`B^3 > 2 n (\\mathrm{SE}(B) / \\mathrm{SE}(1))^4`, where
            :math:`n` is the number of frames, :math:`\\mathrm{SE}(B)` the
            level's standard error, and :math:`\\mathrm{SE}(1)` the naive
            one at :math:`B = 1`. The selected level is the smallest block
            size of the qualifying run that includes the largest block. The
            rule balances the systematic underestimate of the error at small
            blocks against its growing statistical error at large ones.
        ``'sem_err'``
            The first level whose standard error is within one ``sem_err`` of
            each of the next two levels, that is, where the table stops
            rising. A convenience heuristic, not a published rule.

    Returns
    -------
    int or numpy.ndarray
        The selected level index into the table's axis 0: an ``int`` for a
        scalar series, else an ``int64`` array of shape
        ``result.observable_shape``. ``-1`` marks an observable for which no
        level qualifies.

    Raises
    ------
    ValueError
        If `method` is not one of the rules above.

    Warns
    -----
    RuntimeWarning
        When no level qualifies for some observable, which means the series
        is too short for its correlation time or is not stationary. Also,
        under either `method`, when an observable has a selected level but
        the block size the ``'lee'`` rule selects is not below
        ``n_frames / 50``, so that 50 or fewer blocks remain, or that rule
        finds no level: Lee et al. (2011) then consider the estimate
        unreliable. The selected level is still returned.

    Notes
    -----
    With :math:`n_\\mathrm{levels}` the number of rows of the table and
    :math:`n_\\mathrm{obs}` the number of observables, both rules are
    :math:`O(n_\\mathrm{levels} \\, n_\\mathrm{obs})` in time and space.

    Lee et al. write their Eq. (14) as :math:`B_\\mathrm{opt} = (2 n
    n_\\mathrm{corr}^2)^{1/3}`, the optimal block size for :math:`n` frames,
    with :math:`n_\\mathrm{corr}` the statistical inefficiency :math:`s`, the
    plateau of :math:`(\\mathrm{SE}(B) / \\mathrm{SE}(1))^2`; hence the
    fourth power in the rule above. They cite Wolff (2004), whose Eq. (47)
    reads :math:`B = \\tau (2N/\\tau)^{1/3}` for :math:`N` frames, with
    :math:`\\tau` the integrated autocorrelation time, which is
    :math:`s / 2`. The two prescriptions therefore differ by a factor of 4 in
    :math:`B^3`; neither constant is exact, as Wolff's derivation is an
    order-of-magnitude balance. This function implements Lee's form, as
    ``pyblock`` does.

    Lee et al. accept the analysis only when :math:`n / B_\\mathrm{opt} >
    50`, which with Eq. (14) reads :math:`n / s > 500`: the series must hold
    more than about 500 independent samples. That is a property of the
    series, not of the rule, so this function tests it at the ``'lee'``
    level under either `method`. The ``'sem_err'`` rule can stop at a level
    with more than 50 blocks on a series that fails the test, where the
    table levels off before the blocks outgrow the correlation time.

    References
    ----------
    R. M. Lee et al., Phys. Rev. E 83, 066706 (2011), Sec. IV, Eq. (13)-(15).
    U. Wolff, Comput. Phys. Commun. 156, 143 (2004), Eq. (46)-(48).

    Examples
    --------
    A signal that holds each value for 8 frames has a statistical
    inefficiency :math:`s \\approx 8`, so on the plateau
    :math:`(\\mathrm{SE}(B) / \\mathrm{SE}(1))^4 = s^2 \\approx 8^2`. With
    :math:`n = 16384` frames the Lee rule then needs :math:`B^3 > 2 \\cdot
    16384 \\cdot 8^2`, that is :math:`B > 128`, and it selects the first
    power of two above 128:

    >>> import numpy as np
    >>> from polyphys.analyze.stats import block_average, select_plateau
    >>> rng = np.random.default_rng(7)
    >>> x = np.repeat(rng.standard_normal(2048), 8)
    >>> table = block_average(x)
    >>> level = select_plateau(table, method='lee')
    >>> int(table.block_size[level])
    256

    The ``'sem_err'`` rule stops where the standard error first flattens,
    which for this signal is at or just past the plateau length:

    >>> level = select_plateau(table, method='sem_err')
    >>> bool(table.block_size[level] >= 8)
    True

    Holding each of 64 values for 256 frames gives :math:`s \\approx 256`,
    so the series holds only 64 independent samples. The ``'sem_err'`` rule
    stops at the period, where 64 blocks remain, but the Lee block size
    leaves fewer than 50, so both rules warn:

    >>> import warnings
    >>> rng = np.random.default_rng(7)
    >>> x = np.repeat(rng.standard_normal(64), 256)
    >>> table = block_average(x)
    >>> with warnings.catch_warnings(record=True) as caught:
    ...     warnings.simplefilter('always')
    ...     level = select_plateau(table, method='sem_err')
    >>> int(table.block_size[level]), int(table.n_blocks[level])
    (256, 64)
    >>> [type(w.message).__name__ for w in caught]
    ['RuntimeWarning']

    .. versionadded:: 0.5.0
    """
    if method not in ('lee', 'sem_err'):
        raise ValueError(
            f"Unknown plateau method {method!r}; choose 'lee' or 'sem_err'."
        )
    lee_level = _select_lee(result)
    level = lee_level if method == 'lee' else _select_sem_err(result)

    # Lee's acceptance criterion judges the series, so it uses the Lee block
    # size whichever rule selected the level. An observable without a
    # selected level gets the no-plateau warning below instead.
    lee_size = result.block_size[np.maximum(lee_level, 0)]
    too_short = (level >= 0) & (
        (lee_level < 0)
        | (LEE_MIN_BLOCK_RATIO * lee_size >= result.n_frames)
    )
    n_short = int(np.count_nonzero(too_short))
    if n_short:
        warnings.warn(
            f'For {n_short} observable(s) the block size of Lee et al. '
            f'(2011), Eq. (14), is not below n_frames / '
            f'{LEE_MIN_BLOCK_RATIO}, or no block size satisfies it: the '
            'series is too short for its correlation time, and Lee et al. '
            'consider the error estimate unreliable. Collect more frames.',
            RuntimeWarning,
            stacklevel=2,
        )
    n_missing = int(np.count_nonzero(level < 0))
    if n_missing:
        warnings.warn(
            f'No plateau found for {n_missing} observable(s): the series is '
            'too short for its correlation time, or not stationary. The '
            'selected level is -1 and the wrapped estimate NaN.',
            RuntimeWarning,
            stacklevel=2,
        )
    if level.ndim == 0:
        return int(level)
    return level


def _pick_level(
    table: npt.NDArray[np.float64],
    level: int | npt.NDArray[np.int64],
) -> float | npt.NDArray[np.float64]:
    """
    Take one row of a per-level array per observable; NaN where ``level < 0``.
    """
    index = np.asarray(level, dtype=np.int64)
    picked = np.take_along_axis(
        table, np.expand_dims(np.maximum(index, 0), 0), axis=0
    )[0]
    picked = np.where(index < 0, np.nan, picked)
    if picked.ndim == 0:
        return float(picked)
    return picked


def sem_blocked(
    x: npt.ArrayLike,
    *,
    method: PlateauMethod,
    min_blocks: int = 4,
) -> float | npt.NDArray[np.float64]:
    """
    Standard error of a time average over correlated frames.

    Runs :func:`block_average` and :func:`select_plateau` and returns the
    standard error at the selected level. This is the error bar to report
    for the mean of `x` over its frames.

    Parameters
    ----------
    x : array_like
        Time series with frames along axis 0, shape ``(n_frames, ...)``.
    method : {'lee', 'sem_err'}
        The plateau rule; see :func:`select_plateau`.
    min_blocks : int, default 4
        Passed to :func:`block_average`.

    Returns
    -------
    float or numpy.ndarray
        The blocked standard error, in the units of `x`: a ``float`` for a
        scalar series, else a ``float64`` array of shape ``x.shape[1:]``. NaN
        where no plateau was found.

    Raises
    ------
    ValueError
        As :func:`block_average` and :func:`select_plateau`.

    Warns
    -----
    RuntimeWarning
        As :func:`select_plateau`.

    Notes
    -----
    For independent frames the result agrees with
    :func:`polyphys.analyze.measurer.sem` within its own error bar; for
    correlated frames it is larger by :math:`\\sqrt{s}`, where :math:`s` is
    the statistical inefficiency, the number of frames that make one
    independent sample. With :math:`n` the number of frames and
    :math:`n_\\mathrm{obs}` the number of observables (the product of the
    trailing dimensions of `x`), time is :math:`O(n \\, n_\\mathrm{obs})` and
    space one ``float64`` copy of `x`.

    Examples
    --------
    Each value of this signal is held for 8 frames, so :math:`s \\approx 8`
    and the blocked error is about :math:`\\sqrt{8} \\approx 2.8` times the
    naive one:

    >>> import numpy as np
    >>> from polyphys.analyze.stats import sem_blocked
    >>> from polyphys.analyze.measurer import sem
    >>> rng = np.random.default_rng(7)
    >>> x = np.repeat(rng.standard_normal(2048), 8)
    >>> ratio = sem_blocked(x, method='lee') / sem(x)
    >>> bool(2.3 < ratio < 3.4)
    True

    .. versionadded:: 0.5.0
    """
    result = block_average(x, min_blocks=min_blocks)
    level = select_plateau(result, method=method)
    return _pick_level(result.sem, level)


def statistical_inefficiency(
    x: npt.ArrayLike,
    *,
    method: PlateauMethod,
    min_blocks: int = 4,
) -> float | npt.NDArray[np.float64]:
    """
    Frames per independent sample of a time series, from its blocking table.

    Runs :func:`block_average` and :func:`select_plateau` and returns the
    statistical inefficiency at the selected level: the ratio of the plateau
    variance of the mean to the naive one, so ``n_frames / s`` is the
    effective number of independent frames.

    Parameters
    ----------
    x : array_like
        Time series with frames along axis 0, shape ``(n_frames, ...)``.
    method : {'lee', 'sem_err'}
        The plateau rule; see :func:`select_plateau`.
    min_blocks : int, default 4
        Passed to :func:`block_average`.

    Returns
    -------
    float or numpy.ndarray
        Dimensionless, in frames; 1 for independent frames, below 1 for
        anticorrelated ones. A ``float`` for a scalar series, else a
        ``float64`` array of shape ``x.shape[1:]``. NaN where no plateau was
        found.

    Raises
    ------
    ValueError
        As :func:`block_average` and :func:`select_plateau`.

    Warns
    -----
    RuntimeWarning
        As :func:`select_plateau`.

    Notes
    -----
    For a stationary series :math:`x_t` whose normalized autocorrelation
    function at a lag of :math:`t` frames is :math:`\\rho(t)`, with
    :math:`\\rho(0) = 1`, the inefficiency is :math:`s =
    \\sum_{t=-\\infty}^{\\infty} \\rho(t) = 2 \\tau_\\mathrm{int}`, where
    :math:`\\tau_\\mathrm{int}` is the integrated autocorrelation time (Wolff
    2004, Eq. (23) and (25)). For an AR(1) process, :math:`x_t = \\phi
    x_{t-1} + \\epsilon_t` with coefficient :math:`|\\phi| < 1` and
    independent noise :math:`\\epsilon_t`, :math:`\\rho(t) = \\phi^{|t|}` and
    :math:`s = (1 + \\phi) / (1 - \\phi)`. With :math:`n` the number of
    frames and :math:`n_\\mathrm{obs}` the number of observables (the
    product of the trailing dimensions of `x`), time is :math:`O(n \\,
    n_\\mathrm{obs})` and space one ``float64`` copy of `x`.

    Examples
    --------
    Each value of this signal is held for 8 frames, so :math:`s \\approx 8`:

    >>> import numpy as np
    >>> from polyphys.analyze.stats import statistical_inefficiency
    >>> rng = np.random.default_rng(7)
    >>> x = np.repeat(rng.standard_normal(2048), 8)
    >>> s = statistical_inefficiency(x, method='lee')
    >>> bool(5.0 < s < 11.0)
    True

    .. versionadded:: 0.5.0
    """
    result = block_average(x, min_blocks=min_blocks)
    level = select_plateau(result, method=method)
    return _pick_level(result.statistical_inefficiency, level)
