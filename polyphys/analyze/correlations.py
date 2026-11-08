"""
==========================================================
:mod:`polyphys.analyze.correlations`
==========================================================


The :mod:`polyphys.analyze.correlations` module measures how a signal is
correlated with itself: the autocorrelation function of a time series with
its confidence band, the decay model fitted to such a curve, and the cosine
correlation between the bond vectors of a polymer chain, from which a
persistence length follows. Every routine is a pure numpy or scipy kernel on
arrays; reading time series or frames, choosing the lags to fit, and
averaging over frames are the caller's work.

Notation
========

``x``
    A time series, shape ``(n_frames,)`` or ``(n_frames, n_series)``, in the
    caller's unit; the frames are equally spaced in time and run along axis
    0, and each column is an independent series.
``nlags``
    The largest lag :math:`k` at which the autocorrelation is reported, in
    frames; a lag is the number of frames between the two values compared.
``acf``
    The sample autocorrelation function :math:`r_k`, the covariance of the
    series with itself :math:`k` frames later divided by its variance:
    :math:`r_0 = 1` and :math:`|r_k| \\le 1`.
``half_width``
    The half-width of the confidence interval around :math:`r_k` at level
    :math:`1 - \\alpha`, from Bartlett's large-sample variance of
    :math:`r_k`; the interval is ``acf - half_width`` to
    ``acf + half_width``.
*decay model*
    One of the forms of the stretched exponential (Kohlrausch-Williams-Watts)
    function :math:`A \\exp[-(x / \\tau)^\\beta] + B` listed in
    :data:`DECAY_MODELS`, with the *decay time* :math:`\\tau` in the unit of
    ``x``, the *stretching exponent* :math:`\\beta` (1 for a pure
    exponential), the *amplitude* :math:`A`, and the *offset* :math:`B`.
``positions``
    The coordinates of the monomers of one chain in chain order, shape
    ``(..., n_atoms, n_dims)``, unwrapped (no periodic image breaks a bond);
    any leading axes, such as frames, are carried through.
``topology``
    ``'linear'`` or ``'ring'``: whether the chain's last monomer is bonded
    back to its first.
*bond vector*
    :math:`\\mathbf{b}_i = \\mathbf{r}_{i+1} - \\mathbf{r}_i`, and on a ring
    also :math:`\\mathbf{b}_{n-1} = \\mathbf{r}_0 - \\mathbf{r}_{n-1}`;
    ``n_bonds`` is ``n_atoms - 1`` on a linear chain and ``n_atoms`` on a
    ring. :math:`\\hat{\\mathbf{u}}_i` is the unit bond vector.
*bond cosine correlation*
    :math:`C(k) = \\langle \\hat{\\mathbf{u}}_i \\cdot
    \\hat{\\mathbf{u}}_{i+k} \\rangle`, the mean cosine of the angle between
    two bonds :math:`k` bonds apart along the chain, averaged over the pairs
    in one frame; :math:`C(0) = 1`.

Contents
========

Autocorrelation:
:func:`acf` estimates :math:`r_k` for :math:`k = 0, \\dots,` ``nlags`` by
the fast Fourier transform, with Bartlett's confidence band, as an
:class:`AcfResult`; :func:`bartlett_half_width` is that band on its own,
for a stored :math:`r_k` at any level.

Decay fits:
:data:`DECAY_MODELS` names the decay models; :func:`initial_guess` reads a
starting point for each parameter off the data; :func:`fit_decay` fits a
model with :func:`scipy.optimize.curve_fit` and reports the parameters,
their standard errors, and whether the fit converged as a
:class:`DecayFit`, whose :meth:`DecayFit.evaluate` draws the fitted curve.

Bonds:
:func:`bond_vectors` differences the positions along the chain;
:func:`bond_cosine_correlation` averages the cosines between bond pairs at
every lag by the fast Fourier transform.

Conventions
===========

- The autocorrelation uses the ``1/n_frames`` estimator of the
  autocovariance at every lag, which keeps the sequence positive
  semi-definite and is the default of ``statsmodels.tsa.stattools.acf``;
  dividing by ``n_frames - k`` instead is unbiased per lag but can give
  :math:`|r_k| > 1`.
- Bartlett's variance for lag :math:`k` sums the squared autocorrelations at
  the lags below it, :math:`\\mathrm{var}(r_k) = (1 + 2 \\sum_{j=1}^{k-1}
  r_j^2) / n_\\mathrm{frames}` for :math:`k \\ge 2`, with
  :math:`\\mathrm{var}(r_1) = 1 / n_\\mathrm{frames}` and
  :math:`\\mathrm{var}(r_0) = 0`; this is the convention of
  ``statsmodels``, so a band computed here matches one computed there.
- Decay models are written in the decay time :math:`\\tau`, not the rate
  :math:`\\omega = \\tau^{-\\beta}` of some research code; convert with
  :math:`\\tau = \\omega^{-1/\\beta}`.
- A failed fit never raises: :attr:`DecayFit.converged` is false and the
  parameters are NaN.
- On a ring every lag :math:`k` has ``n_bonds`` pairs, going round the ring,
  so :math:`C(k) = C(n_\\mathrm{bonds} - k)` and the lags that carry
  information are :math:`0 \\le k \\le n_\\mathrm{bonds} / 2`. On a linear
  chain lag :math:`k` has ``n_bonds - k`` pairs.
- Complexities are stated per call in ``n_frames`` (:math:`n`),
  ``n_series`` (:math:`m`), ``n_bonds`` (:math:`N`), the number of leading
  frames of a positions array (:math:`F`), and ``n_points`` of a fit.

References
==========

- M. S. Bartlett, "On the theoretical specification and sampling properties
  of autocorrelated time-series", Supplement to the Journal of the Royal
  Statistical Society 8, 27-41 (1946). https://doi.org/10.2307/2983611
- NIST/SEMATECH e-Handbook of Statistical Methods, section 1.3.3.1,
  "Autocorrelation Plot".
  https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm
- G. Williams and D. C. Watts, "Non-symmetrical dielectric relaxation
  behavior arising from a simple empirical decay function", Trans. Faraday
  Soc. 66, 80-85 (1970). https://doi.org/10.1039/TF9706600080
- SciPy, :func:`scipy.optimize.curve_fit`.
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
- M. Rubinstein and R. H. Colby, Polymer Physics (Oxford University Press,
  2003), Chapter 2, on the bond correlation and the persistence length.
- MDAnalysis, :class:`MDAnalysis.analysis.polymer.PersistenceLength`, the
  same bond autocorrelation and exponential fit for linear chains.
  https://docs.mdanalysis.org/stable/documentation_pages/analysis/polymer.html
- MDAnalysis User Guide, "Determining the persistence length of a polymer".
  https://userguide.mdanalysis.org/stable/examples/analysis/polymers_and_membranes/polymer.html

.. versionadded:: 0.6.0
"""

import logging
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
from scipy.optimize import curve_fit
from scipy.stats import norm

__all__ = [
    'AcfResult',
    'DECAY_MODELS',
    'DecayFit',
    'DecayModel',
    'Topology',
    'acf',
    'bartlett_half_width',
    'bond_cosine_correlation',
    'bond_vectors',
    'fit_decay',
    'initial_guess',
]

logger = logging.getLogger(__name__)

Topology = Literal['linear', 'ring']
"""The polymer topologies the bond vectors are defined for."""

FloatArray = npt.NDArray[np.float64]
ModelFunction = Callable[..., FloatArray]


def _check_topology(topology: str) -> None:
    """Raise ``ValueError`` unless `topology` is ``'linear'`` or ``'ring'``."""
    if topology not in ('linear', 'ring'):
        raise ValueError(
            f"topology must be 'linear' or 'ring'; got '{topology}'."
        )


def _fft_length(n: int) -> int:
    """Return the smallest power of two that is at least ``2 n - 1``."""
    return 1 << (2 * n - 1).bit_length()


# ---------------------------------------------------------------------------
# Autocorrelation
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AcfResult:
    """
    The sample autocorrelation of one or more time series.

    Attributes
    ----------
    values : numpy.ndarray
        :math:`r_k` for lags ``0`` to ``nlags``, shape ``(nlags + 1,)`` for a
        one-dimensional series or ``(nlags + 1, n_series)`` for columns.
        ``values[0]`` is 1.
    half_width : numpy.ndarray or None
        The half-width of the confidence interval around each value at
        level ``1 - alpha``, the same shape as `values`; ``half_width[0]``
        is 0. ``None`` when no band was asked for.
    alpha : float or None
        The significance level of the band, or ``None``.

    .. versionadded:: 0.6.0
    """

    values: FloatArray
    half_width: FloatArray | None
    alpha: float | None


def acf(
    x: npt.ArrayLike,
    nlags: int,
    alpha: float | None = 0.05,
) -> AcfResult:
    """
    Estimate the autocorrelation function of a time series by the FFT.

    Parameters
    ----------
    x : array_like
        The time series, shape ``(n_frames,)`` or ``(n_frames, n_series)``,
        with the frames along axis 0 and one series per column; every value
        finite, and no series constant.
    nlags : int
        The largest lag, in frames, from 0 to ``n_frames - 1``.
    alpha : float or None, default 0.05
        The significance level of the confidence band, in ``(0, 1)``; 0.05
        gives a 95 % band. ``None`` skips the band.

    Returns
    -------
    AcfResult
        The autocorrelation at lags 0 to `nlags` and its band.

    Raises
    ------
    ValueError
        If `x` is not one- or two-dimensional or has fewer than two frames,
        holds a non-finite value, or has a constant series; if `nlags` is
        outside ``[0, n_frames - 1]``; or if `alpha` is outside ``(0, 1)``.

    Notes
    -----
    With :math:`\\bar{x}` the mean of the series, the autocovariance is
    :math:`c_k = \\frac{1}{n} \\sum_{t=0}^{n-k-1} (x_t - \\bar{x})
    (x_{t+k} - \\bar{x})` and :math:`r_k = c_k / c_0` (NIST/SEMATECH
    e-Handbook, section 1.3.3.1). The sum over :math:`t` for every
    :math:`k` is one linear correlation, computed as the inverse transform
    of the squared magnitude of the series' Fourier transform after
    zero-padding to at least :math:`2n - 1` points so that no wrap-around
    enters (Wiener-Khinchin).

    The band is Bartlett's (1946) large-sample approximation under the
    hypothesis that the autocorrelation vanishes beyond lag :math:`k - 1`:
    :math:`\\mathrm{var}(r_k) \\approx (1 + 2 \\sum_{j=1}^{k-1} r_j^2) / n`,
    so the band widens with the lag, and the half-width is
    :math:`z_{1-\\alpha/2} \\sqrt{\\mathrm{var}(r_k)}` with :math:`z` the
    standard normal quantile. For a white-noise series it reduces to
    :math:`\\pm z_{1-\\alpha/2} / \\sqrt{n}`. The approximation needs
    :math:`n \\gg` `nlags`.

    Time is :math:`O(m\\, n \\log n)` and space :math:`O(m\\, n)`, with
    :math:`n` the frames and :math:`m` the series.

    Examples
    --------
    The alternating series ``1, -1, 1, -1, ...`` of 8 frames has mean 0,
    :math:`c_k = (-1)^k (8 - k) / 8`, and so :math:`r_1 = -7/8` and
    :math:`r_2 = 6/8`:

    >>> import numpy as np
    >>> from polyphys.analyze.correlations import acf
    >>> x = np.tile([1.0, -1.0], 4)
    >>> result = acf(x, nlags=2)
    >>> np.round(result.values, 6)
    array([ 1.   , -0.875,  0.75 ])

    The 95 % band is 0 at lag 0, :math:`1.96 / \\sqrt{8}` at lag 1, and
    wider at lag 2 because :math:`r_1^2` enters:

    >>> np.round(result.half_width, 4)
    array([0.    , 0.693 , 1.1025])

    .. versionadded:: 0.6.0
    """
    series = np.asarray(x, dtype=np.float64)
    if series.ndim not in (1, 2):
        raise ValueError(
            f'x must be one- or two-dimensional; got {series.ndim}'
            ' dimension(s).'
        )
    n_frames = series.shape[0]
    if n_frames < 2:
        raise ValueError(f'x must have at least two frames; got {n_frames}.')
    if not np.all(np.isfinite(series)):
        raise ValueError('x must be finite.')
    if not 0 <= nlags < n_frames:
        raise ValueError(
            f'nlags must be between 0 and n_frames - 1 = {n_frames - 1};'
            f' got {nlags}.'
        )
    if alpha is not None and not 0 < alpha < 1:
        raise ValueError(f'alpha must be in (0, 1); got {alpha}.')
    if np.any(np.all(series == series[0], axis=0)):
        raise ValueError('x has a constant series; its acf is undefined.')

    centered = series - series.mean(axis=0)
    n_fft = _fft_length(n_frames)
    spectrum = np.fft.rfft(centered, n=n_fft, axis=0)
    autocov = np.fft.irfft(np.abs(spectrum) ** 2, n=n_fft, axis=0)
    autocov = autocov[: nlags + 1] / n_frames
    values: FloatArray = autocov / autocov[0]
    values[0] = 1.0

    half_width: FloatArray | None = None
    if alpha is not None:
        half_width = bartlett_half_width(values, n_frames, alpha)
    return AcfResult(values=values, half_width=half_width, alpha=alpha)


def bartlett_half_width(
    values: npt.ArrayLike,
    n_frames: int,
    alpha: float = 0.05,
) -> FloatArray:
    """
    The half-width of Bartlett's confidence band around an autocorrelation.

    Parameters
    ----------
    values : array_like
        The sample autocorrelation :math:`r_k` for lags 0 to ``nlags``,
        shape ``(nlags + 1,)`` or ``(nlags + 1, n_series)``, as
        :func:`acf` returns it.
    n_frames : int
        The number of frames the autocorrelation was estimated from, at
        least 2.
    alpha : float, default 0.05
        The significance level of the band, in ``(0, 1)``; 0.05 gives a 95 %
        band.

    Returns
    -------
    numpy.ndarray
        The half-width at each lag, the same shape as `values`, ``float64``;
        0 at lag 0. The band is ``values - half_width`` to
        ``values + half_width``.

    Raises
    ------
    ValueError
        If `values` is not one- or two-dimensional, `n_frames` is below 2,
        or `alpha` is outside ``(0, 1)``.

    Notes
    -----
    Bartlett's (1946) large-sample variance under the hypothesis that the
    autocorrelation vanishes beyond lag :math:`k - 1` is
    :math:`\\mathrm{var}(r_k) \\approx (1 + 2 \\sum_{j=1}^{k-1} r_j^2) / n`,
    and the half-width is :math:`z_{1-\\alpha/2} \\sqrt{\\mathrm{var}(r_k)}`.
    The band depends on the stored :math:`r_k`, `n_frames`, and `alpha`
    only, so a pipeline need store the autocorrelation alone and draw the
    band at any level later. Time and space are :math:`O(m\\, k)`, with
    :math:`k` the lags and :math:`m` the series.

    Examples
    --------
    The band of the alternating series of :func:`acf`'s example, rebuilt
    from its autocorrelation:

    >>> import numpy as np
    >>> from polyphys.analyze.correlations import bartlett_half_width
    >>> np.round(bartlett_half_width([1.0, -0.875, 0.75], 8), 4)
    array([0.    , 0.693 , 1.1025])

    .. versionadded:: 0.6.0
    """
    acf_values = np.asarray(values, dtype=np.float64)
    if acf_values.ndim not in (1, 2):
        raise ValueError(
            f'values must be one- or two-dimensional; got {acf_values.ndim}'
            ' dimension(s).'
        )
    if n_frames < 2:
        raise ValueError(f'n_frames must be at least 2; got {n_frames}.')
    if not 0 < alpha < 1:
        raise ValueError(f'alpha must be in (0, 1); got {alpha}.')
    variance = np.full_like(acf_values, 1.0 / n_frames)
    variance[0] = 0.0
    if acf_values.shape[0] > 2:
        variance[2:] *= 1.0 + 2.0 * np.cumsum(acf_values[1:-1] ** 2, axis=0)
    z_score = float(norm.ppf(1.0 - alpha / 2.0))
    half_width: FloatArray = z_score * np.sqrt(variance)
    return half_width


# ---------------------------------------------------------------------------
# Decay fits
# ---------------------------------------------------------------------------


def _exp(x: FloatArray, tau: float) -> FloatArray:
    """:math:`\\exp(-x / \\tau)`."""
    return np.asarray(np.exp(-x / tau), dtype=np.float64)


def _stretched(x: FloatArray, tau: float, beta: float) -> FloatArray:
    """:math:`\\exp[-(x / \\tau)^\\beta]`."""
    return np.asarray(np.exp(-((x / tau) ** beta)), dtype=np.float64)


def _stretched_amp(
    x: FloatArray, tau: float, beta: float, amp: float
) -> FloatArray:
    """:math:`A \\exp[-(x / \\tau)^\\beta]`."""
    return amp * _stretched(x, tau, beta)


def _stretched_amp_offset(
    x: FloatArray, tau: float, beta: float, amp: float, offset: float
) -> FloatArray:
    """:math:`A \\exp[-(x / \\tau)^\\beta] + B`."""
    return amp * _stretched(x, tau, beta) + offset


@dataclass(frozen=True)
class DecayModel:
    """
    A decay model: its function, parameter names, and formula.

    Attributes
    ----------
    function : callable
        ``function(x, *params)``, vectorized over `x`.
    names : tuple of str
        The parameter names in the order `function` takes them, from
        ``'tau'``, ``'beta'``, ``'amp'``, and ``'offset'``.
    formula : str
        The model as a formula, for labels.

    .. versionadded:: 0.6.0
    """

    function: ModelFunction
    names: tuple[str, ...]
    formula: str


DECAY_MODELS: Mapping[str, DecayModel] = {
    'exp': DecayModel(_exp, ('tau',), 'exp(-x/tau)'),
    'stretched': DecayModel(
        _stretched, ('tau', 'beta'), 'exp(-(x/tau)^beta)'
    ),
    'stretched_amp': DecayModel(
        _stretched_amp, ('tau', 'beta', 'amp'), 'amp exp(-(x/tau)^beta)'
    ),
    'stretched_amp_offset': DecayModel(
        _stretched_amp_offset,
        ('tau', 'beta', 'amp', 'offset'),
        'amp exp(-(x/tau)^beta) + offset',
    ),
}
"""
The decay models :func:`fit_decay` accepts, by name.

``'exp'`` is the pure exponential :math:`\\exp(-x/\\tau)`; ``'stretched'``
the Kohlrausch-Williams-Watts form :math:`\\exp[-(x/\\tau)^\\beta]`;
``'stretched_amp'`` scales it by an amplitude :math:`A`; and
``'stretched_amp_offset'`` adds an offset :math:`B`. The parameters are
ordered ``tau``, ``beta``, ``amp``, ``offset``. :math:`\\tau` is in the
unit of ``x`` and positive; :math:`\\beta` is positive, 1 for a pure
exponential.

.. versionadded:: 0.6.0
"""


def _get_model(model: str) -> DecayModel:
    """Return the :class:`DecayModel` named `model` or raise ``ValueError``."""
    try:
        return DECAY_MODELS[model]
    except KeyError:
        names = ', '.join(repr(name) for name in DECAY_MODELS)
        raise ValueError(
            f"model must be one of {names}; got '{model}'."
        ) from None


def _as_curve(
    x: npt.ArrayLike, y: npt.ArrayLike
) -> tuple[FloatArray, FloatArray]:
    """Return `x` and `y` as finite one-dimensional arrays of one length."""
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    if x_arr.ndim != 1 or y_arr.ndim != 1:
        raise ValueError(
            'x and y must be one-dimensional; got'
            f' {x_arr.ndim} and {y_arr.ndim} dimension(s).'
        )
    if x_arr.shape != y_arr.shape:
        raise ValueError(
            f'x and y must have the same length; got {x_arr.size} and'
            f' {y_arr.size}.'
        )
    if not (np.all(np.isfinite(x_arr)) and np.all(np.isfinite(y_arr))):
        raise ValueError('x and y must be finite.')
    return x_arr, y_arr


@dataclass(frozen=True)
class DecayFit:
    """
    The result of fitting a decay model to a curve.

    Attributes
    ----------
    model : str
        The key of the model in :data:`DECAY_MODELS`.
    names : tuple of str
        The parameter names, in the order of `params`.
    params : numpy.ndarray
        The fitted parameters, shape ``(n_params,)``; all NaN when the fit
        did not converge.
    stderr : numpy.ndarray
        The standard error of each parameter, the square root of the
        diagonal of `pcov`; NaN where the covariance could not be estimated
        or the fit did not converge.
    pcov : numpy.ndarray
        The estimated covariance of the parameters, shape
        ``(n_params, n_params)``, as :func:`scipy.optimize.curve_fit`
        returns it; all NaN when the fit did not converge.
    p0 : numpy.ndarray
        The starting point the fit used, shape ``(n_params,)``.
    converged : bool
        Whether :func:`scipy.optimize.curve_fit` returned a solution.

    .. versionadded:: 0.6.0
    """

    model: str
    names: tuple[str, ...]
    params: FloatArray
    stderr: FloatArray
    pcov: FloatArray
    p0: FloatArray
    converged: bool

    def evaluate(self, x: npt.ArrayLike) -> FloatArray:
        """
        Evaluate the fitted model at `x`.

        Parameters
        ----------
        x : array_like
            Points in the unit of the fitted ``x``.

        Returns
        -------
        numpy.ndarray
            The model at `x`, the same shape; all NaN when the fit did not
            converge.

        Examples
        --------
        >>> import numpy as np
        >>> from polyphys.analyze.correlations import fit_decay
        >>> x = np.arange(0.0, 20.0)
        >>> fit = fit_decay(x, np.exp(-x / 4.0), 'exp')
        >>> np.round(fit.evaluate([0.0, 4.0]), 6)
        array([1.      , 0.367879])
        """
        return DECAY_MODELS[self.model].function(
            np.asarray(x, dtype=np.float64), *self.params
        )


def initial_guess(
    x: npt.ArrayLike, y: npt.ArrayLike, model: str
) -> FloatArray:
    """
    Read a starting point for each parameter of a decay model off the data.

    Parameters
    ----------
    x : array_like
        The abscissa, shape ``(n_points,)``, finite, increasing, and
        non-negative, in the caller's unit (lags in frames, or time).
    y : array_like
        The curve at `x`, the same shape, finite.
    model : str
        A key of :data:`DECAY_MODELS`.

    Returns
    -------
    numpy.ndarray
        One value per parameter of the model, in its order.

    Raises
    ------
    ValueError
        If `model` is unknown, or `x` and `y` are not finite
        one-dimensional arrays of one length.

    Notes
    -----
    The offset :math:`B` is the mean of the last tenth of `y` (at least one
    point); the amplitude :math:`A` is ``y[0] - B``, or 1 when that is 0;
    the decay time :math:`\\tau` is the first `x` at which
    :math:`(y - B) / A` is at or below :math:`1/e`, which is exact for a
    pure exponential sampled at :math:`\\tau`, or the last `x` when the
    curve never falls that far; and :math:`\\beta` is 1. A model without an
    offset takes :math:`B = 0`; one without an amplitude takes
    :math:`A = 1`. A :math:`\\tau` of 0 (the curve is already at
    :math:`1/e` at ``x[0] = 0``) is replaced by the first positive `x`, or
    1.

    Time and space are :math:`O(n_\\mathrm{points})`.

    Examples
    --------
    A pure exponential of decay time 4 crosses :math:`1/e` at ``x = 4``:

    >>> import numpy as np
    >>> from polyphys.analyze.correlations import initial_guess
    >>> x = np.arange(0.0, 20.0)
    >>> initial_guess(x, np.exp(-x / 4.0), 'exp')
    array([4.])

    The offset is read off the tail, so it is a starting point, not the
    answer:

    >>> y = 2.0 * np.exp(-x / 4.0) + 0.5
    >>> np.round(initial_guess(x, y, 'stretched_amp_offset'), 2)
    array([4.  , 1.  , 1.98, 0.52])

    .. versionadded:: 0.6.0
    """
    spec = _get_model(model)
    x_arr, y_arr = _as_curve(x, y)
    n_points = x_arr.size
    if n_points == 0:
        raise ValueError('x and y must not be empty.')
    offset = 0.0
    if 'offset' in spec.names:
        offset = float(np.mean(y_arr[-max(1, n_points // 10):]))
    amp = 1.0
    if 'amp' in spec.names:
        amp = float(y_arr[0] - offset)
        if amp == 0.0:
            amp = 1.0
    decayed = (y_arr - offset) / amp
    below = np.flatnonzero(decayed <= np.exp(-1.0))
    tau = float(x_arr[below[0]]) if below.size else float(x_arr[-1])
    if tau <= 0.0:
        positive = np.flatnonzero(x_arr > 0.0)
        tau = float(x_arr[positive[0]]) if positive.size else 1.0
    guesses = {'tau': tau, 'beta': 1.0, 'amp': amp, 'offset': offset}
    return np.array([guesses[name] for name in spec.names], dtype=np.float64)


def fit_decay(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    model: str,
    p0: npt.ArrayLike | None = None,
    **kwargs: Any,
) -> DecayFit:
    """
    Fit a decay model to a curve by nonlinear least squares.

    Parameters
    ----------
    x : array_like
        The abscissa, shape ``(n_points,)``, finite, in the caller's unit
        (lags in frames, or time). Fit only the range meant for the model:
        on an autocorrelation, the lags the model describes; on a ring's
        bond correlation, lags up to ``n_bonds // 2``.
    y : array_like
        The curve at `x`, the same shape, finite.
    model : str
        A key of :data:`DECAY_MODELS`.
    p0 : array_like or None, default None
        The starting point, one value per parameter of the model, in its
        order. ``None`` takes :func:`initial_guess`.
    **kwargs
        Passed to :func:`scipy.optimize.curve_fit`: ``bounds``, ``sigma``,
        ``absolute_sigma``, ``maxfev``, ``method``, and the rest.

    Returns
    -------
    DecayFit
        The parameters, their standard errors, the covariance, the starting
        point used, and whether the fit converged.

    Raises
    ------
    ValueError
        If `model` is unknown; `x` and `y` are not finite one-dimensional
        arrays of one length; there are fewer points than parameters; or
        `p0` has the wrong length. :func:`scipy.optimize.curve_fit` also
        raises :exc:`ValueError` for an invalid option in `kwargs`, such as
        an unknown ``method`` or a `p0` outside ``bounds``.

    Notes
    -----
    :func:`scipy.optimize.curve_fit` minimizes the sum of squared residuals
    from `p0`, by default with the Levenberg-Marquardt method (``'lm'``),
    or with ``'trf'`` when ``bounds`` are given, as its documentation
    states; ``method`` in `kwargs` chooses another. A nonlinear model can
    have several local minima, so `p0` matters, and the one used is
    recorded.
    A fit fails when the model is not finite at the start, when the
    optimizer stops without a solution, or when it does not move from the
    start and the covariance cannot be estimated, which is how the
    ``'lm'`` method signals a Jacobian that is zero or rank-deficient (the
    :func:`scipy.optimize.curve_fit` documentation, under ``pcov``). A
    failed fit is reported with ``converged=False``, NaN parameters, and a
    warning in the :mod:`logging` log, not an exception; invalid inputs or
    options still raise (see Raises). When the fit
    moves but the covariance cannot be estimated (SciPy emits an
    ``OptimizeWarning``), it is converged and the standard errors of the
    affected parameters are NaN. A start at the exact solution of
    noise-free data does not move either, but its covariance is finite, so
    it is converged. The standard
    error is the square root of the covariance diagonal, as in SciPy's
    documentation.

    The models are the stretched exponential of Williams and Watts (1970),
    Trans. Faraday Soc. 66, 80, https://doi.org/10.1039/TF9706600080, with
    :math:`\\beta = 1` for the pure exponential.

    Time is :math:`O(n_\\mathrm{iter}\\, n_\\mathrm{points}\\,
    n_\\mathrm{params})` and space :math:`O(n_\\mathrm{points}\\,
    n_\\mathrm{params})`.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.correlations import fit_decay
    >>> x = np.arange(0.0, 40.0)
    >>> y = 2.0 * np.exp(-((x / 5.0) ** 0.8)) + 0.3
    >>> fit = fit_decay(x, y, 'stretched_amp_offset')
    >>> fit.converged
    True
    >>> fit.names
    ('tau', 'beta', 'amp', 'offset')
    >>> np.round(fit.params, 6)
    array([5. , 0.8, 2. , 0.3])

    .. versionadded:: 0.6.0
    .. versionchanged:: 0.9.0
       A start where the model is not finite, and a fit that does not move
       and has no finite covariance, are failed fits; ``method='trf'`` no
       longer raises on a non-finite start.
    """
    spec = _get_model(model)
    x_arr, y_arr = _as_curve(x, y)
    n_params = len(spec.names)
    if x_arr.size < n_params:
        raise ValueError(
            f"model '{model}' has {n_params} parameter(s) but x has only"
            f' {x_arr.size} point(s).'
        )
    if p0 is None:
        start = initial_guess(x_arr, y_arr, model)
    else:
        start = np.asarray(p0, dtype=np.float64)
        if start.shape != (n_params,):
            raise ValueError(
                f"p0 must have one value per parameter of '{model}',"
                f' {spec.names}; got shape {start.shape}.'
            )
    nan_params = np.full(n_params, np.nan)
    failed = DecayFit(
        model=model,
        names=spec.names,
        params=nan_params,
        stderr=nan_params.copy(),
        pcov=np.full((n_params, n_params), np.nan),
        p0=start,
        converged=False,
    )
    # A start where the model is not finite cannot be fitted: 'lm' stays
    # there and 'trf' raises, so it is reported as a failed fit up front.
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        at_start = spec.function(x_arr, *start)
    if not np.all(np.isfinite(at_start)):
        logger.warning(
            "fit_decay: '%s' is not finite at the start %s; the fit failed.",
            model,
            start,
        )
        return failed
    try:
        params, pcov = curve_fit(
            spec.function, x_arr, y_arr, p0=start, **kwargs
        )
    except RuntimeError as error:
        logger.warning("fit_decay: '%s' did not converge: %s", model, error)
        return failed
    params = np.asarray(params, dtype=np.float64)
    pcov = np.asarray(pcov, dtype=np.float64)
    diagonal = np.diag(pcov)
    if np.array_equal(params, start) and not np.all(np.isfinite(pcov)):
        logger.warning(
            "fit_decay: '%s' did not move from the start %s and its"
            ' covariance could not be estimated; the fit failed.',
            model,
            start,
        )
        return failed
    stderr = np.where(np.isfinite(diagonal), np.sqrt(np.abs(diagonal)), np.nan)
    if not np.all(np.isfinite(diagonal)):
        logger.warning(
            "fit_decay: '%s' converged but its covariance could not be"
            ' estimated.',
            model,
        )
    return DecayFit(
        model=model,
        names=spec.names,
        params=params,
        stderr=np.asarray(stderr, dtype=np.float64),
        pcov=pcov,
        p0=start,
        converged=True,
    )


# ---------------------------------------------------------------------------
# Bonds
# ---------------------------------------------------------------------------


def bond_vectors(
    positions: npt.ArrayLike, topology: Topology
) -> FloatArray:
    """
    Difference the positions of a chain into its bond vectors.

    Parameters
    ----------
    positions : array_like
        The monomer coordinates in chain order, shape
        ``(..., n_atoms, n_dims)``, unwrapped, in the caller's length unit;
        ``n_atoms`` at least 2 on a linear chain and 3 on a ring.
    topology : {'linear', 'ring'}
        Whether the last monomer is bonded back to the first.

    Returns
    -------
    numpy.ndarray
        The bond vectors :math:`\\mathbf{r}_{i+1} - \\mathbf{r}_i`, shape
        ``(..., n_bonds, n_dims)``, with ``n_bonds = n_atoms - 1`` on a
        linear chain and ``n_atoms`` on a ring, whose last bond is
        :math:`\\mathbf{r}_0 - \\mathbf{r}_{n-1}`. Bond lengths are
        ``numpy.linalg.norm(bonds, axis=-1)``.

    Raises
    ------
    ValueError
        If `topology` is unknown, `positions` has fewer than two dimensions,
        or the chain is too short for its topology.

    Notes
    -----
    Positions must be unwrapped: a periodic image jump inside a bond gives
    a bond of the box length. Time and space are
    :math:`O(F\\, n_\\mathrm{atoms}\\, n_\\mathrm{dims})`.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.correlations import bond_vectors
    >>> square = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    >>> bond_vectors(square, 'linear')
    array([[ 1.,  0.],
           [ 0.,  1.],
           [-1.,  0.]])
    >>> bond_vectors(square, 'ring')[-1]
    array([ 0., -1.])

    .. versionadded:: 0.6.0
    """
    _check_topology(topology)
    coords = np.asarray(positions, dtype=np.float64)
    if coords.ndim < 2:
        raise ValueError(
            'positions must have shape (..., n_atoms, n_dims); got'
            f' {coords.ndim} dimension(s).'
        )
    n_atoms = coords.shape[-2]
    minimum = 2 if topology == 'linear' else 3
    if n_atoms < minimum:
        raise ValueError(
            f"a '{topology}' chain needs at least {minimum} atoms; got"
            f' {n_atoms}.'
        )
    if topology == 'linear':
        return np.asarray(coords[..., 1:, :] - coords[..., :-1, :])
    return np.asarray(np.roll(coords, -1, axis=-2) - coords)


def bond_cosine_correlation(
    positions: npt.ArrayLike, topology: Topology
) -> FloatArray:
    """
    Average the cosine between bond vectors at every lag along a chain.

    Parameters
    ----------
    positions : array_like
        The monomer coordinates in chain order, shape
        ``(..., n_atoms, n_dims)``, unwrapped; no two bonded monomers
        coincide.
    topology : {'linear', 'ring'}
        Whether the last monomer is bonded back to the first.

    Returns
    -------
    numpy.ndarray
        :math:`C(k)` for lags ``0`` to ``n_bonds - 1``, shape
        ``(..., n_bonds)``, dimensionless in ``[-1, 1]`` with ``C(0) = 1``.
        Leading axes are kept, so the average over frames is
        ``result.mean(axis=0)`` for a ``(n_frames, n_atoms, n_dims)``
        input.

    Raises
    ------
    ValueError
        If `topology` is unknown, `positions` has a bad shape or a chain too
        short for its topology, or a bond has zero length.

    Notes
    -----
    With :math:`\\hat{\\mathbf{u}}_i` the unit bond vectors, on a linear
    chain :math:`C(k) = \\frac{1}{N - k} \\sum_{i=0}^{N-k-1}
    \\hat{\\mathbf{u}}_i \\cdot \\hat{\\mathbf{u}}_{i+k}` over the
    :math:`N - k` pairs at lag :math:`k`; on a ring :math:`C(k) =
    \\frac{1}{N} \\sum_{i=0}^{N-1} \\hat{\\mathbf{u}}_i \\cdot
    \\hat{\\mathbf{u}}_{(i+k) \\bmod N}` over all :math:`N` pairs, so that
    the result does not depend on which monomer is labeled first and
    :math:`C(k) = C(N - k)`. The sums over :math:`i` are linear (zero-padded)
    or circular autocorrelations, computed by the FFT of each Cartesian
    component and summed over components.

    For a worm-like chain :math:`C(k) = \\exp(-k\\, b / \\ell_p)` with
    :math:`b` the bond length and :math:`\\ell_p` the persistence length
    (Rubinstein and Colby 2003, section 2.3), so fitting ``'exp'`` with
    :func:`fit_decay` on lags up to ``n_bonds // 2`` for a ring gives
    :math:`\\tau = \\ell_p / b`, and :math:`\\ell_p` is :math:`\\tau` times
    the mean bond length. On a linear chain this is the quantity and the
    fit of :class:`MDAnalysis.analysis.polymer.PersistenceLength`
    (``results.bond_autocorrelation``, ``results.lb``, ``results.lp``),
    which averages over chains and frames and does not treat rings.

    Time is :math:`O(F\\, N \\log N)` and space :math:`O(F\\, N)` for
    :math:`F` frames of :math:`N` bonds, against the :math:`O(F\\, N^2)`
    of the pairwise matrix.

    Examples
    --------
    The bonds of a square ring turn by a right angle each, so
    :math:`C(k) = \\cos(\\pi k / 2)`:

    >>> import numpy as np
    >>> from polyphys.analyze.correlations import bond_cosine_correlation
    >>> square = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    >>> np.round(bond_cosine_correlation(square, 'ring'), 12)
    array([ 1.,  0., -1.,  0.])

    A straight rod has every bond parallel:

    >>> rod = np.outer(np.arange(4.0), [1.0, 2.0, 2.0])
    >>> np.round(bond_cosine_correlation(rod, 'linear'), 12)
    array([1., 1., 1.])

    .. versionadded:: 0.6.0
    """
    bonds = bond_vectors(positions, topology)
    lengths = np.linalg.norm(bonds, axis=-1)
    if np.any(lengths == 0.0):
        raise ValueError('every bond must have a positive length.')
    units = bonds / lengths[..., np.newaxis]
    # Put the bond axis last so the FFT runs along it, per component.
    components = np.swapaxes(units, -1, -2)
    n_bonds = components.shape[-1]
    if topology == 'linear':
        n_fft = _fft_length(n_bonds)
        pairs: FloatArray = np.arange(n_bonds, 0, -1, dtype=np.float64)
    else:
        n_fft = n_bonds
        pairs = np.full(n_bonds, float(n_bonds))
    spectrum = np.fft.rfft(components, n=n_fft, axis=-1)
    power = np.sum(np.abs(spectrum) ** 2, axis=-2)
    sums = np.fft.irfft(power, n=n_fft, axis=-1)[..., :n_bonds]
    return np.asarray(sums / pairs, dtype=np.float64)
