# Correlations: autocorrelation, decay fits, and bond cosines

This page explains what {mod}`polyphys.analyze.correlations` computes: the
autocorrelation function of a time series with its confidence band, the
decay model fitted to such a curve, and the cosine correlation between the
bond vectors of a chain, from which a persistence length follows. Every
routine is a pure numpy or scipy kernel on arrays; reading time series or
frames, choosing the lags to fit, and averaging over frames are the caller's
job.

## Notation

$x_t$
: A time series of $n$ equally spaced frames, `x`, shape `(n_frames,)` or
  `(n_frames, n_series)` with one independent series per column.

$k$
: A lag, in frames: the number of frames between the two values compared.
  `nlags` is the largest lag reported.

$r_k$
: The sample autocorrelation at lag $k$, `acf`: the autocovariance of the
  series with itself $k$ frames later, divided by its variance; $r_0 = 1$.

$\tau$, $\beta$, $A$, $B$
: The decay time (in the unit of the fitted $x$), the stretching exponent,
  the amplitude, and the offset of a decay model.

$\mathbf{b}_i$, $\hat{\mathbf{u}}_i$
: The bond vector $\mathbf{r}_{i+1} - \mathbf{r}_i$ of a chain with
  monomers in chain order, and its unit vector. A ring adds the closing
  bond $\mathbf{r}_0 - \mathbf{r}_{N-1}$, so a chain of $n$ monomers has
  $N = n - 1$ bonds when linear and $N = n$ when a ring.

$C(k)$
: The bond cosine correlation, the mean of
  $\hat{\mathbf{u}}_i \cdot \hat{\mathbf{u}}_{i+k}$ over the pairs of bonds
  $k$ bonds apart in one frame; $C(0) = 1$.

## The autocorrelation and its band

{func}`~polyphys.analyze.correlations.acf` estimates $r_k$ for
$k = 0, \dots,$ `nlags`:

```python
import numpy as np
from polyphys.analyze import correlations

result = correlations.acf(gyr_t, nlags=2000, alpha=0.05)
result.values       # r_k, shape (2001,)
result.half_width   # the 95 % band is values ± half_width
```

With $\bar{x}$ the mean of the series, the autocovariance is
$c_k = \frac{1}{n} \sum_{t=0}^{n-k-1} (x_t - \bar{x})(x_{t+k} - \bar{x})$
and $r_k = c_k / c_0$ (NIST/SEMATECH e-Handbook, section 1.3.3.1). The
$1/n$ estimator keeps the sequence positive semi-definite and matches the
default of `statsmodels.tsa.stattools.acf`, which the research code used,
so stored autocorrelations are unchanged; the port agrees with statsmodels
0.15.0 to $3 \times 10^{-16}$. The sums over $t$ for every $k$ are one
linear correlation, computed by the fast Fourier transform after zero
padding to at least $2n - 1$ points so that no wrap-around enters:
$O(n \log n)$ time per series instead of $O(n \cdot$ `nlags`$)$.

The band is Bartlett's (1946) large-sample approximation under the
hypothesis that the autocorrelation vanishes beyond lag $k - 1$,
$\mathrm{var}(r_k) \approx (1 + 2 \sum_{j=1}^{k-1} r_j^2) / n$, with
$\mathrm{var}(r_1) = 1/n$ and $\mathrm{var}(r_0) = 0$, and the half-width is
$z_{1-\alpha/2} \sqrt{\mathrm{var}(r_k)}$. It widens with the lag, reduces to
$\pm z_{1-\alpha/2} / \sqrt{n}$ for white noise, and needs $n \gg$ `nlags`.
This is the statsmodels convention, so a band computed here matches one
computed there; the research `acf_generator` stored $-$`half_width` and
$+$`half_width` as its `acfLowerCi` and `acfUpperCi` columns.

For the uncertainty of a *time average* over correlated frames, which is a
different question, see [Statistics](statistics.md).

## Decay fits

{func}`~polyphys.analyze.correlations.fit_decay` fits one of the models in
{data}`~polyphys.analyze.correlations.DECAY_MODELS` by
{func}`scipy.optimize.curve_fit`:

| Model                  | Formula                            | Parameters                     |
| ---------------------- | ---------------------------------- | ------------------------------ |
| `exp`                  | $\exp(-x/\tau)$                    | `tau`                          |
| `stretched`            | $\exp[-(x/\tau)^\beta]$            | `tau`, `beta`                  |
| `stretched_amp`        | $A \exp[-(x/\tau)^\beta]$          | `tau`, `beta`, `amp`           |
| `stretched_amp_offset` | $A \exp[-(x/\tau)^\beta] + B$      | `tau`, `beta`, `amp`, `offset` |

The stretched form is the Kohlrausch-Williams-Watts function (Williams and
Watts, 1970); $\beta = 1$ is the pure exponential.

```python
lags = np.arange(0, 501)
fit = correlations.fit_decay(lags, result.values[:501], 'stretched_amp_offset')
fit.converged           # False when the fit failed; see below
fit.names, fit.params   # ('tau', 'beta', 'amp', 'offset'), the fitted values
fit.stderr              # sqrt of the covariance diagonal
fit.p0                  # the start the fit used
fit.evaluate(lags)      # the fitted curve
```

Conventions to note:

- Models are written in the decay time $\tau$, in the unit of $x$. The
  research `mono_*` functions came in two parametrizations, $\tau$ and the
  rate $\omega = \tau^{-\beta}$; a stored $\omega$ converts as
  $\tau = \omega^{-1/\beta}$.
- A nonlinear fit starts from a point `p0` and can stop in a wrong minimum
  or not at all. The research `fit_wholes` passed no start, so SciPy began
  at all ones, and `fit_exp_wholes` built one for a single model. Here
  `p0=None` takes {func}`~polyphys.analyze.correlations.initial_guess`,
  which reads the start off the data: $B$ from the mean of the last tenth
  of the curve, $A = y_0 - B$, $\tau$ as the first $x$ at which
  $(y - B)/A$ is at or below $1/e$, and $\beta = 1$. The start used is
  recorded in `fit.p0`, and a caller can pass its own.
- A failed fit never raises: a start where the model is not finite, an
  optimizer that gives up, or a fit that does not move from its start and
  has no finite covariance gives `fit.converged` false,
  the parameters are NaN with the model's length, and a warning goes to the
  `polyphys.analyze.correlations` logger. The research drivers padded four
  zeros whatever the model.
- `fit.stderr` is the square root of the covariance diagonal. The research
  drivers stored the diagonal itself, the variance, under a `-std` name.
- Which model, and which lags, describe an observable are the study's
  choices, not the library's: `fit_decay` fits exactly the `x` it is given.
  Note that lag $k$ of an autocorrelation sits at `x = k`, not `k + 1`.

## Bond cosines and the persistence length

{func}`~polyphys.analyze.correlations.bond_cosine_correlation` averages the
cosine between bonds $k$ bonds apart along a chain, for a single frame or
an array of frames:

```python
C = correlations.bond_cosine_correlation(positions, 'ring')  # (n_frames, N)
C_mean = C.mean(axis=0)
lb = np.linalg.norm(correlations.bond_vectors(positions, 'ring'), axis=-1).mean()
fit = correlations.fit_decay(np.arange(N // 2 + 1), C_mean[: N // 2 + 1], 'exp')
lp = fit.params[0] * lb   # persistence length, in the length unit of positions
```

Positions must be unwrapped: no periodic image may break a bond. On a
linear chain lag $k$ has $N - k$ pairs and
$C(k) = \frac{1}{N-k} \sum_{i=0}^{N-k-1} \hat{\mathbf{u}}_i \cdot
\hat{\mathbf{u}}_{i+k}$. On a ring every lag has $N$ pairs, going round the
ring, $C(k) = \frac{1}{N} \sum_{i=0}^{N-1} \hat{\mathbf{u}}_i \cdot
\hat{\mathbf{u}}_{(i+k) \bmod N}$, so that the result does not depend on
which monomer is labeled first and $C(k) = C(N - k)$: lags $0$ to $N/2$
carry the information, and a fit stops there. The sums are linear or
circular autocorrelations by the FFT of each Cartesian component,
$O(N \log N)$ per frame against the $O(N^2)$ pairwise matrix of the
research `bond_info`.

For a worm-like chain $C(k) = \exp(-k\, b / \ell_p)$ with $b$ the bond
length and $\ell_p$ the persistence length (Rubinstein and Colby, 2003,
section 2.3), so the `exp` fit gives $\tau = \ell_p / b$ in bonds. On a
linear chain this is the quantity and the fit of
{class}`MDAnalysis.analysis.polymer.PersistenceLength`, whose `results.lp`
is in length units; the Phase 4 probe layer checks the two against each
other.

Two differences from the research code to keep in mind when comparing with
stored `bondCosineCorrVecMon` files:

- The research probes summed cosines per lag and divided by the number of
  *open* pairs, $N - k$, on rings as well as chains, so their ring values
  depend on the label origin and rest on one pair at the last lag. The
  cyclic average here has the same expectation for an equilibrated ring but
  differs beyond lag 0 by the pairs that cross the origin.
- `bond_info` raised `IndexError` on every linear chain; only ring probes
  ever ran it. {func}`~polyphys.analyze.correlations.bond_vectors` and the
  cosine kernel take both topologies.

The research `acf_generator`, `fit_wholes`, `fit_exp_wholes`, and
`set_x_property`, which read files, called parsers, and wrote CSV reports,
are not ported: that work belongs to the analysis pipeline.

## References

- M. S. Bartlett, "On the theoretical specification and sampling properties
  of autocorrelated time-series", Supplement to the Journal of the Royal
  Statistical Society 8, 27-41 (1946). <https://doi.org/10.2307/2983611>
- NIST/SEMATECH e-Handbook of Statistical Methods, section 1.3.3.1,
  "Autocorrelation Plot".
  <https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm>
- G. Williams and D. C. Watts, "Non-symmetrical dielectric relaxation
  behavior arising from a simple empirical decay function", Trans. Faraday
  Soc. 66, 80-85 (1970). <https://doi.org/10.1039/TF9706600080>
- SciPy, {func}`scipy.optimize.curve_fit`.
- M. Rubinstein and R. H. Colby, Polymer Physics (Oxford University Press,
  2003), Chapter 2.
- MDAnalysis, {class}`MDAnalysis.analysis.polymer.PersistenceLength`, and
  the User Guide page "Determining the persistence length of a polymer".
  <https://userguide.mdanalysis.org/stable/examples/analysis/polymers_and_membranes/polymer.html>
