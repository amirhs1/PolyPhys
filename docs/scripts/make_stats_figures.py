"""
Regenerate the two figures of ``docs/source/statistics.md``.

Both figures come from seeded synthetic AR(1) series, whose statistical
inefficiency is known analytically, so the plots show the blocking table
against the exact answer:

- ``stats_ar1_table.png``: AR(1) with phi = 0.9 and 2**16 frames, a run long
  enough for a clear plateau.
- ``stats_too_short.png``: AR(1) with phi = 0.99 and 2**14 frames, about 82
  independent samples, where the table runs out of blocks before the plateau.

Run from the repository root, with PolyPhys importable::

    python docs/scripts/make_stats_figures.py

matplotlib is needed only to regenerate these figures; it is not a PolyPhys
dependency, and nothing in the package or the documentation build imports
this script.
"""

import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from polyphys.analyze import stats

matplotlib.use("Agg")

OUT = Path(__file__).resolve().parents[1] / "source" / "_static"


def ar1(phi, n_frames, seed):
    """Return a stationary AR(1) series with unit-variance innovations."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(n_frames)
    x = np.empty(n_frames)
    x[0] = noise[0] / np.sqrt(1.0 - phi ** 2)
    for t in range(1, n_frames):
        x[t] = phi * x[t - 1] + noise[t]
    return x


def ar1_inefficiency(phi, block_size):
    """Return the exact AR(1) inefficiency at each integer block size.

    ``s(B) = 1 + 2 sum_{t=1}^{B-1} (1 - t/B) phi**t``; the derivation is in
    the module docstring of ``polyphys/tests/analyze/test_stats.py`` (R2).
    """
    out = []
    for size in block_size:
        t = np.arange(1, size)
        out.append(1.0 + 2.0 * np.sum((1.0 - t / size) * phi ** t))
    return np.array(out)


def plot(phi, n_frames, seed, filename, title):
    """Plot the SEM and inefficiency tables of one AR(1) series."""
    x = ar1(phi, n_frames, seed)
    table = stats.block_average(x)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lee = stats.select_plateau(table, method="lee")
        sem_err = stats.select_plateau(table, method="sem_err")
    b = table.block_size
    s_inf = (1.0 + phi) / (1.0 - phi)
    se_true = np.sqrt(s_inf / (1.0 - phi ** 2) / n_frames)
    few = table.n_blocks <= stats.LEE_MIN_BLOCK_RATIO

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
    for ax in (ax1, ax2):
        if few.any():
            ax.axvspan(b[few][0] / 1.4, b[-1] * 1.4, color="C3", alpha=0.08,
                       label="50 blocks or fewer")
        ax.set_xscale("log", base=2)
        ax.set_xlabel("block size $B$ (frames)")

    ax1.errorbar(b, table.sem, yerr=table.sem_err, fmt="o-", ms=4,
                 capsize=3, color="C0", label="blocked SEM ± sem_err")
    ax1.axhline(table.sem[0], color="grey", ls=":", label="naive SEM")
    ax1.axhline(se_true, color="k", ls="--", lw=1, label="exact SE")
    ax1.plot(b[lee], table.sem[lee], "^", ms=13, mfc="none", mew=2,
             color="C2", label=f"method='lee', $B$ = {b[lee]}")
    ax1.plot(b[sem_err], table.sem[sem_err], "s", ms=13, mfc="none", mew=2,
             color="C1", label=f"method='sem_err', $B$ = {b[sem_err]}")
    ax1.set_ylabel("standard error of the mean")
    ax1.legend(fontsize=7)

    fine = np.unique(np.round(2 ** np.linspace(0, np.log2(b[-1]), 200)))
    ax2.errorbar(b, table.statistical_inefficiency,
                 yerr=table.inefficiency_err, fmt="o-", ms=4, capsize=3,
                 color="C0", label="estimated $s(B)$ ± err")
    ax2.plot(fine, ar1_inefficiency(phi, fine.astype(int)), "k-", lw=1,
             label="exact $s(B)$")
    ax2.axhline(s_inf, color="k", ls="--", lw=1,
                label=f"exact $s$ = {s_inf:.0f}")
    ax2.plot(b[lee], table.statistical_inefficiency[lee], "^", ms=13,
             mfc="none", mew=2, color="C2")
    ax2.plot(b[sem_err], table.statistical_inefficiency[sem_err], "s",
             ms=13, mfc="none", mew=2, color="C1")
    ax2.set_ylabel("statistical inefficiency $s$")
    ax2.legend(fontsize=7)

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=110)
    plt.close(fig)
    print(f"{filename}: exact SE {se_true:.4g}; lee B={b[lee]} "
          f"SE={table.sem[lee]:.4g} "
          f"s={table.statistical_inefficiency[lee]:.4g};"
          f" sem_err B={b[sem_err]} SE={table.sem[sem_err]:.4g} "
          f"s={table.statistical_inefficiency[sem_err]:.4g}")


if __name__ == "__main__":
    plot(0.9, 2 ** 16, 7, "stats_ar1_table.png",
         r"AR(1), $\phi$ = 0.9, $2^{16}$ frames: a plateau is reached")
    plot(0.99, 2 ** 14, 22, "stats_too_short.png",
         r"AR(1), $\phi$ = 0.99, $2^{14}$ frames: the run is too short")
