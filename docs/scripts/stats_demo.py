"""
Demonstrate polyphys.analyze.stats on synthetic series with known answers.

Seven seeded series (i.i.d. noise; AR(1) with phi = 0.5, 0.9, -0.5, and a
too-short run with phi = 0.99; a block-constant series; a linear drift) and
one three-column call are blocked, both plateau rules are applied, and the
results are compared with the exact inefficiency. Writes, into the output
directory:

- <case>.csv       the full blocking table plus the exact s(B) and the pull
- <case>.png       SEM and inefficiency against block size
- overview.png     the inefficiency curves of all cases on one page
- report.md        one summary row per case, the selected levels, warnings

Run with PolyPhys importable, choosing an output directory outside the
repository::

    python docs/scripts/stats_demo.py --out /path/to/stats-demo | tee log.txt

matplotlib is needed only for this script; it is not a PolyPhys dependency.
"""

import argparse
import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from polyphys.analyze import stats
from polyphys.analyze.measurer import sem

matplotlib.use("Agg")

_parser = argparse.ArgumentParser(description="Demonstrate block averaging.")
_parser.add_argument("--out", type=Path, required=True,
                     help="directory for the tables, plots, and report")
OUT = _parser.parse_args().out
OUT.mkdir(parents=True, exist_ok=True)
print("polyphys.analyze.stats loaded from", stats.__file__)
print("numpy", np.__version__, "matplotlib", matplotlib.__version__)


# --- synthetic inputs with analytic answers ----------------------------------

def ar1(phi, n, seed):
    rng = np.random.default_rng(seed)
    e = rng.standard_normal(n)
    x = np.empty(n)
    x[0] = e[0] / np.sqrt(1 - phi ** 2)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x


def ar1_s(phi, block_size):
    """Inefficiency of AR(1) at block size B (R2 in the test module)."""
    out = []
    for b in block_size:
        t = np.arange(1, b)
        out.append(1 + 2 * np.sum((1 - t / b) * phi ** t))
    return np.array(out)


def block_constant(period, n, seed):
    rng = np.random.default_rng(seed)
    return np.repeat(rng.standard_normal(n // period), period)


def block_constant_s(period, n, block_size):
    """Inefficiency of a block-constant series (R3 in the test module)."""
    eff = np.minimum(block_size, period)
    return eff * (n - 1) / (n - eff)


def case(name, title, x, s_of_b=None, s_inf=None, var_x=None):
    return dict(name=name, title=title, x=x, s_of_b=s_of_b, s_inf=s_inf,
                var_x=var_x)


N16, N14 = 2 ** 16, 2 ** 14
cases = [
    case("01_iid", "i.i.d. N(0,1), n = 2^16",
         np.random.default_rng(4).standard_normal(N16),
         lambda b: np.ones(len(b)), 1.0, 1.0),
    case("02_ar1_phi0.5", "AR(1), phi = 0.5, n = 2^16", ar1(0.5, N16, 6),
         lambda b: ar1_s(0.5, b), 3.0, 1 / (1 - 0.25)),
    case("03_ar1_phi0.9", "AR(1), phi = 0.9, n = 2^16", ar1(0.9, N16, 7),
         lambda b: ar1_s(0.9, b), 19.0, 1 / (1 - 0.81)),
    case("04_ar1_phi-0.5", "AR(1), phi = -0.5 (anticorrelated), n = 2^16",
         ar1(-0.5, N16, 5), lambda b: ar1_s(-0.5, b), 1 / 3, 1 / (1 - 0.25)),
    case("05_block_constant_L8", "block-constant, L = 8, n = 2^14",
         block_constant(8, N14, 10),
         lambda b: block_constant_s(8, N14, b), 8 * (N14 - 1) / (N14 - 8),
         None),
    case("06_ar1_phi0.99_short", "AR(1), phi = 0.99, n = 2^14 (too short?)",
         ar1(0.99, N14, 22), lambda b: ar1_s(0.99, b), 199.0,
         1 / (1 - 0.99 ** 2)),
    case("07_drift", "linear drift + noise, n = 2^14 (not stationary)",
         np.linspace(0, 5, N14)
         + np.random.default_rng(23).standard_normal(N14),
         None, None, None),
]

report_rows = []
overview = []


def select(table, method):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        level = stats.select_plateau(table, method=method)
    return level, [str(w.message) for w in caught]


for c in cases:
    x = c["x"]
    n = len(x)
    table = stats.block_average(x)
    b = table.block_size
    s_analytic = c["s_of_b"](b) if c["s_of_b"] else np.full(len(b), np.nan)
    pull = ((table.statistical_inefficiency - s_analytic)
            / table.inefficiency_err)

    # Table to CSV.
    header = ("level,block_size,n_blocks,mean,sem,sem_err,"
              "inefficiency,inefficiency_err,inefficiency_analytic,pull_sigma")
    rows = np.column_stack([
        np.arange(table.n_levels), b, table.n_blocks, table.mean, table.sem,
        table.sem_err, table.statistical_inefficiency,
        table.inefficiency_err, s_analytic, pull,
    ])
    np.savetxt(OUT / f"{c['name']}.csv", rows, delimiter=",", header=header,
               comments="", fmt=["%d", "%d", "%d"] + ["%.6g"] * 7)

    lee, lee_w = select(table, "lee")
    sme, sme_w = select(table, "sem_err")
    se_true = (np.sqrt(c["s_inf"] * c["var_x"] / n)
               if c["var_x"] is not None else np.nan)

    print(f"\n=== {c['name']}: {c['title']}")
    print(f"{'lvl':>3} {'B':>6} {'blocks':>7} {'sem':>10} {'sem_err':>9} "
          f"{'s':>8} {'s_err':>8} {'s_exact':>8} {'pull':>6}")
    for k in range(table.n_levels):
        mark = ("  <- lee" if k == lee else "") + (
            "  <- sem_err" if k == sme else "")
        print(f"{k:3d} {b[k]:6d} {table.n_blocks[k]:7d} "
              f"{table.sem[k]:10.4g} {table.sem_err[k]:9.2g} "
              f"{table.statistical_inefficiency[k]:8.3f} "
              f"{table.inefficiency_err[k]:8.3f} {s_analytic[k]:8.3f} "
              f"{pull[k]:6.2f}{mark}")
    for w in lee_w:
        print("  warning (lee):", w)
    for w in sme_w:
        print("  warning (sem_err):", w)

    def pick(arr, level):
        return float(arr[level]) if level >= 0 else np.nan

    report_rows.append(dict(
        name=c["name"], title=c["title"], n=n, s_inf=c["s_inf"],
        se_naive=float(sem(x)), se_true=se_true,
        lee_level=lee, lee_B=int(b[lee]) if lee >= 0 else -1,
        lee_se=pick(table.sem, lee),
        lee_s=pick(table.statistical_inefficiency, lee),
        sme_level=sme, sme_B=int(b[sme]) if sme >= 0 else -1,
        sme_se=pick(table.sem, sme),
        sme_s=pick(table.statistical_inefficiency, sme),
        warnings=[f"lee: {w}" for w in lee_w]
        + [f"sem_err: {w}" for w in sme_w],
    ))

    # Per-case figure.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.errorbar(b, table.sem, yerr=table.sem_err, fmt="o-", ms=4,
                 capsize=3, color="C0", label="blocked SEM ± sem_err")
    ax1.axhline(table.sem[0], color="grey", ls=":", label="naive SEM (B = 1)")
    if np.isfinite(se_true):
        ax1.axhline(se_true, color="k", ls="--", lw=1,
                    label="analytic SE, B → ∞")
    few = table.n_blocks <= stats.LEE_MIN_BLOCK_RATIO
    if few.any():
        ax1.axvspan(b[few][0] / 1.4, b[-1] * 1.4, color="C3", alpha=0.08,
                    label="≤ 50 blocks (Lee: unreliable)")
    if lee >= 0:
        ax1.plot(b[lee], table.sem[lee], "^", ms=13, mfc="none", mew=2,
                 color="C2", label=f"lee: B = {b[lee]}")
    if sme >= 0:
        ax1.plot(b[sme], table.sem[sme], "s", ms=13, mfc="none", mew=2,
                 color="C1", label=f"sem_err: B = {b[sme]}")
    ax1.set_xscale("log", base=2)
    ax1.set_xlabel("block size B (frames)")
    ax1.set_ylabel("standard error of the mean")
    ax1.set_title("SEM against block size")
    ax1.legend(fontsize=7, loc="best")

    ax2.errorbar(b, table.statistical_inefficiency,
                 yerr=table.inefficiency_err, fmt="o-", ms=4, capsize=3,
                 color="C0", label="estimated s(B) ± err")
    if c["s_of_b"]:
        bb = 2 ** np.linspace(0, np.log2(b[-1]), 200)
        ax2.plot(bb, c["s_of_b"](np.round(bb).astype(int)), color="k",
                 lw=1, label="analytic s(B)")
        ax2.axhline(c["s_inf"], color="k", ls="--", lw=1,
                    label=f"analytic s, B → ∞ = {c['s_inf']:.3g}")
    if lee >= 0:
        ax2.plot(b[lee], table.statistical_inefficiency[lee], "^", ms=13,
                 mfc="none", mew=2, color="C2", label="lee")
    if sme >= 0:
        ax2.plot(b[sme], table.statistical_inefficiency[sme], "s", ms=13,
                 mfc="none", mew=2, color="C1", label="sem_err")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("block size B (frames)")
    ax2.set_ylabel("statistical inefficiency s")
    ax2.set_title("inefficiency against block size")
    ax2.legend(fontsize=7, loc="best")
    fig.suptitle(c["title"])
    fig.tight_layout()
    fig.savefig(OUT / f"{c['name']}.png", dpi=130)
    plt.close(fig)
    overview.append((c, table, s_analytic, lee, sme))

# 2-D case: three independent observables in one call.
cols = np.column_stack([np.random.default_rng(30).standard_normal(N16),
                        ar1(0.5, N16, 31), ar1(0.9, N16, 32)])
table2 = stats.block_average(cols)
lee2, lee2_w = select(table2, "lee")
sme2, sme2_w = select(table2, "sem_err")
print("\n=== 08_three_columns: columns are iid, AR(1) 0.5, AR(1) 0.9")
print("table shapes: sem", table2.sem.shape, "observable_shape",
      table2.observable_shape)
print("lee level per column    :", lee2, "-> B =", table2.block_size[lee2])
print("sem_err level per column:", sme2, "-> B =", table2.block_size[sme2])
print("sem_blocked(method='lee'):",
      stats.sem_blocked(cols, method="lee"))
print("statistical_inefficiency(method='lee'):",
      stats.statistical_inefficiency(cols, method="lee"),
      "(analytic 1, 3, 19)")
for w in lee2_w + sme2_w:
    print("  warning:", w)

# Overview figure.
fig, axes = plt.subplots(2, 4, figsize=(15, 7))
for ax, (c, table, s_an, lee, sme) in zip(axes.flat, overview):
    b = table.block_size
    ax.errorbar(b, table.statistical_inefficiency,
                yerr=table.inefficiency_err, fmt="o-", ms=3, capsize=2)
    if c["s_of_b"]:
        ax.plot(b, s_an, "k-", lw=1)
    if lee >= 0:
        ax.plot(b[lee], table.statistical_inefficiency[lee], "^", ms=11,
                mfc="none", mew=2, color="C2")
    if sme >= 0:
        ax.plot(b[sme], table.statistical_inefficiency[sme], "s", ms=11,
                mfc="none", mew=2, color="C1")
    ax.set_xscale("log", base=2)
    ax.set_title(c["title"], fontsize=8)
    ax.set_xlabel("B")
    ax.set_ylabel("s(B)")
ax = axes.flat[-1]
for j, lab in enumerate(["iid", "AR(1) 0.5", "AR(1) 0.9"]):
    ax.errorbar(table2.block_size, table2.statistical_inefficiency[:, j],
                yerr=table2.inefficiency_err[:, j], fmt="o-", ms=3,
                capsize=2, label=lab)
    ax.plot(table2.block_size[lee2[j]],
            table2.statistical_inefficiency[lee2[j], j], "^", ms=11,
            mfc="none", mew=2, color="C2")
ax.set_xscale("log", base=2)
ax.set_title("one 2-D call, three columns", fontsize=8)
ax.set_xlabel("B")
ax.legend(fontsize=7)
fig.suptitle("statistical inefficiency against block size; black: analytic;"
             " triangle: lee; square: sem_err")
fig.tight_layout()
fig.savefig(OUT / "overview.png", dpi=110)
plt.close(fig)

# Report.
lines = [
    "# polyphys.analyze.stats demo, synthetic inputs",
    "",
    "SE values are standard errors of the mean of the whole series; "
    "`analytic SE` is sqrt(s var(x) / n) with the B → ∞ inefficiency.",
    "",
    "| case | n | analytic s | naive SE | analytic SE | lee B | lee SE | "
    "lee s | sem_err B | sem_err SE | sem_err s |",
    "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
]
for r in report_rows:
    def f(v, fmt=".4g"):
        return "—" if v is None or (isinstance(v, float) and np.isnan(v)) \
            else format(v, fmt)
    lines.append(
        f"| {r['name']} | {r['n']} | {f(r['s_inf'])} | {f(r['se_naive'])} | "
        f"{f(r['se_true'])} | {r['lee_B'] if r['lee_B'] > 0 else '—'} | "
        f"{f(r['lee_se'])} | {f(r['lee_s'])} | "
        f"{r['sme_B'] if r['sme_B'] > 0 else '—'} | {f(r['sme_se'])} | "
        f"{f(r['sme_s'])} |"
    )
lines += ["", "## Warnings raised", ""]
for r in report_rows:
    for w in r["warnings"]:
        lines.append(f"- {r['name']}, {w}")
lines += ["", "## Files", "",
          "- `<case>.csv`: the full table with analytic s(B) and the pull "
          "(s - s_exact) / s_err",
          "- `<case>.png`: SEM and s against block size",
          "- `overview.png`: all cases",
          "- `log.txt`: printed tables, if the run was piped through `tee`"]
(OUT / "report.md").write_text("\n".join(lines) + "\n")
print("\nwrote", OUT)
