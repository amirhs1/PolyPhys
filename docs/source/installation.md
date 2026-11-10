(installation_target)=
# Installation

From v1.0.0, PolyPhys is published to PyPI, and each release is archived on
Zenodo:

```bash
python -m pip install polyphys
```

With [uv](https://docs.astral.sh/uv/), install it into the active environment,
or add it as a dependency of a uv project:

```bash
uv pip install polyphys
uv add polyphys
```

Releases before v1.0.0 are not on PyPI. To install one, give its tag:

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git@v1.0.0"
```

PolyPhys requires Python 3.11 or newer. Its runtime dependencies are NumPy 2.0
or newer, pandas 3.0 or newer, and SciPy 1.14.1 or newer. The probe layer,
`polyphys.probe`, reads LAMMPS trajectories with
[MDAnalysis](https://www.mdanalysis.org/), which the `mda` extra installs:

```bash
python -m pip install "polyphys[mda]"
```

Without it, `import polyphys.probe` raises `ImportError` with that hint, and
the rest of the package works as before.

## With conda

From v1.0.0, PolyPhys is also on
[conda-forge](https://anaconda.org/conda-forge/polyphys):

```bash
conda install -c conda-forge polyphys
```

The conda package installs the runtime dependencies only. `mdanalysis` is
needed by the probe layer and `matplotlib-base` by the figure helpers of
`polyphys.viz`; install the ones you need alongside PolyPhys:

```bash
conda install -c conda-forge mdanalysis matplotlib-base
```

When they are installed, the package requires MDAnalysis 2.8 and matplotlib
3.8.4 or newer, the floors of the `mda` and `viz` extras.

## Development install

If you are working from a source checkout and want the documentation and
development tools:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[docs,dev,mda,viz]"
```
