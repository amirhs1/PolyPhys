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
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git@v0.9.0"
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

PolyPhys is not on conda-forge yet; a conda-forge package is planned after
v1.0.0 is on PyPI. Until it is listed, create a conda environment with the
dependencies from conda-forge, then install PolyPhys into it with pip, as
[conda's guide](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
advises: conda packages first, then pip. `mdanalysis` is needed only by the
probe layer and `matplotlib-base` only by the figure helpers of
`polyphys.viz`; leave out the ones you do not need:

```bash
conda create -n polyphys -c conda-forge "python>=3.11" "numpy>=2.0" "pandas>=3.0" "scipy>=1.14.1" "mdanalysis>=2.8" "matplotlib-base>=3.8.4"
conda activate polyphys
python -m pip install polyphys
```

## Development install

If you are working from a source checkout and want the documentation and
development tools:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[docs,dev,mda,viz]"
```
