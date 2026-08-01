# Usage

To begin using **PolyPhys**, first install the package as described in the
[Installation](#installation_target) section below.

(installation_target)=
## Installation

PolyPhys is distributed from source and archived on Zenodo. It is not published
to PyPI, so install it directly from the repository:

```bash
python -m pip install "git+https://github.com/amirhs1/poly-phys.git"
```

To install a specific archived release, append the tag:

```bash
python -m pip install "git+https://github.com/amirhs1/poly-phys.git@v0.3.0"
```

PolyPhys requires Python 3.11 or newer. Its runtime dependencies are NumPy and
pandas.

### Development install

If you are working from a source checkout and want the documentation and
development tools:

```bash
git clone https://github.com/amirhs1/poly-phys.git
cd poly-phys
python -m pip install -e ".[docs,dev]"
```

## Quickstart

`polyphys.manage.parser` decodes a simulation filename into typed physical
attributes and derives the artifact's position in the project hierarchy:

```python
>>> from polyphys.manage.parser import SumRuleCyl
>>> artifact = SumRuleCyl(
...     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1.j03',
...     lineage='segment',
...     group='all',
... )
>>> artifact.nmon
200
>>> float(artifact.dcyl)
20.0
>>> artifact.space
'N200D20.0ac2.0'

```

See the [API reference](api.rst) for the full set of modules, and the
[README](https://github.com/amirhs1/poly-phys#readme) for a longer introduction
to the artifact-lineage model.
