# Usage

To begin using **PolyPhys**, first install the package as described in the
[Installation](#installation_target) section below.

(installation_target)=
## Installation

PolyPhys is distributed from source and archived on Zenodo. It is not published
to PyPI, so install it directly from the repository:

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git"
```

To install a specific archived release, append the tag:

```bash
python -m pip install "git+https://github.com/amirhs1/PolyPhys.git@v0.4.0"
```

PolyPhys requires Python 3.11 or newer. Its runtime dependencies are NumPy and
pandas.

### Development install

If you are working from a source checkout and want the documentation and
development tools:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
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
'N200D20ac2'

```

## Lineage names

A parser sets the artifact's name at its own lineage and at each parent
lineage, as the attributes `segment`, `whole`, `ensemble_long`, `ensemble`,
and `space`; the organizer uses them as join keys. Two rules fix their text.

- The long names, `segment`, `whole`, and `ensemble_long`, are cut from the
  parsed name, so they keep its text, and they carry no topology suffix. A
  segment `…ens1.j03.ring` and a whole `…ens1.ring` of one simulation both
  give the `whole` name `…ens1`.
- The short names, `ensemble` and `space`, are built from the parsed values.
  A float that is a whole number prints as an integer, and any other float
  prints as Python's `repr`, the shortest text that reads back as the same
  float. The `space` above is `N200D20ac2` although the name wrote `sig2.0`
  and `r10.5`, while a crowder size of `1.5` stays `ac1.5`.

```python
>>> from polyphys.manage.parser import HnsCub
>>> artifact = HnsCub('N200kbmm2nh12ac2l25epshc1.5nc0ens1.ring', 'whole',
...                   'nucleoid')
>>> artifact.whole, artifact.space
('N200kbmm2nh12ac2l25epshc1.5nc0ens1', 'N200kbmm2nh12ac2epshc1.5')

```

Names written by the research code that preceded PolyPhys can differ from
these in two ways: a `whole` it parsed from a whole file keeps `.ring` or
`.linear`, and a short name can print `ac2.0`. Parsing such a short name at
its lineage and reading the rebuilt name gives the current form.

## Declaring a parser for a new study

A parser class is a declaration: it lists the study's filename keywords, the
order in which each lineage writes them, any fixed values, and the particle
species whose bulk number density and volume fraction to compute.
`ParserBase` does the parsing. Each keyword gets a `KeywordSpec` naming the
attribute that stores its value and the value's type; a `float` keyword may
also convert the value as `scale * value + offset`, here `lcube = 2 l` from
the half-side `l`:

```python
>>> from polyphys.manage.parser import KeywordSpec, ParserBase, Species
>>> class ChainCub(ParserBase):
...     _geometry = 'cubic'
...     _topology = 'linear'
...     _groups = ['bug', 'all']
...     _keywords = {
...         'N': KeywordSpec('nmon', int),
...         'ac': KeywordSpec('dcrowd', float),
...         'nc': KeywordSpec('ncrowd', int),
...         'l': KeywordSpec('lcube', float, scale=2.0),
...         'ens': KeywordSpec('ensemble_id', int),
...         'j': KeywordSpec('segment_id', int),
...     }
...     _patterns = {
...         'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
...         'whole': ('N', 'ac', 'nc', 'l', 'ens'),
...         'ensemble_long': ('N', 'ac', 'nc', 'l'),
...         'ensemble': ('N', 'ac', 'nc'),
...         'space': ('N', 'ac'),
...     }
...     _defaults = {'dmon': 1.0}
...     _species = (Species('m', 'nmon', 'dmon'),
...                 Species('c', 'ncrowd', 'dcrowd'))
...     _project_attributes = {
...         'segment': ['dmon', 'rho_bulk_m', 'rho_bulk_c'],
...         'whole': ['dmon', 'rho_bulk_m', 'rho_bulk_c'],
...         'ensemble_long': ['dmon', 'rho_bulk_m', 'rho_bulk_c'],
...         'ensemble': ['dmon'],
...         'space': ['dmon'],
...     }
>>> artifact = ChainCub('N100ac2.0nc50l10.0ens1.bug', 'whole', 'bug')
>>> artifact.lcube, artifact.ensemble, artifact.space
(20.0, 'N100ac2nc50', 'N100ac2')

```

The declaration is checked when the class is created: a pattern keyword
without a `KeywordSpec`, an attribute that appears twice in one pattern, or a
species whose size or number no long-form name sets raises `TypeError`. A
name that lacks a keyword of its lineage's pattern, or whose value is not a
number, raises `ValueError` naming the keyword, the lineage, and the name.

A value that is a function of parsed values but fits none of the tables,
such as the complementary percentage `pct_two = 100 - pct_one` of
`SumRuleCubHeteroRingBiDisperse`, goes in an override of
`_derived_attributes`, which runs at every lineage after the bulk values.

See the [API reference](api.rst) for the full set of modules, and the
[README](https://github.com/amirhs1/PolyPhys#readme) for a longer introduction
to the artifact-lineage model.
