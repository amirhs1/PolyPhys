# Studies: adding a study

This page shows how to add a study to PolyPhys: the parser class that reads
its file names, the specs that describe its probe and its analysis, the
folder that holds them, and the checks that run on them. It also shows how
to read the tables that the research code before PolyPhys wrote.

## Notation

*study*
: One research project: one simulation model, swept over its parameters,
  such as `HnsCub`, a ring polymer with H-NS proteins among crowders in a
  periodic cube.

*keyword*
: A short label in a file name followed by its value, such as `N200` for a
  chain of 200 monomers or `ens1` for the first replicate.

*parser class*
: A subclass of {class}`~polyphys.manage.parser.ParserBase` that declares a
  study's keywords and reads a name into typed attributes and lineage names.

*spec*
: A record of a study's choices: the
  {class}`~polyphys.manage.projects.ProjectSpec` of its analysis, the
  {class}`~polyphys.manage.projects.ProbeSpec` of its probe, and one
  {class}`~polyphys.analyze.pipeline.AnalysisTable` per particle group.

*species*
: A kind of particle, labeled as in the property names: `Mon` for the
  monomers, `Crd` for the crowders, `Hns` for the H-NS cores.

$a_i$
: The diameter of the particles of species $i$, in the length unit of the
  simulation.

The [Lineage](lineage.md) page defines the lineage levels and the groups.

## The steps

PolyPhys holds no study: adding one changes nothing in `polyphys/`.

1. Declare the study's parser class, as in "Declaring the parser" below.
2. Make a folder for the study, `projects/<study>/`, as a Python package.
   Its `spec.py` builds the study's `ProjectSpec` (as `SPEC`), its
   `ProbeSpec` (as `PROBE`), and a tuple of its `AnalysisTable` records,
   one per particle group (as `ANALYSIS`); its `__init__.py` exports the
   three. The [Probe](probe.md) page, "A study's probe as data", and the
   [Pipeline](pipeline.md) page, "A study's analysis as data", describe the
   last two.
3. Build the specs. Each attribute a `ProjectSpec` names is checked
   against the parser when the spec is built, so a typo raises when
   `spec.py` is imported, not when the analysis runs.
4. Test the study. In this repository, `projects/test_studies.py` checks
   each study's specs against its parser and against the values of the
   research code; add the new study to its parametrizations. CI lints,
   type-checks, and tests `projects/` with the package.
5. Run the study's data through the phases, as the [Workflow](workflow.md)
   page describes.

The maintainer's published studies are in the repository's top-level
[`projects/`](https://github.com/amirhs1/PolyPhys/tree/main/projects)
folder, which is not installed with the package. A study elsewhere uses the
same layout in its own repository. The [Design decisions](decisions.md)
page records why each study is a folder outside the package.

## Declaring the parser

A parser class is a declaration: it lists the study's filename keywords, the
order in which each lineage writes them, any fixed values, the particle
species whose bulk number density and volume fraction to compute, and the
attribute that holds each group's dump period (here fixed in `_defaults`).
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
...     _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
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
...     _defaults = {'dmon': 1.0, 'bdump': 1000, 'adump': 5000}
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
>>> artifact.dump_keyword, artifact.bdump
('bdump', 1000)

```

The declaration is checked when the class is created: a pattern keyword
without a `KeywordSpec`, an attribute that appears twice in one pattern, a
species whose size or number no long-form name sets, or a group without a
dump keyword raises `TypeError`. A
name that lacks a keyword of its lineage's pattern, or whose value is not a
number, raises `ValueError` naming the keyword, the lineage, and the name.

A value that is a function of parsed values but fits none of the tables,
such as the complementary percentage `pct_two = 100 - pct_one` of
`SumRuleCubHeteroRingBiDisperse`, goes in an override of
`_derived_attributes`, which runs at every lineage after the bulk values.

## Declaring the ProjectSpec

A parser says how a study's names read; a `ProjectSpec` from
`polyphys.manage.projects` says how its analysis reads them. It holds what the
parser does not: the polymer's particle group; the step `divisor` to which the
bulk volume fraction of crowders, `phi_bulk_c`, is rounded in the
`phi_bulk_c_round` column; the parser attributes that become columns of the
`space` tables; each density-profile species label mapped to the attribute
holding that species' diameter; and the study's time series. The geometry,
topology, and glob patterns come from the parser:

```python
>>> from polyphys.manage.projects import ProjectSpec
>>> chain = ProjectSpec(
...     parser=ChainCub, group='bug', divisor=0.025,
...     physical_attrs=('space', 'ensemble_long', 'nmon', 'dcrowd'),
...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
...     tseries_properties=('gyrTMon',))
>>> chain.geometry, chain.directions
('cubic', ('r',))
>>> chain.space_pattern, chain.hierarchy
('N*ac*', 'N*')

```

Each attribute a spec names is checked against the attributes the parser sets
at `ensemble_long`, so a typo or a missing declaration raises `ValueError` when
the spec is built. Pass the fields to the organizer, as in
`space_tseries(database, 'gyrTMon', chain.parser, chain.hierarchy,
list(chain.physical_attrs), chain.group, divisor=chain.divisor)`.

A study that probes its trajectories with `polyphys.probe` builds a
`ProbeSpec` beside its `ProjectSpec`: the MDAnalysis selection of each
species, the label of the chain, the chain-size observables to measure,
whether the chain is unwrapped for its asphericity, and the settings of its
spatial histograms, pair distances, and binders; see [Probe](probe.md).
The study's drivers and figure code will join its specs in its folder.

## Reading tables the research code wrote

PolyPhys names a bulk quantity `<quantity>_bulk_<species>`; the research code
it replaces named it `<quantity>_<species>_bulk`, and wrote every parser
attribute as a column of its `*-stamps.csv` files, from where the names carry
into its ensemble-average and `space` tables. `BULK_NAME_MAP` in
`polyphys.manage.organizer` maps each old name to its new one, and
`rename_bulk_columns` applies the map to a table read from such a file,
leaving every other column unchanged:

| Research code      | PolyPhys           |
| ------------------ | ------------------ |
| `phi_m_bulk`       | `phi_bulk_m`       |
| `rho_m_bulk`       | `rho_bulk_m`       |
| `phi_c_bulk`       | `phi_bulk_c`       |
| `rho_c_bulk`       | `rho_bulk_c`       |
| `phi_c_bulk_one`   | `phi_bulk_c_one`   |
| `phi_c_bulk_two`   | `phi_bulk_c_two`   |
| `rho_c_bulk_one`   | `rho_bulk_c_one`   |
| `rho_c_bulk_two`   | `rho_bulk_c_two`   |
| `phi_hns_bulk`     | `phi_bulk_hns`     |
| `rho_hns_bulk`     | `rho_bulk_hns`     |
| `phi_c_bulk_round` | `phi_bulk_c_round` |

```python
>>> import pandas as pd
>>> from polyphys.manage.organizer import rename_bulk_columns
>>> old = pd.DataFrame({'nmon': [20], 'phi_c_bulk': [0.1], 'n_frames': [100]})
>>> list(rename_bulk_columns(old).columns)
['nmon', 'phi_bulk_c', 'n_frames']

```

The parser classes carry no alias for an old name: a name reaches a file only
when a parser lists it in `attributes`, and an alias there would give every
new table both columns.
