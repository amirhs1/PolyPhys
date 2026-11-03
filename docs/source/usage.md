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

## Adding a study

A parser says how a study's names read; a `ProjectSpec` from
`polyphys.manage.projects` says how its analysis reads them. It holds what the
parser does not: the polymer's particle group; the step `divisor` to which the
bulk volume fraction of crowders, `phi_bulk_c`, is rounded in the
`phi_c_bulk_round` column; the parser attributes that become columns of the
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

PolyPhys holds no study: a study is a folder outside the package, and adding
one changes nothing in `polyphys/`. The maintainer's published studies are in
the repository's top-level
[`projects/`](https://github.com/amirhs1/PolyPhys/tree/main/projects) folder,
one package per study, whose `spec.py` builds the `ProjectSpec` and whose
`__init__.py` exports it as `SPEC`; the study's probe, analysis, and figure
code will join it there. The folder is not installed, and CI lints,
type-checks, and tests it with the package. A study elsewhere uses the same
layout in its own repository.

## The density sum rule

The density sum rule of Jeon et al. (2016, eqns 5 to 7) says that in a
crowded polymer solution the volume fractions of the species, each divided by
its particle diameter, sum to the far-field crowder term:
`sum_i phi_i(x) / a_i ~ phi_c(inf) / a_c`, or equivalently
`sum_i a_i**2 rho_i(x) ~ a_c**2 rho_c(inf)` for the number densities. Here
`x` is `r`, the distance from the cylinder axis in a cylinder and from the
polymer's centre of mass in a periodic cube, or `z`, the axial distance from
the centre of mass.

`space_profiles` in `polyphys.analyze.sumrule` builds one term of the sum from a space of
ensemble-averaged profiles: it scales each `Phi` profile by `1 / a_i` and each
`Rho` profile by `a_i**2`, with `a_i` read from the parser attribute that the
study's `species_sizes` names, normalizes the scaled profile by a reference
value, folds a `z` profile about `z = 0`, and adds the state point's
attributes:

```python
monomers = space_profiles(database, 'Phi', 'r', 'Mon', chain.parser,
                          chain.hierarchy, list(chain.physical_attrs), 'all',
                          chain.species_sizes, divisor=chain.divisor)
crowders = space_profiles(database, 'Phi', 'r', 'Crd', chain.parser,
                          chain.hierarchy, list(chain.physical_attrs), 'all',
                          chain.species_sizes, divisor=chain.divisor)
```

The reference value follows the research code unless `normalizer` says
otherwise: the value at `r = 0` for `r` in a cylinder, where the outermost
bins are the crowder layers on the wall and no far field exists, and the
maximum otherwise (`PROFILE_NORMALIZERS`). `normalizer` also takes `'last'`,
a number, or a function of the scaled frame and the file's parser instance,
such as `lambda df, info: info.phi_bulk_c / info.dcrowd` for the nominal bulk
value of a `Phi` profile. `normalize_r` and `normalize_z` are the two steps on
their own. For `r` in a cylinder, `space_profiles` also adds the bin centre in
crowder radii, recentred so that the wall sits at an integer number of crowder
diameters and the layers of different state points overlay.

`sum_rule_profiles` adds the species' scaled profiles state point by state
point and normalizes the sum by its own reference value, so that the sum is
one at the reference and its departure from one measures the rule's failure:

```python
total = sum_rule_profiles({'Mon': monomers, 'Crd': crowders}, 'rPhi',
                          chain.geometry, 'r')
```

The default (`SUM_NORMALIZERS`) is the value at `r = 0` for `r` in a cylinder
and the far-field value otherwise, where the monomer term vanishes and the sum
is the crowder term, which for several crowders is eqn (7) of Jeon et al.
Under cylindrical confinement the radial rule holds only where the wall
layering is weak, which the ratio `kappa = (D - a_c) / a_c` controls (Sadeghi
Isfahani, 2024, chapter 3).

See the [API reference](api.rst) for the full set of modules, and the
[README](https://github.com/amirhs1/PolyPhys#readme) for a longer introduction
to the artifact-lineage model.

## Reading LAMMPS log files

`read_log` reads a LAMMPS log into a `LogFile` record of three DataFrames:
`thermo` holds one row per thermo output step of every run, `run_stats` one
row per run from the summary LAMMPS prints after it, and `wall_time` the total
wall time. A log holds one run per `run` or `minimize` command, so a
production loop of `run` commands gives one run per iteration; the `run`
column, counted from 0 in log order, picks runs:

```python
from polyphys.simulate.log import read_log

log = read_log("log.lammps")
production = log.thermo[log.thermo["run"] >= 2]
timing = log.run_stats[["run", "loop_time", "timesteps_per_s"]]
```

The thermo styles `one`, `custom`, `multi`, and `yaml` are all read.
`WARNING` and `fix halt` lines inside a table are skipped and logged at level
INFO, which `logging.basicConfig(level=logging.INFO)` shows. A log cut off
before a run's `Loop time of` line, as after a walltime kill, or one that
stopped with an `ERROR` raises `BrokenLogError`, so a loop over many logs can
skip it:

```python
from polyphys.simulate.log import BrokenLogError, read_log

for path in paths:
    try:
        log = read_log(path)
    except BrokenLogError as err:
        print(f"skipped {path}: {err}")
```

`parse_log` does the same for text already in memory, such as a member of a
tarball.

## Reading LAMMPS dumps and writing data files

`iter_dump` reads a text dump of `dump atom` or `dump custom`, plain or
gzipped, one frame at a time, so a long trajectory is never in memory whole.
Each `DumpFrame` holds its `timestep`, its `box`, the column names, and the
per-atom `data`, with rows sorted by atom ID unless `sort_by_id=False` is
passed. `positions()` returns Cartesian coordinates from whichever columns
the dump has: `x y z`, or the scaled `xs ys zs` that `dump atom` writes by
default; with `unwrapped=True`, `xu yu zu`, `xsu ysu zsu`, or the wrapped
coordinates moved by the image flags `ix iy iz`. Orthogonal and restricted
triclinic boxes are both handled:

```python
from polyphys.simulate.dump import iter_dump

for frame in iter_dump("whole.bug.lammpstrj.gz"):
    unwrapped = frame.positions(unwrapped=True)
```

A frame cut off by a killed run raises `BrokenDumpError` once the complete
frames before it have been yielded. `read_dump_frame(path, index)` returns a
single frame.

`read_data` and `write_data` read and write LAMMPS data files in the
`atomic`, `bond`, `angle`, `molecular`, and `full` atom styles, as
`DataFile` records of numpy arrays with a `Box`. `with_frame` gives a dump the topology
it lacks: the atom types, masses, bonds, and angles come from a template, and
the positions, box, and image flags from a dump frame, matched by atom ID:

```python
from polyphys.simulate.data import read_data, with_frame, write_data
from polyphys.simulate.dump import read_dump_frame

template = read_data("initial.bug.data")
frame = read_dump_frame("whole.bug.lammpstrj", 0)
write_data(with_frame(template, frame), "whole.bug.data")
```

`write_data` also writes start configurations for LAMMPS. Build a
`DataFile` from arrays, leave `masses` and `velocities` as `None` for the
input script to set, and close a ring with one more bond. Sections of
force-field coefficients are skipped when read and never written, so the
input script sets the force field.

A dump frame and a data file carry the same `Box` from
`polyphys.simulate.core`, which converts scaled coordinates and applies image
flags; `DataFile.unwrapped_positions()` uses it too. The log, dump, and
data-file readers each read a `.gz` file through gzip, and each raises a
subclass of `LammpsFormatError`, itself a `ValueError`, for a broken file, so
one `except LammpsFormatError` skips a broken file of any kind.

## Building start configurations

`polyphys.simulate.build` writes the data file a run starts from: a
bead-spring polymer in a predictable, unphysical shape that the run then
equilibrates. It follows moltemplate's confined-DNA example in four steps: a
curve, beads placed along it, the topology, then the file. The notation is
that of the module documentation: $N$ beads a bond length $b$ apart on a
helix of radius $R$ and pitch $p$, in the Lennard-Jones units of the run,
with $\sigma$ the bead diameter and $\theta_A$ the helical angle,
$\tan\theta_A = 2\pi R / p$.

`helix` places beads on a helix and `plectoneme_ring` on the closed plectoneme
of Neukirch and Starostin: two coaxial helices, half a turn apart, joined by
smooth end loops. Both place each bead at the first point of the curve a
distance $b$ from the one before, so every bond has exactly the length $b$;
the ring also solves its height so that the $N$ chords close, with $R$, $p$,
$b$, and $N$ exactly as given. `resample_curve` does the same along any
curve given as points. `chain_data` adds the IDs, the molecule, the bonds,
and, if asked, the angles, and `check_conformation` tests the file against
the conditions under which LAMMPS warns or stops before `write_data` writes
it:

```python
from polyphys.simulate import Box, write_data
from polyphys.simulate.build import (
    chain_data, check_conformation, helix, plectoneme_ring,
)

positions = helix(200, radius=2.0, pitch=2.5)
box = Box(lo=[-5.0, -5.0, -21.0], hi=[5.0, 5.0, 21.0])
data = chain_data(positions, box)
check_conformation(data, r0=1.5, ghost_cutoff=1.122462 + 0.3,
                   periodic=(False, False, True), dcyl=6.0, dmon=1.0)
write_data(data, "N200D6.data")

ring = chain_data(plectoneme_ring(200, radius=2.0, pitch=2.5),
                  Box(lo=[-5.0, -5.0, -13.0], hi=[5.0, 5.0, 13.0]),
                  ring=True, angles=True)
```

The default bond length is $0.97\,\sigma$, the start bond length Kremer and
Grest used in 1990 and the mean bond length they measured in the
equilibrated melt; the force of the FENE bond with its WCA term is zero at
$0.961\,\sigma$, so a bond of $0.97$ stores $0.04\,\epsilon$. The old
research default, $0.95\,R_0 = 1.425\,\sigma$, is above the length at which
LAMMPS warns `FENE bond too long`, $\sqrt{0.9}\,R_0 = 1.423$, and above the
default ghost-atom cutoff of a WCA run with a skin of $0.3$, $1.42$. In
runs of this model ($K = 30$, $R_0 = 1.5$, $k_B T = \epsilon$, 10 000 steps)
the start bond length decided which integrator survived:

| Start bond length | plain `nve`, dt 0.005 | `nve` + Langevin, dt 0.01 | `nve/limit 0.01` + Langevin |
| ----------------- | --------------------- | ------------------------- | --------------------------- |
| 0.97 to 1.10      | runs, no warning      | runs                      | runs                        |
| 1.15 to 1.35      | runs, warns           | stops from 1.30           | runs                        |
| 1.40 to 1.49      | stops                 | stops                     | runs, warns                 |
| 3.0 and above     | stops at step 0       | stops at step 0           | stops at step 0             |

The radius and the pitch are the user's choice, with no default. Adjacent
turns of a helix are $p\sin\theta_A$ apart and the two strands of a
plectoneme $(p/2)\sin\theta_A$, so a helix needs $p \ge d_{\min} /
\sin\theta_A$ and a ring $p \ge 2 d_{\min} / \sin\theta_A$, about $2.0$ to
$2.1\,\sigma$ for $R \ge 2\sigma$, to keep every pair of beads the minimum
distance $d_{\min}$ apart; `check_conformation` enforces it, with
$d_{\min} = \sigma$ by default, where the WCA energy of a pair is
$\epsilon$. Its other rules are the FENE warning length, the ghost cutoff,
coincident atoms, the box on non-periodic axes, and, in a cylinder of
diameter $D$ built from beads of diameter $a$, the effective radius
$(D - a)/2$. Equilibration stays in the input script: `fix nve/limit` with a
Langevin thermostat, as the research runs do, or the soft pair potential
Kremer and Grest raised slowly to remove overlaps.

`polyphys/tests/data/lammps_build/` holds `in.kg`, a Kremer-Grest input
script that runs a built start under a chosen integrator, and
`check_starts.sh`, which builds the three example starts and runs each
under all of them; the same matrix runs as the `lammps`-marked tests when a
LAMMPS executable is available (see the Testing page).
