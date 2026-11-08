# Usage

To begin using **PolyPhys**, first install the package as described in the
[Installation](#installation_target) section below.

(installation_target)=
## Installation

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

### With conda

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

### Development install

If you are working from a source checkout and want the documentation and
development tools:

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[docs,dev,mda,viz]"
```

## Quickstart

A parser decodes a simulation filename into typed physical attributes and
derives the artifact's position in the project hierarchy. A study declares
its parser as a subclass of {class}`~polyphys.manage.parser.ParserBase`, as
the [Studies](studies.md) page shows; this one reads a linear chain among
crowders in a cylinder:

```python
>>> from polyphys.manage.parser import KeywordSpec, ParserBase, Species
>>> long_keys = ('N', 'r', 'lz', 'ac', 'nc', 'dt', 'bdump', 'adump')
>>> bulk = ['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c']
>>> class ChainCyl(ParserBase):
...     """A linear chain among crowders in a cylinder, periodic along z."""
...     _geometry = 'cylindrical'
...     _topology = 'linear'
...     _groups = ['bug', 'all']
...     _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
...     _keywords = {
...         'N': KeywordSpec('nmon', int),
...         'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
...         'D': KeywordSpec('dcyl', float),
...         'lz': KeywordSpec('lcyl', float, scale=2.0),
...         'ac': KeywordSpec('dcrowd', float),
...         'nc': KeywordSpec('ncrowd', int),
...         'dt': KeywordSpec('dt', float),
...         'bdump': KeywordSpec('bdump', int),
...         'adump': KeywordSpec('adump', int),
...         'ens': KeywordSpec('ensemble_id', int),
...         'j': KeywordSpec('segment_id', int),
...     }
...     _patterns = {
...         'segment': (*long_keys, 'ens', 'j'),
...         'whole': (*long_keys, 'ens'),
...         'ensemble_long': long_keys,
...         'ensemble': ('N', 'D', 'ac', 'nc'),
...         'space': ('N', 'D', 'ac'),
...     }
...     _defaults = {'dmon': 1}
...     _species = (Species('m', 'nmon', 'dmon'),
...                 Species('c', 'ncrowd', 'dcrowd'))
...     _project_attributes = {
...         'segment': bulk, 'whole': bulk, 'ensemble_long': bulk,
...         'ensemble': ['dmon'], 'space': ['dmon']}
>>> artifact = ChainCyl(
...     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1.j03',
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
  segment `…ens1.j03.linear` and a whole `…ens1.linear` of one simulation
  both give the `whole` name `…ens1`.
- The short names, `ensemble` and `space`, are built from the parsed values.
  A float that is a whole number prints as an integer, and any other float
  prints as Python's `repr`, the shortest text that reads back as the same
  float. The `space` above is `N200D20ac2` although the name wrote `ac2.0`
  and `r10.5`, while a crowder size of `1.5` stays `ac1.5`.

```python
>>> artifact = ChainCyl(
...     'N200r10.5lz25ac1.5nc0dt0.005bdump2000adump5000ens1.linear', 'whole',
...     'bug')
>>> artifact.whole, artifact.space
('N200r10.5lz25ac1.5nc0dt0.005bdump2000adump5000ens1', 'N200D20ac1.5')

```

Names written by the research code that preceded PolyPhys can differ from
these in two ways: a `whole` it parsed from a whole file keeps `.ring` or
`.linear`, and a short name can print `ac2.0`. Parsing such a short name at
its lineage and reading the rebuilt name gives the current form.

## Adding a study

A study is a parser class, which says how its file names read, and its
specs, which say how its probe and its analysis read them. The
[Studies](studies.md) page shows how to declare both, where a study's
folder lives, and how to read the tables the research code wrote; the
examples below use the spec `chain` built there.

## The density sum rule

The density sum rule of Jeon et al. (2016, eqns 5 to 7) says that in a
crowded polymer solution the volume fractions of the species, each divided by
its particle diameter, sum to the far-field crowder term:
`sum_i phi_i(x) / a_i ~ phi_c(inf) / a_c`, or equivalently
`sum_i a_i**2 rho_i(x) ~ a_c**2 rho_c(inf)` for the number densities. Here
`x` is `r`, the distance from the cylinder axis in a cylinder and from the
polymer's center of mass in a periodic cube, or `z`, the axial distance from
the center of mass.

`space_profiles` in `polyphys.analyze.sumrule` builds one term of the sum from a space of
ensemble-averaged profiles: it scales each `Phi` profile by `1 / a_i` and each
`Rho` profile by `a_i**2`, with `a_i` read from the parser attribute that the
study's `species_sizes` names, normalizes the scaled profile by a reference
value, folds a `z` profile about `z = 0`, and adds the state point's
attributes. The profiles must be the raw ensemble averages, in their
physical units: the rule relates the species through their scales, so a
profile normalized beforehand, by its sum, its integral, or its maximum
(the `normalize` option of `polyphys.analyze.density`, or by hand), enters
the sum with an arbitrary weight and its `-scale` column loses its unit.
`space_profiles` warns, naming the file, when a profile's `-mean` sums to
one, the signature of a profile normalized by its sum; it cannot tell the
other normalizations from one profile.

`space_profiles_from_spec` reads the parser, the `hierarchy` glob, the
`physical_attrs`, the `species_sizes`, the `divisor`, and the `crowder_attr`
from the study's spec, and takes the profiles from the `'all'` group unless
`group` says otherwise:

```python
monomers = space_profiles_from_spec(database, chain, 'Phi', 'r', 'Mon')
crowders = space_profiles_from_spec(database, chain, 'Phi', 'r', 'Crd')
```

`space_profiles` takes the same values one by one.

The reference value follows the research code unless `normalizer` says
otherwise: the value at `r = 0` for `r` in a cylinder, where the outermost
bins are the crowder layers on the wall and no far field exists, and the
maximum otherwise (`PROFILE_NORMALIZERS`). `normalizer` also takes `'last'`,
a number, or a function of the scaled frame and the file's parser instance,
such as `lambda df, info: info.phi_bulk_c / info.dcrowd` for the nominal bulk
value of a `Phi` profile. `normalize_r` and `normalize_z` are the two steps on
their own. For `r` in a cylinder, `space_profiles` also adds the bin center in
crowder radii, recentered so that the wall sits at an integer number of crowder
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

The ensAvg file's `-sem` is the error of the mean profile, not of the ratio
of a bin to its normalizer. The pipeline (`reduce_space` with
`LocalDensities`) also writes `<prop>-norm-sem` beside the ensAvg columns:
each whole's profile is divided by its own value in the bin the default
normalizer names, the first, the last, or the bin where the ensemble average
peaks, and the standard error is taken over the wholes, which are
independent runs. `<prop>-norm-sem` is the SEM across profiles normalized
separately for each whole. It is not generally an uncertainty estimate for
the normalized ensemble mean: normalization is nonlinear, so the mean of the
normalized wholes and the normalized mean differ. The curve itself stays the
normalized ensemble average, the published estimator.
`space_profiles` carries the column through under the default normalizer,
folding it along `z`, and drops it with a warning under an explicit one.
`sum_rule_profiles` adds the species' columns in quadrature as
`<prop>Sum-norm-sem`. The quadrature sum assumes that the species' errors
are independent, although excluded volume correlates the species within a
whole; it refers to the sum of the separately normalized species, and to
`Sum-norm` only when the species share one normalizer value.

See the [API reference](api.rst) for the full set of modules, and the
[README](https://github.com/amirhs1/PolyPhys#readme) for a longer introduction
to the artifact-lineage model.

## Measures of one frame

{mod}`polyphys.analyze.measurer` holds the small functions the probes and
the analysis build on, each on the positions or values it is given. The
[API reference](api.rst) documents each one; this table says which to reach
for:

| Kind | Functions |
| --- | --- |
| Chain size, on one frame's positions | {func}`~polyphys.analyze.measurer.end_to_end`, {func}`~polyphys.analyze.measurer.transverse_size`, {func}`~polyphys.analyze.measurer.fsd`, {func}`~polyphys.analyze.measurer.max_distance` |
| Distances under periodic boundaries | {func}`~polyphys.analyze.measurer.apply_pbc_orthogonal`, {func}`~polyphys.analyze.measurer.pair_distance` |
| Bins | {func}`~polyphys.analyze.measurer.fixedsize_bins`, {func}`~polyphys.analyze.measurer.create_bin_edge_and_hist` |
| Histograms of one frame | {func}`~polyphys.analyze.measurer.radial_histogram` (spherical $r$), {func}`~polyphys.analyze.measurer.radial_cyl_histogram` (cylindrical $r$), {func}`~polyphys.analyze.measurer.axial_histogram`, {func}`~polyphys.analyze.measurer.azimuth_cyl_histogram`, {func}`~polyphys.analyze.measurer.planar_cartesian_histogram` |
| Volumes of a bead cut by planes or a sphere | {func}`~polyphys.analyze.measurer.spherical_segment`, {func}`~polyphys.analyze.measurer.sphere_sphere_intersection` |
| Statistics of independent samples | {func}`~polyphys.analyze.measurer.simple_stats`, {func}`~polyphys.analyze.measurer.sem` |

`simple_stats` and `sem` treat their samples as independent, as the
replicates of a state point are; for the correlated frames of one run, use
the blocked error of {mod}`polyphys.analyze.stats`, as the
[Statistics](statistics.md) page explains.

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
