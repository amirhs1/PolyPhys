# Quickstart: a synthetic sweep from segments to a space

This page runs PolyPhys on a sweep small enough to check by hand: two state
points, two replicates of each, and two segments of each replicate, each
segment a time series of two frames. It writes the probe files itself,
reduces them down the lineage with
{func}`~polyphys.analyze.pipeline.reduce_space`, and stacks the result into
a `space` table with {func}`~polyphys.manage.organizer.space_tseries`. Every
example on this page runs as a test, so its output is what the code gives.

## Notation

$e$, $w$, $s$, $k$
: The indices of a state point, a replicate (a whole), a segment, and a
  frame within a segment; $e$, $w$, and $s$ start at 1 and $k$ at 0.

$x_{ewsk}$
: The value in frame $k$ of segment $s$ of replicate $w$ of state point
  $e$.

$R$
: The number of replicates of a state point, here 2.

The [Lineage](lineage.md) page defines the lineage levels, phases, and
stages, and the [Workflow](workflow.md) page the steps this page runs.

## The synthetic sweep

The values follow the scheme of the test suite's synthetic lineage trees:

$$
x_{ewsk} = 1000\,e + 100\,w + 10\,s + k,
$$

so each value names its file and frame, and every average below can be
written down by hand. The files are named for the `SumRuleCyl` parser, a
linear chain of 200 monomers in a cylinder among crowders; the two state
points differ only in the number of crowders, `nc1000` and `nc2000`. Each
segment file is named `<whole>.j<NN>.<group>.<property>.npy`, in a directory
per whole under the space's probe directory:

```python
>>> import tempfile
>>> from pathlib import Path
>>> import numpy as np
>>> tmp = tempfile.TemporaryDirectory()
>>> root = Path(tmp.name)
>>> probe = root / 'probes' / 'N200D20ac2-bug-probe'
>>> for e, nc in enumerate((1000, 2000), start=1):
...     state = f'N200epsilon5.0r10.5lz25.0sig2.0nc{nc}dt0.005bdump2000adump5000'
...     for w in (1, 2):
...         whole = probe / f'{state}ens{w}'
...         whole.mkdir(parents=True)
...         for s in (1, 2):
...             x = 1000 * e + 100 * w + 10 * s + np.arange(2.0)
...             np.save(whole / f'{state}ens{w}.j{s:02d}.bug.gyrTMon.npy', x)

```

## The study's specs

A {class}`~polyphys.spec.ProjectSpec` says which parser reads the
names, which particle group the polymer is, which attributes become columns
of the `space` table, and the step to which the crowder volume fraction is
rounded. An {class}`~polyphys.spec.AnalysisTable` says which
properties to reduce: here one time series, `gyrTMon`, read as segments:

```python
>>> from polyphys.manage.parser import SumRuleCyl
>>> from polyphys.spec import AnalysisTable, ProjectSpec, tseries
>>> spec = ProjectSpec(
...     parser=SumRuleCyl, group='bug', divisor=0.025,
...     crowder_attr='phi_bulk_c',
...     physical_attrs=('space', 'ensemble_long', 'ncrowd'),
...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
...     tseries_properties=('gyrTMon',))
>>> table = AnalysisTable(group='bug', lineage='segment', stamps=False,
...                       properties=(tseries('gyrT', 'Mon'),))

```

The [Studies](studies.md) page explains each field.

## Down the lineage

{func}`~polyphys.analyze.pipeline.reduce_space` joins the two segments of
each whole into one series of four frames, collects the two wholes of each
state point into an ensemble, and averages each ensemble over its
replicates:

```python
>>> from polyphys.analyze.pipeline import reduce_space
>>> result = reduce_space(probe, spec, table)
>>> state = 'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000'
>>> average = result['gyrTMon'][state]
>>> average['gyrTMon-mean'].tolist()
[1160.0, 1161.0, 1170.0, 1171.0]
>>> average['gyrTMon-var'].tolist()
[5000.0, 5000.0, 5000.0, 5000.0]
>>> average['gyrTMon-sem'].tolist()
[50.0, 50.0, 50.0, 50.0]

```

By hand, for $e = 1$: whole $w$ is segment 1 then segment 2,
$(1010 + 100w,\ 1011 + 100w,\ 1020 + 100w,\ 1021 + 100w)$. The two wholes
differ by 100 in every frame, so their mean is the value at $w = 1.5$,
$1160, 1161, 1170, 1171$; their sample variance is $(50^2 + 50^2) / (R - 1)
= 5000$; and the standard error of the mean is $\sqrt{5000 / R} = 50$.

`reduce_space` wrote each level to its stage, in the `analysis` phase beside
the `probes` phase:

```python
>>> for stage in sorted((root / 'analysis').iterdir()):
...     print(stage.name, len(list(stage.iterdir())))
N200D20ac2-bug-ens 2
N200D20ac2-bug-ensAvg 2
N200D20ac2-bug-wholeSim 4

```

Four wholes, and two state points at each of the ensemble stages.

## The space table

{func}`~polyphys.analyze.pipeline.space_tseries_from_spec` stacks the
ensemble averages of the two state points into one table and adds the time
and the state point's attributes, with the columns and the rounding step the
spec declares:

```python
>>> from polyphys.analyze.pipeline import space_tseries_from_spec
>>> ens_avg = root / 'analysis' / 'N200D20ac2-bug-ensAvg'
>>> space = space_tseries_from_spec(ens_avg, spec, 'gyrTMon', save=True)
>>> space['ncrowd'].tolist()
[1000, 1000, 1000, 1000, 2000, 2000, 2000, 2000]
>>> space['time'].tolist()
[0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0]
>>> space['gyrTMon-mean'].tolist()
[1160.0, 1161.0, 1170.0, 1171.0, 2160.0, 2161.0, 2170.0, 2171.0]
>>> sorted(set(space['phi_bulk_c_round']))
[0.325, 0.65]
>>> sorted(p.name for p in (root / 'analysis').glob('*-space/*'))
['N200D20ac2-bug-gyrTMon-space.csv']

```

By hand: the `bug` trajectory is dumped every `bdump2000` steps of
`dt0.005`, so one frame is $2000 \times 0.005 = 10$ time units. The bulk
crowder volume fraction of the first state point is 0.3292 (the README's
Quickstart computes it), and of the second, with twice the crowders, twice
that, 0.6584; rounded to the nearest multiple of the spec's `divisor`,
0.025, they are 0.325 and 0.65.

```python
>>> tmp.cleanup()

```

## Next

- The [Workflow](workflow.md) page adds the probe step before and the
  equilibrium values of the wholes after.
- The [Statistics](statistics.md) page explains why the error of one whole's
  time average needs block averaging, where the error across replicates
  above does not.
- The [Studies](studies.md) page shows how to declare a study's own parser
  and specs.
