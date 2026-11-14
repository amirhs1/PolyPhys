# Viz: color lists and matplotlib helpers for figures

This page explains what {mod}`polyphys.viz` offers a figure script: the
color lists the studies' figures use, a way to cut a color map down to part
of its range, a rule for placing an axis's ticks, two legend-handle makers,
a legend mover, and a band shader. Everything here is a thin layer over
matplotlib; the figures themselves, with their study-specific labels,
filters, and layouts, belong to each study, not to the package.

The layer is an optional extra, since the rest of PolyPhys never draws:

```bash
python -m pip install "polyphys[viz]"
```

{mod}`polyphys.viz.colors` imports no plotting library, so its color lists
work on a bare install; its one function,
{func}`~polyphys.viz.colors.truncate_colormap`, imports matplotlib when
called. {mod}`polyphys.viz.helpers` needs matplotlib at import. Both raise
`ImportError` with the install hint when matplotlib is missing, and
`import polyphys.viz` imports neither; import the module you need.

The examples on this page run as doctests in CI, headless, so they select
the Agg backend first; a script that shows figures on screen does not need
that line.

## Colors

{mod}`polyphys.viz.colors` holds five tuples of `#RRGGBB` strings.

`TOL_BRIGHT`, `TOL_VIBRANT`, `TOL_MUTED`
: Paul Tol's qualitative schemes, from *Colour Schemes*, SRON Technical
  Note SRON/EPS/TN/09-002, issue 3.2, 18 August 2021, Figures 1, 3, and 4,
  in the order the note recommends for picking colors in a fixed sequence.
  They are designed to stay distinct for color-blind readers and in print.
  Each ends with the scheme's gray, which the note reserves for bad or
  missing data; take `palette[:-1]` for the data colors alone.

`DUTCH_FIELD`
: A nine-color palette of that name from a HEAVY.AI blog post that could
  not be opened when the module was written; the values are the ones a
  wpDataTables article lists under the same name. The module docstring
  gives both references. Neither page states terms of reuse.

`STUDY_COLORS`
: The maintainer's twelve-color list for the published figures.

```python
>>> from polyphys.viz.colors import TOL_BRIGHT, TOL_MUTED
>>> TOL_BRIGHT[:3]
('#4477AA', '#EE6677', '#228833')
>>> len(TOL_MUTED)
10

```

{func}`polyphys.viz.colors.truncate_colormap` takes a color map, or the
name of a registered one, and the sub-range `[lower, upper]` of its domain
to keep, and returns a linear segmented map through `n` samples of it. The
new map's value at 0 is the parent's at `lower`, and at 1 the parent's at
`upper`, so a sequential map can start past its pale end.

```python
>>> import matplotlib
>>> matplotlib.use("Agg")
>>> from polyphys.viz.colors import truncate_colormap
>>> dark = truncate_colormap("viridis", 0.5, 1.0, n=4)
>>> dark.name
'trunc(viridis,0.5,1.0)'
>>> dark(0.0) == matplotlib.colormaps["viridis"](0.5)
True

```

## Ticks

{func}`polyphys.viz.helpers.set_xticks` and
{func}`polyphys.viz.helpers.set_yticks` place the major ticks of one axis
by one rule. Given `lower`, `upper`, and `step`:

- the ticks are $t_k = \mathrm{lower} + k\,\mathrm{step}$ for
  $k = 0, 1, \dots, \lfloor (\mathrm{upper} - \mathrm{lower}) /
  \mathrm{step} \rfloor$;
- the *tolerance* is $\mathrm{step} / 10$ rounded down to its first
  non-zero digit (so $0.025 \to 0.02$ and $0.05 \to 0.05$), and a tick
  closer to zero than that is written as exactly `0.0`;
- the axis limits are $(\mathrm{lower} - \mathrm{tol}, \mathrm{upper} +
  \mathrm{tol})$, a frame that keeps the end ticks off the axes' edges;
- the labels are the ticks rounded to `decimals` places, or hidden with
  `labels=False`; further keyword arguments, such as `fontsize`, go to the
  labels.

A step that is not positive and smaller than the span, or `lower >= upper`,
raises `ValueError`. The ticks are not snapped to multiples of a power of
ten; choose `lower`, `upper`, and `step` with that in mind.

```python
>>> import matplotlib.pyplot as plt
>>> from polyphys.viz.helpers import set_xticks
>>> fig, ax = plt.subplots()
>>> set_xticks(ax, -0.2, 0.2, 0.1, decimals=1)
array([-0.2, -0.1,  0. ,  0.1,  0.2])
>>> [label.get_text() for label in ax.get_xticklabels()]
['-0.2', '-0.1', '0.0', '0.1', '0.2']
>>> [round(float(limit), 6) for limit in ax.get_xlim()]
[-0.21, 0.21]
>>> plt.close(fig)

```

Here `step / 10 = 0.01`, already one significant digit, so the frame is
$0.01$ on each side (the rounding in the last line only hides float
round-off in the printed limits).

## Legends

A multi-panel figure often needs one legend for its colors and another for
its line styles, whatever artists the plotting call produced. The helpers in
{mod}`polyphys.viz.helpers` make *proxy* handles, artists that draw nothing
on the axes and exist only to sit in a legend:

- {func}`~polyphys.viz.helpers.color_handles` gives one filled patch per
  label, for a legend keyed by color;
- {func}`~polyphys.viz.helpers.line_handles` gives one empty line per
  label, with a line style, width, and marker each given once for all
  labels or once per label, for a legend keyed by style.

{func}`~polyphys.viz.helpers.move_legend_to_axes` removes an existing
legend, such as the one a seaborn facet grid draws beside its panels, and
re-creates it, with the same entries and title, inside the panels you name.
{func}`~polyphys.viz.helpers.band_under_lines` shades a band, such as a
confidence interval, under each line already on an axes, in that line's
color, from the lower and upper edges you pass for each line.

```python
>>> from polyphys.viz.helpers import color_handles, line_handles
>>> from polyphys.viz.helpers import move_legend_to_axes
>>> fig, (left, right) = plt.subplots(ncols=2)
>>> _ = left.plot([0, 1], [0, 1], label="rise")
>>> outer = fig.legend(title="curves")
>>> inner = move_legend_to_axes(outer, [right], ["upper left"])
>>> fig.legends
[]
>>> right.get_legend().get_title().get_text()
'curves'
>>> styles = line_handles(["ring", "linear"], linestyles=["-", "--"],
...                       markers="o")
>>> [(h.get_label(), h.get_linestyle()) for h in styles]
[('ring', '-'), ('linear', '--')]
>>> [h.get_label() for h in color_handles(["bug", "crowder"],
...                                       ["#4477AA", "#EE6677"])]
['bug', 'crowder']
>>> plt.close(fig)

```

## What was not ported

The research code's `visualize` package also held its study tables (labels,
titles, and per-project settings), fifteen figure scripts tied to the
research code's column names, four end-to-end-distribution plotters that
nothing called, and seven helpers and seven palettes with no callers. None
of that is in PolyPhys: the tables and figures are study code, for
each study's folder when its figures are ported, and the rest had no
use. The one bug the ported helpers carried, a legend attribute that
matplotlib 3.9 removed, is fixed with its replacement, which matplotlib 3.7
introduced. The extra requires `matplotlib>=3.8.4`.
