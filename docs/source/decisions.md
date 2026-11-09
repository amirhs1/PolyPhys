# Design decisions

This page records design decisions that shape PolyPhys, with the evidence
behind each and what it gives up.

## `manage` builds on no workflow framework

Decided on 2026-10-31 (amirhs1/PolyPhys#19).

**Decision.** `polyphys.manage` stays a small internal layer that turns the
names of existing simulation artifacts into lineage and feeds
`polyphys.analyze`. It does not build on [signac](https://github.com/glotzerlab/signac) or
another data-space framework, and it is not developed as a standalone data
manager. A study is added by declaring a parser class and a `ProjectSpec`
(see [Studies](studies.md)).

**The question.** signac manages a data space of simulations indexed by
*state points*, the key-value records that describe each job [^adorf2018].
The PolyPhys lineage `segment → whole → ensemble_long → ensemble → space` is
structurally a specialization of such a data space, so a thin signac adapter
was the alternative.

**Evidence.** From the signac 2.4.1 documentation and source:

- A signac project stores its data in a `workspace` directory with one
  subdirectory per job, named by the job id, the hash of the job's state
  point. A state point is unique, so replicas need a replica key.
- `Project.import_from(origin, schema)` brings an existing data space into a
  project by walking the *directories* under `origin` and mapping each to a
  state point, through a schema string such as `a_{a:int}/b_{b:int}` or a
  function of the path. It copies the data into the workspace by default;
  moving instead is possible, and the documentation warns that a failed import
  can then leave an inconsistent state.
- Linked views give the job directories human-readable names without copying.

PolyPhys starts from a different place. Its inputs are LAMMPS output that
was never organized by signac, and the physics is in the *filenames*: one
`ensemble_long` directory holds the files of many wholes, segments, and
properties, and the analysis stages (`ens`, `ensAvg`, `space`) are
directories of aggregates rather than state points. Under signac each whole
would become a job in a hashed directory, the lineage names would survive
only in linked views, and trajectories of tens of megabytes to gigabytes,
kept on the clusters that produced them, would be copied or moved into the
workspace.

**Given up.** signac's state-point queries, linked views, job documents, and
the workflow layer signac-flow. PolyPhys keeps its own small equivalent of a
state point: the attributes a parser reads from a name.

**Gained.** No new dependency, and the lineage, its names, and the files stay
exactly where the simulations wrote them.

**Not done.** The comparison rests on signac's documentation and source; no
prototype mapped a real space onto signac jobs. Revisit this decision if
PolyPhys starts to create simulation workspaces itself rather than read
existing ones.

## Each study is a folder outside the package

Decided on 2026-10-31 and 2026-11-02 (amirhs1/PolyPhys#76).

**Decision.** The package holds mechanisms and no study. `polyphys.spec`
(`polyphys.manage.projects` before 0.9.0) defines `ProjectSpec`, the record of what a study's analysis needs beyond
its parser, and checks it against the parser; PolyPhys holds no table of
studies. Each study is one folder, `projects/<study>/`, which
`pyproject.toml` excludes from the installed package and CI lints,
type-checks, and tests with it. Today the folder holds the study's
`ProjectSpec`, its `ProbeSpec`, and its `AnalysisTable` records, built
beside one another; its drivers and figure code follow. A correction to one
study's data, such as merging two neighboring crowder volume fractions, is code in
that study's folder with its reason, never an option of the library. The
published studies' folders are in this repository; unpublished studies use
the same layout in the research code (see "Where PolyPhys comes from" on
the [front page](index.md)).

**Why.** A study should run from its own repository without editing a table
inside the installed package (amirhs1/PolyPhys#18). A rounding step, a list
of time series, a selection string, or which attributes become columns are
choices of one study, not facts of the library; and the research code
carried the same study table in three copies, one per stage.

**The parsers.** Until 1.0.0 the parser classes of nine studies, such as
`HnsCub` and `TransFociCub`, shipped in `polyphys.manage.parser`. The
maintainer decided on 2026-11-02 that they move to their studies' folders,
and on 2026-11-07 how (amirhs1/PolyPhys#197): the parsers of the four
published studies are in `projects/<study>/parser.py`, those of the
unpublished studies are in the research code, and importing any of
the nine from `polyphys.manage.parser` raises `ImportError` that says where
it went, since a module in the installed package cannot import from
`projects/`. On 2026-11-09 the maintainer took the five unpublished
studies off that list and off the documentation, so their names are not
public; importing one now fails as for any name the module never
defined. The package's own tests and examples use made-up parsers in
`polyphys.tests.synthetic`.

[^adorf2018]: C. S. Adorf, P. M. Dodd, V. Ramasubramani, and S. C. Glotzer,
    "Simple data and workflow management with the signac framework,"
    *Comput. Mater. Sci.* **146**, 220–229 (2018),
    <https://doi.org/10.1016/j.commatsci.2018.01.035>.
