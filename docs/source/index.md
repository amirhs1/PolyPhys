# Welcome to PolyPhys's Documentation!

**PolyPhys** is a Python library for managing and analyzing data from large-scale
coarse-grained molecular dynamics simulations of polymers under macromolecular
crowding, focusing on computational modeling of bacterial chromosome organization.

It decodes physics-encoding simulation filenames into typed attributes, organizes
artifacts along a five-level lineage — `segment`, `whole`, `ensemble_long`,
`ensemble`, `space` — and reduces replicate ensembles into averaged measurements.

Start with [Installation](installation.md), then the [quickstart](quickstart.md),
which runs a small synthetic sweep from its segments to its `space` table; the
[usage](usage.md) page says where each task is explained.

## Where PolyPhys comes from

PolyPhys was migrated from the maintainer's private research repository,
SciOps. These pages, the docstrings, and the tests call it *the research
code*, and its probe functions *the research probe*. It is not public, and
nothing in PolyPhys needs it: where a page says how PolyPhys differs from
the research code, the difference is written out on that page. Some tests
take their expected values from it, written into the test files; the tests
that compare results with its output files need the maintainer's data and
are skipped without it, as the [Testing](testing.md) page explains. The
parsers and specs of the unpublished studies stay in the research code.

```{note}
This project is under active development.
```

```{toctree}
---
maxdepth: 2
caption: Contents:
---

installation
usage
quickstart
lineage
workflow
studies
statistics
clusters
correlations
density
pipeline
probe
simulate
viz
testing
design
decisions
stability
api
```

## Indices and Tables

- [Index](genindex)
- [Module Index](modindex)
- [Search](search)
