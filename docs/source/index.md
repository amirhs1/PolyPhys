# Welcome to PolyPhys's Documentation!

**PolyPhys** is a Python library for managing and analysing data from large-scale
coarse-grained molecular dynamics simulations of polymers under macromolecular
crowding, focusing on computational modeling of bacterial chromosome organization.

It decodes physics-encoding simulation filenames into typed attributes, organizes
artifacts along a five-level lineage — `segment`, `whole`, `ensemble_long`,
`ensemble`, `space` — and reduces replicate ensembles into averaged measurements.

Check out the [usage](usage.md) section for further information, including how to
[install](#installation_target) the project, and the [quickstart](quickstart.md),
which runs a small synthetic sweep from its segments to its `space` table.

```{note}
This project is under active development.
```

```{toctree}
---
maxdepth: 2
caption: Contents:
---

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
testing
design
decisions
api
```

## Indices and Tables

- [Index](genindex)
- [Module Index](modindex)
- [Search](search)
