"""
Analyze Module --- :mod:`polyphys.analyze`
==========================================


The `polyphys.analyze` module provides tools for analyzing molecular simulation
data: geometric and structural measurements of a particle group, simple
statistics, and the binning and histogram machinery used to build spatial
distributions.

Submodules
==========
- `measurer`: Geometric and structural measurements, summary statistics, and
  fixed-size binning with radial, axial, azimuthal, and planar histograms.
- `stats`: Autocorrelation-aware uncertainty of a time average over
  correlated frames, by Flyvbjerg-Petersen block averaging.
- `clusters`: Contact and cluster statistics of a particle group, pair
  separation densities, and the binding statistics of two-patch binders.
- `correlations`: The autocorrelation function of a time series with
  Bartlett's band, decay-model fits, and the bond cosine correlation of a
  chain.
- `sumrule`: The density sum rule on a space's ensemble-averaged density
  profiles: scaled, normalized, folded, and summed over species.
- `density`: The local number density and volume fraction of a species
  from a histogram of its bead centres, with the share of a finite bead
  spread over the bins it reaches.
"""

from polyphys.analyze import (
    clusters,
    correlations,
    density,
    measurer,
    stats,
    sumrule,
)

__all__ = [
    'clusters',
    'correlations',
    'density',
    'measurer',
    'stats',
    'sumrule',
]
