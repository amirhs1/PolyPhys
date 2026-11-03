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
"""

from polyphys.analyze import clusters, correlations, measurer, stats

__all__ = [
    'clusters',
    'correlations',
    'measurer',
    'stats',
]
