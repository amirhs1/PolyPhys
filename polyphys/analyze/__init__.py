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
"""

from polyphys.analyze import measurer, stats

__all__ = [
    'measurer',
    'stats',
]
