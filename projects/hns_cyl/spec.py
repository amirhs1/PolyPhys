"""
The *HnsCyl* study: a ring polymer with H-NS proteins among crowders in a
cylinder that is periodic along its axis. ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.manage.projects.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.analyze.pipeline.AnalysisTable` per particle group.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk`` and ``rho_bulk_hns`` for ``rho_hns_bulk``, and from its
analysis tables (SciOps ``sciops/analyze/analyzer.py``,
``ANALYSIS_DETAILS_NUCLEOID`` and ``ANALYSIS_DETAILS_ALL``,
``'HnsCylWhole'``), with each artifact listed once.

Provisional: whether the study has molecular dynamics data is not yet
known.
"""
from polyphys.analyze.pipeline import (
    AnalysisTable,
    Autocorrelation,
    LocalDensities,
    edge,
    histogram,
    histogram_2d,
    matrix_tseries,
    tseries,
)
import numpy as np

from polyphys.manage.parser import HnsCyl, ParserBase
from polyphys.manage.projects import ProbeSpec, ProjectSpec

SPEC = ProjectSpec(
    parser=HnsCyl,
    group='nucleoid',
    divisor=0.01,
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'dcyl', 'nmon',
                    'nhns', 'dcrowd', 'phi_bulk_c', 'rho_bulk_hns'),
    species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd', 'Hns': 'dhns'},
    tseries_properties=('asphericityTMon', 'fsdTMon', 'shapeTMon',
                        'gyrTMon', 'transSizeTMon', 'bondLengthVecMon',
                        'nBoundTHnsPatch', 'nFreeTHnsPatch',
                        'nEngagedTHnsPatch', 'nFreeTHnsCore',
                        'nBridgeTHnsCore', 'nDangleTHnsCore',
                        'nCisTHnsCore', 'nTransTHnsCore'),
    provisional=True,
)


def _r_bin(info: ParserBase) -> float:
    """A tenth of the smallest of the monomer, H-NS core, and crowder
    diameters: the bin size along r. The research probe left the H-NS core
    out of the minimum here, unlike ``hns_cub_all``; the maintainer called
    that a mistake on 2026-11-04 and both studies now agree."""
    return 0.1 * min(float(getattr(info, 'dmon')),
                     float(getattr(info, 'dhns')),
                     float(getattr(info, 'dcrowd')))


def _l_bin(info: ParserBase) -> float:
    """A tenth of the smaller of the monomer and crowder diameters: the bin
    size along x and y (SciOps ``hns_cyl_all``, ``bin_edges``, ``lEdge``)."""
    return 0.1 * min(float(getattr(info, 'dmon')),
                     float(getattr(info, 'dcrowd')))


def _z_bin(info: ParserBase) -> float:
    """Half the smaller of the monomer and crowder diameters: the bin size
    along z (SciOps ``hns_cyl_all``, ``bin_edges``)."""
    return 0.5 * min(float(getattr(info, 'dmon')),
                     float(getattr(info, 'dcrowd')))


# The probe of the polymer group: the selections, the chain-size observables
# of its per-frame loop, and the unwrap convention of its asphericity (SciOps
# ``sciops/probe/prober.py``, lines 4268-4319); and the histograms of the all
# group (``hns_cyl_all``, lines 4444-4500): its selections, bin sizes, and
# species, the H-NS cores as 'Hns'. The block length and the recentring are
# the maintainer's decisions of 2026-11-04.
PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'HnsPatch': 'type 2', 'Crd': 'resid 0',
                'Hns': 'type 3'},
    chain='Mon',
    observables=('asphericityT', 'fsdT', 'shapeT', 'gyrT', 'transSizeT',
                 'principalT'),
    unwrap=False,
    bin_sizes={'r': _r_bin, 'z': _z_bin, 'theta': np.pi / 36, 'x': _l_bin,
               'y': _l_bin},
    histogram_species=('Crd', 'Mon', 'Hns'),
    block_frames=100,
    recentre='half',
)

# The largest lag of the autocorrelation functions, in frames (SciOps
# ``analyze_measures``, ``nlags``).
NLAGS = 7000

ANALYSIS = (
    AnalysisTable(group='nucleoid', lineage='whole', stamps=True, properties=(
        matrix_tseries('principalT', 'Mon'),
        # Per-whole vectors reduced as time series.
        tseries('loopLengthHist', 'Mon'),
        tseries('bondCosineCorrVec', 'Mon'),
        tseries('bondLengthVec', 'Mon'),
        tseries('nBoundT', 'HnsPatch'),
        tseries('nFreeT', 'HnsPatch'),
        tseries('nEngagedT', 'HnsPatch'),
        tseries('nFreeT', 'HnsCore'),
        tseries('nBridgeT', 'HnsCore'),
        tseries('nDangleT', 'HnsCore'),
        tseries('nCisT', 'HnsCore'),
        tseries('nTransT', 'HnsCore'),
        tseries('fsdT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('gyrT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('transSizeT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('shapeT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('asphericityT', 'Mon', post_hook=Autocorrelation(NLAGS)),
    )),
    AnalysisTable(group='all', lineage='segment', stamps=False, properties=(
        *(histogram(direction, species, post_hook=LocalDensities())
          for direction in ('r', 'z')
          for species in ('Crd', 'Mon', 'Hns')),
        *(histogram('theta', species) for species in ('Crd', 'Mon', 'Hns')),
        *(histogram_2d(plane, species)
          for plane in ('xy', 'xz', 'yz')
          for species in ('Crd', 'Mon', 'Hns')),
        edge('x'),
        edge('y'),
        edge('z'),
    )),
)
