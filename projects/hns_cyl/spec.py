"""
The *HnsCyl* study: a ring polymer with H-NS proteins among crowders in a
cylinder that is periodic along its axis. ``SPEC`` is its
:class:`~polyphys.spec.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.spec.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.spec.AnalysisTable` per particle group.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk`` and ``rho_bulk_hns`` for ``rho_hns_bulk``, and from its
analysis tables (SciOps ``sciops/analyze/analyzer.py``,
``ANALYSIS_DETAILS_NUCLEOID`` and ``ANALYSIS_DETAILS_ALL``,
``'HnsCylWhole'``), with each artifact listed once.

Provisional: whether the study has molecular dynamics data is not yet
known.
"""
from polyphys.analyze.pipeline import Autocorrelation, LocalDensities
import numpy as np

from polyphys.manage.parser import ParserBase
from polyphys.spec import (
    AnalysisTable,
    BinderSpec,
    ProbeSpec,
    ProjectSpec,
    edge,
    histogram,
    histogram_2d,
    matrix_tseries,
    tseries,
)

from .parser import HnsCyl

SPEC = ProjectSpec(
    parser=HnsCyl,
    group='nucleoid',
    divisor=0.01,
    crowder_attr='phi_bulk_c',
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
# ``sciops/probe/prober.py``, lines 4268-4319); the histograms of the all
# group (``hns_cyl_all``, lines 4444-4500): its selections, bin sizes, and
# species, the H-NS cores as 'Hns', with the block length and the recentering
# the maintainer's decisions of 2026-11-04;
# patches (lines 4283-4383, ``cis_threshold = 4`` and the distance matrices
# saved as ``distMatTMonHnsPatch``). The cutoff is the published binding
# criterion, r <= r_cut,mp = 0.685 (thesis, Section 5.3 and Table 5.1),
# which the maintainer chose on 2026-11-04 over the research probe's
# ``round(2**(1/6) * (dmon + dhns_patch) / 2, 3) = 0.661`` with a strict
# comparison; the probe guide, "Choosing the cutoff", gives both.
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
    recenter='half',
    binders=BinderSpec(patches='HnsPatch', binders='HnsCore', cutoff=0.685,
                       inclusive=True, cis_threshold=4,
                       store_distances=True),
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
