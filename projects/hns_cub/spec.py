"""
The *HnsCub* study: a ring polymer with H-NS proteins among crowders in a
periodic cube. ``SPEC`` is its
:class:`~polyphys.spec.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.spec.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.spec.AnalysisTable` per particle group.

Values from the research code's registry, with the PolyPhys attribute names
``phi_bulk_c`` for ``phi_c_bulk`` and ``rho_bulk_hns`` for ``rho_hns_bulk``,
and from its analysis tables, with each artifact listed once.
"""
from polyphys.analyze.pipeline import (
    Autocorrelation,
    LocalDensities,
    TimeAverage,
)
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

from .parser import HnsCub

SPEC = ProjectSpec(
    parser=HnsCub,
    group='nucleoid',
    divisor=0.01,
    crowder_attr='phi_bulk_c',
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'bend_mm',
                    'eps_hc', 'nmon', 'nhns', 'dcrowd', 'phi_bulk_c',
                    'rho_bulk_hns'),
    species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd', 'Hns': 'dhns'},
    tseries_properties=('asphericityTMon', 'gyrTMon', 'shapeTMon',
                        'bondLengthVecMon', 'nBoundTHnsPatch',
                        'nFreeTHnsPatch', 'nEngagedTHnsPatch',
                        'nFreeTHnsCore', 'nBridgeTHnsCore',
                        'nDangleTHnsCore', 'nCisTHnsCore', 'nTransTHnsCore'),
)


def _r_bin(info: ParserBase) -> float:
    """A tenth of the smallest of the monomer, H-NS core, and crowder
    diameters: the bin size along r (the research code's
    ``hns_cub_all``, ``bin_edges``)."""
    return 0.1 * min(float(getattr(info, 'dmon')),
                     float(getattr(info, 'dhns')),
                     float(getattr(info, 'dcrowd')))


def _l_bin(info: ParserBase) -> float:
    """Half the smaller of the monomer and crowder diameters: the bin size
    along x, y, and z (the research code's ``hns_cub_all``,
    ``bin_edges``, ``lEdge``)."""
    return 0.5 * min(float(getattr(info, 'dmon')),
                     float(getattr(info, 'dcrowd')))


# The probe of the polymer group: the selections, the chain-size observables
# of its per-frame loop, and the unwrap convention of its asphericity (the
# research code's ``hns_cub_nucleoid``); the histograms of the all group
# (``hns_cub_all``): its selections, bin sizes, and species, the H-NS cores
# as 'Hns', with the block length and the recentering of the probe guide
# ("The spatial histograms" and "Recentering a frame"); and the binding of
# the H-NS patches (``hns_cub_nucleoid``, ``cis_threshold = 4`` and the
# distance matrices saved as ``distMatTMonHnsPatch``). The cutoff is the
# published binding criterion, r <= r_cut,mp = 0.685 (thesis, Section 5.3
# and Table 5.1), which the maintainer chose on 2026-11-04 over the research
# probe's ``round(2**(1/6) * (dmon + dhns_patch) / 2, 3) = 0.661`` with a
# strict comparison; the probe guide, "Choosing the cutoff", gives both.
# The cis threshold is 3, the maintainer's definition (thesis, chapter 5: a
# bridge with at most two monomers between its sites, at most 3 bonds, is
# cis), which the maintainer chose on 2026-11-08 over the research probe's 4;
# the probe guide, "Choosing the cis threshold", gives both and why.
PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'HnsPatch': 'type 2', 'Crd': 'resid 0',
                'Hns': 'type 3'},
    chain='Mon',
    observables=('asphericityT', 'gyrT', 'shapeT', 'principalT'),
    unwrap=True,
    bin_sizes={'r': _r_bin, 'x': _l_bin, 'y': _l_bin, 'z': _l_bin},
    histogram_species=('Crd', 'Mon', 'Hns'),
    block_frames=100,
    recenter='half',
    binders=BinderSpec(patches='HnsPatch', binders='HnsCore', cutoff=0.685,
                       inclusive=True, cis_threshold=3,
                       store_distances=True),
)

# The largest lag of the autocorrelation functions, in frames (the research
# code's ``analyze_measures``, ``nlags``).
NLAGS = 7000

ANALYSIS = (
    AnalysisTable(group='nucleoid', lineage='whole', stamps=True, properties=(
        matrix_tseries('principalT', 'Mon'),
        # The per-frame bond and cluster-size histograms of the H-NS cores,
        # averaged over the frames; a cluster histogram's column 0 is the
        # count of empty clusters.
        matrix_tseries('bondsHistT', 'HnsCore',
                       post_hook=TimeAverage('bondsHist')),
        matrix_tseries('bondsHistDangleT', 'HnsCore',
                       post_hook=TimeAverage('bondsHistDangle')),
        matrix_tseries('bondsHistBridgeT', 'HnsCore',
                       post_hook=TimeAverage('bondsHistBridge')),
        matrix_tseries('bondsHistDirDepT', 'HnsCore',
                       post_hook=TimeAverage('bondsHistDirDep')),
        matrix_tseries('clustersHistT', 'HnsCore',
                       post_hook=TimeAverage('clustersHist',
                                             drop_first_column=True)),
        matrix_tseries('clustersHistDangleT', 'HnsCore',
                       post_hook=TimeAverage('clustersHistDangle',
                                             drop_first_column=True)),
        matrix_tseries('clustersHistBridgeT', 'HnsCore',
                       post_hook=TimeAverage('clustersHistBridge',
                                             drop_first_column=True)),
        matrix_tseries('clustersHistDirDepT', 'HnsCore',
                       post_hook=TimeAverage('clustersHistDirDep',
                                             drop_first_column=True)),
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
        tseries('gyrT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('shapeT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('asphericityT', 'Mon', post_hook=Autocorrelation(NLAGS)),
    )),
    AnalysisTable(group='all', lineage='segment', stamps=False, properties=(
        histogram('r', 'Crd', post_hook=LocalDensities()),
        histogram('r', 'Mon', post_hook=LocalDensities()),
        histogram('r', 'Hns', post_hook=LocalDensities()),
        *(histogram_2d(plane, species)
          for plane in ('xy', 'xz', 'yz')
          for species in ('Crd', 'Mon', 'Hns')),
        edge('x'),
        edge('y'),
        edge('z'),
    )),
)
