"""
The *TransFociCub* study: a ring of small and large monomers among crowders in
a periodic cube. ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.manage.projects.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.analyze.pipeline.AnalysisTable` per particle group.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk``, and from its analysis tables (SciOps
``sciops/analyze/analyzer.py``, ``ANALYSIS_DETAILS_NUCLEOID`` and
``ANALYSIS_DETAILS_ALL``, ``'TransFociCubWhole'``), with each artifact
listed once.
"""
from polyphys.analyze.pipeline import (
    AnalysisTable,
    Autocorrelation,
    FociPairDistances,
    LocalDensities,
    TimeAverage,
    edge,
    half_contour,
    histogram,
    histogram_2d,
    matrix_tseries,
    tseries,
)
from polyphys.manage.parser import ParserBase, TransFociCub
from polyphys.manage.projects import ProbeSpec, ProjectSpec

SPEC = ProjectSpec(
    parser=TransFociCub,
    group='bug',
    divisor=0.025,
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'nmon_small',
                    'nmon_large', 'dmon_large', 'dcrowd', 'phi_bulk_c'),
    species_sizes={'Mon': 'dmon_small', 'Crd': 'dcrowd',
                   'Foci': 'dmon_large'},
    tseries_properties=('asphericityTMon', 'gyrTMon', 'shapeTMon'),
)

# The probe of the polymer group: the selections, the chain-size observables
# of its per-frame loop, and the unwrap convention of its asphericity (SciOps
# ``sciops/probe/prober.py``, lines 1708-1734).
PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'Foci': 'type 2'},
    chain='Mon',
    observables=('asphericityT', 'gyrT', 'shapeT', 'principalT'),
    unwrap=True,
)

# The largest lag of the autocorrelation functions, in frames (SciOps
# ``analyze_measures``, ``nlags``).
NLAGS = 7000


def _n_monomers(info: ParserBase) -> int:
    """The monomers on the ring: the small ones and the large ones."""
    return int(getattr(info, 'nmon_small') + getattr(info, 'nmon_large'))


def _bin_size(info: ParserBase) -> float:
    """Half the small-monomer diameter (SciOps ``whole_dist_mat_foci``)."""
    return 0.5 * float(getattr(info, 'dmon_small'))


def _r_max(info: ParserBase) -> float:
    """Half the ring contour at the FENE maximum ``1.5 sigma_ls``, with
    ``sigma_ls`` the large-small collision diameter: the largest bond."""
    sigma_ls = 0.5 * (float(getattr(info, 'dmon_large'))
                      + float(getattr(info, 'dmon_small')))
    return half_contour(_n_monomers(info), 1.5 * sigma_ls)


ANALYSIS = (
    AnalysisTable(group='bug', lineage='whole', stamps=True, properties=(
        matrix_tseries('bondsHistT', 'Foci',
                       post_hook=TimeAverage('bondsHist')),
        matrix_tseries('clustersHistT', 'Foci',
                       post_hook=TimeAverage('clustersHist',
                                             drop_first_column=True)),
        matrix_tseries('distMatT', 'Foci',
                       post_hook=FociPairDistances(
                           bin_size=_bin_size, r_max=_r_max,
                           n_monomers=_n_monomers)),
        matrix_tseries('principalT', 'Mon'),
        tseries('gyrT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('shapeT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('asphericityT', 'Mon', post_hook=Autocorrelation(NLAGS)),
    )),
    AnalysisTable(group='all', lineage='segment', stamps=False, properties=(
        histogram('r', 'Crd', post_hook=LocalDensities()),
        histogram('r', 'Mon', post_hook=LocalDensities()),
        histogram('r', 'Foci', post_hook=LocalDensities()),
        histogram('r', 'Dna'),
        *(histogram_2d(plane, species)
          for plane in ('xy', 'xz', 'yz')
          for species in ('Crd', 'Mon', 'Dna', 'Foci')),
        edge('x'),
        edge('y'),
        edge('z'),
    )),
)
