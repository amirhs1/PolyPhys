"""
The *TransFociCyl* study: a ring of small and large monomers among crowders in
a cylinder that is periodic along its axis. ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.manage.projects.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.analyze.pipeline.AnalysisTable` per particle group.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk``, and from its analysis tables (SciOps
``sciops/analyze/analyzer.py``, ``ANALYSIS_DETAILS_NUCLEOID`` and
``ANALYSIS_DETAILS_ALL``, ``'TransFociCylWhole'``), with each artifact
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
import numpy as np

from polyphys.manage.parser import ParserBase, TransFociCyl
from polyphys.manage.projects import ProbeSpec, ProjectSpec

SPEC = ProjectSpec(
    parser=TransFociCyl,
    group='bug',
    divisor=0.025,
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'nmon_small',
                    'nmon_large', 'dmon_large', 'dcyl', 'dcrowd',
                    'phi_bulk_c'),
    species_sizes={'Mon': 'dmon_small', 'Crd': 'dcrowd',
                   'Foci': 'dmon_large'},
    # 'transSizeTMon' is not in the research registry's list, although the
    # research probe writes it and the analysis table below reads it; the
    # maintainer added it on 2026-11-03.
    tseries_properties=('asphericityTMon', 'fsdTMon', 'gyrTMon',
                        'shapeTMon', 'transSizeTMon'),
)


def _r_bin(info: ParserBase) -> float:
    """A tenth of the smaller of the small-monomer and crowder diameters:
    the bin size along r, x, and y (SciOps ``trans_foci_cyl_all``,
    ``bin_edges``)."""
    return 0.1 * min(float(getattr(info, 'dmon_small')),
                     float(getattr(info, 'dcrowd')))


def _z_bin(info: ParserBase) -> float:
    """Half the smaller of the small-monomer and crowder diameters: the bin
    size along z (SciOps ``trans_foci_cyl_all``, ``bin_edges``)."""
    return 0.5 * min(float(getattr(info, 'dmon_small')),
                     float(getattr(info, 'dcrowd')))


def _contact_cutoff(info: ParserBase) -> float:
    """The large-monomer diameter plus the crowder diameter: two foci are in
    direct contact when a crowder no longer fits between them (SciOps
    ``trans_foci_cyl_bug``, ``cluster_cutoff``)."""
    return float(getattr(info, 'dmon_large')) + float(getattr(info, 'dcrowd'))


# The probe of the polymer group: the selections, the chain-size observables
# of its per-frame loop, and the unwrap convention of its asphericity (SciOps
# ``sciops/probe/prober.py``, lines 968-998); the foci pair distances of the
# same loop; and the histograms of the all group (``trans_foci_cyl_all``,
# lines 1084-1140): its selections, bin sizes, and species. The research
# all-group probe selected the crowders, the small monomers, and the foci;
# 'Mon', the whole ring, is added because the analysis table below reads
# its histograms. The block length and the recentering are the maintainer's
# decisions of 2026-11-04.
PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'Foci': 'type 2', 'Crd': 'resid 0',
                'Dna': 'type 1'},
    chain='Mon',
    observables=('asphericityT', 'fsdT', 'gyrT', 'shapeT', 'transSizeT',
                 'principalT'),
    unwrap=False,
    bin_sizes={'r': _r_bin, 'z': _z_bin, 'theta': np.pi / 36, 'x': _r_bin,
               'y': _r_bin},
    histogram_species=('Crd', 'Mon', 'Dna', 'Foci'),
    block_frames=100,
    recentre='half',
    pair_species='Foci',
    contact_cutoff=_contact_cutoff,
    minimum_image=True,
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
        tseries('fsdT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('gyrT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('transSizeT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('shapeT', 'Mon', post_hook=Autocorrelation(NLAGS)),
        tseries('asphericityT', 'Mon', post_hook=Autocorrelation(NLAGS)),
    )),
    AnalysisTable(group='all', lineage='segment', stamps=False, properties=(
        *(histogram(direction, species, post_hook=LocalDensities())
          for direction in ('r', 'z')
          for species in ('Crd', 'Mon', 'Foci')),
        histogram('r', 'Dna'),
        histogram('z', 'Dna'),
        *(histogram('theta', species)
          for species in ('Crd', 'Mon', 'Dna', 'Foci')),
        *(histogram_2d(plane, species)
          for plane in ('xy', 'xz', 'yz')
          for species in ('Crd', 'Mon', 'Dna', 'Foci')),
        edge('x'),
        edge('y'),
        edge('z'),
    )),
)
