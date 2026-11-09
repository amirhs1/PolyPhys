"""
The *TransFociCub* study: a ring of small and large monomers among crowders in
a periodic cube. ``SPEC`` is its
:class:`~polyphys.spec.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.spec.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.spec.AnalysisTable` per particle group.

Values from the research code's registry, with the PolyPhys attribute names
``phi_bulk_c`` for ``phi_c_bulk``, and from its analysis tables, with each
artifact listed once.
"""
from polyphys.analyze.pipeline import (
    Autocorrelation,
    FociPairDistances,
    LocalDensities,
    TimeAverage,
    half_contour,
)
from polyphys.manage.parser import ParserBase
from polyphys.spec import (
    AnalysisTable,
    ProbeSpec,
    ProjectSpec,
    edge,
    histogram,
    histogram_2d,
    matrix_tseries,
    tseries,
)

from .parser import TransFociCub

SPEC = ProjectSpec(
    parser=TransFociCub,
    group='bug',
    divisor=0.025,
    crowder_attr='phi_bulk_c',
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'nmon_small',
                    'nmon_large', 'dmon_large', 'dcrowd', 'phi_bulk_c'),
    species_sizes={'Mon': 'dmon_small', 'Crd': 'dcrowd',
                   'Foci': 'dmon_large'},
    tseries_properties=('asphericityTMon', 'gyrTMon', 'shapeTMon'),
)


def _r_bin(info: ParserBase) -> float:
    """A tenth of the smaller of the small-monomer and crowder diameters:
    the bin size along r (the research code's ``trans_foci_cub_all``,
    ``bin_edges``)."""
    return 0.1 * min(float(getattr(info, 'dmon_small')),
                     float(getattr(info, 'dcrowd')))


def _l_bin(info: ParserBase) -> float:
    """Half the smaller of the small-monomer and crowder diameters: the bin
    size along x, y, and z (the research code's ``trans_foci_cub_all``,
    ``bin_edges``)."""
    return 0.5 * min(float(getattr(info, 'dmon_small')),
                     float(getattr(info, 'dcrowd')))


def _contact_cutoff(info: ParserBase) -> float:
    """The large-monomer diameter plus the crowder diameter: two foci are in
    direct contact when a crowder no longer fits between them (the research
    code's ``trans_foci_cub_bug``, ``cluster_cutoff``)."""
    return float(getattr(info, 'dmon_large')) + float(getattr(info, 'dcrowd'))


# The probe of the polymer group: the selections, the chain-size observables
# of its per-frame loop, and the unwrap convention of its asphericity (the
# research code's ``trans_foci_cub_bug``); the foci pair distances of the
# same loop; and the histograms of the all group (``trans_foci_cub_all``):
# its selections, bin sizes, and species. The research all-group probe
# selected the crowders, the small monomers, and the foci; 'Mon', the whole
# ring, is added because the analysis table below reads its histograms. The
# block length and the recentering are the probe guide's ("The spatial
# histograms" and "Recentering a frame").
PROBE = ProbeSpec(
    project=SPEC,
    selections={'Mon': 'resid 1', 'Foci': 'type 2', 'Crd': 'resid 0',
                'Dna': 'type 1'},
    chain='Mon',
    observables=('asphericityT', 'gyrT', 'shapeT', 'principalT'),
    unwrap=True,
    bin_sizes={'r': _r_bin, 'x': _l_bin, 'y': _l_bin, 'z': _l_bin},
    histogram_species=('Crd', 'Mon', 'Dna', 'Foci'),
    block_frames=100,
    recenter='half',
    pair_species='Foci',
    contact_cutoff=_contact_cutoff,
    minimum_image=True,
)

# The largest lag of the autocorrelation functions, in frames (the research
# code's ``analyze_measures``, ``nlags``).
NLAGS = 7000


def _n_monomers(info: ParserBase) -> int:
    """The monomers on the ring: the small ones and the large ones."""
    return int(getattr(info, 'nmon_small') + getattr(info, 'nmon_large'))


def _bin_size(info: ParserBase) -> float:
    """Half the small-monomer diameter (the research code's
    ``whole_dist_mat_foci``)."""
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
