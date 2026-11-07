"""
The parser of the *HnsCub* study: a ring polymer with H-NS proteins among
crowders in a cube.

.. versionchanged:: 1.0.0
   Moved here from :mod:`polyphys.manage.parser`, which ships no study.
"""
from polyphys.manage.parser import KeywordSpec, ParserBase, Species

__all__ = ["HnsCub"]


class HnsCub(ParserBase):
    """
    Parse an artifact name of the *HnsCub* project: a ring polymer with H-NS
    proteins among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: N#kbmm#nh#ac#l#epshc#nc#ens#.j#.ring
    - ``whole``: N#kbmm#nh#ac#l#epshc#nc#ens#.ring
    - ``ensemble_long``: N#kbmm#nh#ac#l#epshc#nc#
    - ``ensemble``: N#kbmm#nh#ac#epshc#nc#
    - ``space``: N#kbmm#nh#ac#epshc#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'nucleoid', 'all'}
        Particle group: ``'nucleoid'`` or ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon : int
        Monomer diameter, fixed at 1; not in the name.
    nmon : int
        Number of monomers; keyword ``N``.
    dhns : int
        Diameter of an H-NS protein, fixed at 1; not in the name.
    dhns_patch : float
        Diameter of the patch at each pole of an H-NS protein, fixed at
        0.178; not in the name.
    nhns : int
        Number of H-NS proteins; keyword ``nh``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    bend_mm : float
        Bending rigidity of the DNA monomers; keyword ``kbmm``.
    eps_hm : int
        Strength of the attractive LJ interaction between the H-NS poles and
        the monomers, fixed at 29; not in the name.
    eps_hc : float
        Value of the keyword ``epshc``, named for the strength of the LJ
        interaction between the H-NS cores and the crowders. The project's
        LAMMPS inputs do not read it: they fix that strength at 1 with a
        purely repulsive cutoff.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step, fixed at 0.005; not in the name.
    ndump : int
        Dump period of the ``nucleoid`` trajectory, fixed at 5000 time steps;
        not in the name.
    adump : int
        Dump period of the ``all`` trajectory, fixed at 10000 time steps; not
        in the name.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of monomers.
    rho_bulk_hns, phi_bulk_hns : float
        Bulk number density and volume fraction of H-NS proteins.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`~polyphys.manage.parser.ParserBase`.

    The simulation design values that the name does not carry, such as
    the interaction strengths, the bonds, and the thermostat, are listed
    in :doc:`/design`.

    Examples
    --------
    >>> artifact = HnsCub(
    ... 'N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring',
    ... 'whole',
    ... 'nucleoid'
    ... )
    >>> print(artifact.nhns)
    8
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']
    _dump_keywords = {'nucleoid': 'ndump', 'all': 'adump'}
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'kbmm': KeywordSpec('bend_mm', float),
        'nh': KeywordSpec('nhns', int),
        'ac': KeywordSpec('dcrowd', float),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'epshc': KeywordSpec('eps_hc', float),
        'nc': KeywordSpec('ncrowd', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'kbmm', 'nh', 'ac', 'l', 'epshc', 'nc', 'ens', 'j'),
        'whole': ('N', 'kbmm', 'nh', 'ac', 'l', 'epshc', 'nc', 'ens'),
        'ensemble_long': ('N', 'kbmm', 'nh', 'ac', 'l', 'epshc', 'nc'),
        'ensemble': ('N', 'kbmm', 'nh', 'ac', 'epshc', 'nc'),
        'space': ('N', 'kbmm', 'nh', 'ac', 'epshc'),
    }
    _defaults = {'dmon': 1, 'dhns': 1, 'dhns_patch': 0.178, 'dt': 0.005,
                 'ndump': 5000, 'adump': 10000, 'eps_hm': 29}
    _species = (Species('m', 'nmon', 'dmon'), Species('hns', 'nhns', 'dhns'),
                Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m',
                    'phi_bulk_c', 'rho_bulk_c', 'phi_bulk_hns',
                    'rho_bulk_hns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'whole': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                  'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt', 'ndump',
                  'adump', 'eps_hm'],
        'ensemble_long': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c', 'phi_bulk_hns',
                          'rho_bulk_hns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'ensemble': ['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'space': ['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']
    }
