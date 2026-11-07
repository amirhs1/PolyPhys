"""
The parser of the *TransFociCub* study: a ring of small and large monomers
among crowders in a cube.

.. versionchanged:: 1.0.0
   Moved here from :mod:`polyphys.manage.parser`, which ships no study.
"""
from polyphys.manage.parser import KeywordSpec, ParserBase, Species

__all__ = ["TransFociCub"]


class TransFociCub(ParserBase):
    """
    Parse an artifact name of the *TransFociCub* project: a ring of small and
    large monomers among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
    - ``whole``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
    - ``ensemble_long``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#ac#nc#
    - ``space``: ns#nl#al#ac#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the ring, ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon_small : int
        Diameter of a small monomer, fixed at 1; not in the name.
    nmon_small : int
        Number of small monomers; keyword ``ns``.
    mmon_small : int
        Mass of a small monomer, fixed at 1; not in the name.
    dmon_large : float
        Diameter of a large monomer; keyword ``al``.
    nmon_large : int
        Number of large monomers; keyword ``nl``.
    mmon_large : float
        Mass of a large monomer; keyword ``ml``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    mcrowd : float
        Crowder mass, :math:`m_c = d_c^3`, at ``segment``, ``whole``, and
        ``ensemble_long``; ``-1``, a placeholder, at ``ensemble`` and
        ``space``.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m_small, phi_bulk_m_small : float
        Bulk number density and volume fraction of small monomers.
    rho_bulk_m_large, phi_bulk_m_large : float
        Bulk number density and volume fraction of large monomers.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of all monomers: the sums of
        the small and large values.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`~polyphys.manage.parser.ParserBase`.

    The mass density is uniform across all species, so a species whose mass
    is not in the name has the mass :math:`m_i = d_i^3`, where :math:`d_i` is
    its diameter.

    The simulation design values that the name does not carry, such as
    the interaction strengths, the bonds, and the thermostat, are listed
    in :doc:`/design`.

    Examples
    --------
    >>> artifact = TransFociCub(
    ... 'al5.0nl5ml125.0ns400ac1.0nc0l100.0dt0.005bdump2000adump5000.ring',
    ... 'ensemble_long',
    ... 'all'
    ... )
    >>> print(artifact.ncrowd)
    0
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
    _keywords = {
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ml': KeywordSpec('mmon_large', float),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                    'adump', 'ens', 'j'),
        'whole': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                  'adump', 'ens'),
        'ensemble_long': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt',
                          'bdump', 'adump'),
        'ensemble': ('ns', 'nl', 'al', 'ac', 'nc'),
        'space': ('ns', 'nl', 'al', 'ac'),
    }
    _defaults = {'dmon_small': 1, 'mmon_small': 1, 'mcrowd': -1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c', 'ncrowd', 'dcrowd'))
    _totals = {'m': ('m_small', 'm_large')}
    _masses = {'mcrowd': 'dcrowd'}
    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                    'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                  'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'rho_bulk_m_large',
                          'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                          'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }
