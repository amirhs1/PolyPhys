"""
The parser of the *TransFociCyl* study: a ring of small and large monomers
among crowders in a cylinder.

.. versionchanged:: 1.0.0
   Moved here from :mod:`polyphys.manage.parser`, which ships no study.
"""
from polyphys.manage.parser import KeywordSpec, ParserBase, Species

__all__ = ["TransFociCyl"]


class TransFociCyl(ParserBase):
    """
    Parse an artifact name of the *TransFociCyl* project: a ring of small and
    large monomers among crowders in a cylinder that is periodic along its
    axis.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.j#.ring
    - ``whole``: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.ring
    - ``ensemble_long``: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#D#ac#nc#
    - ``space``: ns#nl#al#D#ac#

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
    lcyl : float
        Length of the cylinder, :math:`Z = 2 l_z`, from the half-length
        :math:`l_z` (keyword ``lz``).
    dcyl : float
        Size of the cylinder, :math:`D = 2r - 1`, from the radius :math:`r`
        (keyword ``r``) in the long names; :math:`D` itself (keyword ``D``)
        in the short ones.
    epsilon_small : float
        Strength of the wall-small-monomer LJ interaction; keyword ``epss``.
    epsilon_large : float
        Strength of the wall-large-monomer LJ interaction; keyword ``epsl``.
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
    - The bulk values exist only at ``segment``, ``whole``, and
      ``ensemble_long``. The lineage names are documented in
      :class:`~polyphys.manage.parser.ParserBase`.

    - The mass density is uniform across all species, so a species whose
      mass is not in the name has the mass :math:`m_i = d_i^3`, where
      :math:`d_i` is its diameter.

    - The cylindrical wall is a layer of beads of size :math:`a_w = 1` whose
      centers lie on the cylinder of radius :math:`r` set in LAMMPS (keyword
      ``r``), so the cylinder size is :math:`D = 2r - a_w`, stored as `dcyl`.
      The keyword ``lz`` is the half-length of the cylinder, so its length is
      :math:`Z = 2 l_z`, stored as `lcyl`. The centers of particles of
      species :math:`\\alpha` with size :math:`a_\\alpha` are confined to a
      cylinder of size :math:`D_{\\mathrm{eff},\\alpha} = D - a_\\alpha`, and
      the bulk number density and volume fraction of that species both use
      the effective volume :math:`V_{\\mathrm{eff},\\alpha} = \\pi Z
      D_{\\mathrm{eff},\\alpha}^2 / 4`.

    - The simulation design values that the name does not carry, such as
      the interaction strengths, the bonds, and the thermostat, are
      listed in :doc:`/design`.

    Examples
    --------
    >>> artifact = TransFociCyl(
    ... 'ns500nl5al3.0D20.0ac2.0nc0',
    ... 'ensemble',
    ... 'bug'
    ... )
    >>> print(artifact.dcyl)
    20.0
    """
    _geometry = 'cylindrical'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
    _keywords = {
        'epss': KeywordSpec('epsilon_small', float),
        'epsl': KeywordSpec('epsilon_large', float),
        'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
        'D': KeywordSpec('dcyl', float),
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ml': KeywordSpec('mmon_large', float),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'lz': KeywordSpec('lcyl', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('epss', 'epsl', 'r', 'al', 'nl', 'ml', 'ns', 'ac', 'nc',
                    'lz', 'dt', 'bdump', 'adump', 'ens', 'j'),
        'whole': ('epss', 'epsl', 'r', 'al', 'nl', 'ml', 'ns', 'ac', 'nc',
                  'lz', 'dt', 'bdump', 'adump', 'ens'),
        'ensemble_long': ('epss', 'epsl', 'r', 'al', 'nl', 'ml', 'ns', 'ac',
                          'nc', 'lz', 'dt', 'bdump', 'adump'),
        'ensemble': ('ns', 'nl', 'al', 'D', 'ac', 'nc'),
        'space': ('ns', 'nl', 'al', 'D', 'ac'),
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
