"""
The parser of the *HnsCyl* study: a ring polymer with H-NS proteins among
crowders in a cylinder.

.. versionchanged:: 1.0.0
   Moved here from :mod:`polyphys.manage.parser`, which ships no study.
"""
from polyphys.manage.parser import KeywordSpec, ParserBase, Species

__all__ = ["HnsCyl"]


class HnsCyl(ParserBase):
    """
    Parse an artifact name of the *HnsCyl* project: a ring polymer with H-NS
    proteins among crowders in a cylinder that is periodic along its axis.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.j#.ring
    - ``whole``: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.ring
    - ``ensemble_long``: N#kbmm#r#nh#ac#lz#epshc#nc#
    - ``ensemble``: N#D#nh#ac#epshc#nc#
    - ``space``: N#D#nh#ac#epshc#

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
    lcyl : float
        Length of the cylinder, :math:`Z = 2 l_z`, from the half-length
        :math:`l_z` (keyword ``lz``).
    dcyl : float
        Size of the cylinder, :math:`D = 2r - 1`, from the radius :math:`r`
        (keyword ``r``) in the long names; :math:`D` itself (keyword ``D``)
        in the short ones.
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

    The cylindrical wall is a layer of beads of size :math:`a_w = 1` whose
    centers lie on the cylinder of radius :math:`r` set in LAMMPS (keyword
    ``r``), so the cylinder size is :math:`D = 2r - a_w`, stored as `dcyl`.
    The keyword ``lz`` is the half-length of the cylinder, so its length is
    :math:`Z = 2 l_z`, stored as `lcyl`. The centers of particles of
    species :math:`\\alpha` with size :math:`a_\\alpha` are confined to a
    cylinder of size :math:`D_{\\mathrm{eff},\\alpha} = D - a_\\alpha`, and
    the bulk number density and volume fraction of that species both use
    the effective volume :math:`V_{\\mathrm{eff},\\alpha} = \\pi Z
    D_{\\mathrm{eff},\\alpha}^2 / 4`.

    The simulation design values that the name does not carry, such as
    the interaction strengths, the bonds, and the thermostat, are listed
    in :doc:`/design`.

    Examples
    --------
    >>> artifact = HnsCyl(
    ... 'N200D20.0nh16ac1.0epshc1.0nc0',
    ... 'space',
    ... 'nucleoid'
    ... )
    >>> print(artifact.eps_hc)
    1.0
    """
    _geometry = 'cylindrical'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']
    _dump_keywords = {'nucleoid': 'ndump', 'all': 'adump'}
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'kbmm': KeywordSpec('bend_mm', float),
        'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
        'D': KeywordSpec('dcyl', float),
        'nh': KeywordSpec('nhns', int),
        'ac': KeywordSpec('dcrowd', float),
        'lz': KeywordSpec('lcyl', float, scale=2.0),
        'epshc': KeywordSpec('eps_hc', float),
        'nc': KeywordSpec('ncrowd', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'kbmm', 'r', 'nh', 'ac', 'lz', 'epshc', 'nc', 'ens',
                    'j'),
        'whole': ('N', 'kbmm', 'r', 'nh', 'ac', 'lz', 'epshc', 'nc', 'ens'),
        'ensemble_long': ('N', 'kbmm', 'r', 'nh', 'ac', 'lz', 'epshc', 'nc'),
        'ensemble': ('N', 'D', 'nh', 'ac', 'epshc', 'nc'),
        'space': ('N', 'D', 'nh', 'ac', 'epshc'),
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
