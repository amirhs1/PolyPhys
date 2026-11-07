"""\
==========================================================
:mod:`polyphys.tests.synthetic`
==========================================================

Synthetic parsers for the package's own tests and examples. PolyPhys ships
no study parser (see the "Studies" page); these classes declare made-up
studies, so that the tests and the examples of the library do not depend on
any study's naming scheme.

- :class:`ChainCyl`: a linear chain among crowders in a cylinder that is
  periodic along its axis.
- :class:`RingCub`: a ring of small and large monomers among crowders in a
  periodic cube.
- :class:`RingCubTwoCrowders`: a ring of small and large monomers among
  crowders of two sizes in a periodic cube.
- :class:`RingCubBinders`: a ring of monomers with binders and crowders in
  a periodic cube.

In each name, ``#`` stands for a value. The cylinder size is
:math:`D = 2r - 1` from the radius :math:`r` (keyword ``r``) of a wall of
beads of size 1, the cylinder length :math:`Z = 2 l_z` from the half-length
:math:`l_z` (keyword ``lz``), and the cube side :math:`L = 2l` from the
half-side :math:`l` (keyword ``l``).

Examples
========
>>> from polyphys.tests.synthetic import ChainCyl
>>> artifact = ChainCyl('N200D10.0ac2.0nc0', 'ensemble', 'all')
>>> artifact.nmon, artifact.dcyl, artifact.dcrowd
(200, 10.0, 2.0)
"""
from polyphys.manage.parser import KeywordSpec, ParserBase, Species

__all__ = ['ChainCyl', 'RingCub', 'RingCubBinders', 'RingCubTwoCrowders']

_MONOMER_BULK = ['phi_bulk_m', 'rho_bulk_m']
_RING_BULK = ['phi_bulk_m_small', 'rho_bulk_m_small', 'phi_bulk_m_large',
              'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m']


class ChainCyl(ParserBase):
    """
    Parse a name of the synthetic *ChainCyl* study: a linear chain among
    crowders in a cylinder that is periodic along its axis.

    - ``segment``: N#r#lz#ac#nc#dt#bdump#adump#ens#.j#
    - ``whole``: N#r#lz#ac#nc#dt#bdump#adump#ens#
    - ``ensemble_long``: N#r#lz#ac#nc#dt#bdump#adump#
    - ``ensemble``: N#D#ac#nc#
    - ``space``: N#D#ac#

    The monomer diameter `dmon` is 1 and not in the name. The bulk number
    densities and volume fractions of monomers and crowders exist at
    ``segment``, ``whole``, and ``ensemble_long``.

    Examples
    --------
    >>> artifact = ChainCyl(
    ...     'N200r10.5lz25.0ac2.0nc100dt0.005bdump1000adump5000ens1',
    ...     'whole', 'bug')
    >>> artifact.dcyl, artifact.lcyl, artifact.ensemble
    (20.0, 50.0, 'N200D20ac2nc100')
    """
    _geometry = 'cylindrical'
    _topology = 'linear'
    _groups = ['bug', 'all']
    _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
        'D': KeywordSpec('dcyl', float),
        'lz': KeywordSpec('lcyl', float, scale=2.0),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'r', 'lz', 'ac', 'nc', 'dt', 'bdump', 'adump',
                    'ens', 'j'),
        'whole': ('N', 'r', 'lz', 'ac', 'nc', 'dt', 'bdump', 'adump', 'ens'),
        'ensemble_long': ('N', 'r', 'lz', 'ac', 'nc', 'dt', 'bdump',
                          'adump'),
        'ensemble': ('N', 'D', 'ac', 'nc'),
        'space': ('N', 'D', 'ac'),
    }
    _defaults = {'dmon': 1}
    _species = (Species('m', 'nmon', 'dmon'), Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['dmon', *_MONOMER_BULK, 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon', *_MONOMER_BULK, 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon', *_MONOMER_BULK, 'phi_bulk_c',
                          'rho_bulk_c'],
        'ensemble': ['dmon'],
        'space': ['dmon'],
    }


class RingCub(ParserBase):
    """
    Parse a name of the synthetic *RingCub* study: a ring of small and large
    monomers among crowders in a periodic cube.

    - ``segment``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#
    - ``whole``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#
    - ``ensemble_long``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#ac#nc#
    - ``space``: ns#nl#al#ac#

    The small monomers have diameter and mass 1, and the crowders a mass
    set by their diameter; none is in the name.

    Examples
    --------
    >>> artifact = RingCub(
    ...     'al2.0nl2ml8.0ns20ac1.0nc100l10.0dt0.005bdump1000adump5000ens1',
    ...     'whole', 'bug')
    >>> artifact.nmon_small, artifact.nmon_large, artifact.lcube
    (20, 2, 20.0)
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
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', *_RING_BULK,
                    'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', *_RING_BULK,
                  'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small', 'mcrowd', *_RING_BULK,
                          'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd'],
    }


class RingCubTwoCrowders(ParserBase):
    """
    Parse a name of the synthetic *RingCubTwoCrowders* study: a ring of
    small and large monomers among crowders of two sizes in a periodic cube.

    - ``segment``: al#nl#ns#ac#ab#nc#nb#l#dt#bdump#adump#ens#.j#
    - ``whole``: al#nl#ns#ac#ab#nc#nb#l#dt#bdump#adump#ens#
    - ``ensemble_long``: al#nl#ns#ac#ab#nc#nb#l#dt#bdump#adump#
    - ``ensemble``: al#nl#ns#ac#ab#nc#nb#
    - ``space``: al#nl#ns#ac#ab#

    The crowders of type one have diameter `dcrowd_one` (keyword ``ac``)
    and number `ncrowd_one` (``nc``); those of type two `dcrowd_two`
    (``ab``) and `ncrowd_two` (``nb``). `phi_bulk_c` sums the two.

    Examples
    --------
    >>> artifact = RingCubTwoCrowders(
    ...     'al2.0nl2ns20ac1.0ab2.0nc100nb10l10.0dt0.005bdump1000adump5000',
    ...     'ensemble_long', 'all')
    >>> artifact.dcrowd_one, artifact.dcrowd_two, artifact.ncrowd_two
    (1.0, 2.0, 10)
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
    _keywords = {
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd_one', float),
        'ab': KeywordSpec('dcrowd_two', float),
        'nc': KeywordSpec('ncrowd_one', int),
        'nb': KeywordSpec('ncrowd_two', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('al', 'nl', 'ns', 'ac', 'ab', 'nc', 'nb', 'l', 'dt',
                    'bdump', 'adump', 'ens', 'j'),
        'whole': ('al', 'nl', 'ns', 'ac', 'ab', 'nc', 'nb', 'l', 'dt',
                  'bdump', 'adump', 'ens'),
        'ensemble_long': ('al', 'nl', 'ns', 'ac', 'ab', 'nc', 'nb', 'l',
                          'dt', 'bdump', 'adump'),
        'ensemble': ('al', 'nl', 'ns', 'ac', 'ab', 'nc', 'nb'),
        'space': ('al', 'nl', 'ns', 'ac', 'ab'),
    }
    _defaults = {'dmon_small': 1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c_one', 'ncrowd_one', 'dcrowd_one'),
                Species('c_two', 'ncrowd_two', 'dcrowd_two'))
    _totals = {'m': ('m_small', 'm_large'), 'c': ('c_one', 'c_two')}
    _project_attributes = {
        'segment': ['dmon_small', *_RING_BULK, 'phi_bulk_c_one',
                    'rho_bulk_c_one', 'phi_bulk_c_two', 'rho_bulk_c_two',
                    'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', *_RING_BULK, 'phi_bulk_c_one',
                  'rho_bulk_c_one', 'phi_bulk_c_two', 'rho_bulk_c_two',
                  'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', *_RING_BULK, 'phi_bulk_c_one',
                          'rho_bulk_c_one', 'phi_bulk_c_two',
                          'rho_bulk_c_two', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small'],
        'space': ['dmon_small'],
    }


class RingCubBinders(ParserBase):
    """
    Parse a name of the synthetic *RingCubBinders* study: a ring of
    monomers with binders and crowders in a periodic cube, whose particle
    groups are ``nucleoid`` (the ring and the binders) and ``all``.

    - ``segment``: N#nh#ac#nc#l#ens#.j#
    - ``whole``: N#nh#ac#nc#l#ens#
    - ``ensemble_long``: N#nh#ac#nc#l#
    - ``ensemble``: N#nh#ac#nc#
    - ``space``: N#nh#ac#

    The monomers and the binder cores have diameter 1, and the time step
    and dump periods are fixed; none is in the name.

    Examples
    --------
    >>> artifact = RingCubBinders('N200nh30ac2.0nc0l25.0ens1', 'whole',
    ...                           'nucleoid')
    >>> artifact.nmon, artifact.nhns, artifact.lcube
    (200, 30, 50.0)
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']
    _dump_keywords = {'nucleoid': 'ndump', 'all': 'adump'}
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'nh': KeywordSpec('nhns', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'nh', 'ac', 'nc', 'l', 'ens', 'j'),
        'whole': ('N', 'nh', 'ac', 'nc', 'l', 'ens'),
        'ensemble_long': ('N', 'nh', 'ac', 'nc', 'l'),
        'ensemble': ('N', 'nh', 'ac', 'nc'),
        'space': ('N', 'nh', 'ac'),
    }
    _defaults = {'dmon': 1, 'dhns': 1, 'dt': 0.005, 'ndump': 5000,
                 'adump': 10000}
    _species = (Species('m', 'nmon', 'dmon'), Species('hns', 'nhns', 'dhns'),
                Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['dmon', 'dhns', *_MONOMER_BULK, 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump'],
        'whole': ['dmon', 'dhns', *_MONOMER_BULK, 'phi_bulk_c', 'rho_bulk_c',
                  'phi_bulk_hns', 'rho_bulk_hns', 'dt', 'ndump', 'adump'],
        'ensemble_long': ['dmon', 'dhns', *_MONOMER_BULK, 'phi_bulk_c',
                          'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                          'ndump', 'adump'],
        'ensemble': ['dmon', 'dhns', 'dt', 'ndump', 'adump'],
        'space': ['dmon', 'dhns', 'dt', 'ndump', 'adump'],
    }
