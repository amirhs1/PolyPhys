"""
The *HnsCub* study: a ring polymer with H-NS proteins among crowders in a
periodic cube. ``SPEC`` is its :class:`~polyphys.manage.projects.ProjectSpec`.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk`` and ``rho_bulk_hns`` for ``rho_hns_bulk``.
"""
from polyphys.manage.parser import HnsCub
from polyphys.manage.projects import ProjectSpec

SPEC = ProjectSpec(
    parser=HnsCub,
    group='nucleoid',
    divisor=0.01,
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
