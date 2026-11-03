"""
The *HnsCyl* study: a ring polymer with H-NS proteins among crowders in a
cylinder that is periodic along its axis. ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec`.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk`` and ``rho_bulk_hns`` for ``rho_hns_bulk``.

Provisional: whether the study has molecular dynamics data is not yet
known.
"""
from polyphys.manage.parser import HnsCyl
from polyphys.manage.projects import ProjectSpec

SPEC = ProjectSpec(
    parser=HnsCyl,
    group='nucleoid',
    divisor=0.01,
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'dcyl', 'nmon',
                    'nhns', 'dcrowd', 'phi_bulk_c', 'rho_bulk_hns'),
    species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd', 'Hns': 'dhns'},
    tseries_properties=('asphericityTMon', 'fsdTMon', 'shapeTMon',
                        'gyrTMon', 'transSizeTMon', 'bondLengthVecMon',
                        'nBoundTHnsPatch', 'nFreeTHnsPatch',
                        'nEngagedTHnsPatch', 'nFreeTHnsCore',
                        'nBridgeTHnsCore', 'nDangleTHnsCore',
                        'nCisTHnsCore', 'nTransTHnsCore'),
    provisional=True,
)
