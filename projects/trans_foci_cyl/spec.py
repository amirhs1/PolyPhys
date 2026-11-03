"""
The *TransFociCyl* study: a ring of small and large monomers among crowders in
a cylinder that is periodic along its axis. ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec`.

Values from the research code's registry (SciOps ``sciops/api.py``,
``PROJECTS_DETAILS``), with the PolyPhys attribute names ``phi_bulk_c`` for
``phi_c_bulk``.
"""
from polyphys.manage.parser import TransFociCyl
from polyphys.manage.projects import ProjectSpec

SPEC = ProjectSpec(
    parser=TransFociCyl,
    group='bug',
    divisor=0.025,
    physical_attrs=('space', 'ensemble_long', 'ensemble', 'nmon_small',
                    'nmon_large', 'dmon_large', 'dcyl', 'dcrowd',
                    'phi_bulk_c'),
    species_sizes={'Mon': 'dmon_small', 'Crd': 'dcrowd',
                   'Foci': 'dmon_large'},
    tseries_properties=('asphericityTMon', 'fsdTMon', 'gyrTMon',
                        'shapeTMon'),
)
