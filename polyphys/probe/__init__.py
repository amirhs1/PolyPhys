"""
The probe layer: per-frame observables read from LAMMPS trajectories with
MDAnalysis, which the ``mda`` extra installs
(``python -m pip install 'polyphys[mda]'``). The rest of PolyPhys does not
need it, so this package is imported on demand and raises :exc:`ImportError`
with that hint when MDAnalysis is missing.

.. versionadded:: 0.7.0
"""
try:
    import MDAnalysis  # noqa: F401
except ImportError as error:  # pragma: no cover - depends on the install
    raise ImportError(
        "polyphys.probe needs MDAnalysis; install it with "
        "python -m pip install 'polyphys[mda]'.") from error

from .binders import BinderStatistics
from .chain_size import ChainSize
from .core import make_universe, write_stamps

__all__ = ['BinderStatistics', 'ChainSize', 'make_universe', 'write_stamps']
