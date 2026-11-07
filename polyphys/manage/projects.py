"""\
==========================================================
:mod:`polyphys.manage.projects`
==========================================================

.. deprecated:: 0.9.0
   The records of this module moved to :mod:`polyphys.spec`, which imports
   neither MDAnalysis nor :mod:`polyphys.probe`. Importing
   :class:`~polyphys.spec.ProjectSpec`, :class:`~polyphys.spec.ProbeSpec`,
   :class:`~polyphys.spec.BinderSpec`, :data:`~polyphys.spec.PerWhole`,
   :data:`~polyphys.spec.RecenterT`,
   :data:`~polyphys.spec.CHAIN_SIZE_OBSERVABLES`, or
   :data:`~polyphys.spec.BINDER_SERIES` from here still works and warns with
   :exc:`DeprecationWarning`; this module is removed in 1.0.0.

Examples
========
>>> import warnings
>>> from polyphys import spec
>>> with warnings.catch_warnings(record=True) as caught:
...     warnings.simplefilter('always')
...     from polyphys.manage.projects import ProjectSpec
>>> ProjectSpec is spec.ProjectSpec, caught[0].category.__name__
(True, 'DeprecationWarning')
"""
import warnings

from polyphys import spec as _spec

# The names that moved to polyphys.spec in 0.9.0.
_MOVED_TO_SPEC = frozenset({
    'BINDER_SERIES', 'CHAIN_SIZE_OBSERVABLES', 'BinderSpec', 'PerWhole',
    'ProbeSpec', 'ProjectSpec', 'RecenterT'})


def __getattr__(name: str) -> object:
    """
    Serve a name that moved to :mod:`polyphys.spec`, with a warning.

    Raises
    ------
    AttributeError
        If `name` is not one of the moved names.
    """
    if name in _MOVED_TO_SPEC:
        warnings.warn(
            f"polyphys.manage.projects.{name} moved to polyphys.spec in "
            f"0.9.0; import it from there. The old path is removed in 1.0.0.",
            DeprecationWarning, stacklevel=2)
        return getattr(_spec, name)
    raise AttributeError(
        f"module 'polyphys.manage.projects' has no attribute {name!r}")
