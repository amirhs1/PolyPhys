"""
The *TransFociCub* study: ``SPEC`` is its
:class:`~polyphys.spec.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.spec.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.spec.AnalysisTable` per particle group, and
``TransFociCub`` its parser.
"""
from .parser import TransFociCub
from .spec import ANALYSIS, PROBE, SPEC

__all__ = ["ANALYSIS", "PROBE", "SPEC", "TransFociCub"]
