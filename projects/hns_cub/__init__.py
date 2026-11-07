"""
The *HnsCub* study: ``SPEC`` is its
:class:`~polyphys.spec.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.spec.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.spec.AnalysisTable` per particle group, and
``HnsCub`` its parser.
"""
from .parser import HnsCub
from .spec import ANALYSIS, PROBE, SPEC

__all__ = ["ANALYSIS", "PROBE", "SPEC", "HnsCub"]
