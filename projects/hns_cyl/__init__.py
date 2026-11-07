"""
The *HnsCyl* study: ``SPEC`` is its
:class:`~polyphys.spec.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.spec.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.spec.AnalysisTable` per particle group, and
``HnsCyl`` its parser.
"""
from .parser import HnsCyl
from .spec import ANALYSIS, PROBE, SPEC

__all__ = ["ANALYSIS", "PROBE", "SPEC", "HnsCyl"]
