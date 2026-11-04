"""
The *TransFociCub* study: ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec`, ``PROBE`` its
:class:`~polyphys.manage.projects.ProbeSpec`, and ``ANALYSIS`` its
:class:`~polyphys.analyze.pipeline.AnalysisTable` per particle group.
"""
from .spec import ANALYSIS, PROBE, SPEC

__all__ = ["ANALYSIS", "PROBE", "SPEC"]
