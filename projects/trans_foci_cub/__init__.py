"""
The *TransFociCub* study: ``SPEC`` is its
:class:`~polyphys.manage.projects.ProjectSpec` and ``ANALYSIS`` its
:class:`~polyphys.analyze.pipeline.AnalysisTable` per particle group.
"""
from .spec import ANALYSIS, SPEC

__all__ = ["ANALYSIS", "SPEC"]
