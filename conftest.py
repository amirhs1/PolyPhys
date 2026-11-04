"""
Collection rules for the repository root.

``polyphys/probe`` imports MDAnalysis at module level, which the ``mda``
extra installs and the CI matrix jobs do not. With ``--doctest-modules``
pytest imports every module under ``polyphys/``, so that folder is left out
of collection when MDAnalysis is absent; its tests call
``pytest.importorskip("MDAnalysis")`` as well. This file is not part of the
installed package.
"""
import importlib.util
from pathlib import Path

import pytest

PROBE_DIR = Path(__file__).resolve().parent / "polyphys" / "probe"


def pytest_ignore_collect(collection_path: Path,
                          config: pytest.Config) -> bool | None:
    """Skip ``polyphys/probe`` when MDAnalysis is not installed."""
    if importlib.util.find_spec("MDAnalysis") is not None:
        return None
    if collection_path == PROBE_DIR or PROBE_DIR in collection_path.parents:
        return True
    return None
