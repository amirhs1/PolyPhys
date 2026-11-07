"""
Tests for the examples of the documentation pages in ``docs/source/``.

pytest collects doctests from ``polyphys``, ``README.md``, and ``projects``
only, so the ``>>>`` examples of a documentation page run here, one test
per page, through :func:`doctest.testfile`. The pages compute their expected
values by hand in the text beside each example.
"""
import doctest
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

import polyphys


@pytest.fixture
def docs_source() -> Path:
    """The ``docs/source/`` directory of the checkout under test."""
    return Path(__file__).resolve().parents[2] / "docs" / "source"


@pytest.mark.parametrize("page", [
    pytest.param("quickstart.md", id="quickstart"),
    pytest.param("lineage.md", id="lineage"),
    pytest.param("studies.md", id="studies"),
    pytest.param("usage.md", id="usage"),
    pytest.param("viz.md", id="viz"),
])
def test_docs_page_examples_run(docs_source: Path, page: str) -> None:
    """Every example of the page gives the output the page shows."""
    path = docs_source / page

    result = doctest.testfile(str(path), module_relative=False,
                              encoding="utf-8", verbose=False)

    assert result.attempted > 0, f"{page} has no examples"
    assert result.failed == 0, f"{result.failed} example(s) of {page} failed"


def test_mocked_modules_import_with_their_annotations_evaluated():
    """The probe and viz modules import, and every annotation in them
    evaluates, with MDAnalysis and matplotlib replaced by the mocks of the
    docs build."""
    pytest.importorskip("sphinx.ext.autodoc.mock")
    script = textwrap.dedent("""\
        import importlib
        import inspect
        import pkgutil

        from sphinx.ext.autodoc.mock import mock

        with mock(['MDAnalysis', 'matplotlib']):
            for package in ('polyphys.probe', 'polyphys.viz'):
                names = [package] + [
                    f'{package}.{info.name}' for info in pkgutil.iter_modules(
                        importlib.import_module(package).__path__)]
                for name in names:
                    module = importlib.import_module(name)
                    for obj in vars(module).values():
                        if getattr(obj, '__module__', None) != name:
                            continue
                        members = [obj] + (list(vars(obj).values())
                                           if inspect.isclass(obj) else [])
                        for member in members:
                            if inspect.isfunction(member) or inspect.isclass(
                                    member):
                                inspect.get_annotations(member)
        """)
    # The child imports this checkout's package.
    env = dict(os.environ,
               PYTHONPATH=str(Path(polyphys.__file__).parents[1]))

    completed = subprocess.run(
        [sys.executable, '-c', script], env=env, capture_output=True,
        text=True, timeout=120, check=False)

    assert completed.returncode == 0, completed.stderr
