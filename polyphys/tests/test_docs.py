"""
Tests for the examples of the documentation pages in ``docs/source/``.

pytest collects doctests from ``polyphys``, ``README.md``, and ``projects``
only, so the ``>>>`` examples of a documentation page run here, one test
per page, through :func:`doctest.testfile`. The pages compute their expected
values by hand in the text beside each example.
"""
import doctest
from pathlib import Path

import pytest


@pytest.fixture
def docs_source() -> Path:
    """The ``docs/source/`` directory of the checkout under test."""
    return Path(__file__).resolve().parents[2] / "docs" / "source"


@pytest.mark.parametrize("page", [
    pytest.param("quickstart.md", id="quickstart"),
    pytest.param("lineage.md", id="lineage"),
    pytest.param("studies.md", id="studies"),
    pytest.param("usage.md", id="usage"),
])
def test_docs_page_examples_run(docs_source: Path, page: str) -> None:
    """Every example of the page gives the output the page shows."""
    path = docs_source / page

    result = doctest.testfile(str(path), module_relative=False,
                              encoding="utf-8", verbose=False)

    assert result.attempted > 0, f"{page} has no examples"
    assert result.failed == 0, f"{result.failed} example(s) of {page} failed"
