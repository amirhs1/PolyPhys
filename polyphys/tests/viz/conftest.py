"""
Fixtures for the tests of :mod:`polyphys.viz`.

The viz layer needs matplotlib, which the ``viz`` extra installs. The test
modules of the matplotlib-backed modules call
``pytest.importorskip("matplotlib")`` first, and the fixtures here do the
same, so ``test_colors.py``, whose module imports no plotting library,
runs without the extra. Every figure is drawn with the Agg backend, which
needs no display, and closed by the fixture that made it. All data is tier
A: a few hand-written points.
"""
from typing import Any, Iterator

import pytest


@pytest.fixture
def figure() -> Iterator[Any]:
    """An empty Agg figure, closed after the test."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure()
    yield fig
    plt.close(fig)


@pytest.fixture
def ax(figure: Any) -> Any:
    """One empty axes on the test's figure."""
    return figure.add_subplot()
