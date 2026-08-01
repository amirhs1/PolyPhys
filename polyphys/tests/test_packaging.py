"""
Tests for :mod:`polyphys` distribution metadata.

These tests guard the packaging contract declared in ``pyproject.toml``
against defects that are invisible to the rest of the test suite because
they only surface in an *installed* copy of the package: an entry point
naming a module that does not exist, or a malformed version string.
"""
import re
from importlib.metadata import PackageNotFoundError, distribution

import pytest

from polyphys.__version__ import __version__

#: ``MAJOR.MINOR.PATCH`` with an optional pre-release/build suffix.
SEMVER_RE = re.compile(
    r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)"
    r"(?:-[0-9A-Za-z.-]+)?(?:\+[0-9A-Za-z.-]+)?$"
)


def _polyphys_distribution():
    """
    Return the installed ``polyphys`` distribution, or skip the test.

    Skipping keeps the suite runnable from a bare source checkout, where
    no distribution metadata exists to inspect.
    """
    try:
        return distribution("polyphys")
    except PackageNotFoundError:  # pragma: no cover - depends on install
        pytest.skip("polyphys is not installed; no metadata to inspect")


def test_version_is_valid_semver() -> None:
    """``__version__`` must be a valid semantic version string.

    Guards the release process: ``pyproject.toml`` reads this attribute
    dynamically, so a malformed value produces a broken distribution.
    """
    assert SEMVER_RE.match(__version__), (
        f"__version__ = {__version__!r} is not valid semantic versioning"
    )


def test_declared_console_scripts_are_importable() -> None:
    """Every declared console script must resolve to a real callable.

    Regression test: the distribution previously declared
    ``polyphys = polyphys.cli:main`` while ``polyphys/cli.py`` did not
    exist, so every install shipped a command that raised
    ``ModuleNotFoundError`` when run.
    """
    dist = _polyphys_distribution()
    scripts = [
        ep for ep in dist.entry_points if ep.group == "console_scripts"
    ]
    for entry_point in scripts:
        target = entry_point.load()
        assert callable(target), (
            f"console script {entry_point.name!r} points at "
            f"{entry_point.value!r}, which is not callable"
        )
