"""Tests for the whole-level builders in :mod:`polyphys.manage.organizer`.

Each test writes its input under a directory whose name contains a hyphen,
as every lineage directory name does. Expected whole names come from the
``SumRuleCyl`` parser applied to the bare filename, so each test checks that
the builder's result does not depend on the directory the file sits in.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from polyphys.manage.organizer import (
    whole_from_dist_mat_t,
    whole_from_file,
    whole_from_segments,
)
from polyphys.manage.parser import SumRuleCyl

GROUP = "bug"
WHOLE = (
    "N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000adump"
    "10000ens2"
)


def _hyphenated_dir(tmp_path: Path) -> Path:
    """Create and return a directory named like a lineage directory."""
    directory = tmp_path / "N2000-bug-wholeSim"
    directory.mkdir()
    return directory


def _expected_whole(filename: str, lineage: str) -> str:
    """Return the whole name the parser gives for a bare filename."""
    return getattr(SumRuleCyl(filename, lineage, GROUP), "whole")


def test_whole_from_segments_ignores_hyphenated_directory(tmp_path):
    """Segments under a hyphenated directory join into their parsed whole."""
    directory = _hyphenated_dir(tmp_path)
    paths = []
    for seg_id in (1, 2):
        path = directory / f"{WHOLE}.j0{seg_id}.{GROUP}.gyrTMon.npy"
        np.save(path, np.full(2, float(seg_id)))
        paths.append(str(path))

    wholes = whole_from_segments(
        "gyrTMon", paths, SumRuleCyl, GROUP, "tseries"
    )

    expected = _expected_whole(Path(paths[0]).stem, "segment")
    assert list(wholes) == [expected]
    np.testing.assert_array_equal(wholes[expected], [1.0, 1.0, 2.0, 2.0])


def test_whole_from_file_ignores_hyphenated_directory(tmp_path):
    """A whole under a hyphenated directory loads under its parsed name."""
    path = _hyphenated_dir(tmp_path) / f"{WHOLE}.{GROUP}.gyrTMon.npy"
    np.save(path, np.arange(3.0))

    wholes = whole_from_file([str(path)], SumRuleCyl, GROUP)

    expected = _expected_whole(path.stem, "whole")
    assert list(wholes) == [expected]
    np.testing.assert_array_equal(wholes[expected], np.arange(3.0))


def test_whole_from_dist_mat_t_ignores_hyphenated_directory(tmp_path):
    """A distance matrix under a hyphenated directory keys its parsed name."""
    path = _hyphenated_dir(tmp_path) / f"{WHOLE}.{GROUP}.distMatTFoci.npy"
    np.save(path, np.zeros((2, 2)))
    frame = pd.DataFrame({"x": [0.0]})

    def func(whole_path, whole_info):
        return frame, frame, frame

    freqs, rdfs, tseries = whole_from_dist_mat_t(
        [(str(path),)], SumRuleCyl, GROUP, func
    )

    expected = _expected_whole(path.stem, "whole")
    assert list(freqs) == list(rdfs) == list(tseries) == [expected]
