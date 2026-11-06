"""
Fixtures shared by the probe tests.

Every test here needs MDAnalysis, which the ``mda`` extra installs; the test
modules call ``pytest.importorskip("MDAnalysis")``, and the fixtures import
the probe layer inside their factories so that collecting this folder
without MDAnalysis is harmless.

Fixtures
--------
make_chain_files
    A LAMMPS data file and dump of one linear chain, from an array of
    per-frame positions.
make_chain_universe
    The same files opened with :func:`polyphys.probe.make_universe`.
"""
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import numpy as np
import pytest

from polyphys.simulate.core import Box
from polyphys.simulate.data import DataFile, write_data


@pytest.fixture
def make_chain_files(
    tmp_path: Path, make_lammps_dump: Callable[..., Path]
) -> Callable[..., tuple[Path, Path]]:
    """Factory: write a ``bond``-style data file of one linear chain, with
    masses and a bond between neighbors, and a dump of its frames.

    Takes `positions` of shape ``(n_frames, n_atoms, 3)``, `box` as three
    ``(lo, hi)`` pairs, `boundary` as the dump's three flags, `masses` per
    atom type, `atom_types` per atom, and `name` for both files. Returns
    the data and dump paths.
    """

    def _make(
        positions: Any,
        *,
        box: Sequence[Sequence[float]] = ((0.0, 20.0),) * 3,
        boundary: str = "ff ff pp",
        masses: Sequence[float] = (1.0,),
        atom_types: Optional[Sequence[int]] = None,
        name: str = "chain",
    ) -> tuple[Path, Path]:
        frames = np.asarray(positions, dtype=float)
        n_atoms = frames.shape[1]
        types = np.array(atom_types if atom_types is not None
                         else [1] * n_atoms)
        data = DataFile(
            atom_style="bond",
            box=Box(lo=np.array([lo for lo, _ in box], dtype=float),
                    hi=np.array([hi for _, hi in box], dtype=float)),
            ids=np.arange(1, n_atoms + 1), types=types,
            positions=frames[0], n_atom_types=len(masses),
            masses=np.asarray(masses, dtype=float),
            mol_ids=np.ones(n_atoms, dtype=int),
            bonds=np.array([[1, i, i + 1] for i in range(1, n_atoms)]),
            n_bond_types=1)
        data_path = write_data(data, tmp_path / f"{name}.data")
        dump_path = make_lammps_dump(
            name=f"{name}.lammpstrj", positions=frames, box=box,
            boundary=boundary, atom_types=types.tolist())
        return data_path, dump_path

    return _make


@pytest.fixture
def make_chain_universe(
    make_chain_files: Callable[..., tuple[Path, Path]]
) -> Callable[..., Any]:
    """Factory: the universe of the files ``make_chain_files`` writes, with
    the same arguments."""

    def _make(positions: Any, **kwargs: Any) -> Any:
        from polyphys.probe import make_universe

        data_path, dump_path = make_chain_files(positions, **kwargs)
        return make_universe(data_path, dump_path)

    return _make
