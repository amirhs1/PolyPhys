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
study_spec
    The spec module, ``projects/<study>/spec.py``, of a study.
real_polymer_whole, real_window_whole, real_crowded_whole, \
real_foci_whole, real_hns_whole, real_nh0_whole
    The real wholes of the tier B and tier C tests, one test per whole (see
    "Real frames" below).
real_files
    The files of a real whole: its polymer dump and the topology it is
    opened with, and its ``all`` window and data file when they are tier B.
private_window_files
    The ``all`` window and data file of a crowded whole, tier C.
real_universe
    A real dump opened with :func:`polyphys.probe.make_universe`.
trajectory_arrays
    The positions, box lengths, and time steps of every frame of a
    universe, in one pass.
reference_distances
    Minimum-image distances between two sets of positions, frame by frame,
    from :func:`MDAnalysis.lib.distances.distance_array`.

Real frames
-----------
The tier B files are in ``polyphys/tests/data/probe/<whole>/``: for each of
22 wholes of the four live studies, the last 50 frames of the polymer group
dump (``<whole>.ring.<group>.last50.lammpstrj.gz``), and for the wholes
without crowders and the HnsCyl wholes at ``ac2 nc386``, the ``all`` frames
of the same time window (``<whole>.j<NN>.ring.all.window.lammpstrj.gz``)
and the whole's data file (``<whole>.ring.all.data.gz``). The maintainer
trimmed them from the studies' runs. The HnsCyl runs name their wholes
without ``epshc``; their copies here insert ``epshc1``, the value of their
H-NS core-crowder pair coefficients, so that the parser reads them.

The tier C files, the ``all`` windows and data files of the other crowded
wholes, stay in the maintainer's folder ``<POLYPHYS_DATA_DIR>/fixtures/``
under the runs' own names.

A polymer dump holds only the polymer group, so it cannot be opened with
its whole's full data file. It is opened with the data file of the whole
without crowders of the same family and replicate (:attr:`RealWhole.topology`),
whose atoms are the polymer group: the ring for TransFoci, the ring and the
H-NS proteins for Hns. The research code did the same. The polymer atoms
carry the IDs 1 to n in every system, with the crowders after them.
"""
import importlib
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Callable, Mapping, Optional, Sequence

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


# Real frames


@dataclass(frozen=True)
class RealWhole:
    """One real whole of the tier B and tier C tests.

    Attributes
    ----------
    id : str
        The test id, such as ``'HnsCub-nh30-nc10743-ens1'``.
    study : str
        Its study's folder in ``projects/``, such as ``'hns_cub'``.
    whole : str
        The whole's name without ``.ring``, as its parser reads it.
    topology : str
        The whole whose data file opens the polymer dump: the whole itself
        when it has no crowders, else the whole without crowders of its
        family and replicate.
    group : str
        The polymer group: ``'nucleoid'`` or ``'bug'``.
    n_mon, n_foci, n_hns, n_crowd : int
        The monomers, the large monomers among them, the H-NS proteins, and
        the crowders, from the name (keywords ``N`` or ``ns`` plus ``nl``,
        ``nl``, ``nh``, and ``nc``).
    window_tier : str
        ``'B'`` when the whole's ``all`` window is in the repository,
        ``'C'`` when it is the maintainer's, ``''`` for none.
    n_window_frames : int
        The frames of the ``all`` window: 25 for Hns, 20 for TransFoci.
    n_shared_frames : int
        The window's frames that are also polymer frames, by time step:
        all 25 for Hns, whose ``all`` dump period (10000) is twice the
        ``nucleoid`` one (5000); 10 for TransFoci, whose periods are 5000
        and 2000.
    source : str
        The whole's folder in the maintainer's fixtures folder: `whole`,
        but without ``epshc1`` for HnsCyl.
    """

    id: str
    study: str
    whole: str
    topology: str
    group: str
    n_mon: int
    n_foci: int
    n_hns: int
    n_crowd: int
    window_tier: str
    n_window_frames: int
    n_shared_frames: int
    source: str

    def polymer_counts(self) -> dict[str, int]:
        """The atoms of each species of the study's spec in the polymer
        dump: two patches per H-NS core (Jung, Sadeghi, and Ha, Sci. Rep.
        14, 139 (2024)), the small monomers and the foci of a TransFoci
        ring, and no crowders."""
        if self.group == 'nucleoid':
            return {'Mon': self.n_mon, 'HnsPatch': 2 * self.n_hns,
                    'Hns': self.n_hns, 'Crd': 0}
        return {'Mon': self.n_mon, 'Foci': self.n_foci,
                'Dna': self.n_mon - self.n_foci, 'Crd': 0}

    def window_counts(self) -> dict[str, int]:
        """The atoms of each species of the study's spec in the ``all``
        window: the polymer group plus the crowders."""
        return {**self.polymer_counts(), 'Crd': self.n_crowd}


def _real_wholes() -> list[RealWhole]:
    """The 22 real wholes: per study, two replicates without crowders and
    two at a bulk crowder volume fraction of about 0.36 to 0.38, and for
    the Hns studies the wholes without H-NS; HnsCyl adds crowders of size
    1."""
    hns_cub = 'N200kbmm2nh{nh}ac2l25epshc1nc{nc}ens{e}'
    hns_cyl = 'N200kbmm2r4.5nh{nh}ac{ac}lz75epshc1nc{nc}ens{e}'
    foci_cub = ('al6nl5ml216ns400ac1nc{nc}l{l}dt0.005bdump2000adump5000'
                'ens{e}')
    foci_cyl = ('epss5epsl5r10.5al6nl5ml216ns400ac1nc{nc}lz{l}dt0.005'
                'bdump2000adump5000ens{e}')
    wholes = []
    for e in (1, 2):
        for nh, nc, tier in ((0, 0, 'B'), (30, 0, 'B'), (30, 10743, 'C')):
            name = hns_cub.format(nh=nh, nc=nc, e=e)
            wholes.append(RealWhole(
                f"HnsCub-nh{nh}-nc{nc}-ens{e}", 'hns_cub', name,
                hns_cub.format(nh=nh, nc=0, e=e), 'nucleoid', 200, 0, nh,
                nc, tier, 25, 25, name))
        for nh, ac, nc, tier in ((0, 2, 0, 'B'), (30, 2, 0, 'B'),
                                 (30, 2, 386, 'B'), (30, 1, 4032, 'C')):
            name = hns_cyl.format(nh=nh, ac=ac, nc=nc, e=e)
            wholes.append(RealWhole(
                f"HnsCyl-nh{nh}-ac{ac}-nc{nc}-ens{e}", 'hns_cyl', name,
                hns_cyl.format(nh=nh, ac=2, nc=0, e=e), 'nucleoid', 200, 0,
                nh, nc, tier, 25, 25, name.replace('epshc1', '')))
        for template, study, label, uncrowded, crowded in (
                (foci_cub, 'trans_foci_cub', 'TransFociCub',
                 (0, '46.5'), (112776, '27')),
                (foci_cyl, 'trans_foci_cyl', 'TransFociCyl',
                 (0, '82.5'), (22230, '57'))):
            topology = template.format(nc=uncrowded[0], l=uncrowded[1], e=e)
            for (nc, length), tier in ((uncrowded, 'B'), (crowded, 'C')):
                name = template.format(nc=nc, l=length, e=e)
                wholes.append(RealWhole(
                    f"{label}-nc{nc}-ens{e}", study, name, topology, 'bug',
                    405, 5, 0, nc, tier, 20, 10, name))
    return sorted(wholes, key=lambda whole: whole.id)


def _params(keep: Callable[[RealWhole], bool]) -> list[Any]:
    """The real wholes that `keep` accepts, as pytest params."""
    return [pytest.param(whole, id=whole.id) for whole in _real_wholes()
            if keep(whole)]


@pytest.fixture(params=_params(lambda whole: True))
def real_polymer_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 22 real wholes, whose polymer dumps are all tier B."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.window_tier == 'B'))
def real_window_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each real whole whose ``all`` window and data file are tier B: the
    12 without crowders and the 2 HnsCyl wholes at ``ac2 nc386``."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.window_tier == 'C'))
def real_crowded_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 8 crowded wholes whose ``all`` window and data file are
    tier C."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.n_foci > 0))
def real_foci_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 8 TransFoci wholes."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.n_hns > 0))
def real_hns_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 10 Hns wholes with H-NS proteins."""
    return request.param


@pytest.fixture(params=_params(
    lambda whole: whole.group == 'nucleoid' and whole.n_hns == 0))
def real_nh0_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 4 Hns wholes without H-NS proteins."""
    return request.param


@pytest.fixture(scope='session')
def study_spec() -> Callable[[str], ModuleType]:
    """Factory: the spec module ``projects/<study>/spec.py`` of a study,
    with its ``SPEC`` and ``PROBE``. The study folders are not an installed
    package, so the repository root is put on the import path, as the
    pipeline's tier C test does."""
    root = str(Path(__file__).resolve().parents[3])
    if root not in sys.path:
        sys.path.insert(0, root)

    def _load(study: str) -> ModuleType:
        return importlib.import_module(f"projects.{study}.spec")

    return _load


@dataclass(frozen=True)
class RealFiles:
    """The files of a real whole; `window` and `data` are None when they
    are not in the repository."""

    polymer: Path
    topology: Path
    window: Optional[Path]
    data: Optional[Path]


@pytest.fixture
def real_files() -> Callable[[RealWhole], RealFiles]:
    """Factory: the tier B files of a real whole."""
    root = Path(__file__).resolve().parents[1] / 'data' / 'probe'

    def _files(whole: RealWhole) -> RealFiles:
        folder = root / whole.whole
        windows = sorted(folder.glob(f"{whole.whole}.j*.ring.all.window."
                                     "lammpstrj.gz"))
        data = folder / f"{whole.whole}.ring.all.data.gz"
        return RealFiles(
            polymer=folder / (f"{whole.whole}.ring.{whole.group}.last50."
                              "lammpstrj.gz"),
            topology=root / whole.topology / (f"{whole.topology}.ring.all."
                                              "data.gz"),
            window=windows[0] if windows else None,
            data=data if data.exists() else None)

    return _files


@pytest.fixture
def private_window_files(
    private_data_dir: Path,
) -> Callable[[RealWhole], tuple[Path, Path]]:
    """Factory: the tier C ``all`` window and data file of a crowded whole,
    under ``<POLYPHYS_DATA_DIR>/fixtures/<source>/``. Use it only in a test
    marked ``private_data``."""

    def _files(whole: RealWhole) -> tuple[Path, Path]:
        folder = private_data_dir / 'fixtures' / whole.source
        (window,) = sorted(folder.glob(f"{whole.source}.j*.ring.all.window."
                                       "lammpstrj.gz"))
        return window, folder / f"{whole.source}.ring.all.data.gz"

    return _files


@pytest.fixture
def real_universe() -> Callable[[Path, Path], Any]:
    """Factory: a data file and dump opened with
    :func:`polyphys.probe.make_universe`, without a time step: the tests
    count frames, not time, so the reader's warning about the missing
    ``dt`` is silenced."""

    def _open(topology: Path, trajectory: Path) -> Any:
        from polyphys.probe import make_universe

        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', 'Reader has no dt information', UserWarning)
            return make_universe(topology, trajectory)

    return _open


@pytest.fixture
def trajectory_arrays() -> Callable[[Any], Mapping[str, np.ndarray]]:
    """Factory: every frame of a universe in one pass: ``'positions'``,
    shape ``(n_frames, n_atoms, 3)``, in the reader's frame;
    ``'boxes'``, shape ``(n_frames, 3)``, the box lengths; and
    ``'steps'``, shape ``(n_frames,)``, the LAMMPS time steps."""

    def _arrays(universe: Any) -> Mapping[str, np.ndarray]:
        frames = [(ts.positions.copy(), ts.dimensions[:3].copy(),
                   ts.data['step']) for ts in universe.trajectory]
        return {'positions': np.stack([frame[0] for frame in frames]),
                'boxes': np.stack([frame[1] for frame in frames]),
                'steps': np.array([frame[2] for frame in frames])}

    return _arrays


@pytest.fixture
def reference_distances() -> Callable[..., np.ndarray]:
    """Factory: the distances between the rows of `first` and `second`,
    shapes ``(n_frames, n, 3)`` and ``(n_frames, m, 3)``, frame by frame
    with :func:`MDAnalysis.lib.distances.distance_array` and the minimum
    image along the geometry's periodic directions: all three in a cube;
    only z in a cylinder, whose box lengths across the axis are replaced by
    1e4, far beyond any separation in it, so that only z is imaged.
    Returns shape ``(n_frames, n, m)``."""

    def _distances(first: np.ndarray, second: np.ndarray, boxes: np.ndarray,
                   geometry: str) -> np.ndarray:
        from MDAnalysis.lib.distances import distance_array

        lengths = np.array(boxes, dtype=np.float64)
        if geometry == 'cylindrical':
            lengths[:, :2] = 1e4
        return np.stack([
            distance_array(a, b, box=np.concatenate([box, [90.0] * 3]))
            for a, b, box in zip(first, second, lengths)])

    return _distances
