"""
Fixtures shared by the probe tests.

Every test here needs MDAnalysis, which the ``mda`` extra installs; the test
modules call ``pytest.importorskip("MDAnalysis")``, and the fixtures import
the probe layer inside their factories so that collecting this folder
without MDAnalysis is harmless.

Fixtures
--------
make_chain_universe
    The files of ``make_chain_files``, in ``polyphys/tests/conftest.py``,
    opened with :func:`polyphys.probe.make_universe`.
real_polymer_whole, real_window_whole
    The real wholes of the tier B tests, one test per whole (see "Real
    frames" below).
real_files
    The files of a real whole: its polymer dump and the topology it is
    opened with, and its ``all`` window and data file when they are tier B.
real_universe
    A real dump opened with :func:`polyphys.probe.make_universe`.
trajectory_arrays
    The positions, box lengths, and time steps of every frame of a
    universe, in one pass.

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

The tests here open these files with the library alone. The tests that
probe them with a study's specs, and those that compare them with the
research code's output, are in the polyphys-lab repository
(https://github.com/amirhs1/polyphys-lab), with the studies' parsers.

A polymer dump holds only the polymer group, so it cannot be opened with
its whole's full data file. It is opened with the data file of the whole
without crowders of the same family and replicate (:attr:`RealWhole.topology`),
whose atoms are the polymer group: the ring for TransFoci, the ring and the
H-NS proteins for Hns. The research code did the same. The polymer atoms
carry the IDs 1 to n in every system, with the crowders after them.
"""
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

import numpy as np
import pytest


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
    """One real whole of the tier B tests.

    Attributes
    ----------
    id : str
        The test id, such as ``'HnsCub-nh30-nc10743-ens1'``.
    study : str
        Its study's package name in polyphys-lab, such as ``'hns_cub'``.
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
