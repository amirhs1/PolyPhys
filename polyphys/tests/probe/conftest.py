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
real_polymer_whole, real_window_whole, real_foci_whole, real_hns_whole, \
real_nh0_whole
    The real wholes of the tier B tests, one test per whole (see "Real
    frames" below).
real_parser
    The parser of a real whole's name in one of its groups.
probe_spec
    The test-owned probe spec of a synthetic scheme.
real_files
    The files of a real whole: its polymer dump and the topology it is
    opened with, and its ``all`` window and data file when they are tier B.
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
14 wholes of the three published studies (HnsCub, TransFociCub, and
TransFociCyl), the last 50 frames of the polymer group dump
(``<whole>.ring.<group>.last50.lammpstrj.gz``), and for the wholes without
crowders, the ``all`` frames of the same time window
(``<whole>.j<NN>.ring.all.window.lammpstrj.gz``) and the whole's data file
(``<whole>.ring.all.data.gz``). The maintainer trimmed them from the
studies' runs.

The files keep their real names, and the tests parse them with the
synthetic scheme of the same geometry and species in
``polyphys/tests/synthetic.py``: ``RingCubBinders`` for HnsCub, ``RingCub``
for TransFociCub, and ``RingCyl`` for TransFociCyl. A parser reads
only the keywords of its own pattern, by name and in any order, so it skips
the others (``kbmm`` and ``epshc`` of HnsCub, ``epss`` and ``epsl`` of
TransFociCyl), and a whole carries fewer attributes than its study's parser
gives it. The probe specs are the tests' own, with the values of the studies'
specs in the polyphys-lab repository (https://github.com/amirhs1/polyphys-lab),
which also holds the tests that need a study's parser or the maintainer's
data.

A polymer dump holds only the polymer group, so it cannot be opened with
its whole's full data file. It is opened with the data file of the whole
without crowders of the same family and replicate (:attr:`RealWhole.topology`),
whose atoms are the polymer group: the ring for TransFoci, the ring and the
H-NS proteins for HnsCub. The research code did the same. The polymer atoms
carry the IDs 1 to n in every system, with the crowders after them.
"""
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

import numpy as np
import pytest

from polyphys.manage.parser import ParserBase
from polyphys.spec import BinderSpec, ProbeSpec, ProjectSpec
from polyphys.tests.synthetic import RingCub, RingCubBinders, RingCyl


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
    scheme : str
        The synthetic scheme that parses its name: ``'RingCubBinders'``,
        ``'RingCub'``, or ``'RingCyl'``.
    whole : str
        The whole's name without ``.ring``, as its files carry it.
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
        ``'C'`` when it is the maintainer's.
    n_window_frames : int
        The frames of the ``all`` window: 25 for HnsCub, 20 for TransFoci.
    n_shared_frames : int
        The window's frames that are also polymer frames, by time step:
        all 25 for HnsCub, whose ``all`` dump period (10000) is twice the
        ``nucleoid`` one (5000); 10 for TransFoci, whose periods are 5000
        and 2000.
    """

    id: str
    scheme: str
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

    def polymer_counts(self) -> dict[str, int]:
        """The atoms of each species of the probe spec in the polymer
        dump: two patches per H-NS core (Jung, Sadeghi, and Ha, Sci. Rep.
        14, 139 (2024)), the small monomers and the foci of a TransFoci
        ring, and no crowders."""
        if self.group == 'nucleoid':
            return {'Mon': self.n_mon, 'HnsPatch': 2 * self.n_hns,
                    'Hns': self.n_hns, 'Crd': 0}
        return {'Mon': self.n_mon, 'Foci': self.n_foci,
                'Dna': self.n_mon - self.n_foci, 'Crd': 0}

    def window_counts(self) -> dict[str, int]:
        """The atoms of each species of the probe spec in the ``all``
        window: the polymer group plus the crowders."""
        return {**self.polymer_counts(), 'Crd': self.n_crowd}


def _real_wholes() -> list[RealWhole]:
    """The 14 real wholes: per study, two replicates without crowders and
    two at a bulk crowder volume fraction of about 0.36 to 0.38, and for
    HnsCub the two wholes without H-NS."""
    hns_cub = 'N200kbmm2nh{nh}ac2l25epshc1nc{nc}ens{e}'
    foci_cub = ('al6nl5ml216ns400ac1nc{nc}l{l}dt0.005bdump2000adump5000'
                'ens{e}')
    foci_cyl = ('epss5epsl5r10.5al6nl5ml216ns400ac1nc{nc}lz{l}dt0.005'
                'bdump2000adump5000ens{e}')
    wholes = []
    for e in (1, 2):
        for nh, nc, tier in ((0, 0, 'B'), (30, 0, 'B'), (30, 10743, 'C')):
            wholes.append(RealWhole(
                f"HnsCub-nh{nh}-nc{nc}-ens{e}", 'RingCubBinders',
                hns_cub.format(nh=nh, nc=nc, e=e),
                hns_cub.format(nh=nh, nc=0, e=e), 'nucleoid', 200, 0, nh,
                nc, tier, 25, 25))
        for template, scheme, label, uncrowded, crowded in (
                (foci_cub, 'RingCub', 'TransFociCub',
                 (0, '46.5'), (112776, '27')),
                (foci_cyl, 'RingCyl', 'TransFociCyl',
                 (0, '82.5'), (22230, '57'))):
            topology = template.format(nc=uncrowded[0], l=uncrowded[1], e=e)
            for (nc, length), tier in ((uncrowded, 'B'), (crowded, 'C')):
                wholes.append(RealWhole(
                    f"{label}-nc{nc}-ens{e}", scheme,
                    template.format(nc=nc, l=length, e=e), topology, 'bug',
                    405, 5, 0, nc, tier, 20, 10))
    return sorted(wholes, key=lambda whole: whole.id)


def _params(keep: Callable[[RealWhole], bool]) -> list[Any]:
    """The real wholes that `keep` accepts, as pytest params."""
    return [pytest.param(whole, id=whole.id) for whole in _real_wholes()
            if keep(whole)]


@pytest.fixture(params=_params(lambda whole: True))
def real_polymer_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 14 real wholes, whose polymer dumps are all tier B."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.window_tier == 'B'))
def real_window_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 8 real wholes without crowders, whose ``all`` window and
    data file are tier B."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.n_foci > 0))
def real_foci_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 8 TransFoci wholes."""
    return request.param


@pytest.fixture(params=_params(lambda whole: whole.n_hns > 0))
def real_hns_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 4 HnsCub wholes with H-NS proteins."""
    return request.param


@pytest.fixture(params=_params(
    lambda whole: whole.group == 'nucleoid' and whole.n_hns == 0))
def real_nh0_whole(request: pytest.FixtureRequest) -> RealWhole:
    """Each of the 2 HnsCub wholes without H-NS proteins."""
    return request.param


@pytest.fixture
def real_parser() -> Callable[[RealWhole, str], ParserBase]:
    """Factory: the parser of a real whole's name, with its topology
    suffix, at the ``whole`` lineage in `group`, by the whole's synthetic
    scheme."""
    schemes: dict[str, type[ParserBase]] = {
        'RingCubBinders': RingCubBinders, 'RingCub': RingCub,
        'RingCyl': RingCyl}

    def _parse(whole: RealWhole, group: str) -> ParserBase:
        return schemes[whole.scheme](f"{whole.whole}.ring", 'whole', group)

    return _parse


@pytest.fixture
def probe_spec() -> Callable[[str], ProbeSpec]:
    """Factory: the probe spec of a synthetic scheme. Each takes the
    selections, chain-size observables, unwrap convention, bin sizes,
    histogram species, and binder or pair settings of its study's
    ``PROBE`` in polyphys-lab, ``projects/<study>/python/spec.py``: HnsCub
    for ``RingCubBinders``, TransFociCub for ``RingCub``, and TransFociCyl
    for ``RingCyl``."""

    def _tenth(*attrs: str) -> Callable[[ParserBase], float]:
        """A tenth of the smallest of the diameters `attrs`."""
        return lambda info: 0.1 * min(float(getattr(info, a)) for a in attrs)

    def _half(*attrs: str) -> Callable[[ParserBase], float]:
        """Half the smallest of the diameters `attrs`."""
        return lambda info: 0.5 * min(float(getattr(info, a)) for a in attrs)

    def _contact(info: ParserBase) -> float:
        """The large-monomer diameter plus the crowder diameter."""
        return (float(getattr(info, 'dmon_large'))
                + float(getattr(info, 'dcrowd')))

    def _foci_project(parser: type[ParserBase],
                      tseries: tuple[str, ...]) -> ProjectSpec:
        return ProjectSpec(
            parser=parser, group='bug', divisor=0.025,
            crowder_attr='phi_bulk_c',
            physical_attrs=('space', 'nmon_small', 'nmon_large',
                            'dmon_large', 'dcrowd'),
            species_sizes={'Mon': 'dmon_small', 'Crd': 'dcrowd',
                           'Foci': 'dmon_large'},
            tseries_properties=tseries)

    foci_selections = {'Mon': 'resid 1', 'Foci': 'type 2', 'Crd': 'resid 0',
                       'Dna': 'type 1'}
    foci_species = ('Crd', 'Mon', 'Dna', 'Foci')

    def _ring_cub_binders() -> ProbeSpec:
        project = ProjectSpec(
            parser=RingCubBinders, group='nucleoid', divisor=0.01,
            crowder_attr='phi_bulk_c',
            physical_attrs=('space', 'nmon', 'nhns', 'dcrowd'),
            species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd', 'Hns': 'dhns'},
            tseries_properties=(
                'asphericityTMon', 'gyrTMon', 'shapeTMon', 'nBoundTHnsPatch',
                'nFreeTHnsPatch', 'nEngagedTHnsPatch', 'nFreeTHnsCore',
                'nBridgeTHnsCore', 'nDangleTHnsCore', 'nCisTHnsCore',
                'nTransTHnsCore'))
        return ProbeSpec(
            project,
            {'Mon': 'resid 1', 'HnsPatch': 'type 2', 'Crd': 'resid 0',
             'Hns': 'type 3'},
            'Mon', ('asphericityT', 'gyrT', 'shapeT', 'principalT'),
            unwrap=True,
            bin_sizes={'r': _tenth('dmon', 'dhns', 'dcrowd'),
                       'x': _half('dmon', 'dcrowd'),
                       'y': _half('dmon', 'dcrowd'),
                       'z': _half('dmon', 'dcrowd')},
            histogram_species=('Crd', 'Mon', 'Hns'), block_frames=100,
            recenter='half',
            binders=BinderSpec(patches='HnsPatch', binders='HnsCore',
                               cutoff=0.685, inclusive=True,
                               cis_threshold=3, store_distances=True))

    def _ring_cub() -> ProbeSpec:
        return ProbeSpec(
            _foci_project(RingCub,
                          ('asphericityTMon', 'gyrTMon', 'shapeTMon')),
            foci_selections, 'Mon',
            ('asphericityT', 'gyrT', 'shapeT', 'principalT'), unwrap=True,
            bin_sizes={'r': _tenth('dmon_small', 'dcrowd'),
                       'x': _half('dmon_small', 'dcrowd'),
                       'y': _half('dmon_small', 'dcrowd'),
                       'z': _half('dmon_small', 'dcrowd')},
            histogram_species=foci_species, block_frames=100,
            recenter='half', pair_species='Foci', contact_cutoff=_contact,
            minimum_image=True)

    def _ring_cyl() -> ProbeSpec:
        return ProbeSpec(
            _foci_project(RingCyl,
                          ('asphericityTMon', 'fsdTMon', 'gyrTMon',
                           'shapeTMon', 'transSizeTMon')),
            foci_selections, 'Mon',
            ('asphericityT', 'fsdT', 'gyrT', 'shapeT', 'transSizeT',
             'principalT'), unwrap=False,
            bin_sizes={'r': _tenth('dmon_small', 'dcrowd'),
                       'z': _half('dmon_small', 'dcrowd'),
                       'theta': np.pi / 36,
                       'x': _tenth('dmon_small', 'dcrowd'),
                       'y': _tenth('dmon_small', 'dcrowd')},
            histogram_species=foci_species, block_frames=100,
            recenter='half', pair_species='Foci', contact_cutoff=_contact,
            minimum_image=True)

    builders = {'RingCubBinders': _ring_cub_binders, 'RingCub': _ring_cub,
                'RingCyl': _ring_cyl}

    def _spec(scheme: str) -> ProbeSpec:
        return builders[scheme]()

    return _spec


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
