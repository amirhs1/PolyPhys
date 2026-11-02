"""Tests of :mod:`polyphys.simulate.build`.

Every expected value is analytic, a property of the construction, or a
threshold read from the LAMMPS source and measured with LAMMPS (see the
module documentation); the ``lammps``-marked tests take LAMMPS itself as
the oracle and are skipped without it.
"""

import math
import subprocess
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.simulate.build import (
    ConformationError,
    _plectoneme_curve,
    chain_angles,
    chain_bonds,
    chain_data,
    check_conformation,
    helix,
    plectoneme_ring,
    resample_curve,
)
from polyphys.simulate.core import Box
from polyphys.simulate.data import DataFile, read_data, write_data


@pytest.fixture
def build_dir() -> Path:
    """The folder of the LAMMPS input script ``in.kg``."""
    return Path(__file__).parent.parent / "data" / "lammps_build"


@pytest.fixture
def bond_lengths() -> Callable[..., np.ndarray]:
    """The lengths of the bonds along a chain, the closing one for a ring."""
    def _lengths(positions: np.ndarray, ring: bool = False) -> np.ndarray:
        if ring:
            return np.linalg.norm(np.roll(positions, -1, axis=0) - positions,
                                  axis=1)
        return np.linalg.norm(np.diff(positions, axis=0), axis=1)
    return _lengths


@pytest.fixture
def closest_non_bonded() -> Callable[..., float]:
    """The shortest distance between beads that are not bonded, by brute
    force without periodic images."""
    def _closest(positions: np.ndarray, ring: bool = False) -> float:
        n = positions.shape[0]
        distances = np.linalg.norm(
            positions[:, None, :] - positions[None, :, :], axis=-1)
        distances[np.diag_indices(n)] = np.inf
        first = np.arange(n - 1)
        distances[first, first + 1] = distances[first + 1, first] = np.inf
        if ring:
            distances[0, n - 1] = distances[n - 1, 0] = np.inf
        return float(distances.min())
    return _closest


@pytest.fixture
def make_box() -> Callable[..., Box]:
    """A box centred on the origin, by default a cube of side 10."""
    def _box(half: float = 5.0, half_height: float | None = None,
             **kwargs: Any) -> Box:
        height = half if half_height is None else half_height
        return Box(lo=np.array([-half, -half, -height]),
                   hi=np.array([half, half, height]), **kwargs)
    return _box


@pytest.fixture
def make_record(make_box: Callable[..., Box]) -> Callable[..., DataFile]:
    """A ``bond``-style data file from positions, with the bonds given as
    ``(atom1, atom2)`` pairs of IDs from 1, none by default."""
    def _record(positions: Any, bonds: Any = (), box: Box | None = None,
                ids: Any = None) -> DataFile:
        coordinates = np.array(positions, dtype=np.float64)
        n = coordinates.shape[0]
        pairs = np.array(bonds, dtype=np.int64).reshape(-1, 2)
        return DataFile(
            atom_style="bond",
            box=make_box() if box is None else box,
            ids=np.arange(1, n + 1) if ids is None else np.array(ids),
            types=np.ones(n, dtype=np.int64),
            positions=coordinates,
            n_atom_types=1,
            mol_ids=np.ones(n, dtype=np.int64),
            bonds=np.column_stack([np.ones(len(pairs), dtype=np.int64),
                                   pairs]) if len(pairs) else np.empty(
                (0, 3), dtype=np.int64),
            n_bond_types=1 if len(pairs) else 0,
        )
    return _record


@pytest.fixture
def start_file(tmp_path: Path, make_box: Callable[..., Box]
               ) -> Callable[[str], tuple[Path, str]]:
    """Write a start file of 200 beads (radius 2, pitch 2.5, bond 0.97) of
    one kind, ``helix``, ``ring``, or ``ring_angle``, in a box periodic
    along z only, the configurations of ``make_starts.py``."""
    def _write(kind: str) -> tuple[Path, str]:
        if kind == "helix":
            positions = helix(200, radius=2.0, pitch=2.5)
        else:
            positions = plectoneme_ring(200, radius=2.0, pitch=2.5)
        box = make_box(half=5.0,
                       half_height=0.5 * np.ptp(positions[:, 2]) + 2.0)
        data = chain_data(positions, box, ring=kind != "helix",
                          angles=kind == "ring_angle")
        return write_data(data, tmp_path / f"{kind}.data"), data.atom_style
    return _write


# helix


def test_helix_steps_by_the_chord_root_angle() -> None:
    """A bond of sqrt(2) on a helix of radius 1 and pitch 6 steps by pi/3 and
    rises by 1 per bead, centred on z = 0."""
    # 4 sin^2(pi/6) + (6 / 2 pi)^2 (pi/3)^2 = 1 + 1 = 2 = b^2 (analytic)
    k = np.arange(7)
    expected = np.column_stack(
        [np.cos(k * np.pi / 3), np.sin(k * np.pi / 3), k - 3.0])

    beads = helix(7, radius=1.0, pitch=6.0, bond_length=math.sqrt(2.0))

    # float64 round-off of the root and the trigonometry
    assert_allclose(beads, expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    ("radius", "pitch", "bond"),
    [
        pytest.param(2.0, 2.5, 0.97, id="narrow"),
        pytest.param(10.0, 1.5, 0.97, id="wide"),
        pytest.param(1.0, 0.5, 1.99, id="near-half-a-turn"),
    ],
)
def test_helix_bonds_equal_the_bond_length(
    bond_lengths: Callable[..., np.ndarray], radius: float, pitch: float,
    bond: float,
) -> None:
    """Every bond of a helix has the requested length."""
    beads = helix(50, radius=radius, pitch=pitch, bond_length=bond)

    # exactness up to float64 round-off is the contract
    assert_allclose(bond_lengths(beads), bond, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize(
    ("radius", "pitch"),
    [pytest.param(2.0, 2.5, id="narrow"), pytest.param(10.0, 1.5, id="wide")],
)
def test_helix_lies_on_its_cylinder_centred_on_z(
    radius: float, pitch: float
) -> None:
    """Every bead is the radius from the z axis and the beads are centred on
    z = 0."""
    beads = helix(50, radius=radius, pitch=pitch)

    assert_allclose(np.hypot(beads[:, 0], beads[:, 1]), radius, rtol=1e-12)
    assert beads[:, 2].mean() == pytest.approx(0.0, abs=1e-12)  # round-off


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param({"n_beads": 0}, "'n_beads' must be 1 or more",
                     id="no-beads"),
        pytest.param({"n_beads": 5.0}, "'n_beads' must be an integer",
                     id="float-beads"),
        pytest.param({"radius": 0.0}, "'radius' must be positive",
                     id="zero-radius"),
        pytest.param({"pitch": -1.0}, "'pitch' must be positive",
                     id="negative-pitch"),
        pytest.param({"bond_length": 0.0}, "'bond_length' must be positive",
                     id="zero-bond"),
        # 4 R^2 + (p / 2 pi)^2 pi^2 = 4.0625 < 2.2^2 (analytic)
        pytest.param({"radius": 1.0, "pitch": 0.5, "bond_length": 2.2},
                     "does not fit in half a turn", id="bond-too-long"),
    ],
)
def test_helix_rejects_impossible_parameters(
    kwargs: dict[str, Any], match: str
) -> None:
    """A bead count below 1, a non-positive length, or a bond longer than
    half a turn raises ``ValueError``."""
    arguments: dict[str, Any] = {"n_beads": 5, "radius": 1.0, "pitch": 1.0}
    arguments.update(kwargs)

    with pytest.raises(ValueError, match=match):
        helix(**arguments)


@pytest.mark.parametrize(
    ("build", "ring"),
    [
        pytest.param(lambda: helix(30, radius=2.0, pitch=2.5), False,
                     id="helix"),
        pytest.param(lambda: plectoneme_ring(60, radius=2.0, pitch=2.5),
                     True, id="plectoneme-ring"),
        pytest.param(lambda: resample_curve(
            [[0.0, 0.0, 0.0], [30.0, 0.0, 0.0]], 30), False, id="resample"),
    ],
)
def test_builders_default_to_the_kremer_grest_bond(
    bond_lengths: Callable[..., np.ndarray], build: Callable[[], np.ndarray],
    ring: bool,
) -> None:
    """Without a bond length, every builder uses 0.97, the start bond length
    of Kremer and Grest (1990)."""
    beads = build()

    assert_allclose(bond_lengths(beads, ring=ring), 0.97, rtol=1e-12,
                    atol=0.0)


# resample_curve


def test_resample_curve_spaces_beads_along_a_line() -> None:
    """On a straight segment the beads sit one bond length apart from the
    first point."""
    expected = np.column_stack(
        [0.97 * np.arange(10), np.zeros(10), np.zeros(10)])

    beads = resample_curve([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], 10, 0.97)

    assert_allclose(beads, expected, rtol=1e-12, atol=1e-12)  # round-off


def test_resample_curve_crosses_a_corner_at_the_sphere() -> None:
    """The next bead is where the sphere of one bond length about the last
    bead first meets the polyline, past a corner."""
    # |(1, y) - (0, 0)| = sqrt(2) at y = 1; then one more bond along y
    corner = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 5.0, 0.0]]
    expected = [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0],
                [1.0, 1.0 + math.sqrt(2.0), 0.0]]

    beads = resample_curve(corner, 3, math.sqrt(2.0))

    assert_allclose(beads, expected, rtol=1e-12, atol=1e-12)  # round-off


def test_resample_curve_rejects_a_short_curve() -> None:
    """An open curve that holds fewer beads than asked raises ``ValueError``
    naming how many it holds."""
    # a segment of length 1.5 holds beads at 0 and 1 only
    with pytest.raises(ValueError, match="holds 2 beads 1.0 apart, not 4"):
        resample_curve([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], 4, 1.0)


def test_resample_curve_closes_a_circle_as_a_regular_polygon(
    bond_lengths: Callable[..., np.ndarray]
) -> None:
    """Twelve unit bonds close a circle as a regular 12-gon, whose
    circumradius is 1 / (2 sin(pi / 12))."""
    angles = np.linspace(0.0, 2.0 * np.pi, 2000, endpoint=False)
    circle = np.column_stack(
        [np.cos(angles), np.sin(angles), np.zeros_like(angles)])
    circumradius = 1.0 / (2.0 * math.sin(math.pi / 12.0))  # analytic

    beads = resample_curve(circle, 12, 1.0, closed=True)

    # the polyline's vertices are on the circle and its edges are inside it
    # by 1 - cos(pi / 2000) = 1.2e-6 of the radius
    assert_allclose(np.hypot(beads[:, 0], beads[:, 1]), circumradius,
                    rtol=1e-5)
    assert_allclose(bond_lengths(beads, ring=True), 1.0, rtol=1e-12,
                    atol=0.0)


def test_resample_curve_keeps_a_closed_curve_that_already_fits() -> None:
    """A square of side 4 holds 16 unit bonds as it is: the beads land on
    the corners and the unit marks, unscaled."""
    square = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0],
              [0.0, 4.0, 0.0]]
    marks = np.arange(4, dtype=np.float64)
    expected = np.vstack([
        np.column_stack([marks, np.zeros(4), np.zeros(4)]),
        np.column_stack([np.full(4, 4.0), marks, np.zeros(4)]),
        np.column_stack([4.0 - marks, np.full(4, 4.0), np.zeros(4)]),
        np.column_stack([np.zeros(4), 4.0 - marks, np.zeros(4)]),
    ])

    beads = resample_curve(square, 16, 1.0, closed=True)

    assert_allclose(beads, expected, rtol=0.0, atol=1e-12)  # round-off


@pytest.mark.parametrize(
    ("curve", "n_beads"),
    [
        pytest.param("ellipse", 40, id="ellipse"),
        pytest.param("figure-eight", 60, id="figure-eight-in-3d"),
    ],
)
def test_resample_curve_closes_with_equal_bonds(
    bond_lengths: Callable[..., np.ndarray], curve: str, n_beads: int
) -> None:
    """Every bond of a closed curve, the closing one included, has the
    requested length."""
    t = np.linspace(0.0, 2.0 * np.pi, 5000, endpoint=False)
    if curve == "ellipse":
        points = np.column_stack([3.0 * np.cos(t), 1.5 * np.sin(t), 0 * t])
    else:
        points = np.column_stack(
            [np.cos(t), 0.5 * np.sin(2.0 * t), 0.3 * np.sin(t)])

    beads = resample_curve(points, n_beads, 0.97, closed=True)

    assert beads.shape == (n_beads, 3)
    assert_allclose(bond_lengths(beads, ring=True), 0.97, rtol=1e-12,
                    atol=0.0)


def test_resample_curve_drops_a_repeated_first_point() -> None:
    """A closed curve that repeats its first point last gives the same beads
    as one that does not."""
    square = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0],
              [0.0, 4.0, 0.0]]

    beads = resample_curve(square + square[:1], 16, 1.0, closed=True)

    assert_allclose(beads, resample_curve(square, 16, 1.0, closed=True),
                    rtol=0.0, atol=1e-12)  # round-off


def test_resample_curve_refuses_a_curve_it_cannot_close(
) -> None:
    """A closed curve with teeth smaller than the bond is refused: no scale
    closes it with equal bonds, so it raises rather than returning a ring
    with one odd bond."""
    t = np.linspace(0.0, 2.0 * np.pi, 4000, endpoint=False)
    radius = 2.0 + np.sign(np.sin(12.0 * t)) * np.abs(np.sin(12.0 * t)) ** 0.2
    teeth = np.column_stack([radius * np.cos(t), radius * np.sin(t), 0 * t])

    with pytest.raises(ValueError, match="cannot be closed with equal bonds"):
        resample_curve(teeth, 30, 0.97, closed=True)


@pytest.mark.parametrize(
    ("points", "kwargs", "match"),
    [
        pytest.param([[0.0, 0.0], [1.0, 0.0]], {}, r"shape \(n_points, 3\)",
                     id="not-3d"),
        pytest.param([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]], {},
                     "must be finite", id="nan"),
        pytest.param([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], {},
                     "2 or more distinct points", id="one-point"),
        pytest.param([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                     {"n_beads": 2, "closed": True},
                     "'n_beads' must be 3 or more", id="closed-two-beads"),
    ],
)
def test_resample_curve_rejects_bad_input(
    points: Any, kwargs: dict[str, Any], match: str
) -> None:
    """Points that are not ``(n, 3)``, not finite, or not two distinct, and
    a closed curve with fewer than three beads, raise ``ValueError``."""
    arguments: dict[str, Any] = {"n_beads": 2, "bond_length": 0.5}
    arguments.update(kwargs)

    with pytest.raises(ValueError, match=match):
        resample_curve(points, **arguments)


# plectoneme_ring


@pytest.mark.parametrize(
    ("n_beads", "radius", "pitch"),
    [
        pytest.param(200, 2.0, 2.5, id="docs-ring"),
        pytest.param(60, 2.0, 2.5, id="short"),
        pytest.param(90, 1.0, 2.0, id="narrow"),
    ],
)
def test_plectoneme_ring_bonds_equal_the_bond_length(
    bond_lengths: Callable[..., np.ndarray], n_beads: int, radius: float,
    pitch: float,
) -> None:
    """Every bond of the ring, the closing one included, has the requested
    length, and the ring has the requested number of beads."""
    beads = plectoneme_ring(n_beads, radius=radius, pitch=pitch)

    assert beads.shape == (n_beads, 3)
    assert_allclose(bond_lengths(beads, ring=True), 0.97, rtol=1e-12,
                    atol=0.0)


@pytest.mark.parametrize(
    ("n_beads", "radius", "pitch"),
    [pytest.param(200, 2.0, 2.5, id="docs-ring"),
     pytest.param(90, 1.0, 2.0, id="narrow")],
)
def test_plectoneme_ring_stays_inside_its_radius_centred_on_z(
    n_beads: int, radius: float, pitch: float
) -> None:
    """No bead is farther than the radius from the axis, the strands reach
    it, and the ring is centred on z = 0."""
    beads = plectoneme_ring(n_beads, radius=radius, pitch=pitch)
    radial = np.hypot(beads[:, 0], beads[:, 1])

    assert radial.max() == pytest.approx(radius, rel=1e-12)  # round-off
    assert np.all(radial <= radius * (1.0 + 1e-12))
    assert beads[:, 2].max() == pytest.approx(-beads[:, 2].min(), abs=1e-12)


def test_plectoneme_curve_follows_neukirch_starostin_appendix_a() -> None:
    """The plasmid curve is continuous at its four joints and closes, loop B
    crosses the axis at its midpoint, strand C is strand A rotated by pi
    about the axis and run downward, and the loops' interior points are
    those of Eqs. (A2) and (A4) of Neukirch and Starostin (2008)."""
    strand, theta = 5.0, 1.0
    turn = 1.0 + 2.0 * theta / np.pi
    twist = strand * np.sin(theta)
    joints = np.array([strand, strand + np.pi, 2.0 * strand + np.pi,
                       2.0 * strand + 2.0 * np.pi])
    quarter = np.pi / 4.0
    # Eq. (A2) at s_B = pi/4 and Eq. (A4) at s_D = pi/4, transcribed
    theta_b = theta - quarter * turn
    half_sine = np.sin(theta) * np.sin(quarter) * np.cos(quarter)
    loop_b = [np.sin(twist) * np.cos(quarter) + np.cos(twist) * half_sine,
              -np.cos(twist) * np.cos(quarter) + np.sin(twist) * half_sine,
              strand * np.cos(theta)
              + (np.sin(theta) - np.sin(theta_b)) / turn]
    theta_d = -np.pi - theta + quarter * turn
    loop_d = [half_sine, np.cos(quarter),
              (np.sin(theta_d) - np.sin(theta)) / turn]

    after = _plectoneme_curve(joints, strand, theta)
    before = _plectoneme_curve(joints - 1e-9, strand, theta)
    apex = _plectoneme_curve(np.array([strand + np.pi / 2.0]), strand, theta)
    s = np.linspace(0.0, strand, 7)
    strand_a = _plectoneme_curve(s, strand, theta)
    strand_c = _plectoneme_curve(2.0 * strand + np.pi - s, strand, theta)
    points = _plectoneme_curve(
        np.array([strand + quarter, 2.0 * strand + np.pi + quarter]),
        strand, theta)

    assert_allclose(after, before, rtol=0.0, atol=1e-7)  # 1e-9 x speed
    assert_allclose(after[3], _plectoneme_curve(np.zeros(1), strand, theta)[0],
                    rtol=0.0, atol=1e-12)
    assert_allclose(apex[0], [0.0, 0.0, strand * np.cos(theta)
                              + (1.0 + np.sin(theta)) / turn],
                    rtol=0.0, atol=1e-12)
    assert_allclose(strand_c[:, :2], -strand_a[:, :2], rtol=0.0, atol=1e-12)
    assert_allclose(strand_c[:, 2], strand_a[:, 2], rtol=0.0, atol=1e-12)
    assert_allclose(points, [loop_b, loop_d], rtol=1e-12, atol=1e-12)


def test_plectoneme_ring_starts_at_the_foot_of_the_first_strand() -> None:
    """The first bead is the foot of strand A, at angle zero on the
    cylinder, (0, -R) in Neukirch and Starostin's frame, and the chain
    climbs from there."""
    beads = plectoneme_ring(60, radius=2.0, pitch=2.5)

    assert_allclose(beads[0, :2], [0.0, -2.0], rtol=0.0, atol=1e-12)
    assert beads[1, 2] > beads[0, 2]


def test_plectoneme_ring_passes_the_conformation_check(
    make_box: Callable[..., Box]
) -> None:
    """The ring of the documentation, 200 beads of radius 2 and pitch 2.5,
    starts a Kremer-Grest run: its strands stay a bead diameter apart."""
    beads = plectoneme_ring(200, radius=2.0, pitch=2.5)
    box = make_box(half=5.0, half_height=0.5 * np.ptp(beads[:, 2]) + 2.0)
    data = chain_data(beads, box, ring=True)

    report = check_conformation(data, r0=1.5, ghost_cutoff=1.122462 + 0.3,
                                periodic=(False, False, True),
                                min_distance=1.0)

    assert report.min_pair >= 1.0
    assert report.max_bond == pytest.approx(0.97, rel=1e-12)


def test_plectoneme_ring_with_a_close_pitch_fails_the_pair_check(
    make_box: Callable[..., Box]
) -> None:
    """A pitch of 1.5 at radius 3 puts the two strands about 0.75 apart, so
    the ring fails the non-bonded distance check."""
    # (p / 2) sin(theta_A) with tan(theta_A) = 2 pi R / p gives 0.748; the
    # measured closest pair is 0.7555 (reproducible: this ring)
    beads = plectoneme_ring(1000, radius=3.0, pitch=1.5)
    box = make_box(half=6.0, half_height=0.5 * np.ptp(beads[:, 2]) + 2.0)
    data = chain_data(beads, box, ring=True)

    with pytest.raises(ConformationError,
                       match=r"non-bonded pairs closer than 1: closest 0.75"):
        check_conformation(data, r0=1.5, ghost_cutoff=1.122462 + 0.3,
                           periodic=(False, False, True))


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        # the two loops alone are 15.4 long at radius 2: about 16 beads
        pytest.param({"n_beads": 12}, "'n_beads' must be at least about",
                     id="too-few-for-the-loops"),
        pytest.param({"n_beads": 2}, "'n_beads' must be 3 or more",
                     id="two-beads"),
        pytest.param({"radius": -2.0}, "'radius' must be positive",
                     id="negative-radius"),
        pytest.param({"bond_length": 5.0}, "does not fit in half a turn",
                     id="bond-too-long"),
    ],
)
def test_plectoneme_ring_rejects_impossible_parameters(
    kwargs: dict[str, Any], match: str
) -> None:
    """Too few beads for the two end loops, a non-positive length, or a bond
    longer than half a turn raises ``ValueError``."""
    arguments: dict[str, Any] = {"n_beads": 60, "radius": 2.0, "pitch": 2.5}
    arguments.update(kwargs)

    with pytest.raises(ValueError, match=match):
        plectoneme_ring(**arguments)


# topology


@pytest.mark.parametrize(
    ("ring", "expected"),
    [
        pytest.param(False, [[1, 1, 2], [1, 2, 3], [1, 3, 4]], id="linear"),
        pytest.param(True, [[1, 1, 2], [1, 2, 3], [1, 3, 4], [1, 4, 1]],
                     id="ring"),
    ],
)
def test_chain_bonds_rows(ring: bool, expected: list[list[int]]) -> None:
    """A chain of four beads has the three bonds along it, and a ring the
    bond from the last bead to the first as well."""
    assert_array_equal(chain_bonds(4, ring=ring), expected)


@pytest.mark.parametrize(
    ("ring", "expected"),
    [
        pytest.param(False, [[2, 1, 2, 3], [2, 2, 3, 4]], id="linear"),
        pytest.param(True, [[2, 1, 2, 3], [2, 2, 3, 4], [2, 3, 4, 1],
                            [2, 4, 1, 2]], id="ring"),
    ],
)
def test_chain_angles_rows(ring: bool, expected: list[list[int]]) -> None:
    """A chain of four beads has the two angles with an inner vertex, and a
    ring the two that span the closing bond as well, all of the given type."""
    assert_array_equal(chain_angles(4, ring=ring, angle_type=2), expected)


def test_chain_topology_defaults_to_a_linear_chain_without_angles(
    make_box: Callable[..., Box]
) -> None:
    """Without ``ring`` or ``angles``, the topology is that of a linear
    chain in the ``bond`` style: n - 1 bonds, no closing bond, no angles."""
    positions = np.column_stack([0.97 * np.arange(4), np.zeros(4),
                                 np.zeros(4)])

    data = chain_data(positions, make_box())

    assert_array_equal(chain_bonds(4), [[1, 1, 2], [1, 2, 3], [1, 3, 4]])
    assert_array_equal(chain_angles(4), [[1, 1, 2, 3], [1, 2, 3, 4]])
    assert (data.atom_style, data.n_bonds, data.n_angles) == ("bond", 3, 0)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        pytest.param(lambda: chain_bonds(1), "'n_beads' must be 2 or more",
                     id="one-bead-chain"),
        pytest.param(lambda: chain_bonds(2, ring=True),
                     "'n_beads' must be 3 or more", id="two-bead-ring"),
        pytest.param(lambda: chain_angles(2), "'n_beads' must be 3 or more",
                     id="two-bead-angles"),
        pytest.param(lambda: chain_bonds(3, bond_type=0),
                     "'bond_type' must be 1 or more", id="type-zero"),
    ],
)
def test_chain_topology_rejects_too_few_beads(
    call: Callable[[], np.ndarray], match: str
) -> None:
    """Too few beads for the topology, or a type below 1, raises
    ``ValueError``."""
    with pytest.raises(ValueError, match=match):
        call()


@pytest.mark.parametrize(
    ("angles", "style", "n_angles"),
    [pytest.param(False, "bond", 0, id="bonds"),
     pytest.param(True, "angle", 3, id="angles")],
)
def test_chain_data_builds_the_record(
    make_box: Callable[..., Box], angles: bool, style: str, n_angles: int
) -> None:
    """The record numbers the atoms from 1 in order, puts them in molecule
    1, gives the types, bonds, and angles asked for, and leaves masses and
    velocities to the input script."""
    positions = np.column_stack([0.97 * np.arange(5), np.zeros(5),
                                 np.zeros(5)])

    data = chain_data(positions, make_box(), angles=angles,
                      types=[1, 2, 1, 2, 1])

    assert data.atom_style == style
    assert_array_equal(data.ids, [1, 2, 3, 4, 5])
    assert_array_equal(data.mol_ids, [1, 1, 1, 1, 1])
    assert_array_equal(data.types, [1, 2, 1, 2, 1])
    assert (data.n_atom_types, data.n_bonds, data.n_angles) == (
        2, 4, n_angles)
    assert_array_equal(data.positions, positions)
    assert data.masses is None and data.velocities is None


def test_chain_data_round_trips_through_write_and_read(
    make_box: Callable[..., Box], tmp_path: Path
) -> None:
    """A ring with angles written by ``write_data`` reads back with the same
    positions, bonds, and angles."""
    data = chain_data(helix(6, radius=1.0, pitch=2.0), make_box(), ring=True,
                      angles=True)

    read = read_data(write_data(data, tmp_path / "ring.data"))

    assert read.atom_style == "angle"
    assert_array_equal(read.positions, data.positions)
    assert_array_equal(read.bonds, chain_bonds(6, ring=True))
    assert_array_equal(read.angles, chain_angles(6, ring=True))


@pytest.mark.parametrize(
    ("positions", "types", "match"),
    [
        pytest.param([[0.0, 0.0], [1.0, 0.0]], None,
                     r"shape \(n_beads, 3\)", id="not-3d"),
        pytest.param([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [1, 1, 1],
                     r"'types' must have shape \(2,\)", id="types-shape"),
    ],
)
def test_chain_data_rejects_bad_shapes(
    make_box: Callable[..., Box], positions: Any, types: Any, match: str
) -> None:
    """Positions that are not ``(n, 3)`` or types of another length raise
    ``ValueError``."""
    with pytest.raises(ValueError, match=match):
        chain_data(positions, make_box(), types=types)


# check_conformation


@pytest.mark.parametrize(
    ("length", "passes"),
    [pytest.param(1.4230, True, id="below"),
     pytest.param(1.4231, False, id="above")],
)
def test_check_conformation_bond_warning_threshold(
    make_record: Callable[..., DataFile], length: float, passes: bool
) -> None:
    """A bond just below sqrt(0.9) R0 = 1.4230 passes and one just above
    fails, the length at which LAMMPS warns ``FENE bond too long``."""
    # bond_fene.cpp lines 89 to 94: 1 - (r / R0)^2 < 0.1; R0 = 1.5
    data = make_record([[0.0, 0.0, 0.0], [length, 0.0, 0.0]], bonds=[(1, 2)])
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 2.0,
                                 "periodic": (False, False, False)}

    if passes:
        assert check_conformation(data, **arguments).max_bond == length
    else:
        with pytest.raises(ConformationError,
                           match=r"1 bond at or above sqrt\(0.9\) R0"):
            check_conformation(data, **arguments)


@pytest.mark.parametrize(
    ("length", "passes"),
    [pytest.param(1.4224, True, id="below"),
     pytest.param(1.4226, False, id="above")],
)
def test_check_conformation_ghost_cutoff_threshold(
    make_record: Callable[..., DataFile], length: float, passes: bool
) -> None:
    """A bond just below the ghost cutoff passes and one at or above it
    fails: measured with LAMMPS, a bond of 1.4224 across a periodic face
    runs and one of 1.4226 stops with ``Bad FENE bond``."""
    data = make_record([[0.0, 0.0, 0.0], [length, 0.0, 0.0]], bonds=[(1, 2)])
    arguments: dict[str, Any] = {"r0": 3.0, "ghost_cutoff": 1.122462 + 0.3,
                                 "periodic": (True, True, True)}

    if passes:
        assert check_conformation(data, **arguments).max_bond == length
    else:
        with pytest.raises(ConformationError,
                           match="1 bond at or above the ghost cutoff"):
            check_conformation(data, **arguments)


def test_check_conformation_lets_a_bonded_pair_below_min_distance_pass(
    make_record: Callable[..., DataFile]
) -> None:
    """The non-bonded distance rule does not apply to a bonded pair."""
    data = make_record([[0.0, 0.0, 0.0], [0.95, 0.0, 0.0]], bonds=[(1, 2)])

    report = check_conformation(data, r0=1.5, ghost_cutoff=1.42,
                                periodic=(False, False, False),
                                min_distance=1.0)

    assert report.min_bond == 0.95
    assert math.isnan(report.min_pair)


def test_check_conformation_counts_a_1_3_pair_as_non_bonded(
    make_record: Callable[..., DataFile]
) -> None:
    """Beads two apart along the chain repel like any other pair under
    ``special_bonds fene``, so a 1-3 distance below the minimum fails."""
    # bead 3 one bond from bead 2 and 0.9 from bead 1: x = 0.9^2 / (2 b)
    x = 0.81 / (2.0 * 0.97)
    trimer = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0],
              [x, math.sqrt(0.81 - x * x), 0.0]]
    data = make_record(trimer, bonds=[(1, 2), (2, 3)])

    with pytest.raises(ConformationError,
                       match="1 non-bonded pair closer than 1: closest 0.9"):
        check_conformation(data, r0=1.5, ghost_cutoff=1.42,
                           periodic=(False, False, False), min_distance=1.0)


def test_check_conformation_rejects_coincident_atoms(
    make_record: Callable[..., DataFile]
) -> None:
    """Two atoms at one position fail even when they are bonded, since their
    force is not finite."""
    data = make_record([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], bonds=[(1, 2)])

    with pytest.raises(ConformationError, match="1 pair of atoms closer"):
        check_conformation(data, r0=1.5, ghost_cutoff=1.42,
                           periodic=(False, False, False))


@pytest.mark.parametrize(
    ("periodic", "passes"),
    [pytest.param((False, False, True), False, id="periodic-z"),
     pytest.param((False, False, False), True, id="fixed-z")],
)
def test_check_conformation_uses_the_minimum_image_on_periodic_axes(
    make_record: Callable[..., DataFile], periodic: tuple[bool, bool, bool],
    passes: bool,
) -> None:
    """Two beads 0.5 apart through the z faces of a box of height 10 are an
    overlap on a periodic z axis and 9.5 apart on a fixed one."""
    data = make_record([[0.0, 0.0, -4.8], [0.0, 0.0, 4.7]])
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 1.42,
                                 "periodic": periodic, "min_distance": 1.0}

    if passes:
        assert check_conformation(data, **arguments).min_pair == (
            pytest.approx(9.5))
    else:
        with pytest.raises(ConformationError, match="closest 0.5"):
            check_conformation(data, **arguments)


@pytest.mark.parametrize(
    ("x", "periodic", "passes"),
    [
        pytest.param(5.1, (True, False, False), True, id="periodic-outside"),
        pytest.param(5.1, (False, False, False), False, id="fixed-outside"),
        pytest.param(5.0, (False, False, False), False, id="fixed-on-face"),
        pytest.param(-5.0, (False, False, False), True, id="fixed-on-lo"),
    ],
)
def test_check_conformation_box_rule_depends_on_periodicity(
    make_record: Callable[..., DataFile], x: float,
    periodic: tuple[bool, bool, bool], passes: bool,
) -> None:
    """An atom past the box, or on its upper face, fails on a fixed axis,
    where ``read_data`` loses it, and passes on a periodic one, where
    ``read_data`` remaps it."""
    data = make_record([[x, 0.0, 0.0]])
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 1.42,
                                 "periodic": periodic}

    if passes:
        assert check_conformation(data, **arguments).n_atoms == 1
    else:
        with pytest.raises(ConformationError,
                           match=r"1 atom outside the box.*\(atoms 1\)"):
            check_conformation(data, **arguments)


@pytest.mark.parametrize(
    ("x", "passes"),
    [pytest.param(5.5, True, id="inside-the-tilt"),
     pytest.param(6.5, False, id="beyond-the-tilt")],
)
def test_check_conformation_box_rule_uses_fractions_in_a_triclinic_box(
    make_record: Callable[..., DataFile], make_box: Callable[..., Box],
    x: float, passes: bool,
) -> None:
    """In a box tilted by xy = 2 the atom at (x, 1) is inside when its
    fractional coordinate along A is below 1: (x + 5 - 0.6 * 2) / 10."""
    box = make_box(tilt=np.array([2.0, 0.0, 0.0]), triclinic=True)
    data = make_record([[x, 1.0, 0.0]], box=box)
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 1.42,
                                 "periodic": (False, False, False)}

    if passes:
        assert check_conformation(data, **arguments).n_atoms == 1
    else:
        with pytest.raises(ConformationError, match="outside the box"):
            check_conformation(data, **arguments)


@pytest.mark.parametrize(
    ("position", "passes"),
    [pytest.param([7.49, 3.0, 0.0], True, id="inside-along-x"),
     pytest.param([5.0, 5.51, 0.0], False, id="beyond-along-y")],
)
def test_check_conformation_cylinder_axis_is_the_box_centre(
    make_record: Callable[..., DataFile], position: list[float], passes: bool
) -> None:
    """The cylinder axis runs through the box centre, here (5, 3), not the
    origin: a bead 2.49 from it passes and one 2.51 from it fails."""
    box = Box(lo=np.array([0.0, 0.0, -5.0]), hi=np.array([10.0, 6.0, 5.0]))
    data = make_record([position], box=box)
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 1.42,
                                 "periodic": (False, False, True),
                                 "dcyl": 6.0, "dmon": 1.0}

    if passes:
        assert check_conformation(data, **arguments).max_radial == (
            pytest.approx(2.49))
    else:
        with pytest.raises(ConformationError, match="farthest 2.51"):
            check_conformation(data, **arguments)


@pytest.mark.parametrize(
    ("radial", "passes"),
    [pytest.param(2.49, True, id="inside"),
     pytest.param(2.5, True, id="on-the-limit"),
     pytest.param(2.51, False, id="beyond")],
)
def test_check_conformation_cylinder_limit(
    make_record: Callable[..., DataFile], radial: float, passes: bool
) -> None:
    """A bead within (D - a)/2 = 2.5 of the axis passes for D = 6 and a = 1,
    and one beyond fails: the effective diameter D - a."""
    data = make_record([[radial, 0.0, 0.0]])
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 1.42,
                                 "periodic": (False, False, True),
                                 "dcyl": 6.0, "dmon": 1.0}

    if passes:
        assert check_conformation(data, **arguments).max_radial == radial
    else:
        with pytest.raises(ConformationError,
                           match=r"farther than \(D - a\)/2 = 2.5"):
            check_conformation(data, **arguments)


def test_check_conformation_chunked_minimum_matches_brute_force(
    make_record: Callable[..., DataFile], make_box: Callable[..., Box]
) -> None:
    """The closest non-bonded pair found chunk by chunk equals the one a
    full distance matrix gives, minimum image included: here a pair 0.02
    apart through the z faces, closer than any of 1498 random beads."""
    rng = np.random.default_rng(20261102)
    positions = np.vstack([
        rng.uniform(-10.0, 10.0, size=(1498, 3)),  # two row chunks in all
        [[0.0, 0.0, -9.99], [0.0, 0.0, 9.99]],
    ])
    delta = positions[:, None, :] - positions[None, :, :]
    delta -= 20.0 * np.round(delta / 20.0)
    distances = np.linalg.norm(delta, axis=-1)
    distances[np.diag_indices(1500)] = np.inf
    data = make_record(positions, box=make_box(half=10.0))

    report = check_conformation(data, r0=1.5, ghost_cutoff=1.42,
                                periodic=(True, True, True),
                                min_distance=0.0)

    assert report.min_pair == pytest.approx(distances.min(), rel=1e-12)
    assert report.min_pair == pytest.approx(0.02)  # the pair across the face


def test_check_conformation_reports_every_failure_at_once(
    make_record: Callable[..., DataFile]
) -> None:
    """One error names every failed condition: here a long bond, a close
    non-bonded pair, and an atom outside the box."""
    data = make_record(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 0.6, 0.0], [9.0, 0.0, 0.0]],
        bonds=[(1, 2)])

    with pytest.raises(ConformationError) as caught:
        check_conformation(data, r0=1.5, ghost_cutoff=2.0,
                           periodic=(False, False, False), min_distance=1.0)

    message = str(caught.value)
    assert "1 bond at or above sqrt(0.9) R0" in message
    assert "1 non-bonded pair closer than 1" in message
    assert "1 atom outside the box" in message


def test_check_conformation_report_gives_the_extremes(
    make_record: Callable[..., DataFile]
) -> None:
    """A passing file's report holds the shortest and longest bond, the
    closest non-bonded pair, and the largest distance from the axis."""
    trimer = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0], [0.97, 1.0, 0.0]]
    data = make_record(trimer, bonds=[(1, 2), (2, 3)])

    report = check_conformation(data, r0=1.5, ghost_cutoff=1.42,
                                periodic=(False, False, False),
                                dcyl=6.0, dmon=1.0)

    assert report.n_atoms == 3
    assert (report.min_bond, report.max_bond) == (0.97, 1.0)
    assert report.min_pair == pytest.approx(math.hypot(0.97, 1.0))
    assert report.max_radial == pytest.approx(math.hypot(0.97, 1.0))


def test_check_conformation_rejects_a_periodic_triclinic_box(
    make_record: Callable[..., DataFile], make_box: Callable[..., Box]
) -> None:
    """A triclinic box with a periodic axis raises ``ValueError``."""
    box = make_box(tilt=np.array([1.0, 0.0, 0.0]), triclinic=True)
    data = make_record([[0.0, 0.0, 0.0]], box=box)

    with pytest.raises(ValueError, match="orthogonal box"):
        check_conformation(data, r0=1.5, ghost_cutoff=1.42,
                           periodic=(True, False, False))


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param({"r0": 0.0}, "'r0' must be positive", id="r0"),
        pytest.param({"min_distance": -1.0}, "'min_distance' must be 0",
                     id="negative-min-distance"),
        pytest.param({"periodic": (True, False)}, "three entries",
                     id="two-flags"),
        pytest.param({"dcyl": 6.0}, "given together", id="dcyl-alone"),
        pytest.param({"dcyl": 1.0, "dmon": 1.0}, "'dcyl' must exceed 'dmon'",
                     id="dcyl-not-above-dmon"),
    ],
)
def test_check_conformation_rejects_bad_arguments(
    make_record: Callable[..., DataFile], kwargs: dict[str, Any], match: str
) -> None:
    """An argument out of range raises ``ValueError`` before any check."""
    data = make_record([[0.0, 0.0, 0.0]])
    arguments: dict[str, Any] = {"r0": 1.5, "ghost_cutoff": 1.42,
                                 "periodic": (False, False, False)}
    arguments.update(kwargs)

    with pytest.raises(ValueError, match=match):
        check_conformation(data, **arguments)


# LAMMPS


@pytest.mark.lammps
@pytest.mark.parametrize(
    ("kind", "dcyl"),
    [
        pytest.param("helix", 0.0, id="helix"),
        pytest.param("ring", 0.0, id="ring"),
        pytest.param("ring_angle", 0.0, id="semiflexible-ring"),
        pytest.param("helix", 5.2, id="helix-in-cylinder"),
    ],
)
@pytest.mark.parametrize(
    ("mode", "dt"),
    [
        pytest.param("nve", 0.005, id="nve"),
        pytest.param("nvelang", 0.01, id="langevin"),
        pytest.param("limitlang", 0.005, id="nve-limit"),
    ],
)
def test_lammps_runs_from_a_built_start(
    lammps_executable: str, build_dir: Path, tmp_path: Path,
    start_file: Callable[[str], tuple[Path, str]], kind: str, dcyl: float,
    mode: str, dt: float,
) -> None:
    """LAMMPS runs 10 000 steps of a Kremer-Grest chain from a built start
    under plain nve, Langevin dynamics, and nve/limit, with no error and no
    warning beyond the communication-cutoff heuristic."""
    data, style = start_file(kind)
    log = tmp_path / f"{kind}-{mode}.log"

    subprocess.run(
        [lammps_executable, "-in", str(build_dir / "in.kg"), "-log",
         str(log), "-screen", "none", "-var", "data", str(data), "-var",
         "style", style, "-var", "mode", mode, "-var", "dt", str(dt),
         "-var", "steps", "10000", "-var", "dcyl", str(dcyl)],
        cwd=tmp_path, capture_output=True, text=True, timeout=600,
        check=False)

    lines = log.read_text().splitlines()
    assert any(line.startswith("Loop time") for line in lines)
    assert [line for line in lines if line.startswith("ERROR")] == []
    assert [line for line in lines if line.startswith("WARNING")
            and "Communication cutoff" not in line] == []
