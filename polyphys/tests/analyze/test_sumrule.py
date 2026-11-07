"""Tests for :mod:`polyphys.analyze.sumrule`.

Every expected value is hand arithmetic on the written ``-mean`` columns and
the sizes encoded in the filenames: crowder diameter ``sig1.5``
(``SumRuleCyl``) or ``ac1.5`` (``SumRuleCubHeteroRing``), cylinder radius
``r10.5``, so ``D = 2 * 10.5 - 1 = 20`` and ``D mod a_c = 0.5``. The sum-rule
tests build a small ``analysis`` phase on disk with one ensAvg CSV per state
point. The spec-taking tests take ``phi_bulk_c`` and ``phi_bulk_m`` from the
analytic volume fraction ``2 n a**3 / (3 Z (D - a)**2)`` of ``n`` spheres of
diameter ``a`` in a cylinder of length ``Z = 2 * lz`` and diameter ``D``, or
check a result against ``space_profiles`` called with the spec's fields. Two
tests state properties instead: the ``Phi`` and ``Rho`` forms of
the rule normalize to one curve, and the normalized sum is one wherever the
rule holds.
"""
import logging
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from polyphys.analyze.sumrule import (
    PROFILE_NORMALIZERS,
    SUM_NORMALIZERS,
    normalize_r,
    normalize_z,
    space_profiles,
    space_profiles_from_spec,
    sum_rule_profiles,
)
from polyphys.manage.parser import (
    SumRuleCubHeteroRing,
    SumRuleCubHeteroRingBiDisperse,
    SumRuleCyl,
)
from polyphys.manage.utils import round_to_multiple
from polyphys.spec import ProjectSpec

# A Pandas4Warning marks a deprecated pandas call; fail on it.
pytestmark = pytest.mark.filterwarnings("error::pandas.errors.Pandas4Warning")


@pytest.fixture
def radial_profile():
    """A scaled radial profile of three bins; its first value is 2."""
    return pd.DataFrame({"bin_center": [0.5, 1.5, 2.5],
                         "rPhiCrd-scale": [2.0, 4.0, 1.0]})


@pytest.fixture
def make_profile_space(tmp_path):
    """Factory: write one '<state>-all-<prop>-ensAvg.csv' per state point.

    Takes the property, the bin centers, one '-mean' column per state
    point, the ensemble_long names (default: two SumRuleCyl points), and
    optionally one '-norm-sem' column per state point. Returns the ensAvg
    directory, ending with '/'.
    """
    cyl_states = [
        "N200epsilon5.0r10.5lz25.0sig1.5nc1000dt0.005bdump2000adump5000",
        "N200epsilon5.0r10.5lz25.0sig1.5nc2000dt0.005bdump2000adump5000",
    ]

    def _make(prop, bin_center, means, states=cyl_states, group="all",
              norm_sems=None):
        ens_avg_dir = tmp_path / "analysis" / f"N200-{group}-ensAvg"
        ens_avg_dir.mkdir(parents=True, exist_ok=True)
        if norm_sems is None:
            norm_sems = [None] * len(means)
        for state, mean, norm_sem in zip(states, means, norm_sems):
            frame = pd.DataFrame({
                "bin_center": bin_center,
                f"{prop}-mean": mean,
                f"{prop}-var": np.full(len(mean), 0.1),
                f"{prop}-sem": np.full(len(mean), 0.05),
            })
            if norm_sem is not None:
                frame[f"{prop}-norm-sem"] = norm_sem
            frame.to_csv(ens_avg_dir / f"{state}-{group}-{prop}-ensAvg.csv",
                         index=False)
        return str(ens_avg_dir) + "/"

    return _make


@pytest.fixture
def cyl_sizes():
    """The species-to-size map of a SumRuleCyl study."""
    return {"Mon": "dmon", "Crd": "dcrowd"}


@pytest.mark.parametrize(
    "normalizer, expected",
    [
        pytest.param("first", [1.0, 2.0, 0.5], id="first"),
        pytest.param("max", [0.5, 1.0, 0.25], id="max"),
        pytest.param("last", [2.0, 4.0, 1.0], id="last"),
        pytest.param(8.0, [0.25, 0.5, 0.125], id="number"),
    ],
)
def test_normalize_r_divides_by_the_chosen_reference(radial_profile,
                                                     normalizer, expected):
    """'first', 'max', 'last', or a number is the divisor of '-scale'."""
    out = normalize_r("rPhiCrd", radial_profile, normalizer)
    # Exact: ratios of small floats
    np.testing.assert_allclose(out["rPhiCrd-norm"], expected, rtol=1e-15)


def test_normalize_r_calls_a_function_with_the_frame(radial_profile):
    """A callable normalizer receives the whole frame and sets the value."""
    out = normalize_r("rPhiCrd", radial_profile,
                      lambda df: df["rPhiCrd-scale"].iloc[:2].mean())
    assert out["rPhiCrd-normalizer"].tolist() == [3.0, 3.0, 3.0]
    np.testing.assert_allclose(out["rPhiCrd-norm"], [2 / 3, 4 / 3, 1 / 3],
                               rtol=1e-15)


@pytest.mark.parametrize(
    "normalizer, message",
    [
        pytest.param("mean", "'mean' is an invalid option", id="string"),
        pytest.param([1.0], "must be 'first', 'last', 'max'", id="list"),
    ],
)
def test_normalize_r_rejects_an_unknown_normalizer(radial_profile,
                                                   normalizer, message):
    """A normalizer that is no literal, number, or callable raises."""
    with pytest.raises(ValueError, match=message):
        normalize_r("rPhiCrd", radial_profile, normalizer)


def test_normalize_r_zero_reference_warns_and_zeroes():
    """A zero reference value warns and sets every '-norm' to zero."""
    profile = pd.DataFrame({"rPhiCrd-scale": [0.0, 1.0, 2.0]})
    with pytest.warns(UserWarning, match="normalizer of 'rPhiCrd' is zero"):
        out = normalize_r("rPhiCrd", profile, "first")
    assert out["rPhiCrd-norm"].tolist() == [0.0, 0.0, 0.0]


def test_normalize_leaves_the_input_unchanged(radial_profile):
    """normalize_r and normalize_z return new frames and add no column."""
    before = radial_profile.copy()
    normalize_r("rPhiCrd", radial_profile)
    normalize_z("rPhiCrd", radial_profile, fold=False)
    pd.testing.assert_frame_equal(radial_profile, before)


def test_normalize_z_folds_an_even_count_onto_abs_z():
    """Mirror bins of an even count average into ceil(n/2) rows at |z|."""
    profile = pd.DataFrame({"bin_center": [-1.5, -0.5, 0.5, 1.5],
                            "zRhoMon-scale": [1.0, 3.0, 3.0, 1.0]})
    out = normalize_z("zRhoMon", profile)
    assert out["bin_center"].tolist() == [0.5, 1.5]
    np.testing.assert_allclose(out["zRhoMon-norm"], [1.0, 1 / 3], rtol=1e-15)


def test_normalize_z_keeps_the_center_bin_of_an_odd_count():
    """The center bin of an odd count is kept as it is, first in |z| order."""
    profile = pd.DataFrame({"bin_center": [-1.0, 0.0, 1.0],
                            "zRhoMon-scale": [2.0, 4.0, 2.0]})
    out = normalize_z("zRhoMon", profile)
    assert out["bin_center"].tolist() == [0.0, 1.0]
    np.testing.assert_allclose(out["zRhoMon-norm"], [1.0, 0.5], rtol=1e-15)


def test_normalize_z_averages_mirror_bins_of_an_asymmetric_profile():
    """Folding averages bin i with bin n-1-i: [1, 2, 3, 4] gives [2.5, 2.5]."""
    profile = pd.DataFrame({"bin_center": [-1.5, -0.5, 0.5, 1.5],
                            "zRhoMon-scale": [1.0, 2.0, 3.0, 4.0]})
    out = normalize_z("zRhoMon", profile)
    # (2 + 3) / 2 and (1 + 4) / 2, each over the maximum 4
    np.testing.assert_allclose(out["zRhoMon-norm"], [2.5 / 4, 2.5 / 4],
                               rtol=1e-15)


def test_normalize_z_without_fold_keeps_every_row():
    """fold=False normalizes by the maximum and keeps all n rows."""
    profile = pd.DataFrame({"bin_center": [-1.5, -0.5, 0.5, 1.5],
                            "zRhoMon-scale": [1.0, 2.0, 3.0, 4.0]})
    out = normalize_z("zRhoMon", profile, fold=False)
    assert len(out) == 4
    np.testing.assert_allclose(out["zRhoMon-norm"], [0.25, 0.5, 0.75, 1.0],
                               rtol=1e-15)


def test_normalize_z_folds_var_and_sem_like_any_column():
    """Folding averages '-var' and '-sem' linearly, as documented."""
    profile = pd.DataFrame({"bin_center": [-1.5, -0.5, 0.5, 1.5],
                            "zRhoMon-scale": [1.0, 1.0, 1.0, 1.0],
                            "zRhoMon-var": [1.0, 1.0, 1.0, 1.0],
                            "zRhoMon-sem": [0.2, 0.4, 0.6, 0.8]})
    out = normalize_z("zRhoMon", profile)
    assert out["zRhoMon-var"].tolist() == [1.0, 1.0]
    np.testing.assert_allclose(out["zRhoMon-sem"], [0.5, 0.5], rtol=1e-15)


@pytest.mark.parametrize(
    "centers",
    [
        pytest.param([-2.0, -1.0, 0.5, 1.5], id="shifted-half"),
        pytest.param([-1.0, 0.1, 1.0], id="off-center-middle"),
        pytest.param([1.5, 0.5, -0.5, -1.5], id="decreasing"),
    ],
)
def test_normalize_z_fold_rejects_bins_not_symmetric_about_zero(centers):
    """Folding pairs bin i with bin n-1-i, which is its mirror only when
    the centers increase and are symmetric about z = 0."""
    profile = pd.DataFrame({"bin_center": centers,
                            "zRhoMon-scale": np.ones(len(centers))})
    with pytest.raises(ValueError, match="symmetric about z = 0"):
        normalize_z("zRhoMon", profile)


def test_normalize_z_fold_accepts_round_off_asymmetry():
    """Centers that miss their mirror by float64 round-off still fold."""
    profile = pd.DataFrame({"bin_center": [-1.5, -0.5, 0.5 + 1e-15, 1.5],
                            "zRhoMon-scale": [1.0, 3.0, 3.0, 1.0]})
    out = normalize_z("zRhoMon", profile)
    assert len(out) == 2


@pytest.mark.parametrize(
    "centers,curve,expected",
    [
        pytest.param([-1.5, -0.5, 0.5, 1.5], [0.1, 0.2, 0.3, 0.4],
                     [0.5, 0.5], id="even"),
        pytest.param([-1.0, 0.0, 1.0], [0.25, 0.5, 0.25], [0.5, 0.5],
                     id="odd"),
    ],
)
def test_normalize_z_folds_the_normalized_curve_onto_abs_z(
    centers, curve, expected
):
    """Folding sums mirror bins of ``-scale-normalized_curve``, so the
    folded curve is the distribution of |z| and still sums to 1:
    P(|z|) = P(z) + P(-z) for z != 0 (analytic)."""
    profile = pd.DataFrame({"bin_center": centers,
                            "zRhoMon-scale": np.ones(len(centers)),
                            "zRhoMon-scale-normalized_curve": curve})
    out = normalize_z("zRhoMon", profile)
    # float64 round-off of one addition:
    np.testing.assert_allclose(out["zRhoMon-scale-normalized_curve"],
                               expected, rtol=1e-15)


def test_normalize_z_fold_needs_bin_center():
    """Folding a frame without 'bin_center' raises ValueError."""
    with pytest.raises(ValueError, match="needs a 'bin_center' column"):
        normalize_z("zRhoMon", pd.DataFrame({"zRhoMon-scale": [1.0, 2.0]}))


def test_default_normalizer_tables_follow_the_research_rule():
    """Profiles use 'first' for cylinder r and 'max' otherwise; sums 'last'."""
    assert PROFILE_NORMALIZERS == {("cylindrical", "r"): "first",
                                   ("cylindrical", "z"): "max",
                                   ("cubic", "r"): "max"}
    assert SUM_NORMALIZERS == {("cylindrical", "r"): "first",
                               ("cylindrical", "z"): "last",
                               ("cubic", "r"): "last"}


def test_space_profiles_scales_phi_by_the_species_size(make_profile_space,
                                                       cyl_sizes):
    """A 'Phi' profile is divided by the species diameter from the name."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    assert space["rPhiCrd-scaler"].tolist() == [1.5, 1.5, 1.5]
    assert space["rPhiCrd-sumrule_constant"].tolist() == [1.5, 1.5, 1.5]
    np.testing.assert_allclose(space["rPhiCrd-scale"],
                               [0.2 / 1.5, 0.3 / 1.5, 0.1 / 1.5], rtol=1e-15)


def test_space_profiles_scales_rho_by_the_size_squared(make_profile_space,
                                                       cyl_sizes):
    """A 'Rho' profile is multiplied by the species diameter squared."""
    database = make_profile_space("rRhoCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Rho", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    assert space["rRhoCrd-scaler"].tolist() == [2.25, 2.25, 2.25]
    np.testing.assert_allclose(space["rRhoCrd-scale"],
                               [0.2 * 2.25, 0.3 * 2.25, 0.1 * 2.25],
                               rtol=1e-15)


def test_space_profiles_phi_and_rho_normalize_to_the_same_curve(
    make_profile_space, cyl_sizes
):
    """With rho = 6 phi / (pi a^3), 'Phi' and 'Rho' give one '-norm' curve.

    Property: phi / a and rho a^2 differ by the constant pi / 6, which the
    normalization cancels.
    """
    phi = np.array([0.2, 0.3, 0.1])
    rho = 6 * phi / (np.pi * 1.5 ** 3)
    make_profile_space("rPhiCrd", [0.5, 1.5, 2.5], [phi])
    database = make_profile_space("rRhoCrd", [0.5, 1.5, 2.5], [rho])
    by_phi = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                            "all", cyl_sizes)
    by_rho = space_profiles(database, "Rho", "r", "Crd", SumRuleCyl, "N*", [],
                            "all", cyl_sizes)
    # float64 round-off of a few multiplications
    np.testing.assert_allclose(by_rho["rRhoCrd-norm"], by_phi["rPhiCrd-norm"],
                               rtol=1e-12)


def test_space_profiles_cylinder_r_normalizes_by_the_axis_value(
    make_profile_space, cyl_sizes
):
    """In a cylinder, an r profile is divided by its r = 0 bin by default."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    np.testing.assert_allclose(space["rPhiCrd-normalizer"], [0.2 / 1.5] * 3,
                               rtol=1e-15)
    np.testing.assert_allclose(space["rPhiCrd-norm"], [1.0, 1.5, 0.5],
                               rtol=1e-15)


def test_space_profiles_cube_r_normalizes_by_the_maximum(make_profile_space):
    """In a cube, an r profile is divided by its maximum and not recentered."""
    state = "al5.0nl5ml125.0ns400ac1.5nc100l100.0dt0.005bdump2000adump5000"
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]], states=[state])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCubHeteroRing,
                           "al*", [], "all",
                           {"Mon": "dmon_small", "Crd": "dcrowd",
                            "Foci": "dmon_large"})
    np.testing.assert_allclose(space["rPhiCrd-norm"], [2 / 3, 1.0, 1 / 3],
                               rtol=1e-15)
    np.testing.assert_allclose(space["bin_center-dcrowd"],
                               [1 / 1.5, 3 / 1.5, 5 / 1.5], rtol=1e-15)
    assert not [c for c in space.columns if "recentered" in c]


def test_space_profiles_z_folds_and_normalizes_by_the_maximum(
    make_profile_space, cyl_sizes
):
    """A z profile is folded about z = 0, divided by its maximum, and not
    recentered."""
    database = make_profile_space("zPhiCrd", [-1.5, -0.5, 0.5, 1.5],
                                  [[0.1, 0.3, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "z", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    assert space["bin_center"].tolist() == [0.5, 1.5]
    np.testing.assert_allclose(space["zPhiCrd-norm"], [1.0, 1 / 3], rtol=1e-15)
    np.testing.assert_allclose(space["bin_center-dcrowd"], [1 / 1.5, 3 / 1.5],
                               rtol=1e-15)
    assert not [c for c in space.columns if "recentered" in c]


def test_space_profiles_recenters_cylinder_r_on_the_wall(make_profile_space,
                                                         cyl_sizes):
    """For r in a cylinder, 2b/a_c - (D mod a_c)/a_c and the normalized
    b - (D mod a_c) are added, with D = 20 and a_c = 1.5 from the name."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    # D mod a_c = 20 mod 1.5 = 0.5, so the shift is 1/3 crowder radii
    np.testing.assert_allclose(space["bin_center-dcrowd"],
                               [1 / 1.5, 3 / 1.5, 5 / 1.5], rtol=1e-15)
    np.testing.assert_allclose(space["bin_center-dcrowd-recentered"],
                               [1 / 1.5 - 1 / 3, 3 / 1.5 - 1 / 3,
                                5 / 1.5 - 1 / 3], rtol=1e-14)
    np.testing.assert_allclose(space["bin_center-recentered-norm"],
                               [0.0, 0.5, 1.0], rtol=1e-15, atol=1e-15)


def test_space_profiles_adds_one_column_set_per_crowder(make_profile_space):
    """Two crowders give 'bin_center-<size>' for each, in a cube."""
    state = ("al5.0nl5ns400ac1.0ab2.0pctc40.0nc100nb50l100.0dt0.005"
             "bdump2000adump5000")
    database = make_profile_space("rPhiCrdOne", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]], states=[state])
    space = space_profiles(
        database, "Phi", "r", "CrdOne", SumRuleCubHeteroRingBiDisperse,
        "al*", [], "all",
        {"Dna": "dmon_small", "CrdOne": "dcrowd_one", "CrdTwo": "dcrowd_two",
         "Foci": "dmon_large"},
        crowders=("CrdOne", "CrdTwo"),
    )
    np.testing.assert_allclose(space["bin_center-dcrowd_one"],
                               [1.0, 3.0, 5.0], rtol=1e-15)
    np.testing.assert_allclose(space["bin_center-dcrowd_two"],
                               [0.5, 1.5, 2.5], rtol=1e-15)


def test_space_profiles_two_crowders_in_a_cylinder_recenter_each(
    make_profile_space
):
    """Two crowder labels of one size in a cylinder give equal recentered
    columns, each suffixed by its size attribute.

    Property: the research code's bidisperse copy raised KeyError here.
    """
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(
        database, "Phi", "r", "Crd", SumRuleCyl, "N*", [], "all",
        {"Mon": "dmon", "Crd": "dcrowd", "Big": "dmon"},
        crowders=("Crd", "Big"),
    )
    assert "bin_center-recentered-norm-dcrowd" in space.columns
    assert "bin_center-recentered-norm-dmon" in space.columns
    # a_c = 1: 2b - (20 mod 1) = 2b
    np.testing.assert_allclose(space["bin_center-dmon-recentered"],
                               [1.0, 3.0, 5.0], rtol=1e-15)


def test_space_profiles_concatenates_files_in_name_order(make_profile_space,
                                                         cyl_sizes):
    """Every matching file contributes its rows, sorted by name, with its
    ensemble_long name as a column."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1], [0.4, 0.2, 0.2]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    assert len(space) == 6
    assert space["ensemble_long"].str.extract(r"nc(\d+)")[0].tolist() == (
        ["1000"] * 3 + ["2000"] * 3)
    np.testing.assert_allclose(space["rPhiCrd-norm"],
                               [1.0, 1.5, 0.5, 1.0, 0.5, 0.5], rtol=1e-15)


def test_space_profiles_adds_attributes_and_rounded_phi_c(make_profile_space,
                                                          cyl_sizes):
    """physical_attrs, phi_bulk_c, and phi_bulk_c_round come from the name."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*",
                           ["nmon", "dcyl"], "all", cyl_sizes, divisor=0.05,
                           round_to=2)
    info = SumRuleCyl(
        "N200epsilon5.0r10.5lz25.0sig1.5nc1000dt0.005bdump2000adump5000",
        "ensemble_long", "all")
    assert space["nmon"].tolist() == [200] * 3
    assert space["dcyl"].tolist() == [20.0] * 3
    assert space["phi_bulk_c"].tolist() == [info.phi_bulk_c] * 3
    assert space["phi_bulk_c_round"].tolist() == [
        round_to_multiple(info.phi_bulk_c, 0.05, 2)] * 3


def test_space_profiles_two_argument_normalizer_sees_the_parser(
    make_profile_space, cyl_sizes
):
    """A callable normalizer receives the frame and the file's parser."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes,
                           normalizer=lambda df, info: info.dcrowd)
    assert space["rPhiCrd-normalizer"].tolist() == [1.5, 1.5, 1.5]


def test_space_profiles_records_and_logs_an_override(make_profile_space,
                                                     cyl_sizes, caplog):
    """An explicit normalizer is used, recorded, and logged with the
    default it replaces."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    caplog.set_level(logging.INFO, logger="polyphys.analyze.sumrule")
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes, normalizer="max")
    np.testing.assert_allclose(space["rPhiCrd-normalizer"], [0.3 / 1.5] * 3,
                               rtol=1e-15)
    assert [r.getMessage() for r in caplog.records] == [
        "Normalizing rPhiCrd profiles with 'max' instead of the cylindrical"
        " r default 'first'."]


def test_space_profiles_saves_with_the_direction_in_the_name(
    make_profile_space, cyl_sizes, tmp_path
):
    """is_save writes '<space>-<group>-<dir><quantity><species>-normalized\
Rescaled-space.csv' in the analysis phase's space directory."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes, is_save=True)
    saved = pd.read_csv(tmp_path / "analysis" / "N200-all-space" /
                        "N200-all-rPhiCrd-normalizedRescaled-space.csv")
    pd.testing.assert_frame_equal(saved, space)


@pytest.mark.parametrize(
    "kwargs, error, message",
    [
        pytest.param({"quantity": "Psi"}, ValueError,
                     "'Psi' is an invalid option", id="quantity"),
        pytest.param({"species": "Hns"}, KeyError,
                     "Species 'Hns' is not in 'species_sizes'", id="species"),
        pytest.param({"crowders": ("Crd", "Hns")}, KeyError,
                     r"Crowders \['Hns'\] are not in", id="crowder"),
        pytest.param({"parser": SumRuleCubHeteroRing, "direction": "z"},
                     ValueError, "No density profile along 'z' in a cubic",
                     id="cubic-z"),
        pytest.param({"hierarchy": "M*"}, ValueError,
                     "No 'M\\*-all-rPhiCrd-ensAvg.csv' file", id="no-file"),
    ],
)
def test_space_profiles_rejects_bad_arguments(make_profile_space, cyl_sizes,
                                              kwargs, error, message):
    """An unknown quantity, species, crowder, (geometry, direction) pair, or
    an empty match raises before any frame is built."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    args = dict(input_database=database, quantity="Phi", direction="r",
                species="Crd", parser=SumRuleCyl, hierarchy="N*",
                physical_attrs=[], group="all", species_sizes=cyl_sizes)
    args.update(kwargs)
    with pytest.raises(error, match=message):
        space_profiles(**args)


@pytest.mark.parametrize(
    "means",
    [
        pytest.param([0.2, 0.3, 0.5], id="sums-to-one"),
        pytest.param([0.2, 0.3, 0.5 + 1e-7], id="within-tolerance"),
    ],
)
def test_space_profiles_warns_on_a_profile_that_sums_to_one(
    make_profile_space, cyl_sizes, means
):
    """A '-mean' column summing to one within 1e-6 warns, naming the file,
    and the file is processed regardless."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5], [means])
    with pytest.warns(UserWarning, match="rPhiCrd-ensAvg.csv' sums to one"):
        space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*",
                               [], "all", cyl_sizes)
    assert len(space) == 3


@pytest.mark.parametrize(
    "means",
    [
        pytest.param([0.2, 0.3, 0.1], id="raw"),
        pytest.param([0.2, 0.3, 0.5 + 1e-5], id="outside-tolerance"),
    ],
)
def test_space_profiles_does_not_warn_on_a_raw_profile(make_profile_space,
                                                       cyl_sizes, means):
    """A '-mean' column that does not sum to one within 1e-6 raises no
    warning."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5], [means])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                       "all", cyl_sizes)


def test_space_profiles_carries_norm_sem_through(make_profile_space,
                                                 cyl_sizes):
    """A file's '-norm-sem' column passes through unchanged under the
    default normalizer."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]],
                                  norm_sems=[[0.0, 0.1, 0.2]])
    space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    assert space["rPhiCrd-norm-sem"].tolist() == [0.0, 0.1, 0.2]


def test_space_profiles_folds_norm_sem_along_z(make_profile_space, cyl_sizes):
    """Along z, '-norm-sem' is averaged with its mirror bin like every
    other numeric column."""
    database = make_profile_space("zPhiCrd", [-1.5, -0.5, 0.5, 1.5],
                                  [[0.1, 0.3, 0.3, 0.1]],
                                  norm_sems=[[0.2, 0.4, 0.6, 0.8]])
    space = space_profiles(database, "Phi", "z", "Crd", SumRuleCyl, "N*", [],
                           "all", cyl_sizes)
    # Exact: (0.4 + 0.6) / 2 and (0.2 + 0.8) / 2 in float64
    np.testing.assert_allclose(space["zPhiCrd-norm-sem"], [0.5, 0.5],
                               rtol=1e-15)


def test_space_profiles_drops_norm_sem_under_an_explicit_normalizer(
    make_profile_space, cyl_sizes
):
    """An explicit normalizer drops a file's '-norm-sem', with a warning
    that names the file and both normalizers."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]],
                                  norm_sems=[[0.0, 0.1, 0.2]])
    with pytest.warns(UserWarning,
                      match=r"Dropping 'rPhiCrd-norm-sem' of '.*rPhiCrd-"
                            r"ensAvg.csv'.*'first', not under 'max'"):
        space = space_profiles(database, "Phi", "r", "Crd", SumRuleCyl, "N*",
                               [], "all", cyl_sizes, normalizer="max")
    assert "rPhiCrd-norm-sem" not in space.columns


@pytest.fixture
def make_cyl_spec():
    """Factory: a SumRuleCyl ProjectSpec; keyword arguments override fields."""
    def _make(**overrides):
        fields = {
            "parser": SumRuleCyl, "group": "bug", "divisor": 0.025,
            "crowder_attr": "phi_bulk_c", "physical_attrs": ("nmon", "dcyl"),
            "species_sizes": {"Mon": "dmon", "Crd": "dcrowd"},
            "tseries_properties": ("gyrTMon",),
        }
        return ProjectSpec(**{**fields, **overrides})

    return _make


@pytest.mark.parametrize(
    "kwargs",
    [
        pytest.param({}, id="defaults"),
        pytest.param({"normalizer": "last", "round_to": 2}, id="overrides"),
    ],
)
def test_space_profiles_from_spec_matches_space_profiles_with_spec_fields(
        make_profile_space, make_cyl_spec, kwargs):
    """The spec-taking table equals space_profiles called with the spec's
    fields one by one."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1], [0.4, 0.2, 0.2]])
    spec = make_cyl_spec(divisor=0.05)
    expected = space_profiles(
        database, "Phi", "r", "Crd", spec.parser, spec.hierarchy,
        list(spec.physical_attrs), "all", spec.species_sizes,
        divisor=spec.divisor, crowder_attr=spec.crowder_attr, **kwargs)
    space = space_profiles_from_spec(database, spec, "Phi", "r", "Crd",
                                     **kwargs)
    pd.testing.assert_frame_equal(space, expected)


@pytest.mark.parametrize(
    "divisor, rounded",
    [
        pytest.param(0.001, 0.329, id="step-0.001"),
        pytest.param(0.025, 0.325, id="step-0.025"),
    ],
)
def test_space_profiles_from_spec_rounds_with_the_spec_divisor(
        make_profile_space, make_cyl_spec, divisor, rounded):
    """phi_bulk_c_round rounds phi_bulk_c to a multiple of the spec's
    divisor."""
    # n = 1000 crowders of a = 2 in D = 20, Z = 50: phi_c = 16000 / 48600
    # = 0.32922, which rounds to 0.329 on a grid of 0.001 and to 0.325 on a
    # grid of 0.025.
    database = make_profile_space(
        "rPhiCrd", [0.5, 1.5, 2.5], [[0.2, 0.3, 0.1]],
        states=["N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000"
                "adump5000"])
    space = space_profiles_from_spec(database, make_cyl_spec(divisor=divisor),
                                     "Phi", "r", "Crd")
    assert space["phi_bulk_c_round"].tolist() == [rounded] * 3


def test_space_profiles_from_spec_rounds_the_spec_crowder_attr(
        make_profile_space, make_cyl_spec):
    """The spec's crowder_attr names the attribute that is added and
    rounded."""
    # N = 200 monomers of a = 1 in D = 20, Z = 50: phi_m = 400 / 54150
    # = 0.0073869, which rounds to 0.007 on a grid of 0.001.
    database = make_profile_space(
        "rPhiCrd", [0.5, 1.5, 2.5], [[0.2, 0.3, 0.1]],
        states=["N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000"
                "adump5000"])
    spec = make_cyl_spec(divisor=0.001, crowder_attr="phi_bulk_m")
    space = space_profiles_from_spec(database, spec, "Phi", "r", "Crd")
    # rel: float64 round-off of the parser's arithmetic is far below the
    # five significant figures of the analytic value.
    assert space["phi_bulk_m"].tolist() == pytest.approx([400 / 54150] * 3,
                                                         rel=1e-12)
    assert space["phi_bulk_c_round"].tolist() == [0.007] * 3


@pytest.mark.parametrize(
    "written, kwargs",
    [
        pytest.param("all", {}, id="default-all"),
        pytest.param("bug", {"group": "bug"}, id="override-bug"),
    ],
)
def test_space_profiles_from_spec_reads_the_group_profiles(
        make_profile_space, make_cyl_spec, written, kwargs):
    """The profiles come from the 'all' group, not spec.group, unless group
    names another."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]], group=written)
    space = space_profiles_from_spec(database, make_cyl_spec(), "Phi", "r",
                                     "Crd", **kwargs)
    # Phi / a_c with a_c = 1.5; rtol: float64 round-off of one division.
    np.testing.assert_allclose(space["rPhiCrd-scale"],
                               np.array([0.2, 0.3, 0.1]) / 1.5, rtol=1e-15)


def test_space_profiles_from_spec_saves_from_a_path(make_profile_space,
                                                    make_cyl_spec, tmp_path):
    """A path-like database with save=True writes the table that
    space_profiles names."""
    database = make_profile_space("rPhiCrd", [0.5, 1.5, 2.5],
                                  [[0.2, 0.3, 0.1]])
    space = space_profiles_from_spec(Path(database), make_cyl_spec(), "Phi",
                                     "r", "Crd", save=True)
    saved = pd.read_csv(tmp_path / "analysis" / "N200-all-space" /
                        "N200-all-rPhiCrd-normalizedRescaled-space.csv")
    pd.testing.assert_frame_equal(saved, space)


@pytest.fixture
def species_profiles():
    """Factory: the scaled profiles of two species over two state points.

    Takes the monomer and crowder '-scale' values of each state point
    (three bins each), optionally their '-norm-sem' values, and returns
    the mapping sum_rule_profiles takes.
    """

    def _make(mon, crd, bin_center=(0.5, 1.5, 2.5), mon_sem=None,
              crd_sem=None):
        keys = {"ensemble_long": ["s1"] * 3 + ["s2"] * 3,
                "bin_center": list(bin_center) * 2}
        profiles = {
            "Mon": pd.DataFrame({**keys, "rPhiMon-scale": np.ravel(mon)}),
            "Crd": pd.DataFrame({**keys, "rPhiCrd-scale": np.ravel(crd)}),
        }
        if mon_sem is not None:
            profiles["Mon"]["rPhiMon-norm-sem"] = np.ravel(mon_sem)
        if crd_sem is not None:
            profiles["Crd"]["rPhiCrd-norm-sem"] = np.ravel(crd_sem)
        return profiles

    return _make


def test_sum_rule_profiles_adds_norm_sem_in_quadrature(species_profiles):
    """'<prop>Sum-norm-sem' is the square root of the sum of the species'
    squared '-norm-sem', bin by bin."""
    profiles = species_profiles(
        mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
        crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]],
        mon_sem=[[0.3, 0.0, 0.06], [0.0, 0.3, 0.0]],
        crd_sem=[[0.4, 0.0, 0.08], [0.0, 0.4, 0.0]],
    )
    total = sum_rule_profiles(profiles, "rPhi", "cylindrical", "r")
    # 3-4-5 triangles: sqrt(0.09 + 0.16) and sqrt(0.0036 + 0.0064);
    # float64 round-off of the squares
    np.testing.assert_allclose(total["rPhiSum-norm-sem"],
                               [0.5, 0.0, 0.1, 0.0, 0.5, 0.0], rtol=1e-12,
                               atol=1e-15)


def test_sum_rule_profiles_leaves_the_error_out_when_a_species_lacks_it(
    species_profiles
):
    """Without '-norm-sem' on every species, the sum has no error column."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]],
                                mon_sem=[[0.3, 0.0, 0.0], [0.0, 0.3, 0.0]])
    total = sum_rule_profiles(profiles, "rPhi", "cylindrical", "r")
    assert "rPhiSum-norm-sem" not in total.columns


def test_sum_rule_profiles_sums_species_per_bin(species_profiles):
    """'<prop>Sum-scale' is the per-bin sum of every species' '-scale'."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    total = sum_rule_profiles(profiles, "rPhi", "cylindrical", "r")
    np.testing.assert_allclose(total["rPhiSum-scale"],
                               [1.0, 1.0, 1.0, 1.0, 1.0, 2.0], rtol=1e-15)
    assert total["ensemble_long"].tolist() == ["s1"] * 3 + ["s2"] * 3


def test_sum_rule_profiles_far_field_is_the_crowder_term(species_profiles):
    """With 'last', the normalizer is the crowder far-field term where the
    monomer term vanishes: Jeon et al. (2016), eqn (7)."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    total = sum_rule_profiles(profiles, "rPhi", "cubic", "r")
    assert total["rPhiSum-normalizer"].tolist() == [1.0] * 3 + [2.0] * 3
    np.testing.assert_allclose(total["rPhiSum-norm"],
                               [1.0, 1.0, 1.0, 0.5, 0.5, 1.0], rtol=1e-15)


def test_sum_rule_profiles_norm_is_one_where_the_rule_holds(species_profiles):
    """Property: when phi_m/a + phi_c/a_c is constant in every state point,
    '<prop>Sum-norm' is identically one whatever the constant."""
    profiles = species_profiles(mon=[[0.7, 0.2, 0.1], [1.5, 0.5, 0.0]],
                                crd=[[0.3, 0.8, 0.9], [0.5, 1.5, 2.0]])
    total = sum_rule_profiles(profiles, "rPhi", "cylindrical", "r")
    assert total["rPhiSum-normalizer"].tolist() == [1.0] * 3 + [2.0] * 3
    np.testing.assert_allclose(total["rPhiSum-norm"], np.ones(6), rtol=1e-15)


def test_sum_rule_profiles_without_by_treats_the_input_as_one(
    species_profiles
):
    """An empty 'by' normalizes the whole input by one reference value."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    total = sum_rule_profiles(profiles, "rPhi", "cubic", "r", by=())
    assert total["rPhiSum-normalizer"].tolist() == [2.0] * 6


def test_sum_rule_profiles_records_and_logs_an_override(species_profiles,
                                                        caplog):
    """An explicit normalizer is used and logged with the default."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    caplog.set_level(logging.INFO, logger="polyphys.analyze.sumrule")
    total = sum_rule_profiles(profiles, "rPhi", "cubic", "r",
                              normalizer=4.0)
    assert total["rPhiSum-normalizer"].tolist() == [4.0] * 6
    assert [r.getMessage() for r in caplog.records] == [
        "Normalizing the rPhi sum with 4.0 instead of the cubic r default "
        "'last'."]


def test_sum_rule_profiles_rejects_a_missing_state_point_key(
        species_profiles):
    """A row with no 'ensemble_long' raises instead of being dropped from
    the normalized sum."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    for frame in profiles.values():
        frame.loc[1, "ensemble_long"] = np.nan

    with pytest.raises(ValueError, match=re.escape(
            "1 row has no value in the 'by' columns ['ensemble_long']")):
        sum_rule_profiles(profiles, "rPhi", "cubic", "r")


def test_sum_rule_profiles_keeps_the_input_row_order(species_profiles):
    """State points whose rows are not contiguous come back in the input
    order, each normalized by its own reference: the result is invariant
    to how the rows are grouped."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    order = [0, 3, 1, 4, 2, 5]
    shuffled = {label: frame.iloc[order].reset_index(drop=True)
                for label, frame in profiles.items()}
    expected = sum_rule_profiles(profiles, "rPhi", "cubic", "r")

    total = sum_rule_profiles(shuffled, "rPhi", "cubic", "r")

    pd.testing.assert_frame_equal(
        total, expected.iloc[order].reset_index(drop=True))


def test_sum_rule_profiles_rejects_mismatched_bins(species_profiles):
    """Species whose 'bin_center' differ in any row raise ValueError."""
    profiles = species_profiles(mon=[[0.8, 0.4, 0.0], [0.6, 0.2, 0.0]],
                                crd=[[0.2, 0.6, 1.0], [0.4, 0.8, 2.0]])
    profiles["Crd"]["bin_center"] += 0.5
    with pytest.raises(ValueError, match="'Crd' profiles do not match"):
        sum_rule_profiles(profiles, "rPhi", "cubic", "r")


@pytest.mark.parametrize(
    "profiles, geometry, direction, message",
    [
        pytest.param({}, "cubic", "r", "'profiles' is empty", id="empty"),
        pytest.param({"Mon": pd.DataFrame()}, "cubic", "z",
                     "No density profile along 'z' in a cubic", id="cubic-z"),
    ],
)
def test_sum_rule_profiles_rejects_bad_arguments(profiles, geometry,
                                                 direction, message):
    """No profiles, or a (geometry, direction) pair without a default,
    raises ValueError."""
    with pytest.raises(ValueError, match=message):
        sum_rule_profiles(profiles, "rPhi", geometry, direction)
