import pytest
from polyphys.manage.parser import TwoMonDepCub


@pytest.mark.parametrize(
    "artifact,lineage,group,expected",
    [
        ("nm2am5.0ac1.0", "space", "bug", {
            "nmon": 2,
            "dmon": 5.0,
            "dcrowd": 1.0,
            "ncrowd": 0,
            "name": "nm2am5.0ac1.0",
            "space": "nm2am5.0ac1.0"
            }),
        ("nm100am1.5ac1.0nc20hl10.0sd2.0dt0.01bdump10adump20tdump100ens0.j00",
         "segment",
         "all", {
            "nmon": 100,
            "dmon": 1.5,
            "dcrowd": 1.0,
            "ncrowd": 20,
            "lcube": 20.0,  # hl=10.0 -> lcube = 2*hl
            "d_sur": 2.0,
            "dt": 0.01,
            "bdump": 10,
            "adump": 20,
            "tdump": 100,
            "ensemble_id": 0,
            "segment_id": 0
            })
    ]
)
def test_twomon_basic_parsing(artifact, lineage, group, expected):
    obj = TwoMonDepCub(artifact, lineage, group)
    for key, val in expected.items():
        assert getattr(obj, key) == val, f"Failed for {key} in {artifact}"


def test_twomon_density_calculation():
    artifact = "nm8am1.0ac1.0nc2hl2.5sd1.0dt0.01bdump10adump20tdump100ens1.j03"
    parser = TwoMonDepCub(artifact, 'segment', 'bug')

    # Expected volume of cube: (2*2.5)^3 = 125
    # rho_bulk_m = 8 / 125 = 0.064
    # phi_bulk_m = 8 * (π/6 * 1^3) / 125 = π/6 * 8 / 125 ≈ 0.0335
    import math
    assert pytest.approx(parser.rho_bulk_m, 0.0001) == 8 / 125
    assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
        (8 * (math.pi / 6)) / 125
    assert pytest.approx(parser.rho_bulk_c, 0.0001) == 2 / 125
    assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
        (2 * (math.pi / 6)) / 125


def test_twomon_invalid_key_warns(caplog):
    artifact = "nm2am5.0XXac1.0"  # Invalid key "XX"
    _ = TwoMonDepCub(artifact, 'space', 'bug')
    assert any("not found in artifact name" in msg for msg in caplog.text)


def test_repr_and_str_output():
    parser = TwoMonDepCub("nm2am5.0ac1.0", 'space', 'bug')
    s = str(parser)
    r = repr(parser)
    assert "Artifact" in s
    assert "geometry" in r.lower()
    assert "TwoMonDepCub" in parser.project_name
