import pytest
import math
import os
from polyphys.manage.parser import TwoMonDepCub


from collections import OrderedDict
from typing import List, Dict
from polyphys.manage.parser import ParserBase


class DummyParser(ParserBase):
    """
    A mock parser class for unit testing ParserBase.
    Implements all required abstract methods and provides mock data.
    """
    _geometry = 'cube'
    _topology = 'none'
    _groups = ['bug', 'all']
    _genealogy_attributes = {
        'segment': OrderedDict({
            'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'lcube': 'l', 'd_sur': 'sd', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'tdump': 'tdump', 'ensemble_id': 'ens',
            'segment_id': 'j'
        }),
        'whole': OrderedDict({
            'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'lcube': 'l', 'd_sur': 'sd', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'tdump': 'tdump', 'ensemble_id': 'ens'
        }),
        'ensemble_long': OrderedDict({
            'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'lcube': 'l', 'd_sur': 'sd', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'tdump': 'tdump'
        }),
        'ensemble': OrderedDict({
            'nmon': 'nm', 'dmon': 'am', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'd_sur': 'sd'
        }),
        'space': OrderedDict({
            'nmon': 'nm', 'dmon': 'am', 'dcrowd': 'ac', 'ncrowd': 'nc'
        }),
    }
    _project_attributes = {
        'segment': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': [],
        'space': []
    }

    def _initiate_attributes(self) -> None:
        """
        Stub method for initiating project attributes.
        """
        pass

    def _parse_name(self) -> None:
        """
        Simulates parsing a name by hardcoding test attributes.
        These mimic attributes extracted from a filename.
        """
        self.dmon = 5.0
        self.nmon = 2
        self.dcrowd = 1.0
        self.ncrowd = 100
        self.lcube = 20.0
        self.d_sur = 0.5
        self.dt = 0.01
        self.bdump = "bd1"
        self.adump = "ad1"
        self.tdump = "td1"
        self.ensemble_id = 10
        self.segment_id = 1

    def _dependant_attributes(self) -> None:
        """
        Stub method for computing dependent system attributes.
        """
        pass

class TestTwoMonDepCub:
    """
    Tests for the TwoMonDepCub parser class. Covers filename parsing, filepath
    handling, attribute computation, error cases, and representation methods.
    """

    test_cases = [
        (
            "segment",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100"
            "tdump500ens3.j01-bug",
            {
                "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                "adump": 100, "tdump": 500, "ensemble_id": 3, "segment_id": 1,
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500ens3.j01",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "whole",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500ens3",
            {
                "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                "adump": 100, "tdump": 500, "ensemble_id": 3,
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500ens3", "project_name": "TwoMonDepCub"
            }
        ),
        (
            "ensemble_long",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500",
            {
                "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                "adump": 100, "tdump": 500,
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500", "project_name": "TwoMonDepCub"
            }
        ),
        (
            "ensemble",
            "bug",
            "nm2am5.0ac1.0nc100sd2.0",
            {
                "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                "d_sur": 2.0, "name": "nm2am5.0ac1.0nc100sd2.0",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "space",
            "bug",
            "nm2am5.0ac1.0nc100",
            {
                "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                "name": "nm2am5.0ac1.0nc100", "project_name": "TwoMonDepCub"
            }
        ),
    ]

    @pytest.mark.parametrize("lineage, group, artifact, expected,", test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        """
        Test filename parsing and attribute extraction for each lineage.
        """
        parser = TwoMonDepCub(artifact, lineage, group)

        assert parser.lineage == lineage
        assert parser.group == group
        assert parser.filename == artifact
        assert parser.filepath == "N/A"
        assert parser.name == expected["name"]
        assert parser.project_name == expected["project_name"]

        for key, val in expected.items():
            if key != "name":
                assert getattr(parser, key) == val, f"{key} mismatch"

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_file_path_parsing(
        self,
        make_temp_file,
        lineage,
        group,
        artifact,
        expected
    ):
        """
        Test filepath handling using a real temporary file for each artifact.
        """
        new_artifact = artifact + '.txt'
        full_path = make_temp_file(new_artifact)
        parser = TwoMonDepCub(full_path, lineage, group)

        assert parser.filename == new_artifact
        assert parser.filepath in full_path
        assert os.path.samefile(parser.filepath, os.path.dirname(full_path))
        assert parser.name == expected["name"]
        assert parser.project_name == expected["project_name"]

        for key, val in expected.items():
            if key != "name":
                assert getattr(parser, key) == val, f"{key} mismatch"

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_computed_attributes(self):
        """
        Test computed density and volume fraction attributes.
        """
        artifact = "nm8am1.0ac1.0nc2hl2.5sd1.0dt0.01bdump10adump20tdump100ens1.j03"
        parser = TwoMonDepCub(artifact, 'segment', 'bug')

        assert parser.lcube == 5.0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == 8 / 125
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == (8 * (math.pi / 6)) / 125
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == 2 / 125
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == (2 * (math.pi / 6)) / 125

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_invalid_group(self):
        """
        Ensure that an invalid group raises a ValueError.
        """
        with pytest.raises(ValueError, match="group"):
            TwoMonDepCub("nm2am5.0ac1.0", 'space', 'invalid')

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_invalid_lineage(self):
        """
        Ensure that an invalid lineage raises a ValueError.
        """
        with pytest.raises(ValueError, match="lineage"):
            TwoMonDepCub("nm2am5.0ac1.0", 'wrong', 'bug')

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_invalid_key_warns(self, caplog):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        _ = TwoMonDepCub("nm2am5.0XXac1.0", 'space', 'bug')
        assert any("not found" in msg for msg in caplog.text)

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_str_and_repr_output(self):
        """
        Check that __str__ and __repr__ methods contain useful info.
        """
        parser = TwoMonDepCub("nm2am5.0ac1.0", 'space', 'bug')
        assert "Artifact" in str(parser)
        assert "geometry" in repr(parser).lower()
        assert "TwoMonDepCub" in parser.project_name
