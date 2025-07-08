import pytest
import os
import warnings
from abc import ABC, abstractmethod
from typing import List, ClassVar, Optional, Dict, Tuple, Any
from polyphys.manage.types import GeometryT, TopologyT, GroupT
from polyphys.manage.parser import ParserBase, TwoMonDepCub
from polyphys.manage.utils import (
    number_density_cube,
    number_density_cylinder,
    volume_fraction_cube,
    volume_fraction_cylinder
    )


def test_parserbase_cannot_be_instantiated():
    """
    Ensure ParserBase cannot be instantiated directly due to abstract methods.
    """
    with pytest.raises(TypeError):
        ParserBase("artifact", "ensemble", "bug")


class ParserBaseTestTemplate(ABC):
    """
    A generic test template for testing subclasses of ParserBase.

    Subclasses must define:
    - `parser_class`
    - `test_cases` (list of tuples)
    - `valid_args` (pytest fixture)
    """
    parser_class: ClassVar = None  # must override
    test_cases: ClassVar[List[Tuple[str, str, str, Dict[str, Any]]]] = []
    geometry: ClassVar[Optional[GeometryT]] = None
    topology: ClassVar[Optional[TopologyT]] = None
    groups: ClassVar[Optional[List[GroupT]]] = None

    @abstractmethod
    @pytest.fixture
    def valid_args(self):
        """
        Subclasses must override this fixture to provide valid parser init
        args.
        """

    def test_empty_artifact(self, valid_args):
        """
        Ensure that an empty artifact name raises a ValueError.
        """
        with pytest.raises(ValueError,
                           match="'artifact' cannot be an empty string."):
            self.parser_class(
                "",
                valid_args['lineage'],
                valid_args['group'])

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        """
        Test filename parsing and attribute extraction for each lineage.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            parser = self.parser_class(artifact, lineage, group)
        assert parser.geometry == self.geometry
        assert parser.groups == self.groups
        assert parser.topology == self.topology
        assert parser.physical_attributes == expected["physical_attributes"]
        assert parser.filename == artifact
        assert parser.filepath == "N/A"
        assert parser.group == group
        assert parser.ispath is False
        assert parser.ext == "N/A"
        assert parser.lineage == lineage
        assert parser.name == expected["name"]
        assert parser.project_name == expected["project_name"]
        assert parser.lineage_genealogy == expected["lineage_genealogy"]
        for key, val in expected['attributes'].items():
            if key != "name":
                assert getattr(parser, key) == val, f"{key} mismatch"

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        """
        Test filepath handling using a real temporary file for each artifact.
        """
        ext = '.txt'
        full_path = make_temp_file(artifact + ext)
        with warnings.catch_warnings():
            parser = self.parser_class(full_path, lineage, group, ispath=True)
        assert parser.filename == artifact + ext
        assert parser.filepath == os.path.dirname(full_path)
        assert parser.ispath is True
        assert parser.ext == ext
        assert parser.name == expected["name"]

    @abstractmethod
    def test_computed_attributes(self, valid_args):
        """
        Subclasses must override this method to test computed density and
        volume fraction attributes.
        """

    def test_invalid_group(self, valid_args):
        """
        Ensure that an invalid group raises a ValueError.
        """
        with pytest.raises(ValueError, match="is an invalid option."):
            self.parser_class(valid_args['artifact'], valid_args['lineage'],
                              "invalid")

    def test_invalid_lineage(self, valid_args):
        """
        Ensure that an invalid lineage raises a ValueError.
        """
        with pytest.raises(ValueError, match="is an invalid option."):
            self.parser_class(valid_args['artifact'], "invalid",
                              valid_args['group'])

    @abstractmethod
    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """

    @abstractmethod
    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """

    def test_str_and_repr_output(self, valid_args):
        """
        Check that __str__ and __repr__ methods contain useful info.
        """
        parser = self.parser_class(**valid_args)

        s = str(parser)
        r = repr(parser)

        # Check that key fields are present in __str__
        assert "Artifact:" in s
        assert f"Name: '{parser.filename}'" in s
        assert f"Geometry: '{parser.geometry}'" in s
        assert f"Group: '{parser.group}'" in s
        assert f"Lineage: '{parser.lineage}'" in s
        assert f"Topology: '{parser.topology}'" in s
        assert f"Project: '{parser.project_name}'" in s

        # Check __repr__ content
        assert parser.filename in r
        assert parser.geometry in r
        assert parser.group in r
        assert parser.lineage in r
        assert parser.topology in r
        assert parser.project_name in r


class TestTwoMonDepCub(ParserBaseTestTemplate):
    """
    Tests for the TwoMonDepCub parser class. Covers filename parsing, filepath
    handling, attribute computation, error cases, and representation methods.
    """
    parser_class = TwoMonDepCub
    geometry = "cubic"
    topology = "atomic"
    groups = ["bug", "all"]
    test_cases = [
        (
            "segment",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100"
            "tdump500ens3.j01-bug",
            {
                "lineage_genealogy": ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                "physical_attributes": ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                    "adump": 100, "tdump": 500, "ensemble_id": 3,
                    "segment_id": 1
                },
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
                "lineage_genealogy": ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                "physical_attributes": ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                    "adump": 100, "tdump": 500, "ensemble_id": 3
                },
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500ens3",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "ensemble_long",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500",
            {
                "lineage_genealogy": ['ensemble_long', 'ensemble', 'space'],
                "physical_attributes": ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                    "adump": 100, "tdump": 500
                },
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "ensemble",
            "bug",
            "nm2am5.0ac1.0nc100sd2.0",
            {
                "lineage_genealogy": ['ensemble', 'space'],
                "physical_attributes": [],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "d_sur": 2.0
                },
                "name": "nm2am5.0ac1.0nc100sd2.0",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "space",
            "bug",
            "nm2am5.0ac1.0nc100",
            {
                "lineage_genealogy": ['space'],
                "physical_attributes": [],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100
                },
                "name": "nm2am5.0ac1.0nc100", "project_name": "TwoMonDepCub"
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            "artifact": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100"
            "tdump500",
            "lineage": "ensemble_long",
            "group": "bug"
        }

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        """
        Test computed density and volume fraction attributes.
        """
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == \
            number_density_cube(parser.nmon, parser.dmon, parser.lcube)
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
            volume_fraction_cube(parser.nmon, parser.dmon, parser.lcube)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube)

    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match="has no attribute"):
            self.parser_class(
                "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100xx500",
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match="has no attribute"):
            self.parser_class(
                "5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100dt500",
                valid_args['lineage'],
                valid_args['group'])


class TestTwoMonDepCub(ParserBaseTestTemplate):
    """
    Tests for the TwoMonDepCub parser class. Covers filename parsing, filepath
    handling, attribute computation, error cases, and representation methods.
    """
    parser_class = TwoMonDepCub
    geometry = "cubic"
    topology = "atomic"
    groups = ["bug", "all"]
    test_cases = [
        (
            "segment",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100"
            "tdump500ens3.j01-bug",
            {
                "lineage_genealogy": ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                "physical_attributes": ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                    "adump": 100, "tdump": 500, "ensemble_id": 3,
                    "segment_id": 1
                },
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
                "lineage_genealogy": ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                "physical_attributes": ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                    "adump": 100, "tdump": 500, "ensemble_id": 3
                },
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500ens3",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "ensemble_long",
            "bug",
            "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500",
            {
                "lineage_genealogy": ['ensemble_long', 'ensemble', 'space'],
                "physical_attributes": ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "lcube": 20.0, "d_sur": 2.0, "dt": 0.01, "bdump": 50,
                    "adump": 100, "tdump": 500
                },
                "name": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01"
                "bdump50adump100tdump500",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "ensemble",
            "bug",
            "nm2am5.0ac1.0nc100sd2.0",
            {
                "lineage_genealogy": ['ensemble', 'space'],
                "physical_attributes": [],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100,
                    "d_sur": 2.0
                },
                "name": "nm2am5.0ac1.0nc100sd2.0",
                "project_name": "TwoMonDepCub"
            }
        ),
        (
            "space",
            "bug",
            "nm2am5.0ac1.0nc100",
            {
                "lineage_genealogy": ['space'],
                "physical_attributes": [],
                "attributes": {
                    "dmon": 5.0, "nmon": 2, "dcrowd": 1.0, "ncrowd": 100
                },
                "name": "nm2am5.0ac1.0nc100", "project_name": "TwoMonDepCub"
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            "artifact": "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100"
            "tdump500",
            "lineage": "ensemble_long",
            "group": "bug"
        }

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize("lineage, group, artifact, expected", test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        """
        Test computed density and volume fraction attributes.
        """
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == \
            number_density_cube(parser.nmon, parser.dmon, parser.lcube)
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
            volume_fraction_cube(parser.nmon, parser.dmon, parser.lcube)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube)

    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match="has no attribute"):
            self.parser_class(
                "am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100xx500",
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match="has no attribute"):
            self.parser_class(
                "5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100dt500",
                valid_args['lineage'],
                valid_args['group'])