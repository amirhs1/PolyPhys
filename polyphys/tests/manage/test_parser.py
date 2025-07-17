import pytest
import os
import warnings
from abc import ABC, abstractmethod
from typing import List, ClassVar, Optional, Dict, Tuple, Any
from polyphys.manage.types import GeometryT, TopologyT, GroupT
from polyphys.manage.parser import (
    ParserBase, TwoMonDepCub, SumRuleCyl, SumRuleCubHeteroRing,
    SumRuleCubHeteroLinear, TransFociCub, TransFociCyl, HnsCub, HnsCyl
    )
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
        ParserBase('artifact', 'ensemble', 'bug')


class ParserBaseTestTemplate(ABC):
    """
    A generic test template for testing subclasses of `'ParserBase'`.

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

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        """
        Test filename parsing and attribute extraction for each lineage.
        """
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            parser = self.parser_class(artifact, lineage, group)
        assert parser.geometry == self.geometry
        assert parser.groups == self.groups
        assert parser.topology == self.topology
        assert parser.physical_attributes == expected['physical_attributes']
        assert parser.filename == artifact
        assert parser.filepath == 'N/A'
        assert parser.group == group
        assert parser.ispath is False
        assert parser.ext == 'N/A'
        assert parser.lineage == lineage
        assert parser.name == expected['name']
        assert parser.project_name == expected['project_name']
        assert parser.lineage_genealogy == expected['lineage_genealogy']
        for key, val in expected['attributes'].items():
            if key != 'name':
                assert getattr(parser, key) == val, f"{key} mismatch"

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
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
        assert parser.name == expected['name']

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
        with pytest.raises(ValueError, match='is an invalid option.'):
            self.parser_class(valid_args['artifact'], valid_args['lineage'],
                              'invalid')

    def test_invalid_lineage(self, valid_args):
        """
        Ensure that an invalid lineage raises a ValueError.
        """
        with pytest.raises(ValueError, match='is an invalid option.'):
            self.parser_class(valid_args['artifact'], 'invalid',
                              valid_args['group'])

    @abstractmethod
    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger an attribute
        error.
        """

    @abstractmethod
    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure artifact with no heading keyword (a heading digit) in artifact
        name trigger an attribute error.
        """

    @abstractmethod
    def test_single_dot_in_keyword(self, valid_args):
        """
        Ensure passing `.` as a keyword value in artifact name trigger riases
        an attribute error.
        """

    @abstractmethod
    def test_trailing_dot_in_keyword(self, valid_args):
        """
        Ensure passing a keyword value with traling dot, e.g. `10.` in
        artifact name handle correctly.
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
    parser_class = TwoMonDepCub
    geometry = 'cubic'
    topology = 'atomic'
    groups = ['bug', 'all']
    test_cases = [
        (
            'segment',
            'bug',
            'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
            'tdump500ens3.j01.bug',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                'attributes': {
                    'dmon': 5.0, 'nmon': 2, 'dcrowd': 1.0, 'ncrowd': 100,
                    'lcube': 20.0, 'd_sur': 2.0, 'dt': 0.01, 'bdump': 50,
                    'adump': 100, 'tdump': 500, 'ensemble_id': 3,
                    'segment_id': 1
                },
                'name': 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01'
                'bdump50adump100tdump500ens3.j01',
                'project_name': 'TwoMonDepCub'
            }
        ),
        (
            'whole',
            'bug',
            'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500'
            'ens3.bug',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                'attributes': {
                    'dmon': 5.0, 'nmon': 2, 'dcrowd': 1.0, 'ncrowd': 100,
                    'lcube': 20.0, 'd_sur': 2.0, 'dt': 0.01, 'bdump': 50,
                    'adump': 100, 'tdump': 500, 'ensemble_id': 3
                },
                'name': 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01'
                'bdump50adump100tdump500ens3',
                'project_name': 'TwoMonDepCub'
            }
        ),
        (
            'ensemble_long',
            'bug',
            'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': ['phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                'attributes': {
                    'dmon': 5.0, 'nmon': 2, 'dcrowd': 1.0, 'ncrowd': 100,
                    'lcube': 20.0, 'd_sur': 2.0, 'dt': 0.01, 'bdump': 50,
                    'adump': 100, 'tdump': 500
                },
                'name': 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01'
                'bdump50adump100tdump500',
                'project_name': 'TwoMonDepCub'
            }
        ),
        (
            'ensemble',
            'bug',
            'nm2am5.0ac1.0nc100sd2.0',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': [],
                'attributes': {
                    'dmon': 5.0, 'nmon': 2, 'dcrowd': 1.0, 'ncrowd': 100,
                    'd_sur': 2.0
                },
                'name': 'nm2am5.0ac1.0nc100sd2.0',
                'project_name': 'TwoMonDepCub'
            }
        ),
        (
            'space',
            'bug',
            'nm2am5.0ac1.0nc100',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': [],
                'attributes': {
                    'dmon': 5.0, 'nmon': 2, 'dcrowd': 1.0, 'ncrowd': 100
                },
                'name': 'nm2am5.0ac1.0nc100', 'project_name': 'TwoMonDepCub'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
            'tdump500',
            'lineage': 'ensemble_long',
            'group': 'bug'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == \
            number_density_cube(parser.nmon, parser.dmon, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
            volume_fraction_cube(parser.nmon, parser.dmon, parser.lcube,
                                 pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                 pbc=True)

    def test_invalid_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100xx500',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'am.nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'am5.0nm2ac1.nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.dcrowd == 1.0
        assert parser.ncrowd == 100


class TestSumRuleCyl(ParserBaseTestTemplate):
    parser_class = SumRuleCyl
    geometry = 'cylindrical'
    topology = 'linear'
    groups = ['bug', 'all']
    test_cases = [
        (
            'segment',
            'bug',
            'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000adump'
            '10000ens2.j02.bug',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': ['dmon', 'phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                'attributes': {
                    'nmon': 2000, 'epsilon': 5.0, 'dcyl': 20, 'lcyl': 600.0,
                    'dcrowd': 1.0, 'ncrowd': 117000, 'dt': 0.005,
                    'bdump': 5000, 'adump': 10000, 'ensemble_id': 2,
                    'segment_id': 2
                },
                'name': 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                'bdump5000adump10000ens2.j02',
                'project_name': 'SumRuleCyl'
            }
        ),
        (
            'whole',
            'bug',
            'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000adump'
            '10000ens2.bug',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': ['dmon', 'phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                'attributes': {
                    'nmon': 2000, 'epsilon': 5.0, 'dcyl': 20, 'lcyl': 600.0,
                    'dcrowd': 1.0, 'ncrowd': 117000, 'dt': 0.005,
                    'bdump': 5000, 'adump': 10000, 'ensemble_id': 2
                },
                'name': 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                'bdump5000adump10000ens2',
                'project_name': 'SumRuleCyl'
            }
        ),
        (
            'ensemble_long',
            'bug',
            'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000'
            'adump10000',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': ['dmon', 'phi_bulk_m', 'rho_bulk_m',
                                        'phi_bulk_c', 'rho_bulk_c'],
                'attributes': {
                    'nmon': 2000, 'epsilon': 5.0, 'dcyl': 20, 'lcyl': 600.0,
                    'dcrowd': 1.0, 'ncrowd': 117000, 'dt': 0.005,
                    'bdump': 5000, 'adump': 10000
                },
                'name': 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                'bdump5000adump10000',
                'project_name': 'SumRuleCyl'
            }
        ),
        (
            'ensemble',
            'bug',
            'N2000D20.0ac1.0nc117000',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon'],
                'attributes': {
                    'dmon': 1.0, 'nmon': 2000, 'dcrowd': 1.0, 'ncrowd': 117000,
                    'dcyl': 20
                },
                'name': 'N2000D20.0ac1.0nc117000',
                'project_name': 'SumRuleCyl'
            }
        ),
        (
            'space',
            'bug',
            'N2000D20.0ac1.0',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon'],
                'attributes': {
                    'dmon': 1.0, 'nmon': 2000, 'dcrowd': 1.0, 'dcyl': 20
                },
                'name': 'N2000D20.0ac1.0', 'project_name': 'SumRuleCyl'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
            'bdump5000adump10000',
            'lineage': 'ensemble_long',
            'group': 'bug'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        """
        Test computed density and volume fraction attributes.
        """
        parser = self.parser_class(**valid_args)

        assert parser.dcyl > 0
        assert parser.lcyl > 0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == \
            number_density_cylinder(parser.nmon, parser.dmon, parser.lcyl,
                                    parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
            volume_fraction_cylinder(parser.nmon, parser.dmon, parser.lcyl,
                                     parser.dcyl, pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cylinder(parser.ncrowd, parser.dcrowd, parser.lcyl,
                                    parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cylinder(parser.ncrowd, parser.dcrowd, parser.lcyl,
                                     parser.dcyl, pbc=True)

    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'N2000epsilon5.0r10.5lz300.0sig1.0nc117000XX0.005'
                'bdump5000adump10000',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                'bdump5000adump10000',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'N.epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                'bdump5000adump10000',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'N2000epsilon5.r10.5lz300.0sig1.0nc117000dt0.005'
            'bdump5000adump10000',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.epsilon == 5.0
        assert parser.dcyl == 20.0


class TestSumRuleCubHeteroRing(ParserBaseTestTemplate):
    parser_class = SumRuleCubHeteroRing
    geometry = 'cubic'
    topology = 'ring'
    groups = ['bug', 'all']
    test_cases = [
        (
            'segment',
            'all',
            'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4.j10.'
            'ring.all',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'mmon_large': 125.0,
                    'nmon_small': 400, 'dcrowd': 3.0, 'ncrowd': 7735,
                    'lcube': 90.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000, 'ensemble_id': 4, 'segment_id': 10
                },
                'name': 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                'ens4.j10.ring',
                'project_name': 'SumRuleCubHeteroRing'
            }
        ),
        (
            'whole',
            'all',
            'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4.ring'
            '.all',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'mmon_large': 125.0,
                    'nmon_small': 400, 'dcrowd': 3.0, 'ncrowd': 7735,
                    'lcube': 90.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000, 'ensemble_id': 4
                },
                'name': 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                'ens4.ring',
                'project_name': 'SumRuleCubHeteroRing'
            }
        ),
        (
            'ensemble_long',
            'all',
            'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'mmon_large': 125.0,
                    'nmon_small': 400, 'dcrowd': 3.0, 'ncrowd': 7735,
                    'lcube': 90.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000,
                },
                'name': 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000'
                'adump5000',
                'project_name': 'SumRuleCubHeteroRing'
            }
        ),
        (
            'ensemble',
            'all',
            'ns400nl5al5ac3nc7735',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'nmon_small': 400,
                    'dcrowd': 3.0, 'ncrowd': 7735
                },
                'name': 'ns400nl5al5ac3nc7735',
                'project_name': 'SumRuleCubHeteroRing'
            }
        ),
        (
            'space',
            'all',
            'ns400nl5al5ac3',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'nmon_small': 400,
                    'dcrowd': 3.0
                },
                'name': 'ns400nl5al5ac3',
                'project_name': 'SumRuleCubHeteroRing'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
            'bdump5000adump10000',
            'lineage': 'ensemble_long',
            'group': 'bug'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m_large, 0.0001) == \
            number_density_cube(parser.nmon_large, parser.dmon_large,
                                parser.lcube, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_large, 0.0001) == \
            volume_fraction_cube(parser.nmon_large, parser.dmon_large,
                                 parser.lcube, pbc=True)
        assert pytest.approx(parser.rho_bulk_m_small, 0.0001) == \
            number_density_cube(parser.nmon_small, parser.dmon_small,
                                parser.lcube, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_small, 0.0001) == \
            volume_fraction_cube(parser.nmon_small, parser.dmon_small,
                                 parser.lcube, pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                 pbc=True)

    def test_invalid_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'al5nl5ml125ns400ac3XX7735l45dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'al.nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'al5nl5ml125.ns400ac3nc7735l45dt0.005bdump2000adump5000',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.mmon_large == 125.0
        assert parser.nmon_small == 400


class TestSumRuleCubHeteroLinear(ParserBaseTestTemplate):
    parser_class = SumRuleCubHeteroLinear
    geometry = 'cubic'
    topology = 'linear'
    groups = ['bug', 'all']
    test_cases = [
        (
            'segment',
            'all',
            'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000ens1.j11.'
            'linear.all',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 6.0, 'nmon_large': 5, 'mmon_large': 216.0,
                    'nmon_small': 400, 'dcrowd': 3.0, 'ncrowd': 7735,
                    'lcube': 90.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000, 'ensemble_id': 1, 'segment_id': 11
                },
                'name': 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                'ens1.j11.linear',
                'project_name': 'SumRuleCubHeteroLinear'
            }
        ),
        (
            'whole',
            'all',
            'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000ens1.linear'
            '.all',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 6.0, 'nmon_large': 5, 'mmon_large': 216.0,
                    'nmon_small': 400, 'dcrowd': 3.0, 'ncrowd': 7735,
                    'lcube': 90.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000, 'ensemble_id': 1
                },
                'name': 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                'ens1.linear',
                'project_name': 'SumRuleCubHeteroLinear'
            }
        ),
        (
            'ensemble_long',
            'all',
            'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 6.0, 'nmon_large': 5, 'mmon_large': 216.0,
                    'nmon_small': 400, 'dcrowd': 3.0, 'ncrowd': 7735,
                    'lcube': 90.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000
                },
                'name': 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000'
                'adump5000',
                'project_name': 'SumRuleCubHeteroLinear'
            }
        ),
        (
            'ensemble',
            'all',
            'ns400nl5al6ac3nc7735',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dmon_large': 6.0, 'nmon_large': 5, 'nmon_small': 400,
                    'dcrowd': 3.0, 'ncrowd': 7735
                },
                'name': 'ns400nl5al6ac3nc7735',
                'project_name': 'SumRuleCubHeteroLinear'
            }
        ),
        (
            'space',
            'all',
            'ns400nl5al6ac3',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dmon_large': 6.0, 'nmon_large': 5, 'nmon_small': 400,
                    'dcrowd': 3.0
                },
                'name': 'ns400nl5al6ac3',
                'project_name': 'SumRuleCubHeteroLinear'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
            'bdump5000adump10000',
            'lineage': 'ensemble_long',
            'group': 'bug'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m_large, 0.0001) == \
            number_density_cube(parser.nmon_large, parser.dmon_large,
                                parser.lcube, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_large, 0.0001) == \
            volume_fraction_cube(parser.nmon_large, parser.dmon_large,
                                 parser.lcube, pbc=True)
        assert pytest.approx(parser.rho_bulk_m_small, 0.0001) == \
            number_density_cube(parser.nmon_small, parser.dmon_small,
                                parser.lcube, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_small, 0.0001) == \
            volume_fraction_cube(parser.nmon_small, parser.dmon_small,
                                 parser.lcube, pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                 pbc=True)

    def test_invalid_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'al6nl5ml216XX400ac3nc7735l45dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'al.nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'al6nl5ml216.ns400ac3nc7735l45dt0.005bdump2000adump5000',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.mmon_large == 216.0
        assert parser.nmon_small == 400


class TestTransFociCyl(ParserBaseTestTemplate):
    parser_class = TransFociCyl
    geometry = 'cylindrical'
    topology = 'ring'
    groups = ['bug', 'all']
    test_cases = [
        (
            'segment',
            'bug',
            'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
            'bdump2000adump5000ens8.j05.ring.bug',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'epsilon_small': 5.0, 'epsilon_large': 5.0, 'dcyl': 20.0,
                    'dmon_large': 3.0, 'nmon_large': 5, 'mmon_large': 27.0,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 8340,
                    'lcyl': 139.0, 'dt': 0.005, 'bdump': 2000, 'adump': 5000,
                    'ensemble_id': 8, 'segment_id': 5
                },
                'name': 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                'bdump2000adump5000ens8.j05.ring',
                'project_name': 'TransFociCyl'
            }
        ),
        (
            'whole',
            'bug',
            'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
            'bdump2000adump5000ens8.ring.bug',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'epsilon_small': 5.0, 'epsilon_large': 5.0, 'dcyl': 20.0,
                    'dmon_large': 3.0, 'nmon_large': 5, 'mmon_large': 27.0,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 8340,
                    'lcyl': 139.0, 'dt': 0.005, 'bdump': 2000, 'adump': 5000,
                    'ensemble_id': 8
                },
                'name': 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                'bdump2000adump5000ens8.ring',
                'project_name': 'TransFociCyl'
            }
        ),
        (
            'ensemble_long',
            'bug',
            'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
            'bdump2000adump5000',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'epsilon_small': 5.0, 'epsilon_large': 5.0, 'dcyl': 20.0,
                    'dmon_large': 3.0, 'nmon_large': 5, 'mmon_large': 27.0,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 8340,
                    'lcyl': 139.0, 'dt': 0.005, 'bdump': 2000, 'adump': 5000
                },
                'name': 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                'bdump2000adump5000',
                'project_name': 'TransFociCyl'
            }
        ),
        (
            'ensemble',
            'bug',
            'ns400nl5al3D20.0ac1nc8340',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dcyl': 20.0, 'dmon_large': 3.0, 'nmon_large': 5,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 8340
                },
                'name': 'ns400nl5al3D20.0ac1nc8340',
                'project_name': 'TransFociCyl'
            }
        ),
        (
            'space',
            'bug',
            'ns400nl5al3D20.0ac1',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dcyl': 20.0, 'dmon_large': 3.0, 'nmon_large': 5,
                    'nmon_small': 400, 'dcrowd': 1.0
                },
                'name': 'ns400nl5al3D20.0ac1', 'project_name': 'TransFociCyl'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
            'bdump2000adump5000',
            'lineage': 'ensemble_long',
            'group': 'bug'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        """
        Test computed density and volume fraction attributes.
        """
        parser = self.parser_class(**valid_args)

        assert parser.dcyl > 0
        assert parser.lcyl > 0
        assert pytest.approx(parser.rho_bulk_m_large, 0.0001) == \
            number_density_cylinder(parser.nmon_large, parser.dmon_large,
                                    parser.lcyl, parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_large, 0.0001) == \
            volume_fraction_cylinder(parser.nmon_large, parser.dmon_large,
                                     parser.lcyl,
                                     parser.dcyl, pbc=True)
        assert pytest.approx(parser.rho_bulk_m_small, 0.0001) == \
            number_density_cylinder(parser.nmon_small, parser.dmon_small,
                                    parser.lcyl, parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_small, 0.0001) == \
            volume_fraction_cylinder(parser.nmon_small, parser.dmon_small,
                                     parser.lcyl, parser.dcyl, pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cylinder(parser.ncrowd, parser.dcrowd, parser.lcyl,
                                    parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cylinder(parser.ncrowd, parser.dcrowd, parser.lcyl,
                                     parser.dcyl, pbc=True)

    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'epss5epsl5r10.5al3nl5xx27ns400ac1nc8340lz69.5dt0.005'
                'bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                'bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'epss.epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                'bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'epss5epsl5r10.5al3.nl5ml27.ns400ac1nc8340lz69.5dt0.005'
            'bdump2000adump5000',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.mmon_large == 27.0
        assert parser.nmon_small == 400


class TestTransFociCub(ParserBaseTestTemplate):
    parser_class = TransFociCub
    geometry = 'cubic'
    topology = 'ring'
    groups = ['bug', 'all']
    test_cases = [
        (
            'segment',
            'all',
            'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000ens5.j13'
            '.ring.all',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'mmon_large': 125.0,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 154699,
                    'lcube': 60.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000, 'ensemble_id': 5, 'segment_id': 13
                },
                'name': 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000'
                'adump5000ens5.j13.ring',
                'project_name': 'TransFociCub'
            }
        ),
        (
            'whole',
            'all',
            'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000ens5'
            '.ring.all',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'mmon_large': 125.0,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 154699,
                    'lcube': 60.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000, 'ensemble_id': 5
                },
                'name': 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000'
                'adump5000ens5.ring',
                'project_name': 'TransFociCub'
            }
        ),
        (
            'ensemble_long',
            'all',
            'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': [
                    'dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'
                    ],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'mmon_large': 125.0,
                    'nmon_small': 400, 'dcrowd': 1.0, 'ncrowd': 154699,
                    'lcube': 60.0, 'dt': 0.005, 'bdump': 2000,
                    'adump': 5000
                },
                'name': 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000'
                'adump5000',
                'project_name': 'TransFociCub'
            }
        ),
        (
            'ensemble',
            'all',
            'ns400nl5al5ac1nc154699',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'nmon_small': 400,
                    'dcrowd': 1.0, 'ncrowd': 154699
                },
                'name': 'ns400nl5al5ac1nc154699',
                'project_name': 'TransFociCub'
            }
        ),
        (
            'space',
            'all',
            'ns400nl5al5ac1',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon_small', 'mmon_small', 'mcrowd'],
                'attributes': {
                    'dmon_large': 5.0, 'nmon_large': 5, 'nmon_small': 400,
                    'dcrowd': 1.0
                },
                'name': 'ns400nl5al5ac1',
                'project_name': 'TransFociCub'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000'
            'adump5000',
            'lineage': 'ensemble_long',
            'group': 'bug'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m_large, 0.0001) == \
            number_density_cube(parser.nmon_large, parser.dmon_large,
                                parser.lcube, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_large, 0.0001) == \
            volume_fraction_cube(parser.nmon_large, parser.dmon_large,
                                 parser.lcube, pbc=True)
        assert pytest.approx(parser.rho_bulk_m_small, 0.0001) == \
            number_density_cube(parser.nmon_small, parser.dmon_small,
                                parser.lcube, pbc=True)
        assert pytest.approx(parser.phi_bulk_m_small, 0.0001) == \
            volume_fraction_cube(parser.nmon_small, parser.dmon_small,
                                 parser.lcube, pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                 pbc=True)

    def test_invalid_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'al5nl5ml125ns400ac1nc154699XX30dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'al5nl.ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'al5nl5ml125.ns400ac1nc154699l30dt0.005bdump2000adump5000',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.mmon_large == 125.0
        assert parser.nmon_small == 400


class TestHnsCub(ParserBaseTestTemplate):
    parser_class = HnsCub
    geometry = 'cubic'
    topology = 'ring'
    groups = ['nucleoid', 'all']
    test_cases = [
        (
            'segment',
            'nucleoid',
            'N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': [
                    'dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump', 'eps_hm'
                    ],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'nhns': 12, 'dcrowd': 2.0,
                    'lcube': 50, 'eps_hc': 1.0, 'ncrowd': 11937,
                    'ensemble_id': 2, 'segment_id': 5
                },
                'name': 'N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring',
                'project_name': 'HnsCub'
            }
        ),
        (
            'whole',
            'nucleoid',
            'N200kbmm2nh12ac2l25epshc1nc11937ens2.ring.nucleoid',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': [
                    'dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump', 'eps_hm'
                    ],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'nhns': 12, 'dcrowd': 2.0,
                    'lcube': 50, 'eps_hc': 1.0, 'ncrowd': 11937,
                    'ensemble_id': 2
                },
                'name': 'N200kbmm2nh12ac2l25epshc1nc11937ens2.ring',
                'project_name': 'HnsCub'
            }
        ),
        (
            'ensemble_long',
            'nucleoid',
            'N200kbmm2nh12ac2l25epshc1nc11937',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': [
                    'dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump', 'eps_hm'
                    ],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'nhns': 12, 'dcrowd': 2.0,
                    'lcube': 50, 'eps_hc': 1.0, 'ncrowd': 11937
                },
                'name': 'N200kbmm2nh12ac2l25epshc1nc11937',
                'project_name': 'HnsCub'
            }
        ),
        (
            'ensemble',
            'nucleoid',
            'N200kbmm2nh12ac2.0epshc1nc11937',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon', 'dhns', 'dt', 'ndump', 'adump',
                                        'eps_hm'],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'nhns': 12, 'dcrowd': 2.0,
                    'eps_hc': 1.0, 'ncrowd': 11937
                },
                'name': 'N200kbmm2nh12ac2.0epshc1nc11937',
                'project_name': 'HnsCub'
            }
        ),
        (
            'space',
            'nucleoid',
            'N200kbmm2nh12ac2.0epshc1',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon', 'dhns', 'dt', 'ndump', 'adump',
                                        'eps_hm'],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'nhns': 12, 'dcrowd': 2.0,
                    'eps_hc': 1.0
                },
                'name': 'N200kbmm2nh12ac2.0epshc1',
                'project_name': 'HnsCub'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'N200kbmm2nh12ac2l25epshc1nc11937',
            'lineage': 'ensemble_long',
            'group': 'nucleoid'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        parser = self.parser_class(**valid_args)

        assert parser.lcube > 0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == \
            number_density_cube(parser.nmon, parser.dmon, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
            volume_fraction_cube(parser.nmon, parser.dmon, parser.lcube,
                                 pbc=True)
        assert pytest.approx(parser.rho_bulk_hns, 0.0001) == \
            number_density_cube(parser.nhns, parser.dhns, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_hns, 0.0001) == \
            volume_fraction_cube(parser.nhns, parser.dhns, parser.lcube,
                                 pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cube(parser.ncrowd, parser.dcrowd, parser.lcube,
                                 pbc=True)

    def test_invalid_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'N200kbmm2XX12ac2l25epshc1nc11937',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '200kbmm2nh12ac2l25epshc1nc11937',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'N.kbmm2nh12ac2l25epshc1nc11937',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'N200kbmm2nh12ac2l25epshc1.nc11937',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.eps_hc == 1.0
        assert parser.ncrowd == 11937


class TestHnsCyl(ParserBaseTestTemplate):
    parser_class = HnsCyl
    geometry = 'cylindrical'
    topology = 'ring'
    groups = ['nucleoid', 'all']
    test_cases = [
        (
            'segment',
            'all',
            'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring.all',
            {
                'lineage_genealogy': ['segment', 'whole', 'ensemble_long',
                                      'ensemble', 'space'],
                'physical_attributes': [
                    'dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump', 'eps_hm'
                    ],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'dcyl': 8, 'nhns': 4,
                    'dcrowd': 1.0, 'lcyl': 150, 'eps_hc': 1.0, 'ncrowd': 1152,
                    'ensemble_id': 1, 'segment_id': 1
                },
                'name': 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring',
                'project_name': 'HnsCyl'
            }
        ),
        (
            'whole',
            'all',
            'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring.all',
            {
                'lineage_genealogy': ['whole', 'ensemble_long', 'ensemble',
                                      'space'],
                'physical_attributes': [
                    'dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump', 'eps_hm'
                    ],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'dcyl': 8, 'nhns': 4,
                    'dcrowd': 1.0, 'lcyl': 150, 'eps_hc': 1.0, 'ncrowd': 1152,
                    'ensemble_id': 1
                },
                'name': 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring',
                'project_name': 'HnsCyl'
            }
        ),
        (
            'ensemble_long',
            'all',
            'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
            {
                'lineage_genealogy': ['ensemble_long', 'ensemble', 'space'],
                'physical_attributes': [
                    'dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                    'ndump', 'adump', 'eps_hm'
                    ],
                'attributes': {
                    'nmon': 200, 'bend_mm': 2.0, 'dcyl': 8, 'nhns': 4,
                    'dcrowd': 1.0, 'lcyl': 150, 'eps_hc': 1.0, 'ncrowd': 1152
                },
                'name': 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
                'project_name': 'HnsCyl'
            }
        ),
        (
            'ensemble',
            'all',
            'N200D8nh4ac1epshc1.0nc1152',
            {
                'lineage_genealogy': ['ensemble', 'space'],
                'physical_attributes': ['dmon', 'dhns', 'dt', 'ndump', 'adump',
                                        'eps_hm'],
                'attributes': {
                    'nmon': 200, 'dcyl': 8, 'nhns': 4, 'dcrowd': 1.0,
                    'eps_hc': 1.0, 'ncrowd': 1152
                },
                'name': 'N200D8nh4ac1epshc1.0nc1152',
                'project_name': 'HnsCyl'
            }
        ),
        (
            'space',
            'all',
            'N200D8nh4ac1epshc1.0',
            {
                'lineage_genealogy': ['space'],
                'physical_attributes': ['dmon', 'dhns', 'dt', 'ndump', 'adump',
                                        'eps_hm'],
                'attributes': {
                    'nmon': 200, 'dcyl': 8, 'nhns': 4, 'dcrowd': 1.0,
                    'eps_hc': 1.0
                },
                'name': 'N200D8nh4ac1epshc1.0', 'project_name': 'HnsCyl'
            }
        ),
    ]

    @pytest.fixture
    def valid_args(self):
        return {
            'artifact': 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
            'lineage': 'ensemble_long',
            'group': 'nucleoid'
        }

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_filename_parsing(self, lineage, group, artifact, expected):
        super().test_filename_parsing(lineage, group, artifact, expected)

    @pytest.mark.parametrize('lineage, group, artifact, expected', test_cases)
    def test_file_path_parsing(self, make_temp_file, lineage, group, artifact,
                               expected):
        super().test_file_path_parsing(make_temp_file, lineage, group,
                                       artifact, expected)

    def test_computed_attributes(self, valid_args):
        """
        Test computed density and volume fraction attributes.
        """
        parser = self.parser_class(**valid_args)

        assert parser.dcyl > 0
        assert parser.lcyl > 0
        assert pytest.approx(parser.rho_bulk_m, 0.0001) == \
            number_density_cylinder(parser.nmon, parser.dmon,
                                    parser.lcyl, parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_m, 0.0001) == \
            volume_fraction_cylinder(parser.nmon, parser.dmon,
                                     parser.lcyl,
                                     parser.dcyl, pbc=True)
        assert pytest.approx(parser.rho_bulk_hns, 0.0001) == \
            number_density_cylinder(parser.nhns, parser.dhns,
                                    parser.lcyl, parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_hns, 0.0001) == \
            volume_fraction_cylinder(parser.nhns, parser.dhns,
                                     parser.lcyl, parser.dcyl, pbc=True)
        assert pytest.approx(parser.rho_bulk_c, 0.0001) == \
            number_density_cylinder(parser.ncrowd, parser.dcrowd, parser.lcyl,
                                    parser.dcyl, pbc=True)
        assert pytest.approx(parser.phi_bulk_c, 0.0001) == \
            volume_fraction_cylinder(parser.ncrowd, parser.dcrowd, parser.lcyl,
                                     parser.dcyl, pbc=True)

    def test_invalid_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'N200kbmm2r4.5nh4ac1lz75XX1.0nc1152',
                valid_args['lineage'],
                valid_args['group'])

    def test_invalid_heading_keyword(self, valid_args):
        """
        Ensure unknown keyword patterns in artifact name trigger a warning.
        """
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                '200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
                valid_args['lineage'],
                valid_args['group'])

    def test_single_dot_in_keyword(self, valid_args):
        with pytest.raises(AttributeError, match='has no attribute'):
            self.parser_class(
                'N.kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
                valid_args['lineage'],
                valid_args['group'])

    def test_trailing_dot_in_keyword(self, valid_args):
        parser = self.parser_class(
            'N200kbmm2r4.5nh4ac1lz75epshc1.nc1152',
            valid_args['lineage'],
            valid_args['group'])
        assert parser.eps_hc == 1.0
        assert parser.ncrowd == 1152
