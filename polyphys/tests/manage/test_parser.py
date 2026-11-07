"""Tests for :mod:`polyphys.manage.parser`.

The parsers here are the synthetic ones of :mod:`polyphys.tests.synthetic`;
each published study tests its own parser in ``projects/<study>/``.
"""
import importlib
import os
import re
from collections import OrderedDict
from operator import attrgetter
from pathlib import Path
from typing import Any, NamedTuple

import pytest

from polyphys.manage.parser import KeywordSpec, ParserBase, Species
from polyphys.manage.utils import (
    number_density_cube,
    number_density_cylinder,
    volume_fraction_cube,
    volume_fraction_cylinder,
)
from polyphys.tests.synthetic import (
    ChainCyl,
    RingCub,
    RingCubBinders,
    RingCubTwoCrowders,
)


class ParseCase(NamedTuple):
    """A hand-written artifact name and what its parser must give."""
    parser_class: type[ParserBase]
    lineage: Any
    group: Any
    artifact: str
    name: str
    values: dict[str, Any]
    physical: list[str]


class ValidCase(NamedTuple):
    """A valid ``ensemble_long`` artifact name for one parser class."""
    parser_class: type[ParserBase]
    artifact: str
    lineage: Any
    group: Any


@pytest.fixture(params=[
    pytest.param(
        ParseCase(
            ChainCyl, 'segment', 'bug',
            artifact='N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000'
                     'adump10000ens2.j02.bug',
            name='N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000'
                 'ens2.j02',
            values=dict(nmon=2000, dcyl=20, lcyl=600.0, dcrowd=1.0,
                        ncrowd=117000, dt=0.005, bdump=5000, adump=10000,
                        ensemble_id=2, segment_id=2),
            physical=['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c']),
        id='ChainCyl-segment'),
    pytest.param(
        ParseCase(
            ChainCyl, 'whole', 'bug',
            artifact='N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000'
                     'adump10000ens2.bug',
            name='N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000'
                 'ens2',
            values=dict(nmon=2000, dcyl=20, lcyl=600.0, dcrowd=1.0,
                        ncrowd=117000, dt=0.005, bdump=5000, adump=10000,
                        ensemble_id=2),
            physical=['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c']),
        id='ChainCyl-whole'),
    pytest.param(
        ParseCase(
            ChainCyl, 'ensemble_long', 'bug',
            artifact='N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000'
                     'adump10000',
            name='N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
            values=dict(nmon=2000, dcyl=20, lcyl=600.0, dcrowd=1.0,
                        ncrowd=117000, dt=0.005, bdump=5000, adump=10000),
            physical=['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c']),
        id='ChainCyl-ensemble_long'),
    pytest.param(
        ParseCase(
            ChainCyl, 'ensemble', 'bug',
            artifact='N2000D20.0ac1.0nc117000',
            name='N2000D20.0ac1.0nc117000',
            values=dict(dmon=1.0, nmon=2000, dcrowd=1.0, ncrowd=117000,
                        dcyl=20),
            physical=['dmon']),
        id='ChainCyl-ensemble'),
    pytest.param(
        ParseCase(
            ChainCyl, 'space', 'bug',
            artifact='N2000D20.0ac1.0',
            name='N2000D20.0ac1.0',
            values=dict(dmon=1.0, nmon=2000, dcrowd=1.0, dcyl=20),
            physical=['dmon']),
        id='ChainCyl-space'),
    pytest.param(
        ParseCase(
            RingCub, 'segment', 'all',
            artifact='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens4.j10.ring.all',
            name='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4.'
                 'j10.ring',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=4,
                        segment_id=10),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='RingCub-segment'),
    pytest.param(
        ParseCase(
            RingCub, 'whole', 'all',
            artifact='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens4.ring.all',
            name='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4.'
                 'ring',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=4),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='RingCub-whole'),
    pytest.param(
        ParseCase(
            RingCub, 'ensemble_long', 'all',
            artifact='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
            name='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='RingCub-ensemble_long'),
    pytest.param(
        ParseCase(
            RingCub, 'ensemble', 'all',
            artifact='ns400nl5al5ac3nc7735',
            name='ns400nl5al5ac3nc7735',
            values=dict(dmon_large=5.0, nmon_large=5, nmon_small=400,
                        dcrowd=3.0, ncrowd=7735),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='RingCub-ensemble'),
    pytest.param(
        ParseCase(
            RingCub, 'space', 'all',
            artifact='ns400nl5al5ac3',
            name='ns400nl5al5ac3',
            values=dict(dmon_large=5.0, nmon_large=5, nmon_small=400,
                        dcrowd=3.0),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='RingCub-space'),
])
def parse_case(request: pytest.FixtureRequest) -> ParseCase:
    """
    One hand-written name per parser class and lineage, with the values its
    parser must give: the digits in the name, converted as the class
    docstring states (``lcube = 2 l``, ``lcyl = 2 lz``, ``dcyl = 2 r - 1``).
    """
    case: ParseCase = request.param
    return case


@pytest.fixture
def parsed(parse_case: ParseCase) -> ParserBase:
    """The parser built from the name of `parse_case`."""
    return parse_case.parser_class(parse_case.artifact, parse_case.lineage,
                                   parse_case.group)


@pytest.fixture(params=[
    pytest.param(
        ValidCase(ChainCyl,
                  'N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                  'ensemble_long', 'bug'),
        id='ChainCyl'),
    pytest.param(
        ValidCase(RingCub,
                  'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                  'ensemble_long', 'bug'),
        id='RingCub'),
])
def valid_case(request: pytest.FixtureRequest) -> ValidCase:
    """A valid ``ensemble_long`` name for each parser class."""
    case: ValidCase = request.param
    return case


@pytest.fixture(params=[
    pytest.param((parser_class, lineage),
                 id=f"{parser_class.__name__}-{lineage}")
    for parser_class in (ChainCyl, RingCub, RingCubTwoCrowders,
                         RingCubBinders)
    for lineage in ('segment', 'whole', 'ensemble_long', 'ensemble', 'space')
])
def declared_case(
    request: pytest.FixtureRequest
) -> tuple[type[ParserBase], Any, str, dict[str, Any]]:
    """
    A name built from a class's declared pattern for one lineage, with 10.5
    as the value of every float keyword and 2 as the value of every integer
    keyword, and the values the parser must store for them. With these
    values every volume fraction stays below 1 in both geometries.
    """
    parser_class, lineage = request.param
    specs = [(keyword, parser_class._keywords[keyword])
             for keyword in parser_class._patterns[lineage]]
    name = ''.join(f"{keyword}{'10.5' if spec.kind is float else '2'}"
                   for keyword, spec in specs)
    expected = {spec.attr: (spec.scale * 10.5 + spec.offset
                            if spec.kind is float else 2)
                for _, spec in specs}
    return parser_class, lineage, name, expected


@pytest.fixture
def chain_declaration() -> dict[str, Any]:
    """The class attributes of a valid parser for a chain in a cube."""
    return {
        '_geometry': 'cubic',
        '_topology': 'linear',
        '_groups': ['bug', 'all'],
        '_dump_keywords': {'bug': 'bdump', 'all': 'adump'},
        '_keywords': {
            'N': KeywordSpec('nmon', int),
            'ac': KeywordSpec('dcrowd', float),
            'nc': KeywordSpec('ncrowd', int),
            'l': KeywordSpec('lcube', float, scale=2.0),
            'ens': KeywordSpec('ensemble_id', int),
            'j': KeywordSpec('segment_id', int),
        },
        '_patterns': {
            'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
            'whole': ('N', 'ac', 'nc', 'l', 'ens'),
            'ensemble_long': ('N', 'ac', 'nc', 'l'),
            'ensemble': ('N', 'ac', 'nc'),
            'space': ('N', 'ac'),
        },
        '_defaults': {'dmon': 1.0, 'bdump': 1000, 'adump': 5000},
        '_species': (Species('m', 'nmon', 'dmon'),
                     Species('c', 'ncrowd', 'dcrowd')),
        '_project_attributes': {
            'segment': ['dmon'],
            'whole': ['dmon'],
            'ensemble_long': ['dmon'],
            'ensemble': ['dmon'],
            'space': ['dmon'],
        },
    }


# Parsing hand-written names


def test_parser_sets_values_from_name(parse_case, parsed):
    """Each lineage's name parses to its digits, converted as documented."""
    actual = attrgetter(*parse_case.values)(parsed)
    # Values are parsed from decimal text; equal to float64 round-off.
    assert actual == pytest.approx(tuple(parse_case.values.values()),
                                   rel=1e-12, abs=1e-12)


def test_parser_reduces_artifact_to_lineage_name(parse_case, parsed):
    """The name drops the group and topology suffixes after the lineage."""
    assert parsed.name == parse_case.name


def test_parser_lists_physical_attributes(parse_case, parsed):
    """Each lineage lists the fixed and computed attributes it sets."""
    assert parsed.physical_attributes == parse_case.physical


def test_parser_names_its_project(parse_case, parsed):
    """The project name is the parser class's name."""
    assert parsed.project_name == parse_case.parser_class.__name__


def test_parser_without_path_reports_no_path(parse_case, parsed):
    """A parser built from a name keeps it and has no filepath or ext."""
    assert (parsed.filename, parsed.filepath, parsed.ext, parsed.ispath,
            parsed.lineage, parsed.group) == (
        parse_case.artifact, None, None, False, parse_case.lineage,
        parse_case.group)


def test_from_name_reports_no_path():
    """A parser built with from_name has no filepath or ext."""
    parser = ChainCyl.from_name('N200D10.0ac2.0', 'space', 'bug')

    assert (parser.filepath, parser.ext, parser.ispath) == (None, None, False)


def test_parser_reads_name_from_file_path(parse_case, make_temp_file):
    """A parser built from a file path parses the file's stem."""
    full_path = make_temp_file(parse_case.artifact + '.txt')

    parser = parse_case.parser_class(full_path, parse_case.lineage,
                                     parse_case.group, ispath=True)

    assert (parser.filename, parser.filepath, parser.ext, parser.ispath,
            parser.name) == (
        parse_case.artifact + '.txt', Path(os.path.dirname(full_path)),
        '.txt', True, parse_case.name)


def test_parser_reports_empty_ext_for_file_without_suffix(make_temp_file):
    """A file path with no suffix gives an empty ext, not None."""
    full_path = make_temp_file('N200D10ac2')

    parser = ChainCyl(full_path, 'space', 'bug', ispath=True)

    assert (parser.ext, parser.filepath, parser.name) == (
        '', Path(os.path.dirname(full_path)), 'N200D10ac2')


@pytest.mark.parametrize('filename, name, ext, dcrowd', [
    pytest.param('N20D10ac1.5', 'N20D10ac1.5', '', 1.5, id='trailing-float'),
    pytest.param('N20D10ac1.5.csv', 'N20D10ac1.5', '.csv', 1.5,
                 id='trailing-float-with-extension'),
    pytest.param('N20D10ac1.', 'N20D10ac1.', '', 1.0, id='trailing-dot'),
])
def test_parser_strips_only_an_extension_with_a_letter(
        make_temp_file, filename, name, ext, dcrowd):
    """With ispath, a suffix made of digits, or a lone dot, stays part of the
    name instead of being dropped as an extension."""
    full_path = make_temp_file(filename)

    parser = ChainCyl(full_path, 'space', 'bug', ispath=True)

    assert (parser.name, parser.ext, parser.dcrowd) == (name, ext, dcrowd)


@pytest.mark.parametrize('parser_class, lineage, group, artifact, expected', [
    pytest.param(ChainCyl, 'ensemble', 'bug', 'N2000D20.0ac1.0nc117000',
                 ['nmon', 'dcyl', 'dcrowd', 'ncrowd'],
                 id='ChainCyl-ensemble'),
])
def test_parser_lists_lineage_attributes_in_pattern_order(
        parser_class, lineage, group, artifact, expected):
    """The parsed attributes are listed in the order the name gives them."""
    parser = parser_class(artifact, lineage, group)

    assert parser.lineage_attributes == expected


def test_name_through_unknown_keyword_raises_value_error():
    """Cutting a name at a keyword it lacks raises ValueError naming it."""
    parser = ChainCyl('N2000D20.0ac1.0', 'space', 'bug')

    with pytest.raises(ValueError,
                       match="'zz' keyword not found in 'N2000D20.0ac1.0'"):
        parser._name_through('zz')


def test_parser_lists_each_attribute_once(parsed):
    """No attribute appears twice in `attributes`."""
    assert len(parsed.attributes) == len(set(parsed.attributes))


def test_parser_sets_every_listed_attribute(parsed):
    """Every attribute in `attributes` is set on the instance."""
    assert set(parsed.attributes) <= vars(parsed).keys()


@pytest.mark.parametrize('lineage, artifact, expected', [
    pytest.param('segment',
                 'N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000'
                 'ens2.j02.bug',
                 ['segment', 'whole', 'ensemble_long', 'ensemble', 'space'],
                 id='segment'),
    pytest.param('whole',
                 'N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000'
                 'ens2.bug',
                 ['whole', 'ensemble_long', 'ensemble', 'space'],
                 id='whole'),
    pytest.param('ensemble_long',
                 'N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                 ['ensemble_long', 'ensemble', 'space'], id='ensemble_long'),
    pytest.param('ensemble', 'N2000D20.0ac1.0nc117000',
                 ['ensemble', 'space'], id='ensemble'),
    pytest.param('space', 'N2000D20.0ac1.0', ['space'], id='space'),
])
def test_lineage_genealogy_lists_lineage_and_parents(lineage, artifact,
                                                     expected):
    """A lineage's genealogy is itself followed by its parent lineages."""
    parser = ChainCyl(artifact, lineage, 'bug')

    assert parser.lineage_genealogy == expected


@pytest.mark.parametrize(
    'parser_class, artifact, group, geometry, topology, groups', [
        pytest.param(ChainCyl, 'N2000D20.0ac1.0', 'bug', 'cylindrical',
                     'linear', ['bug', 'all'], id='ChainCyl'),
        pytest.param(RingCub, 'ns400nl5al5ac3', 'all', 'cubic', 'ring',
                     ['bug', 'all'], id='RingCub'),
    ])
def test_parser_reports_geometry_topology_and_groups(
        parser_class, artifact, group, geometry, topology, groups):
    """Each parser reports its project's geometry, topology, and groups."""
    parser = parser_class(artifact, 'space', group)

    assert (parser.geometry, parser.topology, parser.groups) == (
        geometry, topology, groups)


@pytest.mark.parametrize('parser_class, artifact, group, dump_keywords', [
    pytest.param(ChainCyl, 'N2000D20.0ac1.0', 'bug',
                 {'bug': 'bdump', 'all': 'adump'}, id='ChainCyl'),
    pytest.param(RingCub, 'ns400nl5al5ac3', 'all',
                 {'bug': 'bdump', 'all': 'adump'}, id='RingCub'),
])
def test_parser_maps_each_group_to_its_dump_keyword(
        parser_class, artifact, group, dump_keywords):
    """Each group's dump keyword names the period of its own trajectory."""
    parser = parser_class(artifact, 'space', group)

    assert parser.dump_keywords == dump_keywords


@pytest.mark.parametrize('parser_class, artifact, group, keyword, period', [
    pytest.param(ChainCyl,
                 'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1',
                 'bug', 'bdump', 2000, id='ChainCyl-bug'),
    pytest.param(ChainCyl,
                 'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1',
                 'all', 'adump', 5000, id='ChainCyl-all'),
])
def test_dump_keyword_reads_the_group_dump_period(
        parser_class, artifact, group, keyword, period):
    """dump_keyword names the attribute holding the group's dump period."""
    parser = parser_class(artifact, 'whole', group)

    assert (parser.dump_keyword, getattr(parser, parser.dump_keyword)) == (
        keyword, period)


@pytest.mark.parametrize('parser_class, artifact, group, attrs, expected', [
    pytest.param(RingCub,
                 'al5nl5ml125.ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', ('mmon_large', 'nmon_small'), (125.0, 400),
                 id='RingCub'),
])
def test_parser_reads_value_with_trailing_dot(parser_class, artifact, group,
                                              attrs, expected):
    """A value written with a trailing dot, such as ``1.``, parses."""
    parser = parser_class(artifact, 'ensemble_long', group)

    # float('1.') is exact.
    assert attrgetter(*attrs)(parser) == pytest.approx(expected, rel=1e-12)


def test_declared_pattern_round_trips(declared_case):
    """
    A name built from a declared pattern parses back to the declared types
    and conversions, so float keywords keep non-integer values (property).
    """
    parser_class, lineage, name, expected = declared_case

    actual = attrgetter(*expected)(
        parser_class(name, lineage, parser_class._groups[0]))

    # One multiplication and one addition in float64.
    assert actual == pytest.approx(tuple(expected.values()), rel=1e-12)
    assert tuple(map(type, actual)) == tuple(map(type, expected.values()))


# Error paths of parsing


@pytest.mark.parametrize('parser_class, artifact, group, keyword, attr', [
    pytest.param(ChainCyl,
                 'N2000r10.5lz300.0ac1.0nc117000XX0.005bdump5000adump10000',
                 'bug', 'dt', 'dt', id='ChainCyl-replaced'),
    pytest.param(RingCub,
                 'al5nl5ml125ns400ac3XX7735l45dt0.005bdump2000adump5000',
                 'bug', 'nc', 'ncrowd', id='RingCub-replaced'),
    pytest.param(ChainCyl,
                 '2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                 'bug', 'N', 'nmon', id='ChainCyl-heading-digit'),
    pytest.param(RingCub,
                 '5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'al', 'dmon_large', id='RingCub-heading-digit'),
    pytest.param(RingCubTwoCrowders,
                 'al4nl5ns400ac1ab5nc139229XX1671l45dt0.005bdump2000'
                 'adump5000',
                 'bug', 'nb', 'ncrowd_two', id='RingCubTwoCrowders-replaced'),
])
def test_missing_keyword_raises_value_error(parser_class, artifact, group,
                                            keyword, attr):
    """A name without a pattern keyword raises ValueError naming it."""
    message = (f"keyword '{keyword}' for attribute '{attr}' not found in "
               f"'{artifact}' (lineage 'ensemble_long')")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


@pytest.mark.parametrize('parser_class, artifact, group, keyword', [
    pytest.param(ChainCyl,
                 'N.r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                 'bug', 'N', id='ChainCyl'),
    pytest.param(RingCub,
                 'al.nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'al', id='RingCub'),
])
def test_non_numeric_value_raises_value_error(parser_class, artifact, group,
                                              keyword):
    """A keyword whose value is not a number raises ValueError naming it."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               "(lineage 'ensemble_long'): '.' is not a number")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


@pytest.mark.parametrize('parser_class, artifact, group, keyword, token', [
    pytest.param(ChainCyl,
                 'N2.5r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                 'bug', 'N', '2.5', id='ChainCyl-N'),
])
def test_non_integral_integer_raises_value_error(parser_class, artifact,
                                                 group, keyword, token):
    """A non-integral value of an integer keyword raises ValueError."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               f"(lineage 'ensemble_long'): '{token}' is not an integer")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


def test_empty_artifact_raises_value_error(valid_case):
    """An empty artifact name raises ValueError."""
    with pytest.raises(ValueError,
                       match="'artifact' cannot be an empty string."):
        valid_case.parser_class('', valid_case.lineage, valid_case.group)


def test_invalid_group_raises_value_error(valid_case):
    """A group the project does not have raises ValueError."""
    with pytest.raises(ValueError, match='is an invalid option.'):
        valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                'invalid')


def test_invalid_lineage_raises_value_error(valid_case):
    """A lineage outside the five lineages raises ValueError."""
    with pytest.raises(ValueError, match='is an invalid option.'):
        valid_case.parser_class(valid_case.artifact, 'invalid',
                                valid_case.group)


# Text representations


def test_str_summarizes_artifact(valid_case):
    """`str` lists the name, geometry, group, lineage, topology, project."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)

    assert str(parser) == (
        f"Artifact:\n"
        f"    Name: '{valid_case.artifact}',\n"
        f"    Geometry: '{parser.geometry}',\n"
        f"    Group: '{valid_case.group}',\n"
        f"    Lineage: '{valid_case.lineage}',\n"
        f"    Topology: '{parser.topology}',\n"
        f"    Project: '{valid_case.parser_class.__name__}'")


def test_repr_names_artifact(valid_case):
    """`repr` names the file, geometry, group, lineage, topology, project."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)

    assert repr(parser) == (
        f"Artifact(filename='{valid_case.artifact}', "
        f"geometry='{parser.geometry}', group='{valid_case.group}', "
        f"lineage='{valid_case.lineage}', topology='{parser.topology}', "
        f"project='{valid_case.parser_class.__name__}')")


# Bulk values


@pytest.mark.parametrize('parser_class, artifact, suffix, number, size', [
    pytest.param(RingCub,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'm_small', 'nmon_small', 'dmon_small', id='RingCub-m_small'),
    pytest.param(RingCub,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'm_large', 'nmon_large', 'dmon_large', id='RingCub-m_large'),
    pytest.param(RingCub,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'c', 'ncrowd', 'dcrowd', id='RingCub-c'),
    pytest.param(RingCubTwoCrowders,
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000',
                 'm_small', 'nmon_small', 'dmon_small',
                 id='RingCubTwoCrowders-m_small'),
    pytest.param(RingCubTwoCrowders,
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000',
                 'm_large', 'nmon_large', 'dmon_large',
                 id='RingCubTwoCrowders-m_large'),
    pytest.param(RingCubTwoCrowders,
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000',
                 'c_one', 'ncrowd_one', 'dcrowd_one',
                 id='RingCubTwoCrowders-c_one'),
    pytest.param(RingCubTwoCrowders,
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000',
                 'c_two', 'ncrowd_two', 'dcrowd_two',
                 id='RingCubTwoCrowders-c_two'),
])
def test_cube_bulk_values_match_utils(parser_class, artifact, suffix, number,
                                      size):
    """Cubic bulk values are the cube functions of the parsed values."""
    parser = parser_class(artifact, 'ensemble_long', parser_class._groups[0])
    n, d, lcube = attrgetter(number, size, 'lcube')(parser)

    actual = attrgetter(f'rho_bulk_{suffix}', f'phi_bulk_{suffix}')(parser)

    # The same functions on the same inputs; float64 round-off only.
    assert actual == pytest.approx(
        (number_density_cube(n, d, lcube, pbc=True),
         volume_fraction_cube(n, d, lcube, pbc=True)), rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact, suffix, number, size', [
    pytest.param(ChainCyl,
                 'N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                 'm', 'nmon', 'dmon', id='ChainCyl-m'),
    pytest.param(ChainCyl,
                 'N2000r10.5lz300.0ac1.0nc117000dt0.005bdump5000adump10000',
                 'c', 'ncrowd', 'dcrowd', id='ChainCyl-c'),
])
def test_cylinder_bulk_values_match_utils(parser_class, artifact, suffix,
                                          number, size):
    """Cylindrical bulk values are the cylinder functions of the values."""
    parser = parser_class(artifact, 'ensemble_long', parser_class._groups[0])
    n, d, lcyl, dcyl = attrgetter(number, size, 'lcyl', 'dcyl')(parser)

    actual = attrgetter(f'rho_bulk_{suffix}', f'phi_bulk_{suffix}')(parser)

    # The same functions on the same inputs; float64 round-off only.
    assert actual == pytest.approx(
        (number_density_cylinder(n, d, lcyl, dcyl, pbc=True),
         volume_fraction_cylinder(n, d, lcyl, dcyl, pbc=True)), rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact', [
    pytest.param(RingCub,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 id='RingCub'),
    pytest.param(RingCubTwoCrowders,
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000',
                 id='RingCubTwoCrowders'),
])
def test_monomer_bulk_values_are_sums_over_sizes(parser_class, artifact):
    """All-monomer bulk values are the sums of the small and large ones."""
    parser = parser_class(artifact, 'ensemble_long', 'bug')
    small = attrgetter('rho_bulk_m_small', 'phi_bulk_m_small')(parser)
    large = attrgetter('rho_bulk_m_large', 'phi_bulk_m_large')(parser)

    actual = attrgetter('rho_bulk_m', 'phi_bulk_m')(parser)

    # One float64 addition each.
    assert actual == pytest.approx(
        (small[0] + large[0], small[1] + large[1]), rel=1e-12)


def test_crowder_bulk_values_are_sums_over_species():
    """All-crowder bulk values are the sums of the two species' values."""
    parser = RingCubTwoCrowders(
        'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000adump5000',
        'ensemble_long', 'bug')
    one = attrgetter('rho_bulk_c_one', 'phi_bulk_c_one')(parser)
    two = attrgetter('rho_bulk_c_two', 'phi_bulk_c_two')(parser)

    actual = attrgetter('rho_bulk_c', 'phi_bulk_c')(parser)

    # One float64 addition each.
    assert actual == pytest.approx((one[0] + two[0], one[1] + two[1]),
                                   rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact', [
    pytest.param(RingCub,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 id='RingCub'),
])
def test_crowder_mass_is_size_cubed(parser_class, artifact):
    """At a long lineage the crowder mass is its size cubed, m = d^3."""
    parser = parser_class(artifact, 'ensemble_long', 'bug')

    # One float64 power.
    assert parser.mcrowd == pytest.approx(parser.dcrowd ** 3, rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact', [
    pytest.param(RingCub, 'ns400nl5al5ac3nc7735', id='RingCub'),
])
def test_crowder_mass_is_placeholder_at_ensemble(parser_class, artifact):
    """At the ensemble lineage the crowder mass keeps the placeholder -1."""
    parser = parser_class(artifact, 'ensemble', 'bug')

    assert parser.mcrowd == -1


# Lineage-name rules


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, ensemble, space',
    [
        pytest.param(ChainCyl, 'all', 'segment',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1.j03',
                     'N200D20ac2nc1000', 'N200D20ac2',
                     id='ChainCyl-converted-dcyl'),
    ],
)
def test_short_name_prints_whole_number_float_as_integer(
        parser_class, group, lineage, artifact, ensemble, space):
    """A whole-number float, zero included, prints as an integer."""
    parser = parser_class(artifact, lineage, group)

    assert (parser.ensemble, parser.space) == (ensemble, space)


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, ensemble, space',
    [
        pytest.param(ChainCyl, 'all', 'segment',
                     'N200r10.75lz25.0ac1.000001nc1000dt0.005bdump2000'
                     'adump5000ens1.j03',
                     'N200D20.5ac1.000001nc1000', 'N200D20.5ac1.000001',
                     id='ChainCyl-converted-dcyl'),
    ],
)
def test_short_name_prints_fractional_float_as_repr(
        parser_class, group, lineage, artifact, ensemble, space):
    """A fractional float prints in full, as repr, in short names."""
    parser = parser_class(artifact, lineage, group)

    assert (parser.ensemble, parser.space) == (ensemble, space)


@pytest.mark.parametrize('token, expected', [
    pytest.param('4.1', 7.2, id='r4.1'),
    pytest.param('3.3', 5.6, id='r3.3'),
    pytest.param('5.5', 10.0, id='r5.5'),
])
def test_keyword_spec_converts_decimal_token_exactly(token, expected):
    """The conversion 2r - 1 of a decimal token gives the float nearest the
    exact decimal result, not one carrying binary round-off."""
    spec = KeywordSpec('dcyl', float, scale=2.0, offset=-1.0)

    value = spec.parse(token)

    assert value == expected


@pytest.mark.parametrize('parser_class, group, artifact, ensemble', [
    pytest.param(ChainCyl, 'bug',
                 'N20r4.1lz10ac1nc100dt0.005bdump1000adump5000',
                 'N20D7.2ac1nc100', id='ChainCyl-r4.1'),
    pytest.param(ChainCyl, 'bug',
                 'N20r5.5lz10ac1nc100dt0.005bdump1000adump5000',
                 'N20D10ac1nc100', id='ChainCyl-r5.5'),
])
def test_short_name_prints_cylinder_size_from_r_as_decimal(
        parser_class, group, artifact, ensemble):
    """The size D = 2r - 1 converted from a decimal r prints as the decimal
    result in the ensemble name, so it matches the files named with it."""
    parser = parser_class(artifact, 'ensemble_long', group)

    assert parser.ensemble == ensemble


@pytest.mark.parametrize('short_lineage', [
    pytest.param('ensemble', id='ensemble'),
    pytest.param('space', id='space'),
])
def test_short_name_reparses_to_itself(valid_case, short_lineage):
    """A short name parsed at its own lineage rebuilds the same name."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)
    name = getattr(parser, short_lineage)

    reparsed = valid_case.parser_class(name, short_lineage, valid_case.group)

    assert getattr(reparsed, short_lineage) == name


@pytest.mark.parametrize('artifact, lineage, group, expected', [
    pytest.param('runs/al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.j20.ring.all.lammpstrj',
                 'segment', 'all',
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.j20.ring', id='segment-trajectory'),
    pytest.param('runs/al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring.bug.data',
                 'whole', 'bug',
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring', id='whole-data'),
    pytest.param('runs/al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring.bug.restart.lammpstrj',
                 'whole', 'bug',
                 'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring', id='whole-restart'),
])
def test_parser_reduces_file_path_to_lineage_name(artifact, lineage, group,
                                                  expected):
    """A project file path reduces to its name before the group keyword."""
    parser = RingCubTwoCrowders(artifact, lineage, group, ispath=True)

    assert parser.name == expected


def test_segment_names_its_parents():
    """A segment gives every parent name, short ones without a trailing .0."""
    parser = RingCubTwoCrowders(
        'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000adump5000'
        'ens1.j20.ring.all', 'segment', 'all')

    assert (parser.segment, parser.whole, parser.ensemble_long,
            parser.ensemble, parser.space) == (
        'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000adump5000ens1.j20',
        'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000adump5000ens1',
        'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000adump5000',
        'al4nl5ns400ac1ab5nc139229nb1671',
        'al4nl5ns400ac1ab5')


@pytest.mark.parametrize(
    'parser_class, group, child_lineage, child_name, parent_lineage',
    [
        pytest.param(ChainCyl, 'all', 'segment',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1.j01', 'whole',
                     id='ChainCyl-segment-whole'),
        pytest.param(ChainCyl, 'all', 'whole',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1', 'ensemble_long',
                     id='ChainCyl-whole-ensemble_long'),
        pytest.param(RingCub, 'all', 'whole',
                     'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens4.ring.all', 'ensemble_long',
                     id='RingCub-whole-ensemble_long'),
        pytest.param(RingCubTwoCrowders, 'all', 'segment',
                     'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                     'adump5000ens1.j20.ring.all', 'whole',
                     id='RingCubTwoCrowders-segment-whole'),
        pytest.param(RingCubTwoCrowders, 'bug', 'whole',
                     'al4nl5ns400ac1ab5nc139229nb1671l45dt0.005bdump2000'
                     'adump5000ens1.ring.bug', 'ensemble_long',
                     id='RingCubTwoCrowders-whole-ensemble_long'),
    ],
)
def test_parser_rebuilt_parent_name_keeps_attributes(
        parser_class, group, child_lineage, child_name, parent_lineage):
    """A rebuilt parent name parses to its child's attribute values."""
    child = parser_class(child_name, child_lineage, group)
    parent = parser_class(getattr(child, parent_lineage), parent_lineage,
                          group)
    values = attrgetter(*parser_class._genealogy_attributes[parent_lineage])

    # Both sides parse the same digits, so they agree to float64 round-off.
    assert values(parent) == pytest.approx(values(child), rel=1e-12,
                                           abs=1e-12)


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, parent_lineage, expected',
    [
        pytest.param(ChainCyl, 'all', 'segment',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1.j01', 'segment',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1.j01', id='ChainCyl-segment'),
        pytest.param(ChainCyl, 'all', 'segment',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1.j01', 'whole',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1', id='ChainCyl-whole'),
        pytest.param(ChainCyl, 'all', 'segment',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000'
                     'ens1.j01', 'ensemble_long',
                     'N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000',
                     id='ChainCyl-ensemble_long'),
    ],
)
def test_parser_rebuilt_parent_name_is_name_prefix(
        parser_class, group, lineage, artifact, parent_lineage, expected):
    """A long-form parent name is the artifact name cut after its last key."""
    parser = parser_class(artifact, lineage, group)

    assert getattr(parser, parent_lineage) == expected


# Declarations


@pytest.mark.parametrize('parser_class, lineage, expected', [
    pytest.param(ChainCyl, 'whole',
                 OrderedDict([('nmon', 'N'), ('dcyl', 'r'), ('lcyl', 'lz'),
                              ('dcrowd', 'ac'), ('ncrowd', 'nc'),
                              ('dt', 'dt'), ('bdump', 'bdump'),
                              ('adump', 'adump'), ('ensemble_id', 'ens')]),
                 id='ChainCyl-whole'),
    pytest.param(ChainCyl, 'ensemble',
                 OrderedDict([('nmon', 'N'), ('dcyl', 'D'),
                              ('dcrowd', 'ac'), ('ncrowd', 'nc')]),
                 id='ChainCyl-ensemble'),
])
def test_genealogy_attributes_follow_pattern_order(parser_class, lineage,
                                                   expected):
    """Each lineage maps attributes to keywords in the pattern's order."""
    table = parser_class._genealogy_attributes[lineage]

    assert (type(table), list(table.items())) == (OrderedDict,
                                                  list(expected.items()))


@pytest.mark.parametrize('parser_class', [
    pytest.param(ParserBase, id='ParserBase'),
    pytest.param(type('PatternlessParser', (ParserBase,), {}),
                 id='subclass-without-patterns'),
])
def test_parser_without_patterns_cannot_be_instantiated(parser_class):
    """A parser class that declares no patterns raises TypeError."""
    with pytest.raises(TypeError, match='declares no lineage patterns'):
        parser_class('N200', 'space', 'bug')


# Class-level reads of the declaration


def test_pattern_returns_declared_keywords():
    """pattern gives a lineage's keywords in the order its name writes them."""
    keywords = ChainCyl.pattern('space')

    # The ChainCyl declaration: 'space': ('N', 'D', 'ac').
    assert keywords == ('N', 'D', 'ac')


@pytest.mark.parametrize('lineage, expected', [
    pytest.param(
        'ensemble_long',
        {'nmon', 'dcyl', 'lcyl', 'dcrowd', 'ncrowd', 'dt', 'bdump', 'adump',
         'dmon', 'segment', 'whole', 'ensemble_long', 'ensemble', 'space',
         'rho_bulk_m', 'phi_bulk_m', 'rho_bulk_c', 'phi_bulk_c'},
        id='ensemble_long-with-bulk-values'),
    pytest.param(
        'space',
        {'nmon', 'dcyl', 'dcrowd', 'dmon', 'segment', 'whole',
         'ensemble_long', 'ensemble', 'space'},
        id='space-without-bulk-values'),
])
def test_attributes_at_lists_declared_attributes(lineage, expected):
    """attributes_at gives the parsed, default, lineage, and bulk names."""
    names = ChainCyl.attributes_at(lineage)

    # From the ChainCyl declaration: pattern attributes, 'dmon' in
    # _defaults, the five lineage names, and at ensemble_long the bulk
    # values of the species 'm' and 'c'.
    assert names == frozenset(expected)


def test_attributes_at_adds_unlisted_bulk_values(chain_declaration):
    """attributes_at adds bulk values that _project_attributes omits."""
    chain_cub = type('ChainCub', (ParserBase,), chain_declaration)

    names = chain_cub.attributes_at('ensemble_long')

    # The declaration lists only 'dmon' as a project attribute; its species
    # 'm' and 'c' still get bulk values at ensemble_long.
    assert {'rho_bulk_m', 'phi_bulk_m', 'rho_bulk_c', 'phi_bulk_c'} <= names


def test_attributes_at_matches_parsed_instance(parse_case, parsed):
    """attributes_at names exactly the public attributes a parser sets."""
    public = {name for name in vars(parsed) if not name.startswith('_')}

    names = parse_case.parser_class.attributes_at(parse_case.lineage)

    assert names == public


@pytest.mark.parametrize('method', [
    pytest.param('pattern', id='pattern'),
    pytest.param('attributes_at', id='attributes_at'),
])
def test_class_reads_need_patterns(method):
    """pattern and attributes_at raise TypeError on a patternless class."""
    with pytest.raises(TypeError, match='declares no lineage patterns'):
        getattr(ParserBase, method)('space')


@pytest.mark.parametrize('method', [
    pytest.param('pattern', id='pattern'),
    pytest.param('attributes_at', id='attributes_at'),
])
def test_class_reads_reject_unknown_lineage(method):
    """pattern and attributes_at raise ValueError on an unknown lineage."""
    with pytest.raises(ValueError, match='invalid option'):
        getattr(ChainCyl, method)('galaxy')


def test_declared_parser_parses_its_names(chain_declaration):
    """A class built from a valid declaration parses names of its grammar."""
    chain_cub = type('ChainCub', (ParserBase,), chain_declaration)

    parser = chain_cub('N100ac2.0nc50l10.0ens1.bug', 'whole', 'bug')

    assert (parser.nmon, parser.dcrowd, parser.ncrowd, parser.lcube,
            parser.ensemble_long, parser.ensemble, parser.space) == (
        100, 2.0, 50, 20.0, 'N100ac2.0nc50l10.0', 'N100ac2nc50',
        'N100ac2')


@pytest.mark.parametrize('override, match', [
    pytest.param({'_geometry': None}, 'must declare', id='no-geometry'),
    pytest.param({'_dump_keywords': None}, "'_dump_keywords'",
                 id='no-dump-keywords'),
    pytest.param({'_dump_keywords': {'bug': 'bdump'}},
                 "_dump_keywords must have one entry per group",
                 id='group-no-dump-keyword'),
    pytest.param(
        {'_dump_keywords': {'bug': 'bdump', 'all': 'adump',
                            'nucleoid': 'ndump'}},
        "_dump_keywords must have one entry per group",
        id='dump-keyword-of-undeclared-group'),
    pytest.param(
        {'_dump_keywords': {'bug': 'tdump', 'all': 'adump'}},
        "the dump keyword 'tdump' of group 'bug' is not set by the "
        "'segment' pattern or '_defaults'",
        id='dump-keyword-not-set'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc', 'l'),
                       'ensemble': ('N', 'ac', 'nc')}},
        "_patterns must have one entry per lineage", id='lineage-no-pattern'),
    pytest.param(
        {'_project_attributes': {'segment': [], 'whole': []}},
        "_project_attributes must have one entry per lineage",
        id='lineage-no-project-attributes'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc', 'l'),
                       'ensemble': ('N', 'ac', 'nc'),
                       'space': ('N', 'xx')}},
        "keyword 'xx' in the 'space' pattern has no KeywordSpec",
        id='undeclared-keyword'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc', 'l'),
                       'ensemble': ('N', 'ac', 'nc'),
                       'space': ('N', 'ac', 'N')}},
        "attribute 'nmon' appears more than once in the 'space' pattern",
        id='duplicate-attribute'),
    pytest.param(
        {'_totals': {'all': ('m', 'x')}},
        "total 'all' sums 'x', which is not a declared species",
        id='total-of-undeclared-species'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc'),
                       'ensemble': ('N', 'ac', 'nc'),
                       'space': ('N', 'ac')}},
        "the bulk values need 'lcube', which the 'ensemble_long' pattern",
        id='bulk-input-not-parsed'),
])
def test_invalid_declaration_raises_type_error_at_class_creation(
        chain_declaration, override, match):
    """An incomplete or inconsistent declaration fails when it is defined."""
    with pytest.raises(TypeError, match=re.escape(match)):
        type('BadParser', (ParserBase,), {**chain_declaration, **override})


# KeywordSpec


@pytest.mark.parametrize('spec, token, expected', [
    pytest.param(KeywordSpec('nmon', int), '200', 200, id='int'),
    pytest.param(KeywordSpec('ncrowd', int), '1000.', 1000,
                 id='int-trailing-dot'),
    pytest.param(KeywordSpec('dcrowd', float), '2.5', 2.5, id='float'),
    pytest.param(KeywordSpec('lcube', float, scale=2.0), '25.0', 50.0,
                 id='half-side'),
    pytest.param(KeywordSpec('dcyl', float, scale=2.0, offset=-1.0), '10.5',
                 20.0, id='wall-radius'),
])
def test_keyword_spec_converts_token(spec, token, expected):
    """A token converts to `scale * value + offset`, typed as declared."""
    actual = spec.parse(token)

    # The conversions are exact in float64 for these tokens.
    assert (actual, type(actual)) == (pytest.approx(expected, rel=1e-12),
                                      type(expected))


@pytest.mark.parametrize('spec, token, match', [
    pytest.param(KeywordSpec('dcrowd', float), '.', "'.' is not a number",
                 id='dot'),
    pytest.param(KeywordSpec('dcrowd', float), '', "'' is not a number",
                 id='empty'),
    pytest.param(KeywordSpec('nmon', int), '1.5', "'1.5' is not an integer",
                 id='non-integral'),
])
def test_keyword_spec_rejects_bad_token(spec, token, match):
    """A token that is not a number, or not whole for int, is rejected."""
    with pytest.raises(ValueError, match=re.escape(match)):
        spec.parse(token)


@pytest.mark.parametrize('kwargs, match', [
    pytest.param({'attr': 'name', 'kind': str}, "must be int or float",
                 id='str-kind'),
    pytest.param({'attr': 'nmon', 'kind': int, 'offset': 1.0},
                 "only a float keyword", id='int-with-offset'),
])
def test_keyword_spec_rejects_bad_declaration(kwargs, match):
    """A KeywordSpec of another type, or an int with a conversion, fails."""
    with pytest.raises(TypeError, match=match):
        KeywordSpec(**kwargs)


# Study parsers that left the module


@pytest.mark.parametrize('name, message', [
    pytest.param('HnsCub', 'projects/hns_cub/parser.py', id='HnsCub'),
    pytest.param('HnsCyl', 'projects/hns_cyl/parser.py', id='HnsCyl'),
    pytest.param('TransFociCub', 'projects/trans_foci_cub/parser.py',
                 id='TransFociCub'),
    pytest.param('TransFociCyl', 'projects/trans_foci_cyl/parser.py',
                 id='TransFociCyl'),
    pytest.param('SumRuleCyl', 'removed from PolyPhys in 1.0.0',
                 id='SumRuleCyl'),
    pytest.param('SumRuleCubHeteroRing', 'removed from PolyPhys in 1.0.0',
                 id='SumRuleCubHeteroRing'),
    pytest.param('SumRuleCubHeteroLinear', 'removed from PolyPhys in 1.0.0',
                 id='SumRuleCubHeteroLinear'),
    pytest.param('SumRuleCubHeteroRingBiDisperse',
                 'removed from PolyPhys in 1.0.0',
                 id='SumRuleCubHeteroRingBiDisperse'),
    pytest.param('TwoMonDepCub', 'removed from PolyPhys in 1.0.0',
                 id='TwoMonDepCub'),
])
def test_moved_study_parser_import_raises_import_error(name, message):
    """Importing a study parser that left the module says where it went."""
    module = importlib.import_module('polyphys.manage.parser')

    with pytest.raises(ImportError, match=re.escape(message)):
        getattr(module, name)


def test_unknown_name_raises_attribute_error():
    """A name the module never defined raises AttributeError, as usual."""
    module = importlib.import_module('polyphys.manage.parser')

    with pytest.raises(AttributeError, match="has no attribute 'NoSuchName'"):
        getattr(module, 'NoSuchName')
