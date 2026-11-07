"""Tests for the *HnsCyl* parser, :mod:`hns_cyl.parser`."""
import os
import re
from collections.abc import Callable
from operator import attrgetter
from pathlib import Path
from typing import Any, NamedTuple, cast

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from hns_cyl.parser import HnsCyl
from polyphys.manage.parser import ParserBase
from polyphys.manage.types import LineageT
from polyphys.tests.strategies import lineage_names, stored_values
from polyphys.manage.utils import (
    number_density_cylinder,
    volume_fraction_cylinder,
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
            HnsCyl, 'segment', 'all',
            artifact='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring.all',
            name='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring',
            values=dict(nmon=200, bend_mm=2.0, dcyl=8, nhns=4, dcrowd=1.0,
                        lcyl=150, eps_hc=1.0, ncrowd=1152, ensemble_id=1,
                        segment_id=1),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-segment'),
    pytest.param(
        ParseCase(
            HnsCyl, 'whole', 'all',
            artifact='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring.all',
            name='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring',
            values=dict(nmon=200, bend_mm=2.0, dcyl=8, nhns=4, dcrowd=1.0,
                        lcyl=150, eps_hc=1.0, ncrowd=1152, ensemble_id=1),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-whole'),
    pytest.param(
        ParseCase(
            HnsCyl, 'ensemble_long', 'all',
            artifact='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
            name='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
            values=dict(nmon=200, bend_mm=2.0, dcyl=8, nhns=4, dcrowd=1.0,
                        lcyl=150, eps_hc=1.0, ncrowd=1152),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-ensemble_long'),
    pytest.param(
        ParseCase(
            HnsCyl, 'ensemble', 'all',
            artifact='N200D8nh4ac1epshc1.0nc1152',
            name='N200D8nh4ac1epshc1.0nc1152',
            values=dict(nmon=200, dcyl=8, nhns=4, dcrowd=1.0, eps_hc=1.0,
                        ncrowd=1152),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-ensemble'),
    pytest.param(
        ParseCase(
            HnsCyl, 'space', 'all',
            artifact='N200D8nh4ac1epshc1.0',
            name='N200D8nh4ac1epshc1.0',
            values=dict(nmon=200, dcyl=8, nhns=4, dcrowd=1.0, eps_hc=1.0),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-space'),
])
def parse_case(request: pytest.FixtureRequest) -> ParseCase:
    """
    One hand-written name per lineage, with the values the parser must give:
    the digits in the name, converted as the class Notes state
    (``lcyl = 2 lz``, ``dcyl = 2 r - 1``).
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
        ValidCase(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
                  'ensemble_long', 'nucleoid'),
        id='HnsCyl'),
])
def valid_case(request: pytest.FixtureRequest) -> ValidCase:
    """A valid ``ensemble_long`` name of the parser class."""
    case: ValidCase = request.param
    return case


@pytest.fixture(params=[
    pytest.param((parser_class, lineage),
                 id=f"{parser_class.__name__}-{lineage}")
    for parser_class in (HnsCyl,)
    for lineage in ('segment', 'whole', 'ensemble_long', 'ensemble', 'space')
])
def declared_case(
    request: pytest.FixtureRequest
) -> tuple[type[ParserBase], Any, str, dict[str, Any]]:
    """
    A name built from the class's declared pattern for one lineage, with 10.5
    as the value of every float keyword and 2 as the value of every integer
    keyword, and the values the parser must store for them. With these
    values every volume fraction stays below 1.
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
def make_temp_file(tmp_path: Path) -> Callable[[str], str]:
    """
    Create an empty temporary file with a given filename and return its full
    path. Useful for testing parsers that expect a file path as input.
    """

    def _make_temp_file(filename: str) -> str:
        f = tmp_path / filename
        f.touch()
        return str(f)

    return _make_temp_file


# Parsing hand-written names


def test_parser_sets_values_from_name(
        parse_case: ParseCase, parsed: ParserBase) -> None:
    """Each lineage's name parses to its digits, converted as documented."""
    actual = attrgetter(*parse_case.values)(parsed)
    # Values are parsed from decimal text; equal to float64 round-off.
    assert actual == pytest.approx(tuple(parse_case.values.values()),
                                   rel=1e-12, abs=1e-12)


def test_parser_reduces_artifact_to_lineage_name(
        parse_case: ParseCase, parsed: ParserBase) -> None:
    """The name drops the group and topology suffixes after the lineage."""
    assert parsed.name == parse_case.name


def test_parser_lists_physical_attributes(
        parse_case: ParseCase, parsed: ParserBase) -> None:
    """Each lineage lists the fixed and computed attributes it sets."""
    assert parsed.physical_attributes == parse_case.physical


def test_parser_names_its_project(
        parse_case: ParseCase, parsed: ParserBase) -> None:
    """The project name is the parser class's name."""
    assert parsed.project_name == parse_case.parser_class.__name__


def test_parser_without_path_reports_no_path(
        parse_case: ParseCase, parsed: ParserBase) -> None:
    """A parser built from a name keeps it and has no filepath or ext."""
    assert (parsed.filename, parsed.filepath, parsed.ext, parsed.ispath,
            parsed.lineage, parsed.group) == (
        parse_case.artifact, None, None, False, parse_case.lineage,
        parse_case.group)


def test_parser_reads_name_from_file_path(
        parse_case: ParseCase, make_temp_file: Callable[[str], str]) -> None:
    """A parser built from a file path parses the file's stem."""
    full_path = make_temp_file(parse_case.artifact + '.txt')

    parser = parse_case.parser_class(full_path, parse_case.lineage,
                                     parse_case.group, ispath=True)

    assert (parser.filename, parser.filepath, parser.ext, parser.ispath,
            parser.name) == (
        parse_case.artifact + '.txt', Path(os.path.dirname(full_path)),
        '.txt', True, parse_case.name)


def test_parser_lists_each_attribute_once(parsed: ParserBase) -> None:
    """No attribute appears twice in `attributes`."""
    assert len(parsed.attributes) == len(set(parsed.attributes))


def test_parser_sets_every_listed_attribute(parsed: ParserBase) -> None:
    """Every attribute in `attributes` is set on the instance."""
    assert set(parsed.attributes) <= vars(parsed).keys()


@pytest.mark.parametrize(
    'parser_class, artifact, group, geometry, topology, groups', [
        pytest.param(HnsCyl, 'N200D8nh4ac1epshc1.0', 'all', 'cylindrical',
                     'ring', ['nucleoid', 'all'], id='HnsCyl'),
    ])
def test_parser_reports_geometry_topology_and_groups(
        parser_class: type[ParserBase], artifact: str, group: str,
        geometry: str, topology: str, groups: list[str]) -> None:
    """The parser reports its project's geometry, topology, and groups."""
    parser = parser_class(artifact, 'space', group)

    assert (parser.geometry, parser.topology, parser.groups) == (
        geometry, topology, groups)


@pytest.mark.parametrize('parser_class, artifact, group, dump_keywords', [
    pytest.param(HnsCyl, 'N200D8nh4ac1epshc1.0', 'all',
                 {'nucleoid': 'ndump', 'all': 'adump'}, id='HnsCyl'),
])
def test_parser_maps_each_group_to_its_dump_keyword(
        parser_class: type[ParserBase], artifact: str, group: str,
        dump_keywords: dict[str, str]) -> None:
    """Each group's dump keyword names the period of its own trajectory."""
    parser = parser_class(artifact, 'space', group)

    assert parser.dump_keywords == dump_keywords


@pytest.mark.parametrize('parser_class, artifact, group, attrs, expected', [
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.nc1152', 'nucleoid',
                 ('eps_hc', 'ncrowd'), (1.0, 1152), id='HnsCyl'),
])
def test_parser_reads_value_with_trailing_dot(
        parser_class: type[ParserBase], artifact: str, group: str,
        attrs: tuple[str, ...], expected: tuple[float, ...]) -> None:
    """A value written with a trailing dot, such as ``1.``, parses."""
    parser = parser_class(artifact, 'ensemble_long', group)

    # float('1.') is exact.
    assert attrgetter(*attrs)(parser) == pytest.approx(expected, rel=1e-12)


def test_declared_pattern_round_trips(
    declared_case: tuple[type[ParserBase], Any, str, dict[str, Any]]
) -> None:
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
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75XX1.0nc1152', 'nucleoid',
                 'epshc', 'eps_hc', id='HnsCyl-replaced'),
    pytest.param(HnsCyl, '200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'nucleoid',
                 'N', 'nmon', id='HnsCyl-heading-digit'),
])
def test_missing_keyword_raises_value_error(
        parser_class: type[ParserBase], artifact: str, group: str,
        keyword: str, attr: str) -> None:
    """A name without a pattern keyword raises ValueError naming it."""
    message = (f"keyword '{keyword}' for attribute '{attr}' not found in "
               f"'{artifact}' (lineage 'ensemble_long')")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


@pytest.mark.parametrize('parser_class, artifact, group, keyword', [
    pytest.param(HnsCyl, 'N.kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'nucleoid',
                 'N', id='HnsCyl'),
])
def test_non_numeric_value_raises_value_error(
        parser_class: type[ParserBase], artifact: str, group: str,
        keyword: str) -> None:
    """A keyword whose value is not a number raises ValueError naming it."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               "(lineage 'ensemble_long'): '.' is not a number")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


def test_empty_artifact_raises_value_error(valid_case: ValidCase) -> None:
    """An empty artifact name raises ValueError."""
    with pytest.raises(ValueError,
                       match="'artifact' cannot be an empty string."):
        valid_case.parser_class('', valid_case.lineage, valid_case.group)


def test_invalid_group_raises_value_error(valid_case: ValidCase) -> None:
    """A group the project does not have raises ValueError."""
    with pytest.raises(ValueError, match='is an invalid option.'):
        valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                'invalid')


def test_invalid_lineage_raises_value_error(valid_case: ValidCase) -> None:
    """A lineage outside the five lineages raises ValueError."""
    with pytest.raises(ValueError, match='is an invalid option.'):
        valid_case.parser_class(valid_case.artifact,
                                cast(LineageT, 'invalid'),
                                valid_case.group)


# Text representations


def test_str_summarizes_artifact(valid_case: ValidCase) -> None:
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


def test_repr_names_artifact(valid_case: ValidCase) -> None:
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
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'm',
                 'nmon', 'dmon', id='HnsCyl-m'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'hns',
                 'nhns', 'dhns', id='HnsCyl-hns'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'c',
                 'ncrowd', 'dcrowd', id='HnsCyl-c'),
])
def test_cylinder_bulk_values_match_utils(
        parser_class: type[ParserBase], artifact: str, suffix: str,
        number: str, size: str) -> None:
    """Cylindrical bulk values are the cylinder functions of the values."""
    parser = parser_class(artifact, 'ensemble_long', parser_class._groups[0])
    n, d, lcyl, dcyl = attrgetter(number, size, 'lcyl', 'dcyl')(parser)

    actual = attrgetter(f'rho_bulk_{suffix}', f'phi_bulk_{suffix}')(parser)

    # The same functions on the same inputs; float64 round-off only.
    assert actual == pytest.approx(
        (number_density_cylinder(n, d, lcyl, dcyl, pbc=True),
         volume_fraction_cylinder(n, d, lcyl, dcyl, pbc=True)), rel=1e-12)


# Lineage-name rules


@pytest.mark.parametrize('parser_class, group, artifact, ensemble', [
    pytest.param(HnsCyl, 'nucleoid', 'N200kbmm2r4.1nh4ac1lz75epshc1.0nc1152',
                 'N200D7.2nh4ac1epshc1nc1152', id='HnsCyl-r4.1'),
])
def test_short_name_prints_cylinder_size_from_r_as_decimal(
        parser_class: type[ParserBase], group: str, artifact: str,
        ensemble: str) -> None:
    """The size D = 2r - 1 converted from a decimal r prints as the decimal
    result in the ensemble name, so it matches the files named with it."""
    parser = parser_class(artifact, 'ensemble_long', group)

    assert getattr(parser, 'ensemble') == ensemble


@pytest.mark.parametrize('short_lineage', [
    pytest.param('ensemble', id='ensemble'),
    pytest.param('space', id='space'),
])
def test_short_name_reparses_to_itself(
        valid_case: ValidCase, short_lineage: LineageT) -> None:
    """A short name parsed at its own lineage rebuilds the same name."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)
    name = getattr(parser, short_lineage)

    reparsed = valid_case.parser_class(name, short_lineage, valid_case.group)

    assert getattr(reparsed, short_lineage) == name


@pytest.mark.parametrize(
    'parser_class, group, child_lineage, child_name, parent_lineage',
    [
        pytest.param(HnsCyl, 'all', 'segment',
                     'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring.all',
                     'whole', id='HnsCyl-segment-whole'),
        pytest.param(HnsCyl, 'all', 'whole',
                     'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring.all',
                     'ensemble_long', id='HnsCyl-whole-ensemble_long'),
    ],
)
def test_parser_rebuilt_parent_name_keeps_attributes(
        parser_class: type[ParserBase], group: str, child_lineage: LineageT,
        child_name: str, parent_lineage: LineageT) -> None:
    """A rebuilt parent name parses to its child's attribute values."""
    child = parser_class(child_name, child_lineage, group)
    parent = parser_class(getattr(child, parent_lineage), parent_lineage,
                          group)
    values = attrgetter(*parser_class._genealogy_attributes[parent_lineage])

    # Both sides parse the same digits, so they agree to float64 round-off.
    assert values(parent) == pytest.approx(values(child), rel=1e-12,
                                           abs=1e-12)


# Declarations


def test_attributes_at_matches_parsed_instance(
        parse_case: ParseCase, parsed: ParserBase) -> None:
    """attributes_at names exactly the public attributes a parser sets."""
    public = {name for name in vars(parsed) if not name.startswith('_')}

    names = parse_case.parser_class.attributes_at(parse_case.lineage)

    assert names == public


# Generated names


@pytest.mark.parametrize('lineage', [
    pytest.param('segment', id='segment'),
    pytest.param('whole', id='whole'),
    pytest.param('ensemble_long', id='ensemble_long'),
    pytest.param('ensemble', id='ensemble'),
    pytest.param('space', id='space'),
])
@settings(max_examples=50, deadline=None)
@given(data=st.data())
def test_generated_name_parses_to_its_values_and_back(
        lineage: LineageT, data: st.DataObject) -> None:
    """A name built from HnsCyl's own pattern parses to its values and
    rebuilds to itself: a round-trip property, with no reference value."""
    name, values = data.draw(lineage_names(HnsCyl, lineage))

    parsed = HnsCyl(name, lineage, HnsCyl._groups[0])

    assert getattr(parsed, lineage) == name
    # The parser converts in exact decimal arithmetic, the expectation in
    # float64: a few ulp apart; atol covers a value converted to zero.
    expected = stored_values(HnsCyl, values)
    assert {attr: getattr(parsed, attr) for attr in expected} == \
        pytest.approx(expected, rel=1e-12, abs=1e-12)
