"""Tests for the *HnsCub* parser, :mod:`hns_cub.parser`."""
import os
import re
from collections import OrderedDict
from collections.abc import Callable
from operator import attrgetter
from pathlib import Path
from typing import Any, NamedTuple, cast

import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from hns_cub.parser import HnsCub
from polyphys.manage.parser import ParserBase
from polyphys.manage.types import LineageT
from polyphys.tests.strategies import lineage_names, stored_values
from polyphys.manage.utils import number_density_cube, volume_fraction_cube

DATA_DIR_ENV = 'POLYPHYS_DATA_DIR'


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
            HnsCub, 'segment', 'nucleoid',
            artifact='N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
            name='N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, lcube=50,
                        eps_hc=1.0, ncrowd=11937, ensemble_id=2, segment_id=5),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCub-segment'),
    pytest.param(
        ParseCase(
            HnsCub, 'whole', 'nucleoid',
            artifact='N200kbmm2nh12ac2l25epshc1nc11937ens2.ring.nucleoid',
            name='N200kbmm2nh12ac2l25epshc1nc11937ens2.ring',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, lcube=50,
                        eps_hc=1.0, ncrowd=11937, ensemble_id=2),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCub-whole'),
    pytest.param(
        ParseCase(
            HnsCub, 'ensemble_long', 'nucleoid',
            artifact='N200kbmm2nh12ac2l25epshc1nc11937',
            name='N200kbmm2nh12ac2l25epshc1nc11937',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, lcube=50,
                        eps_hc=1.0, ncrowd=11937),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCub-ensemble_long'),
    pytest.param(
        ParseCase(
            HnsCub, 'ensemble', 'nucleoid',
            artifact='N200kbmm2nh12ac2.0epshc1nc11937',
            name='N200kbmm2nh12ac2.0epshc1nc11937',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, eps_hc=1.0,
                        ncrowd=11937),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCub-ensemble'),
    pytest.param(
        ParseCase(
            HnsCub, 'space', 'nucleoid',
            artifact='N200kbmm2nh12ac2.0epshc1',
            name='N200kbmm2nh12ac2.0epshc1',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0,
                        eps_hc=1.0),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCub-space'),
])
def parse_case(request: pytest.FixtureRequest) -> ParseCase:
    """
    One hand-written name per lineage, with the values the parser must give:
    the digits in the name, converted as the class Notes state
    (``lcube = 2 l``).
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
        ValidCase(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937',
                  'ensemble_long', 'nucleoid'),
        id='HnsCub'),
])
def valid_case(request: pytest.FixtureRequest) -> ValidCase:
    """A valid ``ensemble_long`` name of the parser class."""
    case: ValidCase = request.param
    return case


@pytest.fixture(params=[
    pytest.param((parser_class, lineage),
                 id=f"{parser_class.__name__}-{lineage}")
    for parser_class in (HnsCub,)
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


@pytest.fixture
def private_data_dir() -> Path:  # pragma: no cover - needs data
    """
    Return the directory named by ``POLYPHYS_DATA_DIR``, the maintainer's
    local full-size data; skip the test when the variable is unset or
    empty, and fail it when the variable names something that is not a
    directory.
    """
    data_dir = os.environ.get(DATA_DIR_ENV, '')
    if data_dir == '':
        pytest.skip(f"{DATA_DIR_ENV} is not set")
    if not Path(data_dir).is_dir():
        pytest.fail(f"{DATA_DIR_ENV}={data_dir!r} is not a directory",
                    pytrace=False)
    return Path(data_dir)


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


@pytest.mark.parametrize('parser_class, lineage, group, artifact, expected', [
    pytest.param(HnsCub, 'space', 'nucleoid', 'N200kbmm2nh12ac2.0epshc1',
                 ['nmon', 'bend_mm', 'nhns', 'dcrowd', 'eps_hc'],
                 id='HnsCub-space'),
])
def test_parser_lists_lineage_attributes_in_pattern_order(
        parser_class: type[ParserBase], lineage: LineageT, group: str,
        artifact: str, expected: list[str]) -> None:
    """The parsed attributes are listed in the order the name gives them."""
    parser = parser_class(artifact, lineage, group)

    assert parser.lineage_attributes == expected


def test_parser_lists_each_attribute_once(parsed: ParserBase) -> None:
    """No attribute appears twice in `attributes`."""
    assert len(parsed.attributes) == len(set(parsed.attributes))


def test_parser_sets_every_listed_attribute(parsed: ParserBase) -> None:
    """Every attribute in `attributes` is set on the instance."""
    assert set(parsed.attributes) <= vars(parsed).keys()


@pytest.mark.parametrize(
    'parser_class, artifact, group, geometry, topology, groups', [
        pytest.param(HnsCub, 'N200kbmm2nh12ac2.0epshc1', 'nucleoid', 'cubic',
                     'ring', ['nucleoid', 'all'], id='HnsCub'),
    ])
def test_parser_reports_geometry_topology_and_groups(
        parser_class: type[ParserBase], artifact: str, group: str,
        geometry: str, topology: str, groups: list[str]) -> None:
    """The parser reports its project's geometry, topology, and groups."""
    parser = parser_class(artifact, 'space', group)

    assert (parser.geometry, parser.topology, parser.groups) == (
        geometry, topology, groups)


@pytest.mark.parametrize('parser_class, artifact, group, dump_keywords', [
    pytest.param(HnsCub, 'N200kbmm2nh12ac2.0epshc1', 'nucleoid',
                 {'nucleoid': 'ndump', 'all': 'adump'}, id='HnsCub'),
])
def test_parser_maps_each_group_to_its_dump_keyword(
        parser_class: type[ParserBase], artifact: str, group: str,
        dump_keywords: dict[str, str]) -> None:
    """Each group's dump keyword names the period of its own trajectory."""
    parser = parser_class(artifact, 'space', group)

    assert parser.dump_keywords == dump_keywords


@pytest.mark.parametrize('parser_class, artifact, group, keyword, period', [
    # HnsCub names carry no dump periods; they are fixed in _defaults.
    pytest.param(HnsCub, 'N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring',
                 'nucleoid', 'ndump', 5000, id='HnsCub-nucleoid'),
    pytest.param(HnsCub, 'N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring',
                 'all', 'adump', 10000, id='HnsCub-all'),
])
def test_dump_keyword_reads_the_group_dump_period(
        parser_class: type[ParserBase], artifact: str, group: str,
        keyword: str, period: int) -> None:
    """dump_keyword names the attribute holding the group's dump period."""
    parser = parser_class(artifact, 'whole', group)

    assert (parser.dump_keyword, getattr(parser, parser.dump_keyword)) == (
        keyword, period)


@pytest.mark.parametrize('parser_class, artifact, group, attrs, expected', [
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1.nc11937', 'nucleoid',
                 ('eps_hc', 'ncrowd'), (1.0, 11937), id='HnsCub'),
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
    pytest.param(HnsCub, 'N200kbmm2XX12ac2l25epshc1nc11937', 'nucleoid',
                 'nh', 'nhns', id='HnsCub-replaced'),
    pytest.param(HnsCub, '200kbmm2nh12ac2l25epshc1nc11937', 'nucleoid',
                 'N', 'nmon', id='HnsCub-heading-digit'),
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
    pytest.param(HnsCub, 'N.kbmm2nh12ac2l25epshc1nc11937', 'nucleoid', 'N',
                 id='HnsCub'),
])
def test_non_numeric_value_raises_value_error(
        parser_class: type[ParserBase], artifact: str, group: str,
        keyword: str) -> None:
    """A keyword whose value is not a number raises ValueError naming it."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               "(lineage 'ensemble_long'): '.' is not a number")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


@pytest.mark.parametrize('parser_class, artifact, group, keyword, token', [
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc1.5', 'nucleoid', 'nc',
                 '1.5', id='HnsCub-nc'),
])
def test_non_integral_integer_raises_value_error(
        parser_class: type[ParserBase], artifact: str, group: str,
        keyword: str, token: str) -> None:
    """A non-integral value of an integer keyword raises ValueError."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               f"(lineage 'ensemble_long'): '{token}' is not an integer")

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
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937', 'm', 'nmon',
                 'dmon', id='HnsCub-m'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937', 'hns', 'nhns',
                 'dhns', id='HnsCub-hns'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937', 'c', 'ncrowd',
                 'dcrowd', id='HnsCub-c'),
])
def test_cube_bulk_values_match_utils(
        parser_class: type[ParserBase], artifact: str, suffix: str,
        number: str, size: str) -> None:
    """Cubic bulk values are the cube functions of the parsed values."""
    parser = parser_class(artifact, 'ensemble_long', parser_class._groups[0])
    n, d, lcube = attrgetter(number, size, 'lcube')(parser)

    actual = attrgetter(f'rho_bulk_{suffix}', f'phi_bulk_{suffix}')(parser)

    # The same functions on the same inputs; float64 round-off only.
    assert actual == pytest.approx(
        (number_density_cube(n, d, lcube, pbc=True),
         volume_fraction_cube(n, d, lcube, pbc=True)), rel=1e-12)


# HnsCub float keywords


def test_hnscub_keeps_fractional_kbmm_and_epshc() -> None:
    """Fractional kbmm and epshc stay floats in every lineage name."""
    parser = HnsCub('N200kbmm2.5nh8ac1.0l25.0epshc1.5nc0ens1.ring', 'whole',
                    'nucleoid')

    assert attrgetter('bend_mm', 'eps_hc', 'whole', 'ensemble_long',
                      'ensemble', 'space')(parser) == (
        2.5, 1.5, 'N200kbmm2.5nh8ac1.0l25.0epshc1.5nc0ens1',
        'N200kbmm2.5nh8ac1.0l25.0epshc1.5nc0',
        'N200kbmm2.5nh8ac1epshc1.5nc0', 'N200kbmm2.5nh8ac1epshc1.5')


@pytest.mark.parametrize('lineage', [
    pytest.param('whole', id='whole'),
    pytest.param('ensemble_long', id='ensemble_long'),
    pytest.param('ensemble', id='ensemble'),
    pytest.param('space', id='space'),
])
@pytest.mark.parametrize('epshc, other_epshc', [
    pytest.param('1.5', '1.0', id='half-step'),
    # Six significant digits, as the '%g' format prints, would merge these.
    pytest.param('1.000001', '1', id='seventh-digit'),
])
def test_hnscub_names_differing_in_epshc_stay_distinct(
        lineage: LineageT, epshc: str, other_epshc: str) -> None:
    """Names that differ only in epshc keep distinct lineage names."""
    parser = HnsCub(f'N200kbmm2.5nh8ac1.0l25.0epshc{epshc}nc0ens1.ring',
                    'whole', 'nucleoid')
    other = HnsCub(f'N200kbmm2.5nh8ac1.0l25.0epshc{other_epshc}nc0ens1.ring',
                   'whole', 'nucleoid')

    assert getattr(parser, lineage) != getattr(other, lineage)


# Lineage-name rules


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, ensemble, space',
    [
        pytest.param(HnsCub, 'nucleoid', 'whole',
                     'N200kbmm2nh12ac2l25epshc1nc0ens1.ring',
                     'N200kbmm2nh12ac2epshc1nc0', 'N200kbmm2nh12ac2epshc1',
                     id='HnsCub-integer-text'),
        pytest.param(HnsCub, 'nucleoid', 'whole',
                     'N200kbmm2.0nh12ac2.0l25.0epshc1.0nc0ens1.ring',
                     'N200kbmm2nh12ac2epshc1nc0', 'N200kbmm2nh12ac2epshc1',
                     id='HnsCub-float-text'),
        pytest.param(HnsCub, 'all', 'whole',
                     'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
                     'N200kbmm2nh12ac2epshc0nc7162', 'N200kbmm2nh12ac2epshc0',
                     id='HnsCub-zero'),
        pytest.param(HnsCub, 'all', 'whole',
                     'N200kbmm2nh12ac2l25epshc0.0nc7162ens1.ring',
                     'N200kbmm2nh12ac2epshc0nc7162', 'N200kbmm2nh12ac2epshc0',
                     id='HnsCub-zero-float-text'),
    ],
)
def test_short_name_prints_whole_number_float_as_integer(
        parser_class: type[ParserBase], group: str, lineage: LineageT,
        artifact: str, ensemble: str, space: str) -> None:
    """A whole-number float, zero included, prints as an integer."""
    parser = parser_class(artifact, lineage, group)

    assert attrgetter('ensemble', 'space')(parser) == (ensemble, space)


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, ensemble, space',
    [
        pytest.param(HnsCub, 'nucleoid', 'whole',
                     'N200kbmm2.5nh8ac1.5l25epshc0.25nc0ens1.ring',
                     'N200kbmm2.5nh8ac1.5epshc0.25nc0',
                     'N200kbmm2.5nh8ac1.5epshc0.25', id='HnsCub'),
    ],
)
def test_short_name_prints_fractional_float_as_repr(
        parser_class: type[ParserBase], group: str, lineage: LineageT,
        artifact: str, ensemble: str, space: str) -> None:
    """A fractional float prints in full, as repr, in short names."""
    parser = parser_class(artifact, lineage, group)

    assert attrgetter('ensemble', 'space')(parser) == (ensemble, space)


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


@pytest.mark.parametrize('parser_class, group, segment, whole, expected', [
    pytest.param(HnsCub, 'all',
                 'N200kbmm2nh12ac2l25epshc0nc7162ens1.j07.ring.all',
                 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring.all',
                 'N200kbmm2nh12ac2l25epshc0nc7162ens1', id='HnsCub-ring'),
])
def test_whole_name_drops_topology_suffix_from_segment_and_whole(
        parser_class: type[ParserBase], group: str, segment: str, whole: str,
        expected: str) -> None:
    """A simulation's segment and whole files give one suffix-free whole."""
    from_segment = parser_class(segment, 'segment', group)
    from_whole = parser_class(whole, 'whole', group)

    assert (getattr(from_segment, 'whole'),
            getattr(from_whole, 'whole')) == (expected, expected)


def _stamp_mismatches(path: Path) -> list[tuple[str, ...]]:
    """
    Compare the names in one SciOps ``HnsCub`` stamp file with the rebuilt
    names, under the mapping that the Usage page, "Lineage names", documents.

    A stamp file is ``<lineage_name>-<group>-stamps.csv``; its columns start
    with ``lineage_name`` and the lineage names, with a ``segment`` column
    only for a segment. The expected rebuilt name is the stamp's ``segment``
    or ``whole`` without ``.ring``, its ``ensemble_long`` unchanged, and its
    ``ensemble`` or ``space`` parsed at that lineage and rebuilt; a short
    key must also parse to the artifact's own values. Returns one
    ``(file, column, stamp, expected, rebuilt)`` tuple per mismatch.
    """
    group = path.name.removesuffix('-stamps.csv').rsplit('-', 1)[1]
    frame = pd.read_csv(path, dtype=str)
    lineage: LineageT = ('segment' if 'segment' in frame.columns
                         else 'whole')
    mismatches: list[tuple[str, ...]] = []
    for row in frame.to_dict('records'):
        parser = HnsCub(row['lineage_name'], lineage, group)
        expected: dict[str, str] = {
            column: row[column].removesuffix('.ring')
            for column in ('segment', 'whole') if column in row}
        expected['ensemble_long'] = row['ensemble_long']
        for column in ('ensemble', 'space'):
            old = HnsCub(row[column], column, group)
            expected[column] = getattr(old, column)
            values = attrgetter(*HnsCub._genealogy_attributes[column])
            if values(old) != values(parser):
                mismatches.append((path.name, f'{column} values',
                                   row[column], str(values(parser)),
                                   str(values(old))))
        mismatches.extend(
            (path.name, column, row[column], name, getattr(parser, column))
            for column, name in expected.items()
            if getattr(parser, column) != name)
    return mismatches


@pytest.fixture
def write_stamp(tmp_path: Path) -> Callable[[str, dict[str, str]], Path]:
    """
    Return a factory that writes a one-row stamp file in the SciOps format,
    ``<lineage_name>-<group>-stamps.csv``, from a mapping of column to value,
    and returns its path.
    """
    def _write(group: str, columns: dict[str, str]) -> Path:
        path = tmp_path / f"{columns['lineage_name']}-{group}-stamps.csv"
        path.write_text(','.join(columns) + '\n'
                        + ','.join(columns.values()) + '\n')
        return path
    return _write


@pytest.mark.parametrize('group, columns, flagged', [
    pytest.param(
        'nucleoid',
        {'lineage_name': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'whole': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'ensemble_long': 'N200kbmm2nh12ac2l25epshc0nc7162',
         'ensemble': 'N200kbmm2nh12ac2.0epshc0nc7162',
         'space': 'N200kbmm2nh12ac2.0epshc0', 'nmon': '200', 'n_frames': '10'},
        [], id='whole-matches'),
    pytest.param(
        'all',
        {'lineage_name': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.j07.ring',
         'segment': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.j07.ring',
         'whole': 'N200kbmm2nh12ac2l25epshc0nc7162ens1',
         'ensemble_long': 'N200kbmm2nh12ac2l25epshc0nc7162',
         'ensemble': 'N200kbmm2nh12ac2.0epshc0nc7162',
         'space': 'N200kbmm2nh12ac2.0epshc0', 'nmon': '200', 'n_frames': '10'},
        [], id='segment-matches'),
    pytest.param(
        'nucleoid',
        {'lineage_name': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'whole': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'ensemble_long': 'N200kbmm2nh12ac2l25epshc0nc7162',
         'ensemble': 'N200kbmm2nh12ac2.0epshc0nc7162',
         'space': 'N200kbmm2nh12ac2.0epshc1', 'nmon': '200', 'n_frames': '10'},
        ['space values', 'space'], id='space-value-differs'),
])
def test_stamp_mismatches_flags_only_departures_from_mapping(
        write_stamp: Callable[[str, dict[str, str]], Path], group: str,
        columns: dict[str, str], flagged: list[str]) -> None:
    """A stamp differing from the rebuilt names only by the mapping passes."""
    path = write_stamp(group, columns)

    mismatches = _stamp_mismatches(path)

    assert [mismatch[1] for mismatch in mismatches] == flagged


def test_hnscub_names_match_existing_stamps_under_mapping(  # pragma: no cover
        private_data_dir: Path) -> None:
    """Rebuilt HnsCub names equal the SciOps stamp names under the mapping."""
    stamps = sorted(
        path for path in private_data_dir.rglob('*-stamps.csv')
        if any(part.startswith('HnsCub')
               for part in path.relative_to(private_data_dir).parts[:-1]))

    mismatches = [mismatch for path in stamps
                  for mismatch in _stamp_mismatches(path)]

    assert stamps, 'no stamp file under a directory named HnsCub*'
    assert mismatches == []


@pytest.mark.parametrize(
    'parser_class, group, child_lineage, child_name, parent_lineage',
    [
        pytest.param(HnsCub, 'nucleoid', 'segment',
                     'N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
                     'whole', id='HnsCub-segment-whole'),
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


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, parent_lineage, expected',
    [
        pytest.param(HnsCub, 'nucleoid', 'segment',
                     'N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
                     'ensemble_long', 'N200kbmm2nh12ac2l25epshc1nc11937',
                     id='HnsCub-ensemble_long'),
    ],
)
def test_parser_rebuilt_parent_name_is_name_prefix(
        parser_class: type[ParserBase], group: str, lineage: LineageT,
        artifact: str, parent_lineage: LineageT, expected: str) -> None:
    """A long-form parent name is the artifact name cut after its last key."""
    parser = parser_class(artifact, lineage, group)

    assert getattr(parser, parent_lineage) == expected


# Declarations


@pytest.mark.parametrize('parser_class, lineage, expected', [
    pytest.param(HnsCub, 'space',
                 OrderedDict([('nmon', 'N'), ('bend_mm', 'kbmm'),
                              ('nhns', 'nh'), ('dcrowd', 'ac'),
                              ('eps_hc', 'epshc')]),
                 id='HnsCub-space'),
])
def test_genealogy_attributes_follow_pattern_order(
        parser_class: type[ParserBase], lineage: LineageT,
        expected: OrderedDict[str, str]) -> None:
    """Each lineage maps attributes to keywords in the pattern's order."""
    table = parser_class._genealogy_attributes[lineage]

    assert (type(table), list(table.items())) == (OrderedDict,
                                                  list(expected.items()))


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
    """A name built from HnsCub's own pattern parses to its values and
    rebuilds to itself: a round-trip property, with no reference value."""
    name, values = data.draw(lineage_names(HnsCub, lineage))

    parsed = HnsCub(name, lineage, HnsCub._groups[0])

    assert getattr(parsed, lineage) == name
    # The parser converts in exact decimal arithmetic, the expectation in
    # float64: a few ulp apart; atol covers a value converted to zero.
    expected = stored_values(HnsCub, values)
    assert {attr: getattr(parsed, attr) for attr in expected} == \
        pytest.approx(expected, rel=1e-12, abs=1e-12)
