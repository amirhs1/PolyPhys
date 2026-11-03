"""Tests for :mod:`polyphys.manage.projects`."""
import dataclasses
from typing import Any, Callable

import pandas as pd
import pytest

from polyphys.manage.organizer import space_tseries
from polyphys.manage.parser import (
    KeywordSpec,
    ParserBase,
    Species,
    SumRuleCyl,
)
from polyphys.manage.projects import ProjectSpec


@pytest.fixture
def make_spec() -> Callable[..., ProjectSpec]:
    """Factory: a valid SumRuleCyl spec, with any field overridden."""

    def _make(**overrides: Any) -> ProjectSpec:
        fields: dict[str, Any] = {
            'parser': SumRuleCyl,
            'group': 'bug',
            'divisor': 0.025,
            'physical_attrs': ('space', 'nmon', 'dcyl'),
            'species_sizes': {'Mon': 'dmon', 'Crd': 'dcrowd'},
            'tseries_properties': ('gyrTMon',),
        }
        fields.update(overrides)
        return ProjectSpec(**fields)

    return _make


@pytest.fixture
def chain_cub() -> type[ParserBase]:
    """A parser for a chain in a cube, declared outside PolyPhys."""
    return type('ChainCub', (ParserBase,), {
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
    })


# ProjectSpec


def test_spec_is_frozen(make_spec):
    """A spec's fields cannot be reassigned."""
    spec = make_spec()

    with pytest.raises(dataclasses.FrozenInstanceError):
        spec.divisor = 0.5  # type: ignore[misc]


def test_spec_copies_its_containers(make_spec):
    """Changing the caller's list or dict later leaves a spec unchanged."""
    attrs = ['space', 'nmon']
    sizes = {'Mon': 'dmon'}
    spec = make_spec(physical_attrs=attrs, species_sizes=sizes)

    attrs.append('dcyl')
    sizes['Crd'] = 'dcrowd'

    assert (spec.physical_attrs, dict(spec.species_sizes)) == (
        ('space', 'nmon'), {'Mon': 'dmon'})


def test_spec_species_sizes_are_read_only(make_spec):
    """A spec's species-to-size map cannot be changed."""
    spec = make_spec()

    with pytest.raises(TypeError):
        spec.species_sizes['Crd'] = 'dmon'  # type: ignore[index]


@pytest.mark.parametrize('overrides, match', [
    pytest.param({'physical_attrs': ('nmon', 'nhns')},
                 r"SumRuleCyl does not set \['nhns'\] at 'ensemble_long'",
                 id='unknown-physical-attribute'),
    pytest.param({'species_sizes': {'Hns': 'dhns'}},
                 r"does not set \['dhns'\]", id='unknown-size-attribute'),
    pytest.param({'group': 'nucleoid'}, "'group' must be the polymer group",
                 id='group-not-in-parser'),
    pytest.param({'group': 'all'}, "'group' must be the polymer group",
                 id='group-all'),
    pytest.param({'divisor': 0.0}, "'divisor' must be a finite positive",
                 id='zero-divisor'),
    pytest.param({'divisor': -0.025}, "'divisor' must be a finite positive",
                 id='negative-divisor'),
    pytest.param({'divisor': float('nan')},
                 "'divisor' must be a finite positive", id='nan-divisor'),
])
def test_spec_rejects_invalid_fields(make_spec, overrides, match):
    """A spec with an unknown attribute, group, or bad divisor raises."""
    with pytest.raises(ValueError, match=match):
        make_spec(**overrides)


@pytest.mark.parametrize('parser', [
    pytest.param(ParserBase, id='abstract-parser'),
    pytest.param(str, id='not-a-parser'),
])
def test_spec_rejects_non_concrete_parser(make_spec, parser):
    """A spec needs a concrete ParserBase subclass."""
    with pytest.raises(TypeError, match='concrete ParserBase subclass'):
        make_spec(parser=parser)


def test_spec_rejects_parser_without_crowder_fraction(make_spec, chain_cub):
    """A spec needs a parser that sets phi_bulk_c at ensemble_long."""
    no_crowders = type('ChainOnly', (chain_cub,), {
        '_species': (Species('m', 'nmon', 'dmon'),)})

    with pytest.raises(ValueError, match=r"does not set \['phi_bulk_c'\]"):
        make_spec(parser=no_crowders, physical_attrs=('nmon',),
                  species_sizes={'Mon': 'dmon'})


def test_spec_for_parser_declared_outside_polyphys(chain_cub):
    """A study's own parser gets a valid spec with patterns from its names."""
    spec = ProjectSpec(
        parser=chain_cub, group='bug', divisor=0.025,
        physical_attrs=('space', 'ensemble_long', 'nmon', 'dcrowd'),
        species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
        tseries_properties=('gyrTMon',))

    # From the declaration: cubic; 'space' is ('N', 'ac'); 'ensemble_long'
    # starts with 'N'.
    assert (spec.name, spec.geometry, spec.topology, spec.directions,
            spec.space_pattern, spec.space_hierarchy, spec.hierarchy) == (
        'ChainCub', 'cubic', 'linear', ('r',), 'N*ac*', 'N*', 'N*')


def test_spec_fields_drive_space_tseries(make_spec, tmp_path):
    """space_tseries rounds phi_bulk_c with the spec's divisor."""
    spec = make_spec(divisor=0.001)
    database = tmp_path / 'N200-bug-ensAvg'
    database.mkdir()
    state = 'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000'
    pd.DataFrame({'gyrTMon-mean': [1.0, 2.0]}).to_csv(
        database / f'{state}-bug-gyrTMon-ensAvg.csv', index=False)

    space = space_tseries(
        str(database) + '/', 'gyrTMon', spec.parser, spec.hierarchy,
        list(spec.physical_attrs), spec.group, divisor=spec.divisor)

    # phi_c = n pi a^3 / 6 / (pi Z (D - a)^2 / 4) = 2 n a^3 / (3 Z (D - a)^2)
    # with n = 1000, a = 2, Z = 2 lz = 50, D = 2 r - 1 = 20:
    # 16000 / 48600 = 0.32922, which rounds to 0.329 on the spec's grid of
    # 0.001 (and to 0.325 on the default grid of 0.025).
    assert space['phi_c_bulk_round'].tolist() == [0.329, 0.329]
