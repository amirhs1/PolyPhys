"""Tests for :mod:`polyphys.spec`."""
import dataclasses
import importlib
import subprocess
import sys
from typing import Any, Callable

import numpy as np
import pytest

from polyphys.manage.parser import (
    KeywordSpec,
    ParserBase,
    Species,
)
from polyphys.tests.synthetic import ChainCyl
from polyphys.analyze import pipeline
from polyphys.analyze.pipeline import (
    Autocorrelation,
    FociPairDistances,
    LocalDensities,
    TimeAverage,
)
from polyphys.spec import (
    AnalysisTable,
    BinderSpec,
    PerWhole,
    ProbeSpec,
    ProjectSpec,
    PropertySpec,
    edge,
    histogram,
    histogram_2d,
    matrix_tseries,
    tseries,
)


@pytest.fixture
def make_spec() -> Callable[..., ProjectSpec]:
    """Factory: a valid ChainCyl spec, with any field overridden."""

    def _make(**overrides: Any) -> ProjectSpec:
        fields: dict[str, Any] = {
            'parser': ChainCyl,
            'group': 'bug',
            'divisor': 0.025,
            'crowder_attr': 'phi_bulk_c',
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


# The module


def test_spec_imports_without_mdanalysis():
    """polyphys.spec imports, and builds a spec, where MDAnalysis cannot be
    imported."""
    script = (
        "import sys\n"
        "sys.modules['MDAnalysis'] = None\n"
        "from polyphys.tests.synthetic import ChainCyl\n"
        "from polyphys.spec import ProbeSpec, ProjectSpec\n"
        "spec = ProjectSpec(ChainCyl, 'bug', 0.025, (), {}, ('gyrTMon',),\n"
        "                   crowder_attr='phi_bulk_c')\n"
        "ProbeSpec(spec, {'Mon': 'resid 1'}, 'Mon', ('gyrT',), False)\n"
        "print('polyphys.probe' in sys.modules)\n")

    result = subprocess.run([sys.executable, '-c', script],
                            capture_output=True, text=True, check=False)

    assert (result.returncode, result.stdout.strip()) == (0, 'False'), \
        result.stderr


def test_per_whole_is_defined_once():
    """The pipeline's PerWhole is the spec module's alias, not a copy."""
    assert pipeline.PerWhole is PerWhole


def test_manage_projects_alias_is_removed():
    """The 0.9.0 alias module polyphys.manage.projects is gone in 1.0.0."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module('polyphys.manage.projects')


@pytest.mark.parametrize('name', [
    'AnalysisTable', 'PropertySpec', 'PostHook', 'tseries', 'histogram',
    'histogram_2d', 'edge', 'matrix_tseries',
])
def test_pipeline_no_longer_serves_names_moved_to_spec(name):
    """A table record or builder that moved to polyphys.spec in 0.9.0 is no
    longer an attribute of polyphys.analyze.pipeline in 1.0.0."""
    assert not hasattr(pipeline, name)


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
                 r"ChainCyl does not set \['nhns'\] at 'ensemble_long'",
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
    pytest.param({'crowder_attr': 'phi_c'}, r"does not set \['phi_c'\]",
                 id='unknown-crowder-attribute'),
    pytest.param({'crowder_attr': ''}, "'crowder_attr' must be a non-empty",
                 id='empty-crowder-attribute'),
])
def test_spec_rejects_invalid_fields(make_spec, overrides, match):
    """A spec with an unknown attribute, group, crowder attribute, or bad
    divisor raises."""
    with pytest.raises(ValueError, match=match):
        make_spec(**overrides)


def test_spec_requires_the_crowder_attribute():
    """A spec names its crowder attribute; the library supplies none."""
    with pytest.raises(TypeError, match='crowder_attr'):
        ProjectSpec(parser=ChainCyl, group='bug', divisor=0.025,
                    physical_attrs=(), species_sizes={},
                    tseries_properties=())


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
        crowder_attr='phi_bulk_c',
        physical_attrs=('space', 'ensemble_long', 'nmon', 'dcrowd'),
        species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
        tseries_properties=('gyrTMon',))

    # From the declaration: cubic; 'space' is ('N', 'ac'); 'ensemble_long'
    # starts with 'N'.
    assert (spec.name, spec.geometry, spec.topology, spec.directions,
            spec.space_pattern, spec.space_hierarchy, spec.hierarchy) == (
        'ChainCub', 'cubic', 'linear', ('r',), 'N*ac*', 'N*', 'N*')


# ProbeSpec


@pytest.fixture
def make_probe(make_spec) -> Callable[..., ProbeSpec]:
    """Factory: a valid probe spec of the ChainCyl spec, with any field
    overridden; the project lists the chain-size time series."""

    def _make(**overrides: Any) -> ProbeSpec:
        fields: dict[str, Any] = {
            'project': make_spec(tseries_properties=(
                'fsdTMon', 'gyrTMon', 'rfloryTMon', 'asphericityTMon')),
            'selections': {'Mon': 'resid 1', 'Crd': 'resid 0'},
            'chain': 'Mon',
            'observables': ('rfloryT', 'fsdT', 'gyrT', 'asphericityT'),
            'unwrap': False,
        }
        fields.update(overrides)
        return ProbeSpec(**fields)

    return _make


def test_probe_reads_parser_and_geometry_from_its_project(make_probe):
    """A probe spec's name, parser, and geometry are its project's."""
    probe = make_probe()

    assert (probe.name, probe.parser, probe.geometry) == (
        'ChainCyl', ChainCyl, 'cylindrical')


def test_probe_is_frozen(make_probe):
    """A probe spec's fields cannot be reassigned."""
    probe = make_probe()

    with pytest.raises(dataclasses.FrozenInstanceError):
        probe.unwrap = True  # type: ignore[misc]


def test_probe_copies_its_containers(make_probe):
    """Changing the caller's dict or list later leaves a probe spec
    unchanged, and its maps are read-only."""
    selections = {'Mon': 'resid 1'}
    observables = ['gyrT']
    bin_sizes = {'r': 0.1}
    probe = make_probe(selections=selections, observables=observables,
                       bin_sizes=bin_sizes)

    selections['Crd'] = 'resid 0'
    observables.append('fsdT')
    bin_sizes['z'] = 0.5

    assert (dict(probe.selections), probe.observables,
            dict(probe.bin_sizes)) == (
        {'Mon': 'resid 1'}, ('gyrT',), {'r': 0.1})
    with pytest.raises(TypeError):
        probe.selections['Crd'] = 'resid 0'  # type: ignore[index]


def test_probe_allows_the_matrix_observable_without_a_time_series(
        make_probe):
    """principalT, a per-frame matrix, needs no entry in the project's
    tseries_properties."""
    probe = make_probe(observables=('gyrT', 'principalT'))

    assert probe.observables == ('gyrT', 'principalT')


@pytest.mark.parametrize('overrides, match', [
    pytest.param({'selections': {}}, "'selections' must map",
                 id='no-selections'),
    pytest.param({'selections': {'Mon': ''}}, "'selections' must map",
                 id='empty-selection'),
    pytest.param({'selections': {'': 'resid 1'}}, "'selections' must map",
                 id='empty-label'),
    pytest.param({'chain': 'Foci'}, "'chain' must be a label in",
                 id='chain-not-selected'),
    pytest.param({'observables': ()}, "'observables' must not be empty",
                 id='no-observables'),
    pytest.param({'observables': ('gyrT', 'radius')},
                 r"names \['radius'\], which are not in",
                 id='unknown-observable'),
    pytest.param({'observables': ('gyrT', 'gyrT')}, "repeats a name",
                 id='repeated-observable'),
    pytest.param({'observables': ('gyrT', 'shapeT')},
                 r"names \['shapeTMon'\], which are not among",
                 id='observable-not-a-time-series'),
    pytest.param({'bin_sizes': {'phi': 0.1}}, "has the direction 'phi'",
                 id='unknown-direction'),
    pytest.param({'bin_sizes': {'r': 0.0}}, "finite positive number",
                 id='zero-bin-size'),
    pytest.param({'bin_sizes': {'r': float('inf')}},
                 "finite positive number", id='infinite-bin-size'),
    pytest.param({'bin_sizes': {'r': '0.1'}}, "finite positive number",
                 id='string-bin-size'),
    pytest.param({'bin_sizes': {'r': True}}, "finite positive number",
                 id='bool-bin-size'),
    pytest.param({'histogram_species': ('Foci',), 'bin_sizes': {'r': 0.1},
                  'block_frames': 100},
                 r"'histogram_species' names \['Foci'\]",
                 id='histogram-species-not-selected'),
    pytest.param({'histogram_species': ('Mon', 'Mon'),
                  'bin_sizes': {'r': 0.1}, 'block_frames': 100},
                 "'histogram_species' repeats a label",
                 id='repeated-histogram-species'),
    pytest.param({'histogram_species': ('Mon',), 'block_frames': 100},
                 "'histogram_species' needs 'bin_sizes'",
                 id='histogram-species-without-bin-sizes'),
    pytest.param({'histogram_species': ('Mon',), 'bin_sizes': {'r': 0.1}},
                 "'histogram_species' needs 'block_frames'",
                 id='histogram-species-without-block-frames'),
    pytest.param({'block_frames': 0}, "'block_frames' must be a positive",
                 id='zero-block-frames'),
    pytest.param({'block_frames': 2.0}, "'block_frames' must be a positive",
                 id='float-block-frames'),
    pytest.param({'recenter': 'center'}, "'recenter' must be",
                 id='unknown-recenter'),
    pytest.param({'recenter': (1.0, 2.0)}, "'recenter' must be",
                 id='short-recenter'),
    pytest.param({'recenter': (1.0, 2.0, float('nan'))}, "'recenter' must be",
                 id='nan-recenter'),
    pytest.param({'pair_species': 'Foci', 'contact_cutoff': 1.5},
                 "'pair_species' must be a label in",
                 id='pair-species-not-selected'),
    pytest.param({'pair_species': 'Mon'},
                 "'pair_species' needs 'contact_cutoff'",
                 id='pair-species-without-cutoff'),
    pytest.param({'contact_cutoff': 0.0}, "'contact_cutoff' must be",
                 id='zero-contact-cutoff'),
])
def test_probe_rejects_invalid_fields(make_probe, overrides, match):
    """A probe spec with a bad selection, chain, observable, bin size,
    histogram species, block length, recentering, pair species, or contact
    cutoff raises, naming the project."""
    with pytest.raises(ValueError, match=match):
        make_probe(**overrides)


def test_probe_rejects_a_non_project(make_probe):
    """A probe spec needs a ProjectSpec."""
    with pytest.raises(TypeError, match="'project' must be a ProjectSpec"):
        make_probe(project=ChainCyl)


def test_probe_accepts_positive_bin_sizes(make_probe):
    """Finite positive bin sizes along known directions are kept in order."""
    probe = make_probe(bin_sizes={'r': 0.1, 'z': 0.5, 'theta': 0.0873})

    assert list(probe.bin_sizes.items()) == [
        ('r', 0.1), ('z', 0.5), ('theta', 0.0873)]


def test_probe_accepts_numpy_scalar_bin_sizes(make_probe):
    """A NumPy integer or float bin size is a finite positive number."""
    probe = make_probe(bin_sizes={'r': np.float32(0.5), 'z': np.int64(2)})

    assert list(probe.bin_sizes.items()) == [('r', 0.5), ('z', 2)]


def test_probe_defaults_leave_histograms_and_pairs_unprobed(make_probe):
    """Without the histogram and pair fields, a probe spec recenters by half
    the box, uses the minimum image, and names no species for them."""
    probe = make_probe()

    assert (probe.histogram_species, probe.block_frames, probe.recenter,
            probe.pair_species, probe.contact_cutoff,
            probe.minimum_image) == ((), None, 'half', None, None, True)


def test_probe_keeps_histogram_and_pair_fields(make_probe, whole):
    """The histogram species, block length, pair species, contact cutoff,
    and minimum-image flag are stored as given, the species as a tuple, and
    the cutoff resolves for a whole."""
    def cutoff(info: ParserBase) -> float:
        return 1.5

    probe = make_probe(bin_sizes={'r': 0.1}, histogram_species=['Crd', 'Mon'],
                       block_frames=100, pair_species='Mon',
                       contact_cutoff=cutoff, minimum_image=False)

    assert (probe.histogram_species, probe.block_frames, probe.pair_species,
            probe.contact_cutoff, probe.minimum_image) == (
        ('Crd', 'Mon'), 100, 'Mon', cutoff, False)
    assert ProbeSpec.resolve(probe.contact_cutoff, whole) == 1.5


@pytest.mark.parametrize('given, stored', [
    pytest.param('none', 'none', id='none'),
    pytest.param('half', 'half', id='half'),
    pytest.param([1, 2.5, -3], (1.0, 2.5, -3.0), id='triple-as-floats'),
])
def test_probe_stores_recenter(make_probe, given, stored):
    """A recentering keyword is kept and a triple becomes a tuple of floats."""
    assert make_probe(recenter=given).recenter == stored


def test_probe_accepts_callable_bin_sizes(make_probe, whole):
    """A bin size may be a function of the whole's parser, kept as given and
    resolved for a whole: a tenth of its crowder diameter 0.5."""
    def bin_size(info: ParserBase) -> float:
        return 0.1 * float(getattr(info, 'dcrowd'))

    probe = make_probe(bin_sizes={'r': bin_size, 'z': 0.5})

    assert probe.bin_sizes['r'] is bin_size
    assert ProbeSpec.resolve(probe.bin_sizes['r'], whole) == 0.05


@pytest.fixture
def whole() -> ParserBase:
    """A ChainCyl whole with dmon 1 and dcrowd 0.5."""
    return ChainCyl(
        'N200r10.5lz50ac0.5nc100dt0.005bdump1adump1ens1',
        'whole', 'all')


@pytest.mark.parametrize('value, expected', [
    pytest.param(2, 2.0, id='number-as-float'),
    pytest.param(lambda info: info.dmon + info.dcrowd, 1.5,
                 id='callable-of-the-parser'),
    pytest.param(lambda info: np.float32(0.5), 0.5, id='numpy-float32'),
    pytest.param(lambda info: np.int64(2), 2.0, id='numpy-int64'),
])
def test_probe_resolve_returns_the_number_for_one_whole(whole, value,
                                                        expected):
    """resolve returns a number as a float and calls a callable with the
    parser."""
    assert ProbeSpec.resolve(value, whole) == expected


@pytest.mark.parametrize('bad', [
    pytest.param(lambda info: 0.0, id='zero'),
    pytest.param(lambda info: -1.0, id='negative'),
    pytest.param(lambda info: float('nan'), id='nan'),
    pytest.param(lambda info: '1', id='string'),
    pytest.param(lambda info: True, id='bool'),
    pytest.param(lambda info: np.bool_(True), id='numpy-bool'),
])
def test_probe_resolve_rejects_a_bad_callable_result(whole, bad):
    """resolve raises, naming the whole, when a callable does not return a
    finite positive number."""
    with pytest.raises(ValueError, match=f"{whole.name} must resolve to"):
        ProbeSpec.resolve(bad, whole)


# BinderSpec


@pytest.fixture
def make_binders() -> Callable[..., BinderSpec]:
    """Factory: a valid binder spec with the H-NS labels and the published
    cutoff, with any field overridden."""

    def _make(**overrides: Any) -> BinderSpec:
        fields: dict[str, Any] = {
            'patches': 'HnsPatch',
            'binders': 'HnsCore',
            'cutoff': 0.685,
        }
        fields.update(overrides)
        return BinderSpec(**fields)

    return _make


@pytest.fixture
def make_hns_probe(make_spec, make_binders) -> Callable[..., ProbeSpec]:
    """Factory: a probe spec whose project lists the eight binding series
    and whose selections include the patches, with any field overridden."""

    def _make(**overrides: Any) -> ProbeSpec:
        fields: dict[str, Any] = {
            'project': make_spec(tseries_properties=(
                'gyrTMon', 'nBoundTHnsPatch', 'nFreeTHnsPatch',
                'nEngagedTHnsPatch', 'nFreeTHnsCore', 'nBridgeTHnsCore',
                'nDangleTHnsCore', 'nCisTHnsCore', 'nTransTHnsCore')),
            'selections': {'Mon': 'resid 1', 'HnsPatch': 'type 2'},
            'chain': 'Mon',
            'observables': ('gyrT',),
            'unwrap': False,
            'binders': make_binders(),
        }
        fields.update(overrides)
        return ProbeSpec(**fields)

    return _make


def test_binders_series_take_the_patch_and_binder_labels(make_binders):
    """The three patch counts are named by the patch label and the five
    binder counts by the binder label, in the research code's order."""
    binders = make_binders(patches='Site', binders='Prot')

    assert list(binders.series.items()) == [
        ('n_patch_monomer_contacts', 'nBoundTSite'),
        ('n_patches_free', 'nFreeTSite'),
        ('n_patches_engaged', 'nEngagedTSite'),
        ('n_binders_free', 'nFreeTProt'),
        ('n_binders_bridging', 'nBridgeTProt'),
        ('n_binders_dangling', 'nDangleTProt'),
        ('n_pairs_cis', 'nCisTProt'),
        ('n_pairs_trans', 'nTransTProt')]


def test_binders_defaults(make_binders):
    """A binder spec defaults to an inclusive comparison, a cis threshold
    of 4, two patches per binder, and no stored distances."""
    binders = make_binders()

    assert (binders.inclusive, binders.cis_threshold,
            binders.patches_per_binder, binders.store_distances) == (
        True, 4, 2, False)


def test_binders_is_frozen(make_binders):
    """A binder spec's fields cannot be reassigned."""
    binders = make_binders()

    with pytest.raises(dataclasses.FrozenInstanceError):
        binders.cutoff = 0.661  # type: ignore[misc]


@pytest.mark.parametrize('overrides, match', [
    pytest.param({'patches': ''}, "'patches' must be a non-empty",
                 id='empty-patches'),
    pytest.param({'binders': ''}, "'binders' must be a non-empty",
                 id='empty-binders'),
    pytest.param({'binders': 'HnsPatch'}, "must be different labels",
                 id='same-labels'),
    pytest.param({'cutoff': 0.0}, "'cutoff' must be a finite positive",
                 id='zero-cutoff'),
    pytest.param({'cutoff': -0.685}, "'cutoff' must be a finite positive",
                 id='negative-cutoff'),
    pytest.param({'cutoff': float('inf')},
                 "'cutoff' must be a finite positive", id='infinite-cutoff'),
    pytest.param({'cis_threshold': -1}, "'cis_threshold' must be",
                 id='negative-cis'),
    pytest.param({'cis_threshold': 2.5}, "'cis_threshold' must be",
                 id='fractional-cis'),
    pytest.param({'patches_per_binder': 3}, "'patches_per_binder' must be 2",
                 id='three-patches'),
])
def test_binders_rejects_invalid_fields(make_binders, overrides, match):
    """An empty or repeated label, a non-positive or infinite cutoff, a bad
    cis threshold, or a valence other than two raises."""
    with pytest.raises(ValueError, match=match):
        make_binders(**overrides)


def test_probe_has_no_binders_by_default(make_probe):
    """A probe spec without a binder spec has binders None."""
    assert make_probe().binders is None


def test_probe_carries_its_binder_spec(make_hns_probe, make_binders):
    """A probe spec keeps the binder spec whose patches it selects and
    whose series its project lists."""
    binders = make_binders(store_distances=True)

    probe = make_hns_probe(binders=binders)

    assert probe.binders == binders


@pytest.mark.parametrize('overrides, match', [
    pytest.param({'selections': {'Mon': 'resid 1'}},
                 "patches of 'binders' must be a label in",
                 id='patches-not-selected'),
    pytest.param({'binders': BinderSpec('Mon', 'HnsCore', cutoff=0.685)},
                 "patches of 'binders' must be a label in",
                 id='patches-are-the-chain'),
])
def test_probe_rejects_binders_outside_its_selections(make_hns_probe,
                                                      overrides, match):
    """The patches of the binder spec must be a selected label other than
    the chain."""
    with pytest.raises(ValueError, match=match):
        make_hns_probe(**overrides)


def test_probe_rejects_binder_series_the_project_lacks(make_hns_probe,
                                                       make_spec):
    """Each of the eight binding series must be one of the project's
    tseries_properties."""
    project = make_spec(tseries_properties=(
        'gyrTMon', 'nBoundTHnsPatch', 'nFreeTHnsPatch', 'nEngagedTHnsPatch',
        'nFreeTHnsCore', 'nBridgeTHnsCore', 'nDangleTHnsCore'))

    match = r"names the series \['nCisTHnsCore', 'nTransTHnsCore'\]"
    with pytest.raises(ValueError, match=match):
        make_hns_probe(project=project)


def test_probe_rejects_a_non_binder_spec(make_hns_probe):
    """binders must be a BinderSpec or None."""
    with pytest.raises(TypeError, match="'binders' must be a BinderSpec"):
        make_hns_probe(binders={'patches': 'HnsPatch'})


# PropertySpec and AnalysisTable


@pytest.mark.parametrize("kwargs, match", [
    pytest.param(dict(name="", species="Mon", relation="tseries",
                      whole_type="vector"), "must not be empty",
                 id="empty-name"),
    pytest.param(dict(name="gyrT", species="Mon", relation="series",
                      whole_type="vector"), "relation must be",
                 id="unknown-relation"),
    pytest.param(dict(name="gyrT", species="Mon", relation="tseries",
                      whole_type="frame"), "whole_type must be",
                 id="unknown-whole-type"),
    pytest.param(dict(name="xEdge", species="", relation="bin_edge",
                      whole_type="vector"), "do not match",
                 id="edge-relation-vector-whole"),
    pytest.param(dict(name="xEdge", species="", relation="tseries",
                      whole_type="bin_edge"), "do not match",
                 id="tseries-relation-edge-whole"),
    pytest.param(dict(name="gyrT", species="Mon", relation="tseries",
                      whole_type="vector", edge_prop="gyrEdgeMon"),
                 "edge_prop is only valid", id="edges-on-tseries"),
    pytest.param(dict(name="xyHist", species="Crd", relation="histogram",
                      whole_type="matrix", edge_prop="xyEdgeCrd"),
                 "edge_prop is only valid", id="edges-on-matrix"),
])
def test_property_spec_rejects_inconsistent_kinds(kwargs, match):
    """A property whose relation, shape, and edges disagree is rejected."""
    with pytest.raises(ValueError, match=match):
        PropertySpec(**kwargs)


def test_constructors_build_the_kinds_the_studies_use():
    """Each named constructor sets the property, relation, shape, and edges
    of its kind."""
    assert (tseries("gyrT", "Mon").prop, tseries("gyrT", "Mon").relation,
            tseries("gyrT", "Mon").whole_type) == ("gyrTMon", "tseries",
                                                   "vector")
    hist = histogram("z", "Crd")
    assert (hist.prop, hist.relation, hist.edge_prop) == (
        "zHistCrd", "histogram", "zEdgeCrd")
    assert (histogram_2d("xz", "Dna").prop,
            histogram_2d("xz", "Dna").whole_type) == ("xzHistDna", "matrix")
    assert (edge("y").prop, edge("y").relation, edge("y").has_edges) == (
        "yEdge", "bin_edge", False)
    assert matrix_tseries("principalT", "Mon").whole_type == "matrix"


def test_analysis_table_rejects_bad_lineage_and_repeated_properties():
    """A table needs a known lineage and lists each artifact once."""
    with pytest.raises(ValueError, match="lineage must be"):
        AnalysisTable("bug", "ensemble", False, ())
    with pytest.raises(ValueError, match=r"more than once: \['gyrTMon'\]"):
        AnalysisTable("bug", "whole", False,
                      (tseries("gyrT", "Mon"), tseries("gyrT", "Mon")))


@pytest.mark.parametrize("properties, repeated", [
    pytest.param(
        (tseries("gyrT", "Mon", post_hook=TimeAverage("size")),
         tseries("fsdT", "Mon", post_hook=TimeAverage("size"))),
        "sizeMon", id="two-hooks-one-output"),
    pytest.param(
        (tseries("gyrT", "Mon", post_hook=TimeAverage("fsdT")),
         tseries("fsdT", "Mon")),
        "fsdTMon", id="hook-output-is-a-listed-property"),
    pytest.param(
        (histogram("r", "Crd", post_hook=LocalDensities()),
         tseries("rRho", "Crd")),
        "rRhoCrd", id="density-output-is-a-listed-property"),
    pytest.param(
        (matrix_tseries("distMatT", "Foci",
                        post_hook=FociPairDistances(1.0, 3.0, 10)),
         matrix_tseries("pairDistT", "Foci")),
        "pairDistTFoci", id="foci-output-is-a-listed-property"),
])
def test_analysis_table_rejects_repeated_output_names(properties, repeated):
    """A table whose declared outputs name one artifact twice is rejected,
    so no output overwrites another at run time."""
    match = rf"Outputs named more than once: \['{repeated}'\]"
    with pytest.raises(ValueError, match=match):
        AnalysisTable("bug", "whole", False, properties)


def test_analysis_table_accepts_hooks_that_pass_their_property_through():
    """A hook that passes its own property through under its name, beside
    its derived outputs, is not a repeat."""
    table = AnalysisTable("all", "segment", False, (
        tseries("gyrT", "Mon", post_hook=Autocorrelation(nlags=2)),
        histogram("r", "Crd", post_hook=LocalDensities()),
    ))

    assert [prop.prop for prop in table.properties] == ["gyrTMon", "rHistCrd"]


def test_analysis_table_rejects_local_densities_on_a_non_histogram():
    """LocalDensities on a property that is not a histogram fails when the
    table is built, not when it runs."""
    with pytest.raises(ValueError, match="needs a histogram"):
        AnalysisTable("bug", "whole", False,
                      (tseries("gyrT", "Mon", post_hook=LocalDensities()),))
