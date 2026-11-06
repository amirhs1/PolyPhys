"""Tests for the study specs in this folder.

The expected values come from the research code's registry, SciOps
``sciops/api.py``, ``PROJECTS_DETAILS``, with the PolyPhys names
``phi_bulk_c`` for ``phi_c_bulk`` and ``rho_bulk_hns`` for
``rho_hns_bulk``.
"""
import importlib
from fnmatch import fnmatch
from types import SimpleNamespace
from typing import cast

import numpy as np
import pytest

from polyphys.analyze.pipeline import (
    AnalysisTable,
    FociPairDistances,
    LocalDensities,
)
from polyphys.manage.parser import (
    HnsCub,
    HnsCyl,
    ParserBase,
    TransFociCub,
    TransFociCyl,
)
from polyphys.manage.projects import BinderSpec, ProbeSpec, ProjectSpec


def _spec(module: str) -> ProjectSpec:
    """Return the ``SPEC`` of the study module `module`."""
    spec: ProjectSpec = importlib.import_module(module).SPEC
    return spec


@pytest.mark.parametrize('module, name', [
    pytest.param('trans_foci_cyl', 'TransFociCyl', id='TransFociCyl'),
    pytest.param('trans_foci_cub', 'TransFociCub', id='TransFociCub'),
    pytest.param('hns_cub', 'HnsCub', id='HnsCub'),
    pytest.param('hns_cyl', 'HnsCyl', id='HnsCyl'),
])
def test_study_spec_uses_its_parser(module: str, name: str) -> None:
    """Each study module's spec is built on that study's parser class."""
    assert _spec(module).name == name


@pytest.mark.parametrize('module, expected', [
    pytest.param('trans_foci_cyl', False, id='TransFociCyl'),
    pytest.param('trans_foci_cub', False, id='TransFociCub'),
    pytest.param('hns_cub', False, id='HnsCub'),
    pytest.param('hns_cyl', True, id='HnsCyl'),
])
def test_study_marks_provisional(module: str, expected: bool) -> None:
    """Only HnsCyl is provisional."""
    assert _spec(module).provisional is expected


@pytest.mark.parametrize('module', [
    pytest.param('trans_foci_cyl', id='TransFociCyl'),
    pytest.param('trans_foci_cub', id='TransFociCub'),
    pytest.param('hns_cub', id='HnsCub'),
    pytest.param('hns_cyl', id='HnsCyl'),
])
def test_study_names_only_attributes_its_parser_sets(module: str) -> None:
    """Every attribute a spec names is set at ensemble_long by its parser."""
    spec = _spec(module)
    named = {*spec.physical_attrs, *spec.species_sizes.values(), 'phi_bulk_c'}

    available = spec.parser.attributes_at('ensemble_long')

    assert named <= available


@pytest.mark.parametrize('module, group, divisor', [
    # SciOps api.py, 'group' and 'divisor' of each entry.
    pytest.param('trans_foci_cyl', 'bug', 0.025, id='TransFociCyl'),
    pytest.param('trans_foci_cub', 'bug', 0.025, id='TransFociCub'),
    pytest.param('hns_cub', 'nucleoid', 0.01, id='HnsCub'),
    pytest.param('hns_cyl', 'nucleoid', 0.01, id='HnsCyl'),
])
def test_study_group_and_divisor(module: str, group: str,
                                 divisor: float) -> None:
    """Each study has its polymer group and rounding step."""
    spec = _spec(module)

    # Literals compared with ==: the stored float is the literal itself.
    assert (spec.group, spec.divisor) == (group, divisor)


@pytest.mark.parametrize('module, geometry, topology, directions', [
    # SciOps api.py, 'geometry', 'topology', and 'directions'.
    pytest.param('trans_foci_cyl', 'cylindrical', 'ring', ('r', 'z'),
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub', 'cubic', 'ring', ('r',),
                 id='TransFociCub'),
    pytest.param('hns_cub', 'cubic', 'ring', ('r',), id='HnsCub'),
    pytest.param('hns_cyl', 'cylindrical', 'ring', ('r', 'z'), id='HnsCyl'),
])
def test_study_geometry_follows_parser(module: str, geometry: str,
                                       topology: str,
                                       directions: tuple[str, ...]) -> None:
    """Geometry, topology, and profile directions match the old registry."""
    spec = _spec(module)

    assert (spec.geometry, spec.topology, spec.directions) == (
        geometry, topology, directions)


@pytest.mark.parametrize('module, space_pattern, space_hierarchy', [
    # SciOps api.py, 'space_pat' and 'space_hierarchy'.
    pytest.param('trans_foci_cyl', 'ns*nl*al*D*ac*', 'ns*',
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub', 'ns*nl*al*ac*', 'ns*', id='TransFociCub'),
    pytest.param('hns_cub', 'N*kbmm*nh*ac*epshc*', 'N*', id='HnsCub'),
    pytest.param('hns_cyl', 'N*D*nh*ac*epshc*', 'N*', id='HnsCyl'),
])
def test_study_space_patterns(module: str, space_pattern: str,
                              space_hierarchy: str) -> None:
    """The space glob patterns derived from the parser match the old ones."""
    spec = _spec(module)

    assert (spec.space_pattern, spec.space_hierarchy) == (
        space_pattern, space_hierarchy)


@pytest.mark.parametrize('module, hierarchy', [
    # SciOps api.py, 'hierarchy'; TransFociCyl differs (see the next test).
    pytest.param('trans_foci_cub', 'al*', id='TransFociCub'),
    pytest.param('hns_cub', 'N*', id='HnsCub'),
    pytest.param('hns_cyl', 'N*', id='HnsCyl'),
])
def test_study_hierarchy(module: str, hierarchy: str) -> None:
    """The ensemble_long glob prefix derived from the parser matches."""
    assert _spec(module).hierarchy == hierarchy


def test_trans_foci_cyl_hierarchy_matches_its_names() -> None:
    """TransFociCyl's derived prefix epss* matches its ensemble_long names."""
    # A hand-written TransFociCyl ensemble_long name from test_parser.py; the
    # research code's prefix was 'eps*', which matches the same names.
    parser = TransFociCyl(
        'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
        'bdump2000adump5000', 'ensemble_long', 'bug')
    name = getattr(parser, 'ensemble_long')

    hierarchy = _spec('trans_foci_cyl').hierarchy

    assert (hierarchy, fnmatch(name, hierarchy)) == ('epss*', True)


@pytest.mark.parametrize('module, expected', [
    # SciOps api.py, 'attributes', with PolyPhys names.
    pytest.param('trans_foci_cyl',
                 ('space', 'ensemble_long', 'ensemble', 'nmon_small',
                  'nmon_large', 'dmon_large', 'dcyl', 'dcrowd', 'phi_bulk_c'),
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub',
                 ('space', 'ensemble_long', 'ensemble', 'nmon_small',
                  'nmon_large', 'dmon_large', 'dcrowd', 'phi_bulk_c'),
                 id='TransFociCub'),
    pytest.param('hns_cub',
                 ('space', 'ensemble_long', 'ensemble', 'bend_mm', 'eps_hc',
                  'nmon', 'nhns', 'dcrowd', 'phi_bulk_c', 'rho_bulk_hns'),
                 id='HnsCub'),
    pytest.param('hns_cyl',
                 ('space', 'ensemble_long', 'ensemble', 'dcyl', 'nmon',
                  'nhns', 'dcrowd', 'phi_bulk_c', 'rho_bulk_hns'),
                 id='HnsCyl'),
])
def test_study_physical_attrs(module: str, expected: tuple[str, ...]) -> None:
    """Each study adds its attributes as columns, in order."""
    assert _spec(module).physical_attrs == expected


@pytest.mark.parametrize('module, expected', [
    # SciOps api.py, 'rhosPhisNormalizedScaled', in its order.
    pytest.param('trans_foci_cyl',
                 [('Mon', 'dmon_small'), ('Crd', 'dcrowd'),
                  ('Foci', 'dmon_large')],
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub',
                 [('Mon', 'dmon_small'), ('Crd', 'dcrowd'),
                  ('Foci', 'dmon_large')],
                 id='TransFociCub'),
    pytest.param('hns_cub', [('Mon', 'dmon'), ('Crd', 'dcrowd'),
                             ('Hns', 'dhns')],
                 id='HnsCub'),
    pytest.param('hns_cyl', [('Mon', 'dmon'), ('Crd', 'dcrowd'),
                             ('Hns', 'dhns')],
                 id='HnsCyl'),
])
def test_study_species_sizes(module: str,
                             expected: list[tuple[str, str]]) -> None:
    """Each study maps its species labels to size attributes, in order."""
    assert list(_spec(module).species_sizes.items()) == expected


@pytest.mark.parametrize('module, expected', [
    # SciOps api.py, 'time_varying_props', in its order. TransFociCyl adds
    # 'transSizeTMon', which the registry omitted although the research
    # probe writes it and the analysis table reads it (maintainer,
    # 2026-11-03).
    pytest.param('trans_foci_cyl',
                 ('asphericityTMon', 'fsdTMon', 'gyrTMon', 'shapeTMon',
                  'transSizeTMon'),
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub',
                 ('asphericityTMon', 'gyrTMon', 'shapeTMon'),
                 id='TransFociCub'),
    pytest.param('hns_cub',
                 ('asphericityTMon', 'gyrTMon', 'shapeTMon',
                  'bondLengthVecMon', 'nBoundTHnsPatch', 'nFreeTHnsPatch',
                  'nEngagedTHnsPatch', 'nFreeTHnsCore', 'nBridgeTHnsCore',
                  'nDangleTHnsCore', 'nCisTHnsCore', 'nTransTHnsCore'),
                 id='HnsCub'),
    pytest.param('hns_cyl',
                 ('asphericityTMon', 'fsdTMon', 'shapeTMon', 'gyrTMon',
                  'transSizeTMon', 'bondLengthVecMon', 'nBoundTHnsPatch',
                  'nFreeTHnsPatch', 'nEngagedTHnsPatch', 'nFreeTHnsCore',
                  'nBridgeTHnsCore', 'nDangleTHnsCore', 'nCisTHnsCore',
                  'nTransTHnsCore'),
                 id='HnsCyl'),
])
def test_study_tseries_properties(module: str,
                                  expected: tuple[str, ...]) -> None:
    """Each study lists its time series, in order."""
    assert _spec(module).tseries_properties == expected


# The analysis tables


STUDIES = [
    pytest.param('trans_foci_cyl', id='TransFociCyl'),
    pytest.param('trans_foci_cub', id='TransFociCub'),
    pytest.param('hns_cub', id='HnsCub'),
    pytest.param('hns_cyl', id='HnsCyl'),
]


def _analysis(module: str) -> tuple[AnalysisTable, ...]:
    """Return the ``ANALYSIS`` tables of the study module `module`."""
    tables: tuple[AnalysisTable, ...] = importlib.import_module(
        module).ANALYSIS
    return tables


def _properties(module: str, group: str) -> set[str]:
    """Return the property names of the study's table for `group`."""
    (table,) = [t for t in _analysis(module) if t.group == group]
    return {prop.prop for prop in table.properties}


@pytest.mark.parametrize('module', STUDIES)
def test_study_has_one_table_per_group(module: str) -> None:
    """Each study has a table for its polymer group, as wholes with stamps,
    and one for 'all', as segments without stamps (SciOps ``is_segment``
    and ``has_stamp``)."""
    spec = _spec(module)
    polymer, everything = _analysis(module)

    assert (polymer.group, polymer.lineage, polymer.stamps) == (
        spec.group, 'whole', True)
    assert (everything.group, everything.lineage, everything.stamps) == (
        'all', 'segment', False)


@pytest.mark.parametrize('module', STUDIES)
def test_study_density_species_have_sizes(module: str) -> None:
    """Every species whose density profiles a table derives has a size
    attribute in the spec, which the hook reads."""
    spec = _spec(module)

    for table in _analysis(module):
        for prop in table.properties:
            if isinstance(prop.post_hook, LocalDensities):
                assert prop.species in spec.species_sizes


def _planar(species: tuple[str, ...]) -> set[str]:
    """The planar histograms and edges of the 'all' group."""
    return ({f'{plane}Hist{sp}' for plane in ('xy', 'xz', 'yz')
             for sp in species} | {'xEdge', 'yEdge', 'zEdge'})


@pytest.mark.parametrize('module, group, expected', [
    # SciOps analyzer.py, ANALYSIS_DETAILS_NUCLEOID and ANALYSIS_DETAILS_ALL,
    # the '<Study>Whole' entries, with each artifact listed once.
    pytest.param('trans_foci_cub', 'bug',
                 {'bondsHistTFoci', 'clustersHistTFoci', 'distMatTFoci',
                  'principalTMon', 'gyrTMon', 'shapeTMon',
                  'asphericityTMon'},
                 id='TransFociCub-bug'),
    pytest.param('trans_foci_cub', 'all',
                 {'rHistCrd', 'rHistMon', 'rHistFoci', 'rHistDna'}
                 | _planar(('Crd', 'Mon', 'Dna', 'Foci')),
                 id='TransFociCub-all'),
    pytest.param('trans_foci_cyl', 'bug',
                 {'bondsHistTFoci', 'clustersHistTFoci', 'distMatTFoci',
                  'principalTMon', 'fsdTMon', 'gyrTMon', 'transSizeTMon',
                  'shapeTMon', 'asphericityTMon'},
                 id='TransFociCyl-bug'),
    pytest.param('trans_foci_cyl', 'all',
                 {f'{d}Hist{sp}' for d in ('r', 'z')
                  for sp in ('Crd', 'Mon', 'Foci', 'Dna')}
                 | {f'thetaHist{sp}' for sp in ('Crd', 'Mon', 'Dna', 'Foci')}
                 | _planar(('Crd', 'Mon', 'Dna', 'Foci')),
                 id='TransFociCyl-all'),
    pytest.param('hns_cub', 'nucleoid',
                 {'principalTMon', 'bondsHistTHnsCore',
                  'bondsHistDangleTHnsCore', 'bondsHistBridgeTHnsCore',
                  'bondsHistDirDepTHnsCore', 'clustersHistTHnsCore',
                  'clustersHistDangleTHnsCore', 'clustersHistBridgeTHnsCore',
                  'clustersHistDirDepTHnsCore', 'loopLengthHistMon',
                  'bondCosineCorrVecMon', 'bondLengthVecMon',
                  'nBoundTHnsPatch', 'nFreeTHnsPatch', 'nEngagedTHnsPatch',
                  'nFreeTHnsCore', 'nBridgeTHnsCore', 'nDangleTHnsCore',
                  'nCisTHnsCore', 'nTransTHnsCore', 'gyrTMon', 'shapeTMon',
                  'asphericityTMon'},
                 id='HnsCub-nucleoid'),
    pytest.param('hns_cub', 'all',
                 {'rHistCrd', 'rHistMon', 'rHistHns'}
                 | _planar(('Crd', 'Mon', 'Hns')),
                 id='HnsCub-all'),
    pytest.param('hns_cyl', 'nucleoid',
                 {'principalTMon', 'loopLengthHistMon',
                  'bondCosineCorrVecMon', 'bondLengthVecMon',
                  'nBoundTHnsPatch', 'nFreeTHnsPatch', 'nEngagedTHnsPatch',
                  'nFreeTHnsCore', 'nBridgeTHnsCore', 'nDangleTHnsCore',
                  'nCisTHnsCore', 'nTransTHnsCore', 'fsdTMon', 'gyrTMon',
                  'transSizeTMon', 'shapeTMon', 'asphericityTMon'},
                 id='HnsCyl-nucleoid'),
    pytest.param('hns_cyl', 'all',
                 {f'{d}Hist{sp}' for d in ('r', 'z', 'theta')
                  for sp in ('Crd', 'Mon', 'Hns')}
                 | _planar(('Crd', 'Mon', 'Hns')),
                 id='HnsCyl-all'),
])
def test_study_properties_match_the_research_tables(
        module: str, group: str, expected: set[str]) -> None:
    """Each table reduces the artifacts the research code's tables named."""
    assert _properties(module, group) == expected


@pytest.mark.parametrize('module', [
    pytest.param('trans_foci_cyl', id='TransFociCyl'),
    pytest.param('trans_foci_cub', id='TransFociCub'),
])
def test_foci_pair_distance_bound_is_the_half_contour(module: str) -> None:
    """The foci pair histograms of the TransFoci studies bin at half the
    small-monomer diameter up to half the ring contour at the FENE
    maximum: 405 monomers, bonds up to ``1.5 (a_l + a_s) / 2``."""
    (table,) = [t for t in _analysis(module) if t.group == 'bug']
    (prop,) = [p for p in table.properties if p.prop == 'distMatTFoci']
    hook = prop.post_hook
    assert isinstance(hook, FociPairDistances)
    info = TransFociCub(
        'al3nl5ml27ns400ac1nc8340l50dt0.005bdump2000adump5000ens1',
        'whole', 'bug')

    assert callable(hook.n_monomers) and hook.n_monomers(info) == 405
    assert callable(hook.bin_size) and hook.bin_size(info) == 0.5
    # 202 bonds along the shorter arc, each at most 1.5 * (3 + 1) / 2 = 3.
    assert callable(hook.r_max) and hook.r_max(info) == 202 * 3.0


# The probe specs


def _probe(module: str) -> ProbeSpec:
    """Return the ``PROBE`` of the study module `module`."""
    probe: ProbeSpec = importlib.import_module(module).PROBE
    return probe


@pytest.mark.parametrize('module', STUDIES)
def test_study_probe_is_built_on_its_spec(module: str) -> None:
    """Each study's probe spec reads its parser from that study's SPEC."""
    assert _probe(module).project is _spec(module)


@pytest.mark.parametrize('module, selections', [
    # SciOps prober.py: the select_atoms calls of each study's polymer-group
    # probe, 'resid 1' for the chain and 'type 2' for the foci or the H-NS
    # patches, then those of its all-group probe, 'resid 0' for the
    # crowders, 'type 1' for the small monomers, and 'type 3' for the H-NS
    # cores.
    pytest.param('trans_foci_cyl', [('Mon', 'resid 1'), ('Foci', 'type 2'),
                                    ('Crd', 'resid 0'), ('Dna', 'type 1')],
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub', [('Mon', 'resid 1'), ('Foci', 'type 2'),
                                    ('Crd', 'resid 0'), ('Dna', 'type 1')],
                 id='TransFociCub'),
    pytest.param('hns_cub', [('Mon', 'resid 1'), ('HnsPatch', 'type 2'),
                             ('Crd', 'resid 0'), ('Hns', 'type 3')],
                 id='HnsCub'),
    pytest.param('hns_cyl', [('Mon', 'resid 1'), ('HnsPatch', 'type 2'),
                             ('Crd', 'resid 0'), ('Hns', 'type 3')],
                 id='HnsCyl'),
])
def test_study_probe_selections(module: str,
                                selections: list[tuple[str, str]]) -> None:
    """Each study selects its chain, its second species, and the species of
    its all-group histograms as the research probes did, with the chain
    labeled Mon."""
    probe = _probe(module)

    assert (list(probe.selections.items()), probe.chain) == (
        selections, 'Mon')


@pytest.mark.parametrize('module, observables', [
    # SciOps prober.py: the chain-size measures of each study's per-frame
    # loop, in the order of the study's time series, with the principal
    # axes last.
    pytest.param('trans_foci_cyl',
                 ('asphericityT', 'fsdT', 'gyrT', 'shapeT', 'transSizeT',
                  'principalT'),
                 id='TransFociCyl'),
    pytest.param('trans_foci_cub',
                 ('asphericityT', 'gyrT', 'shapeT', 'principalT'),
                 id='TransFociCub'),
    pytest.param('hns_cub', ('asphericityT', 'gyrT', 'shapeT', 'principalT'),
                 id='HnsCub'),
    pytest.param('hns_cyl',
                 ('asphericityT', 'fsdT', 'shapeT', 'gyrT', 'transSizeT',
                  'principalT'),
                 id='HnsCyl'),
])
def test_study_probe_observables(module: str,
                                 observables: tuple[str, ...]) -> None:
    """Each study measures the chain-size observables of its research
    probe."""
    assert _probe(module).observables == observables


@pytest.mark.parametrize('module, unwrap', [
    # SciOps prober.py: the unwrap argument of bug.asphericity in each
    # study's probe; the cylinders are periodic along z only.
    pytest.param('trans_foci_cyl', False, id='TransFociCyl'),
    pytest.param('trans_foci_cub', True, id='TransFociCub'),
    pytest.param('hns_cub', True, id='HnsCub'),
    pytest.param('hns_cyl', False, id='HnsCyl'),
])
def test_study_probe_unwraps_in_a_cube_only(module: str,
                                            unwrap: bool) -> None:
    """The cubic studies unwrap the chain for its asphericity; the
    cylindrical ones do not."""
    assert _probe(module).unwrap is unwrap


@pytest.mark.parametrize('module, species, keys', [
    # SciOps prober.py: the species each study's all-group probe binned, as
    # its analysis table reads them, and the directions of its geometry.
    pytest.param('trans_foci_cyl', ('Crd', 'Mon', 'Dna', 'Foci'),
                 ['r', 'z', 'theta', 'x', 'y'], id='TransFociCyl'),
    pytest.param('trans_foci_cub', ('Crd', 'Mon', 'Dna', 'Foci'),
                 ['r', 'x', 'y', 'z'], id='TransFociCub'),
    pytest.param('hns_cub', ('Crd', 'Mon', 'Hns'), ['r', 'x', 'y', 'z'],
                 id='HnsCub'),
    pytest.param('hns_cyl', ('Crd', 'Mon', 'Hns'),
                 ['r', 'z', 'theta', 'x', 'y'], id='HnsCyl'),
])
def test_study_probe_histogram_fields(module: str, species: tuple[str, ...],
                                      keys: list[str]) -> None:
    """Each study bins the species its analysis table reads, with a bin size
    along every direction its geometry needs, blocks of 100 frames, and the
    frame recentered by half the box."""
    probe = _probe(module)

    assert probe.histogram_species == species
    assert list(probe.bin_sizes) == keys
    assert (probe.block_frames, probe.recenter) == (100, 'half')


@pytest.mark.parametrize('module, info, r, cartesian, theta', [
    # SciOps prober.py, the bin_edges of each study's all-group probe: a tenth
    # of the smallest diameter along r, half of it along the Cartesian axes of
    # a cube or along z of a cylinder, pi/36 along theta; for the Hns studies
    # the radial minimum includes the H-NS core, as the maintainer settled on
    # 2026-11-04. The wholes: small monomers of 1 with crowders of 1.5 or 3
    # (the minimum is the monomer), and a crowder of 0.5 (the minimum is
    # the crowder).
    pytest.param('trans_foci_cyl', TransFociCyl(
        'epss5epsl5r5al2nl3ml8ns7ac1.5nc0lz10dt0.005bdump1adump1ens1',
        'whole', 'all'), 0.1, 0.5, np.pi / 36, id='TransFociCyl'),
    pytest.param('trans_foci_cub', TransFociCub(
        'al3nl5ml27ns400ac0.5nc8340l50dt0.005bdump2000adump5000ens1',
        'whole', 'all'), 0.05, 0.25, None, id='TransFociCub'),
    pytest.param('hns_cub', HnsCub(
        'N200kbmm2.5nh8ac0.5l25.0epshc1.5nc0ens1.ring', 'whole', 'all'),
        0.05, 0.25, None, id='HnsCub'),
    pytest.param('hns_cyl', HnsCyl(
        'N200kbmm2.0r4.5nh8ac3.0lz50epshc1.0nc0ens1', 'whole', 'all'),
        0.1, 0.5, np.pi / 36, id='HnsCyl'),
])
def test_study_probe_bin_sizes_follow_the_whole(
        module: str, info: ParserBase, r: float, cartesian: float,
        theta: float | None) -> None:
    """Each study's bin sizes resolve for a whole to a tenth of its smallest
    diameter along r and the transverse axes of a cylinder, half of it along
    the other Cartesian axes, and pi/36 along theta."""
    probe = _probe(module)

    sizes = {key: ProbeSpec.resolve(value, info)
             for key, value in probe.bin_sizes.items()}

    assert sizes['r'] == pytest.approx(r)
    if probe.geometry == 'cylindrical':
        assert sizes['x'] == sizes['y'] == pytest.approx(r)
        assert sizes['z'] == pytest.approx(cartesian)
        assert sizes['theta'] == pytest.approx(theta)
    else:
        assert sizes['x'] == sizes['y'] == sizes['z'] == pytest.approx(
            cartesian)


@pytest.mark.parametrize('module', [
    pytest.param('hns_cub', id='HnsCub'),
    pytest.param('hns_cyl', id='HnsCyl'),
])
def test_hns_radial_bin_includes_the_core_diameter(module: str) -> None:
    """An H-NS core smaller than the monomer and the crowder sets the radial
    bin size of both Hns studies, which the research code did for the cube
    only; the Cartesian bins ignore the core as the research code did."""
    # The parsers fix dhns at 1, so a stand-in with the three diameters
    # tells the formulas apart: dhns = 0.5 is the minimum for r alone.
    info = cast(ParserBase, SimpleNamespace(name='stand-in', dmon=1.0,
                                            dhns=0.5, dcrowd=2.0))
    probe = _probe(module)

    assert ProbeSpec.resolve(probe.bin_sizes['r'], info) == pytest.approx(
        0.05)
    assert ProbeSpec.resolve(probe.bin_sizes['x'], info) == pytest.approx(
        0.1 if probe.geometry == 'cylindrical' else 0.5)


@pytest.mark.parametrize('module, pair, cutoff', [
    # SciOps prober.py: the foci of the TransFoci studies, with
    # cluster_cutoff = dmon_large + dcrowd; the Hns studies measure no pairs.
    pytest.param('trans_foci_cyl', 'Foci', 3.5, id='TransFociCyl'),
    pytest.param('trans_foci_cub', 'Foci', 3.5, id='TransFociCub'),
    pytest.param('hns_cub', None, None, id='HnsCub'),
    pytest.param('hns_cyl', None, None, id='HnsCyl'),
])
def test_study_probe_pair_fields(module: str, pair: str | None,
                                 cutoff: float | None) -> None:
    """The TransFoci studies probe the foci pairs with the research cutoff,
    the large-monomer plus the crowder diameter, 2 + 1.5 here, and the
    minimum image; the Hns studies probe no pairs."""
    probe = _probe(module)
    info = TransFociCyl(
        'epss5epsl5r5al2nl3ml8ns7ac1.5nc0lz10dt0.005bdump1adump1ens1',
        'whole', 'bug')

    assert probe.pair_species == pair
    assert probe.minimum_image is True
    if cutoff is None:
        assert probe.contact_cutoff is None
    else:
        assert probe.contact_cutoff is not None
        assert ProbeSpec.resolve(probe.contact_cutoff, info) == pytest.approx(
            cutoff)


@pytest.mark.parametrize('module, binders', [
    pytest.param('trans_foci_cyl', None, id='TransFociCyl'),
    pytest.param('trans_foci_cub', None, id='TransFociCub'),
    # SciOps prober.py, the H-NS binding blocks: cis_threshold = 4 and the
    # distance matrices saved per frame. The cutoff and the inclusive
    # comparison are the thesis's binding criterion, r <= r_cut,mp = 0.685
    # (Section 5.3, Table 5.1), chosen by the maintainer on 2026-11-04 over
    # the research probe's 0.661 with a strict comparison.
    pytest.param('hns_cub',
                 BinderSpec('HnsPatch', 'HnsCore', cutoff=0.685,
                            inclusive=True, cis_threshold=4,
                            patches_per_binder=2, store_distances=True),
                 id='HnsCub'),
    pytest.param('hns_cyl',
                 BinderSpec('HnsPatch', 'HnsCore', cutoff=0.685,
                            inclusive=True, cis_threshold=4,
                            patches_per_binder=2, store_distances=True),
                 id='HnsCyl'),
])
def test_study_probe_binders(module: str,
                             binders: BinderSpec | None) -> None:
    """The H-NS studies count their two-patch binders at the published
    cutoff with the research cis threshold and keep the distance matrices;
    the foci studies have no binders."""
    assert _probe(module).binders == binders
