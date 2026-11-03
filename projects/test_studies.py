"""Tests for the study specs in this folder.

The expected values come from the research code's registry, SciOps
``sciops/api.py``, ``PROJECTS_DETAILS``, with the PolyPhys names
``phi_bulk_c`` for ``phi_c_bulk`` and ``rho_bulk_hns`` for
``rho_hns_bulk``.
"""
import importlib
from fnmatch import fnmatch

import pytest

from polyphys.manage.parser import TransFociCyl
from polyphys.manage.projects import ProjectSpec


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
    # SciOps api.py, 'time_varying_props', in its order.
    pytest.param('trans_foci_cyl',
                 ('asphericityTMon', 'fsdTMon', 'gyrTMon', 'shapeTMon'),
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
