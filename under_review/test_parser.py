# test_parser.py

import os
import re
import pytest

from polyphys.manage.parser import ParserBase, TwoMonDepCub
from polyphys.manage.utils import invalid_keyword
from polyphys.analyze.measurer import (
    number_density_cube,
    volume_fraction_cube
)

# --------------------------------------------------------------------------------
# 1) Fixtures / Helper Strings
# --------------------------------------------------------------------------------

@pytest.fixture
def simple_space_name():
    """
    A “space”‐lineage filename that encodes:
      - nmon (nm) = 2
      - dmon (am) = 5.0
      - dcrowd (ac) = 1.0
      - ncrowd (nc) = 3
    """
    return "nm2am5.0ac1.0nc3"


@pytest.fixture
def full_segment_name():
    """
    A “segment”‐lineage filename that encodes:
      - dmon   (am) = 5.0
      - nmon   (nm) = 2
      - dcrowd (ac) = 1.0
      - ncrowd (nc) = 3
      - hl     (lcube half‐side) = 10   → lcube = 2×10 = 20
      - sd     (d_sur) = 2.0
      - dt     (dt) = 0.5
      - bdump  (bdump) = 10
      - adump  (adump) = 20
      - tdump  (tdump) = 5
      - ens    (ensemble_id) = 1
      - j      (segment_id) = 05 → becomes int(5)
    """
    return (
        "am5.0nm2ac1.0nc3hl10sd2.0dt0.5bdump10adump20tdump5ens1j05"
    )


# --------------------------------------------------------------------------------
# 2) Tests for ParserBase (abstract‐class behaviors)
# --------------------------------------------------------------------------------
class TestParserBase:
    """
    Tests for ParserBase abstract class. These tests ensure that the base class
    behaves correctly, even though it is not meant to be instantiated directly.
    """

    def test_parserbase_initialization(self):
        """
        ParserBase should raise NotImplementedError if you try to instantiate it.
        """
        with pytest.raises(NotImplementedError):
            ParserBase("dummy", lineage="space", group="all")

    def test_parserbase_str_and_repr_raise_if_uninitialized():
        """
        ParserBase is abstract; you cannot instantiate it directly. Instead, check
        that its properties complain if subclass has not set the class‐vars.
        """
        class DummyParser(ParserBase):
            # Does not initialize _geometry, _topology, _groups, etc.
            _geometry = None
            _topology = None
            _groups = None
            _genealogy_attributes = None
            _project_attributes = None

            def _initiate_attributes(self):
                pass
            def _parse_name(self):
                pass
            def _dependant_attributes(self):
                pass

    # Instantiating DummyParser should fail when accessing geometry/groups/topology
    dp = DummyParser("foo", lineage="space", group="all")
    with pytest.raises(AttributeError):
        _ = dp.geometry
    with pytest.raises(AttributeError):
        _ = dp.groups
    with pytest.raises(AttributeError):
        _ = dp.topology
    with pytest.raises(AttributeError):
        _ = dp.genealogy_attributes
    with pytest.raises(AttributeError):
        _ = dp.project_attributes

    # __str__ and __repr__ do not rely on subclass‐specific geometry/topology/groups
    # They will print placeholders (None). We check they at least contain the filename.
    s = str(dp)
    assert "Name: 'foo'" in s
    r = repr(dp)
    assert "Artifact('foo'" in r


def test_invalid_lineage_or_group_triggers_invalid_keyword():
    """
    Passing an invalid lineage or group to any ParserBase subclass should
    raise the same exception that `invalid_keyword` raises (ValueError).
    """
    # TwoMonDepCub defines _groups = ['bug','all'] and valid lineages in ParserBase._lineages
    with pytest.raises(ValueError):
        _ = TwoMonDepCub("nm2am5.0ac1.0nc3", lineage="INVALID_LINEAGE", group="bug")

    with pytest.raises(ValueError):
        _ = TwoMonDepCub("nm2am5.0ac1.0nc3", lineage="space", group="NOT_A_GROUP")


# --------------------------------------------------------------------------------
# 3) Tests for TwoMonDepCub (concrete‐class behaviors)
# --------------------------------------------------------------------------------

def test_space_lineage_parsing(simple_space_name):
    """
    When lineage='space', only the keywords 'nm', 'am', 'ac', 'nc' should be parsed.
    Dependent‐attributes (rho_bulk_*, phi_bulk_*) are NOT computed because space has no
    project_attributes in that lineage (they are only for segment/whole/ensemble_long).
    """
    artifact = simple_space_name
    parser = TwoMonDepCub(artifact, lineage="space", group="bug")

    # Check that filepath/filename were parsed correctly
    assert parser.filename == artifact
    assert parser.filepath == "N/A"

    # Check lineage and group
    assert parser.lineage == "space"
    assert parser.group == "bug"

    # The “name” (unique name parsed by _find_name) for a “space” lineage:
    #    _find_name → name = filename.split("-")[0]
    # Since there is no '-' in `artifact`, name == artifact
    assert parser.name == simple_space_name

    # Extracted (parsed) attributes
    #  nmon (nm → int), dmon (am → float), dcrowd (ac → float), ncrowd (nc → int)
    assert isinstance(parser.nmon, int)
    assert isinstance(parser.dmon, float)
    assert isinstance(parser.dcrowd, float)
    assert isinstance(parser.ncrowd, int)

    assert parser.nmon == 2
    assert pytest.approx(parser.dmon, rel=1e-8) == 5.0
    assert pytest.approx(parser.dcrowd, rel=1e-8) == 1.0
    assert parser.ncrowd == 3

    # Because lineage='space' has no project_attributes for bulk quantities, these
    # should NOT exist (AttributeError if accessed)
    with pytest.raises(AttributeError):
        _ = parser.phi_bulk_m
    with pytest.raises(AttributeError):
        _ = parser.rho_bulk_m


def test_segment_lineage_parsing_and_dependent_attributes(full_segment_name, monkeypatch):
    """
    When lineage='segment', we should parse all the keywords and then compute
    dependent attributes (rho_bulk_* and phi_bulk_* via the measurer module).
    """
    artifact = full_segment_name
    # Make sure the measurer functions give deterministic results for our test:
    #    - number_density_cube(nmon, dmon, lcube) → some float
    #    - volume_fraction_cube(nmon, dmon, lcube) → some float
    # We’ll monkeypatch them to simple lambdas so that we can predict outputs.
    monkeypatch.setattr(
        "polyphys.manage.parser.number_density_cube",
        lambda n, d, L: 100.0 * n + d + L
    )
    monkeypatch.setattr(
        "polyphys.manage.parser.volume_fraction_cube",
        lambda n, d, L: 10.0 * n + d + 0.5 * L
    )

    parser = TwoMonDepCub(artifact, lineage="segment", group="bug")

    # 1) Name logic: for lineage='segment', group='bug', name = filename.split(".bug")[0].split("-")[0]
    #    Our filename “am5.0nm2ac1.0nc3hl10sd2.0dt0.5bdump10adump20tdump5ens1j05”
    #    has no ".bug" and no “-”, so name == full string:
    assert parser.name == full_segment_name

    # 2) Parsed lineage‐attributes (all int/float as specified):
    assert parser.dmon == pytest.approx(5.0)
    assert parser.nmon == 2
    assert parser.dcrowd == pytest.approx(1.0)
    assert parser.ncrowd == 3

    # hl=10 in filename means lcube = 2 * 10 == 20
    assert parser.lcube == pytest.approx(20.0)

    assert parser.d_sur == pytest.approx(2.0)
    assert parser.dt == pytest.approx(0.5)
    assert parser.bdump == 10
    assert parser.adump == 20
    assert parser.tdump == 5

    assert parser.ensemble_id == 1
    # 'j05' → int(float("05")) == 5
    assert parser.segment_id == 5

    # 3) Dependent attributes (monomers)
    #    Because we monkeypatched number_density_cube to 100*n + d + L, and
    #    volume_fraction_cube to (10*n + d + 0.5*L), we compute:
    #      rho_bulk_m = 100*2 + 5.0 + 20 = 225.0
    #      phi_bulk_m = 10*2 + 5.0 + 0.5*20 = 20 + 5 + 10 = 35.0
    assert parser.rho_bulk_m == pytest.approx(225.0)
    assert parser.phi_bulk_m == pytest.approx(35.0)

    # 4) Dependent attributes (crowders)
    #      rho_bulk_c = 100*3 + 1.0 + 20  = 321.0
    #      phi_bulk_c = 10*3 + 1.0 + 0.5*20 = 30 + 1 + 10 = 41.0
    assert parser.rho_bulk_c == pytest.approx(321.0)
    assert parser.phi_bulk_c == pytest.approx(41.0)

    # 5) Parent‐lineage strings (via _set_parents):
    #    The object should now have four new attributes: space, ensemble, ensemble_long, whole, segment
    #    In each string, we concatenate short‐form keys+values. For example:
    #      - For lineage='segment', space includes "am5.0nm2ac1.0nc3", etc.
    #    We do a quick substring check here rather than full equality.
    assert hasattr(parser, "space")
    assert hasattr(parser, "ensemble")
    assert hasattr(parser, "ensemble_long")
    assert hasattr(parser, "whole")
    assert hasattr(parser, "segment")
    # e.g. parser.space should contain "am" + "5.0" and "nm" + "2"
    assert "am5.0" in parser.space
    assert "nm2" in parser.space


def test_filepath_and_filename_extraction_from_path(tmp_path):
    """
    If you pass a full file path (with slashes), ParserBase should split it out so that:
      - parser.filepath == the directory
      - parser.filename == the basename
    """
    # Create a dummy filename in a temporary directory
    p = tmp_path / "subdir"
    p.mkdir()
    fname = p / "nm2am5.0ac1.0nc3.txt"
    fname.write_text("dummy")  # content doesn’t matter

    # Use the POSIX‐style path (pytest/tmp_path is a pathlib.Path; convert to str)
    artifact_path = str(fname)
    # We still pretend it’s a “space” lineage with group="bug"
    parser = TwoMonDepCub(artifact_path, lineage="space", group="bug")
    assert os.path.sep in parser.filepath  # e.g. "/tmp/pytest-.../subdir"
    assert parser.filename == "nm2am5.0ac1.0nc3.txt"

    # The “name” logic: for lineage='space', group='bug', name = filename.split("-")[0]
    # Here filename has no “-”, so name == entire basename
    assert parser.name == "nm2am5.0ac1.0nc3.txt".split("-")[0]


def test_str_and_repr_include_expected_substrings(full_segment_name):
    """
    Check that __str__ and __repr__ print out geometry, group, lineage, topology,
    and filename/project_name in a reasonable way.
    """
    parser = TwoMonDepCub(full_segment_name, lineage="segment", group="bug")
    txt = str(parser)
    # Should mention "Artifact:", the filename, the geometry (cubic), the group, the lineage, and the topology (atomic)
    assert "Artifact:" in txt
    assert f"Name: '{parser.filename}'" in txt
    assert "Geometry: 'cubic'" in txt
    assert "Group: 'bug'" in txt
    assert "Lineage: 'segment'" in txt
    assert "Topology: 'atomic'" in txt

    rp = repr(parser)
    # repr should begin with "Artifact('filename' in geometry 'cubic' ..."
    assert "Artifact(" in rp
    assert parser.filename in rp
    assert "geometry 'cubic'" in rp
    assert "group 'bug'" in rp
    assert "lineage 'segment'" in rp
    assert "topology 'atomic'" in rp


def test_partial_name_with_hyphen_and_group_suffix():
    """
    When a filename contains “.all” or “.bug” (for segment/whole lineages) *and* has hyphens,
    the name extraction should split on the first occurrence of “.group” or hyphen.
    """
    # Example: “prefix‐A‐B.bug” – name should be everything before “.bug”
    artifact = "someprefix-data-b.N5.0-bdump10.bug"
    parser = TwoMonDepCub(artifact, lineage="segment", group="bug")

    # _find_name will do: filename.split(".bug")[0].split("-")[0]
    # so name = “someprefix” (the portion before the first hyphen)
    assert parser.name == "someprefix"
