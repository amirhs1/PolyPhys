"""Tests for :mod:`polyphys.manage.utils`.

The string and filename expectations are worked out by hand from each
function's stated rule. The bulk-quantity expectations are closed forms: the
number density ``n / v_avail``, with ``v_avail = (l - d)^3`` in a cube and
``pi (l - d) (D - d)^2 / 4`` in a cylinder (``l`` in place of ``l - d``
under periodic boundaries), and the volume fraction ``rho pi d^3 / 6``. They
are compared at ``rel=1e-12``, far above the float64 round-off of the few
operations behind them (about 1e-16).
"""
import gzip
import re
from pathlib import PosixPath

import numpy as np
import pytest

from polyphys.manage.utils import (
    read_camel_case,
    to_float_if_possible,
    split_alphanumeric,
    normalize_exts,
    sort_filenames,
    sort_filenames_str,
    openany,
    openany_context,
    round_down_first_non_zero,
    round_to_multiple,
    invalid_keyword,
    number_density_cube,
    volume_fraction_cube,
    number_density_cylinder,
    volume_fraction_cylinder
)


@pytest.fixture
def make_text_file(tmp_path):
    """Return a factory that writes text to a file in ``tmp_path`` and
    returns its path; a name ending in ``.gz`` is written gzipped."""
    def _make(name, content):
        path = tmp_path / name
        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "wt") as handle:
            handle.write(content)
        return str(path)
    return _make


@pytest.mark.parametrize(
    "string,expected",
    [
        pytest.param("camelCase", ["camel", "Case"], id="lower-camel"),
        pytest.param(
            "CamelCaseString", ["Camel", "Case", "String"], id="upper-camel"
        ),
        pytest.param("CAMELCase", ["CAMEL", "Case"], id="leading-acronym"),
        pytest.param("CamelCASE", ["Camel", "CASE"], id="trailing-acronym"),
        pytest.param("camel", ["camel"], id="one-lower-word"),
        pytest.param("CAMEL", ["CAMEL"], id="one-acronym"),
    ],
)
def test_read_camel_case_splits_at_case_changes(string, expected):
    """A camel-case string splits into its words, and a run of capitals is
    kept as one word."""
    assert read_camel_case(string) == expected


@pytest.mark.parametrize(
    "string,expected",
    [
        pytest.param("3.14", 3.14, id="decimal"),
        pytest.param("42", 42.0, id="integer"),
        pytest.param("abc", "abc", id="not-a-number"),
    ],
)
def test_to_float_if_possible_converts_numbers_only(string, expected):
    """A numeric string becomes a float and any other string comes back
    unchanged."""
    assert to_float_if_possible(string) == expected


@pytest.mark.parametrize(
    "string,expected",
    [
        pytest.param(
            "file20.5name10", ["file", 20.5, "name", 10], id="float-and-int"
        ),
        pytest.param("abc123def", ["abc", 123, "def"], id="int-inside"),
        pytest.param("123.45file", [123.45, "file"], id="leading-float"),
    ],
)
def test_split_alphanumeric_separates_numbers_from_letters(string, expected):
    """Runs of digits become ints or floats and the letters between them stay
    strings, in order."""
    assert split_alphanumeric(string) == expected


@pytest.mark.parametrize(
    "exts,expected",
    [
        pytest.param("npy", (".npy",), id="str-without-dot"),
        pytest.param(".trj", (".trj",), id="str-with-dot"),
        pytest.param(
            ("trj", ".lammpstrj"), (".trj", ".lammpstrj"), id="tuple-mixed"
        ),
    ],
)
def test_normalize_exts_gives_dot_prefixed_tuple(exts, expected):
    """An extension or a tuple of them becomes a tuple of extensions, each
    with exactly one leading dot."""
    assert normalize_exts(exts) == expected


@pytest.mark.parametrize(
    "filenames,formats,expected",
    [
        pytest.param(
            ["file2.trj", "file1.trj", "file1.data", "file2.data",
             "file3.lammpstrj", "file3.data"],
            ["data", ("trj", "lammpstrj")],
            [
                (PosixPath("file1.data"), PosixPath("file1.trj")),
                (PosixPath("file2.data"), PosixPath("file2.trj")),
                (PosixPath("file3.data"), PosixPath("file3.lammpstrj")),
            ],
            id="equivalent-extensions",
        ),
        pytest.param(
            [], ["data", ("trj", "lammpstrj")], [], id="empty-input"
        ),
        pytest.param(
            ["file1.data", "file2.data", "file3.data"],
            ["data"],
            [
                (PosixPath("file1.data"),),
                (PosixPath("file2.data"),),
                (PosixPath("file3.data"),),
            ],
            id="one-format",
        ),
        pytest.param(
            ["file1.data", "file1.trj", "file1.npy",
             "file2.npy", "file2.trj", "file2.data"],
            ["data", "trj", "npy"],
            [
                (PosixPath("file1.data"), PosixPath("file1.trj"),
                 PosixPath("file1.npy")),
                (PosixPath("file2.data"), PosixPath("file2.trj"),
                 PosixPath("file2.npy")),
            ],
            id="slot-order-follows-formats",
        ),
        pytest.param(
            ["sample1.0.tar.gz", "sample1.0.gz",
             "sample2.0.tar.gz", "sample2.0.gz"],
            ["tar.gz", "gz"],
            [
                (PosixPath("sample1.0.tar.gz"), PosixPath("sample1.0.gz")),
                (PosixPath("sample2.0.tar.gz"), PosixPath("sample2.0.gz")),
            ],
            id="compound-extension-distinct",
        ),
        pytest.param(
            ["a.data", "b.DATA"],
            ["DATA"],
            [(PosixPath("b.DATA"),)],
            id="case-sensitive",
        ),
    ],
)
def test_sort_filenames_pairs_by_basename_in_natural_order(
    filenames, formats, expected
):
    """Files sharing a basename are grouped one per format, in the order of
    the formats, and the groups are sorted naturally by basename."""
    assert sort_filenames(filenames, formats) == expected


def test_sort_filenames_drops_incomplete_basename_with_warning():
    """A basename missing one of its formats is left out, with a warning that
    names the missing counterpart."""
    files = ["file2.trj", "file1.trj", "file1.data", "file3.lammpstrj",
             "file3.data"]
    expected_message = (
        "'file2.trj' does not have matching 'file2.data' so it is ignored"
    )
    with pytest.warns(UserWarning, match=re.escape(expected_message)):
        result = sort_filenames(files, ["data", ("trj", "lammpstrj")])
    assert result == [
        (PosixPath("file1.data"), PosixPath("file1.trj")),
        (PosixPath("file3.data"), PosixPath("file3.lammpstrj")),
    ]


def test_sort_filenames_missing_format_group_gives_empty_with_warning():
    """When no file matches one format group, every basename is incomplete,
    so the result is empty and each file draws a warning."""
    files = ["file1.data", "file2.data", "file3.data"]
    with pytest.warns(UserWarning) as record:
        result = sort_filenames(files, ["data", ("trj", "lammpstrj")])
    assert result == []
    assert [str(warning.message) for warning in record] == [
        "'file1.data' does not have matching 'file1.trj' so it is ignored",
        "'file2.data' does not have matching 'file2.trj' so it is ignored",
        "'file3.data' does not have matching 'file3.trj' so it is ignored",
    ]


def test_sort_filenames_rejects_two_files_in_one_slot():
    """Two equivalent extensions for one basename fill the same slot twice,
    which is rejected."""
    files = ["file1.trj", "file1.data", "file1.lammpstrj"]
    with pytest.raises(
        ValueError, match="Duplicate files for basename 'file1'"
    ):
        sort_filenames(files, ["data", ("trj", "lammpstrj")])


def test_sort_filenames_rejects_extension_in_two_format_groups():
    """An extension listed in two format groups is ambiguous and
    rejected."""
    formats = ["data", ("trj", "data")]
    files = ["file1.trj", "file1.data"]
    expected = (
        "Extension '.data' appears in multiple format groups: "
        "('.data',) and ('.trj', '.data')"
    )
    with pytest.raises(ValueError, match=re.escape(expected)):
        sort_filenames(files, formats)


@pytest.mark.parametrize(
    "filenames,formats,expected",
    [
        pytest.param(
            ["file2.trj", "file1.trj", "file1.data", "file2.data",
             "file3.lammpstrj", "file3.data"],
            ["data", ("trj", "lammpstrj")],
            [("file1.data", "file1.trj"),
             ("file2.data", "file2.trj"),
             ("file3.data", "file3.lammpstrj")],
            id="equivalent-extensions",
        ),
        pytest.param(
            [], ["data", ("trj", "lammpstrj")], [], id="empty-input"
        ),
        pytest.param(
            ["file1.data", "file2.data", "file3.data"],
            ["data"],
            [("file1.data",), ("file2.data",), ("file3.data",)],
            id="one-format",
        ),
        pytest.param(
            ["file1.data", "file2.data", "file3.data"],
            ["data", ("trj", "lammpstrj")],
            [],
            id="missing-format-group",
        ),
    ],
)
def test_sort_filenames_str_zips_naturally_sorted_groups(
    filenames, formats, expected
):
    """Each format group is sorted naturally and the groups are zipped, so a
    group with no files gives an empty result."""
    assert sort_filenames_str(filenames, formats) == expected


@pytest.mark.parametrize(
    "name,content",
    [
        pytest.param("plain.txt", "test content", id="plain"),
        pytest.param("packed.txt.gz", "gzipped content", id="gzipped"),
    ],
)
def test_openany_reads_plain_and_gzipped_files(make_text_file, name, content):
    """A plain file and a gzipped one both read back as their text."""
    path = make_text_file(name, content)
    with openany(path) as handle:
        text = handle.read()
    assert text == content


@pytest.mark.parametrize(
    "name,content",
    [
        pytest.param("plain.txt", "test content", id="plain"),
        pytest.param("packed.txt.gz", "gzipped content", id="gzipped"),
    ],
)
def test_openany_context_reads_plain_and_gzipped_files(
    make_text_file, name, content
):
    """Inside the context, a plain file and a gzipped one both read back as
    their text."""
    path = make_text_file(name, content)
    with openany_context(path) as handle:
        text = handle.read()
    assert text == content


def test_openany_context_lets_body_error_propagate(make_text_file, capsys):
    """An error raised in the with-block propagates unchanged, unprinted."""
    path = make_text_file("plain.txt", "test content")
    with pytest.raises(RuntimeError, match="raised in the body"):
        with openany_context(path):
            raise RuntimeError("raised in the body")
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize(
    "num,expected",
    [
        pytest.param(1234, 1000, id="integer"),
        pytest.param(0.0456, 0.04, id="fraction"),
        pytest.param(0, 0, id="zero"),
    ],
)
def test_round_down_first_non_zero_keeps_leading_digit(num, expected):
    """A number rounds down to its leading digit times its power of ten, and
    zero stays zero."""
    result = round_down_first_non_zero(num)
    # float64 round-off left after rounding; atol for the zero case:
    assert result == pytest.approx(expected, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize(
    "dividend,divider,round_to,expected",
    [
        pytest.param(52, 10, 0, 50.0, id="integers-down"),
        pytest.param(49.5, 0.1, 1, 49.5, id="already-a-multiple"),
        pytest.param(7, 3, 1, 6.0, id="integers-down-to-6"),
        pytest.param(0.0874, 0.025, 3, 0.075, id="phi-below-midpoint"),
        pytest.param(0.0876, 0.025, 3, 0.1, id="phi-above-midpoint"),
        pytest.param(0.1001, 0.025, 3, 0.1, id="phi-just-above-multiple"),
        pytest.param(0.11, 0.025, 3, 0.1, id="label-below-value"),
    ],
)
def test_round_to_multiple_rounds_to_nearest_multiple(
    dividend, divider, round_to, expected
):
    """The result is the multiple of divider nearest to dividend, which can
    lie below dividend."""
    result = round_to_multiple(dividend, divider, round_to)
    # float64 round-off left after rounding to round_to decimals:
    assert result == pytest.approx(expected, rel=0, abs=1e-12)


@pytest.mark.parametrize(
    "dividend,divider,round_to,expected",
    [
        pytest.param(2.5, 1.0, 0, 2.0, id="tie-down-to-even"),
        pytest.param(3.5, 1.0, 0, 4.0, id="tie-up-to-even"),
        pytest.param(0.0625, 0.025, 3, 0.05, id="phi-tie-to-even"),
    ],
)
def test_round_to_multiple_breaks_ties_to_even(
    dividend, divider, round_to, expected
):
    """A dividend exactly halfway between two multiples goes to the even
    multiple, as Python's round does."""
    result = round_to_multiple(dividend, divider, round_to)
    # float64 round-off left after rounding to round_to decimals:
    assert result == pytest.approx(expected, rel=0, abs=1e-12)


def test_invalid_keyword_accepts_valid_keyword():
    """A keyword in the valid list passes without raising."""
    assert invalid_keyword("opt1", ["opt1", "opt2"]) is None


def test_invalid_keyword_raises_with_custom_message():
    """A keyword outside the valid list raises with the given message
    appended to the quoted keyword."""
    with pytest.raises(ValueError, match=re.escape("'invalid' is wrong.")):
        invalid_keyword("invalid", ["opt1", "opt2"], message=" is wrong.")


@pytest.mark.parametrize(
    "pbc,expected",
    [
        pytest.param(False, 1000 / 9**3, id="walls"),
        pytest.param(True, 1000 / 10**3, id="periodic"),
    ],
)
def test_number_density_cube_divides_by_available_volume(pbc, expected):
    """1000 particles of diameter 1 in a cube of side 10 have the density
    1000 / 9^3 with walls and 1000 / 10^3 under periodic boundaries."""
    result = number_density_cube(1000, 1.0, 10.0, pbc=pbc)
    # float64 round-off of a cube and a division:
    assert result == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "pbc,expected",
    [
        pytest.param(False, 1000 / 9**3 * np.pi / 6, id="walls"),
        pytest.param(True, 1000 / 10**3 * np.pi / 6, id="periodic"),
    ],
)
def test_volume_fraction_cube_is_density_times_bead_volume(pbc, expected):
    """The volume fraction is the number density times pi d^3 / 6."""
    result = volume_fraction_cube(1000, 1.0, 10.0, pbc=pbc)
    # float64 round-off of a few products and divisions:
    assert result == pytest.approx(expected, rel=1e-12)


def test_volume_fraction_cube_allows_nominal_phi_above_one():
    """1000 unit beads in a cube of side 5 with walls give the nominal
    fraction 1000 / 4^3 * pi / 6, above one, instead of raising."""
    result = volume_fraction_cube(n_atom=1000, d_atom=1.0, l_cube=5.0)
    # float64 round-off of a few products and divisions:
    assert result == pytest.approx(1000 / 4**3 * np.pi / 6, rel=1e-12)


@pytest.mark.parametrize(
    "d_atom,l_cube,pbc",
    [
        pytest.param(5.0, 5.0, False, id="as-wide-as-box"),
        pytest.param(6.0, 5.0, False, id="wider-than-box"),
    ],
)
def test_number_density_cube_rejects_box_without_room(d_atom, l_cube, pbc):
    """A particle that leaves its center no room in a walled cube raises
    ValueError instead of returning a meaningless density."""
    with pytest.raises(ValueError, match="no room for particle centers"):
        number_density_cube(10, d_atom, l_cube, pbc=pbc)


@pytest.mark.parametrize(
    "pbc,expected",
    [
        pytest.param(False, 100 / (np.pi * 9 * 4**2 / 4), id="walls"),
        pytest.param(True, 100 / (np.pi * 10 * 4**2 / 4), id="periodic"),
    ],
)
def test_number_density_cylinder_divides_by_available_volume(pbc, expected):
    """100 unit beads in a cylinder of length 10 and diameter 5 have the
    density 100 / (pi 9 4^2 / 4) with end walls and 100 / (pi 10 4^2 / 4)
    under periodic boundaries along the axis."""
    result = number_density_cylinder(100, 1.0, 10.0, 5.0, pbc=pbc)
    # float64 round-off of a few products and a division:
    assert result == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "pbc,expected",
    [
        pytest.param(False, 100 / 216, id="walls"),
        pytest.param(True, 100 / 240, id="periodic"),
    ],
)
def test_volume_fraction_cylinder_is_density_times_bead_volume(pbc, expected):
    """The volume fraction is the number density times pi d^3 / 6, which
    cancels pi: 100 / (9 * 4^2 * 6 / 4) with end walls and
    100 / (10 * 4^2 * 6 / 4) under periodic boundaries."""
    result = volume_fraction_cylinder(100, 1.0, 10.0, 5.0, pbc=pbc)
    # float64 round-off of a few products and divisions:
    assert result == pytest.approx(expected, rel=1e-12)


def test_volume_fraction_cylinder_allows_nominal_phi_above_one():
    """500 unit beads in a cylinder of length 5 and diameter 3 with end walls
    give the nominal fraction 500 / (4 * 2^2 * 6 / 4) = 500 / 24, above one,
    instead of raising."""
    result = volume_fraction_cylinder(
        n_atom=500, d_atom=1.0, l_cyl=5.0, d_cyl=3.0)
    # float64 round-off of a few products and divisions:
    assert result == pytest.approx(500 / 24, rel=1e-12)


@pytest.mark.parametrize(
    "d_atom,l_cyl,d_cyl,pbc",
    [
        pytest.param(5.0, 20.0, 5.0, True, id="as-wide-as-cylinder"),
        pytest.param(6.0, 20.0, 5.0, True, id="wider-than-cylinder"),
        pytest.param(5.0, 20.0, 5.0, False, id="as-wide-with-end-walls"),
        pytest.param(1.0, 1.0, 5.0, False, id="as-long-as-cylinder"),
    ],
)
def test_number_density_cylinder_rejects_cylinder_without_room(
        d_atom, l_cyl, d_cyl, pbc):
    """A particle that leaves its center no room across or along the
    cylinder raises ValueError instead of dividing by zero or returning a
    density from a squared negative diameter."""
    with pytest.raises(ValueError, match="no room for particle centers"):
        number_density_cylinder(10, d_atom, l_cyl, d_cyl, pbc=pbc)
