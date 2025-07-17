import tempfile
import gzip
import os
import pytest
import numpy as np

from polyphys.manage.utils import (
    read_camel_case,
    to_float_if_possible,
    split_alphanumeric,
    sort_filenames,
    openany,
    openany_context,
    round_down_first_non_zero,
    round_up_nearest,
    invalid_keyword,
    number_density_cube,
    volume_fraction_cube,
    number_density_cylinder,
    volume_fraction_cylinder
)


def test_read_camel_case():
    assert read_camel_case("camelCase") == ['camel', 'Case']
    assert read_camel_case("CamelCaseString") \
        == ['Camel', 'Case', 'String']
    assert read_camel_case("CAMELCase") == ['CAMEL', 'Case']
    assert read_camel_case("CamelCASE") == ['Camel', 'CASE']
    assert read_camel_case("camel") == ['camel']
    assert read_camel_case("CAMEL") == ['CAMEL']


def test_to_float_if_possible():
    assert to_float_if_possible("3.14") == 3.14
    assert to_float_if_possible("42") == 42.0
    assert to_float_if_possible("abc") == "abc"


def test_split_alphanumeric():
    assert split_alphanumeric("file20.5name10") == \
        ['file', 20.5, 'name', 10]
    assert split_alphanumeric("abc123def") == ['abc', 123, 'def']
    assert split_alphanumeric("123.45file") == [123.45, 'file']


@pytest.fixture
def sample_formats():
    return ['data', ('trj', 'lammpstrj')]


@pytest.fixture
def sample_input_files():
    return ['file2.trj', 'file1.trj',
            'file1.data', 'file2.data',
            'file3.lammpstrj', 'file3.data']


@pytest.fixture
def sample_formats_single_format():
    return ['data']


@pytest.fixture
def sample_input_files_single_format():
    return ['file1.data', 'file2.data', 'file3.data']


def test_sort_filenames_standard(sample_input_files, sample_formats):
    result = sort_filenames(sample_input_files, sample_formats)
    assert result == [('file1.data', 'file1.trj'),
                      ('file2.data', 'file2.trj'),
                      ('file3.data', 'file3.lammpstrj')]


def test_sort_filenames_empty_input(sample_formats):  # Only needs formats
    result = sort_filenames([], sample_formats)
    assert result == []


def test_sort_filenames_single_file_type(
        sample_input_files_single_format, sample_formats_single_format):
    result = sort_filenames(
        sample_input_files_single_format, sample_formats_single_format)
    assert result == \
        [('file1.data',), ('file2.data',), ('file3.data',)]


def test_sort_filenames_single_file_type_mutiple_formats(
        sample_input_files_single_format, sample_formats):
    result = sort_filenames(
        sample_input_files_single_format, sample_formats)
    assert result == []  # No 'trj' files → empty zip


@pytest.fixture
def sample_files():
    # Create a plain text file
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write("test content")
        txt_path = tmp.name

    # Create a gzipped text file
    with tempfile.NamedTemporaryFile(delete=False, suffix='.gz') as gz_tmp:
        with gzip.open(gz_tmp.name, 'wt') as f:
            f.write("gzipped content")
        gz_path = gz_tmp.name

    yield txt_path, gz_path

    # Teardown: remove files after test
    os.remove(txt_path)
    os.remove(gz_path)


def test_openany(sample_files):
    txt_path, gz_path = sample_files

    with openany(txt_path) as f:
        assert f.read() == "test content"

    with openany(gz_path) as f:
        assert f.read() == "gzipped content"


def test_openany_context(sample_files):
    txt_path, gz_path = sample_files

    with openany_context(txt_path) as f:
        assert f.read() == "test content"

    with openany_context(gz_path) as f:
        assert f.read() == "gzipped content"


def test_round_down_first_non_zero():
    assert round_down_first_non_zero(1234) == 1000
    assert round_down_first_non_zero(0.0456) == 0.04
    assert round_down_first_non_zero(0) == 0


def test_round_up_nearest():
    assert round_up_nearest(52, 10, 0) == 50
    assert round_up_nearest(49.5, 0.1, 1) == 49.5
    assert np.isclose(round_up_nearest(7, 3, 1), 6.0, atol=0.1)


def test_invalid_keyword_valid():
    invalid_keyword("opt1", ["opt1", "opt2"])


def test_invalid_keyword_raises():
    with pytest.raises(ValueError) as exc:
        invalid_keyword(
            "invalid", ["opt1", "opt2"], message=" is wrong.")
    assert "'invalid' is wrong." in str(exc.value)


def test_number_density_cube():
    assert number_density_cube(1000, 1.0, 10.0) == \
        pytest.approx(1.37174, rel=1e-5)
    assert number_density_cube(1000, 1.0, 10.0, pbc=True) == \
        pytest.approx(1.0, rel=1e-5)


def test_volume_fraction_cube():
    assert volume_fraction_cube(1000, 1.0, 10.0) == \
        pytest.approx(0.7182428741954303, rel=1e-5)  # or lower rel
    assert volume_fraction_cube(1000, 1.0, 10.0, pbc=True) == \
        pytest.approx(0.5235987755982988, rel=1e-5)


def test_volume_fraction_cube_phi_greater_than_one():
    with pytest.raises(ValueError, match="Volume fraction exceeds 1.0"):
        volume_fraction_cube(n_atom=1000, d_atom=1.0, l_cube=5.0)
        # Large particles in a small box → phi > 1


def test_number_density_cylinder():
    assert number_density_cylinder(100, 1.0, 10.0, 5.0) == \
        pytest.approx(0.88419, rel=1e-5)
    assert number_density_cylinder(100, 1.0, 10.0, 5.0, pbc=True) == \
        pytest.approx(0.79577, rel=1e-5)


def test_volume_fraction_cylinder():
    assert volume_fraction_cylinder(100, 1.0, 10.0, 5.0) == \
        pytest.approx(0.462962963, rel=1e-5)  # use 9 digits or so
    assert volume_fraction_cylinder(100, 1.0, 10.0, 5.0, pbc=True) == \
        pytest.approx(0.4166666666666667, rel=1e-5)


def test_volume_fraction_cylinder_phi_greater_than_one():
    with pytest.raises(ValueError, match="Volume fraction exceeds 1.0"):
        volume_fraction_cylinder(n_atom=500, d_atom=1.0, l_cyl=5.0, d_cyl=3.0)
        # Too many large particles → phi > 1
