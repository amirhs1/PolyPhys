import tempfile
import gzip
import os
import pytest
import numpy as np

from polyphys.manage import utilizer


def test_read_camel_case():
    assert utilizer.read_camel_case("camelCase") == ['camel', 'Case']
    assert utilizer.read_camel_case("CamelCaseString") \
        == ['Camel', 'Case', 'String']
    assert utilizer.read_camel_case("CAMELCase") == ['CAMEL', 'Case']
    assert utilizer.read_camel_case("CamelCASE") == ['Camel', 'CASE']
    assert utilizer.read_camel_case("camel") == ['camel']
    assert utilizer.read_camel_case("CAMEL") == ['CAMEL']


def test_to_float_if_possible():
    assert utilizer.to_float_if_possible("3.14") == 3.14
    assert utilizer.to_float_if_possible("42") == 42.0
    assert utilizer.to_float_if_possible("abc") == "abc"


def test_split_alphanumeric():
    assert utilizer.split_alphanumeric("file20.5name10") == \
        ['file', 20.5, 'name', 10]
    assert utilizer.split_alphanumeric("abc123def") == ['abc', 123, 'def']
    assert utilizer.split_alphanumeric("123.45file") == [123.45, 'file']


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
    result = utilizer.sort_filenames(sample_input_files, sample_formats)
    assert result == [('file1.data', 'file1.trj'),
                      ('file2.data', 'file2.trj'),
                      ('file3.data', 'file3.lammpstrj')]


def test_sort_filenames_empty_input(sample_formats):  # Only needs formats
    result = utilizer.sort_filenames([], sample_formats)
    assert result == []


def test_sort_filenames_single_file_type(
        sample_input_files_single_format, sample_formats_single_format):
    result = utilizer.sort_filenames(
        sample_input_files_single_format, sample_formats_single_format)
    assert result == \
        [('file1.data',), ('file2.data',), ('file3.data',)]


def test_sort_filenames_single_file_type_mutiple_formats(
        sample_input_files_single_format, sample_formats):
    result = utilizer.sort_filenames(
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

    with utilizer.openany(txt_path) as f:
        assert f.read() == "test content"

    with utilizer.openany(gz_path) as f:
        assert f.read() == "gzipped content"


def test_openany_context(sample_files):
    txt_path, gz_path = sample_files

    with utilizer.openany_context(txt_path) as f:
        assert f.read() == "test content"

    with utilizer.openany_context(gz_path) as f:
        assert f.read() == "gzipped content"


def test_round_down_first_non_zero():
    assert utilizer.round_down_first_non_zero(1234) == 1000
    assert utilizer.round_down_first_non_zero(0.0456) == 0.04
    assert utilizer.round_down_first_non_zero(0) == 0


def test_round_up_nearest():
    assert utilizer.round_up_nearest(52, 10, 0) == 50
    assert utilizer.round_up_nearest(49.5, 0.1, 1) == 49.5
    assert np.isclose(utilizer.round_up_nearest(7, 3, 1), 6.0, atol=0.1)


def test_invalid_keyword_valid():
    utilizer.invalid_keyword("opt1", ["opt1", "opt2"])


def test_invalid_keyword_raises():
    with pytest.raises(ValueError) as exc:
        utilizer.invalid_keyword(
            "invalid", ["opt1", "opt2"], message=" is wrong.")
    assert "'invalid' is wrong." in str(exc.value)
