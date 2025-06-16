import tempfile
import gzip
import os
import io
import pytest
import numpy as np

from polyphys.manage import utilizer

# -----------------------
# Tests for read_camel_case
# -----------------------
def test_read_camel_case():
    assert utilizer.read_camel_case("camelCase") == ['camel', 'Case']
    assert utilizer.read_camel_case("CamelCaseString") == ['Camel', 'Case', 'String']
    assert utilizer.read_camel_case("CAMELCase") == ['CAMEL', 'Case']
    assert utilizer.read_camel_case("CamelCASE") == ['Camel', 'CASE']
    assert utilizer.read_camel_case("camel") == ['camel']
    assert utilizer.read_camel_case("CAMEL") == ['CAMEL']

# -----------------------
# Tests for to_float_if_possible
# -----------------------
def test_to_float_if_possible():
    assert utilizer.to_float_if_possible("3.14") == 3.14
    assert utilizer.to_float_if_possible("42") == 42.0
    assert utilizer.to_float_if_possible("abc") == "abc"

# -----------------------
# Tests for split_alphanumeric
# -----------------------
def test_split_alphanumeric():
    assert utilizer.split_alphanumeric("file20.5name10") == ['file', 20.5, 'name', 10]
    assert utilizer.split_alphanumeric("abc123def") == ['abc', 123, 'def']
    assert utilizer.split_alphanumeric("123.45file") == [123.45, 'file']

# -----------------------
# Tests for sort_filenames
# -----------------------
def test_sort_filenames():
    input_files = ['file2.trj', 'file1.trj', 'file1.data', 'file2.data']
    formats = ['data', ('trj', 'lammpstrj')]
    result = utilizer.sort_filenames(input_files, formats)
    assert result == [('file1.data', 'file1.trj'), ('file2.data', 'file2.trj')]

# -----------------------
# Tests for openany and openany_context
# -----------------------
def test_openany_and_context():
    # Test uncompressed file
    with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write("test content")
        tmp_path = tmp.name

    with utilizer.openany_context(tmp_path, 'r') as f:
        assert f.read() == "test content"

    f2 = utilizer.openany(tmp_path)
    assert f2.read() == "test content"
    f2.close()
    os.remove(tmp_path)

    # Test gzipped file
    with tempfile.NamedTemporaryFile(delete=False, suffix='.gz') as gz_tmp:
        with gzip.open(gz_tmp.name, 'wt') as f:
            f.write("gzipped content")
        gz_path = gz_tmp.name

    with utilizer.openany_context(gz_path, 'r') as f:
        assert f.read() == "gzipped content"

    f2 = utilizer.openany(gz_path)
    assert f2.read() == "gzipped content"
    f2.close()
    os.remove(gz_path)

# -----------------------
# Tests for round_down_first_non_zero
# -----------------------
def test_round_down_first_non_zero():
    assert utilizer.round_down_first_non_zero(1234) == 1000
    assert utilizer.round_down_first_non_zero(0.0456) == 0.04
    assert utilizer.round_down_first_non_zero(0) == 0

# -----------------------
# Tests for round_up_nearest
# -----------------------
def test_round_up_nearest():
    assert utilizer.round_up_nearest(52, 10, 0) == 50
    assert utilizer.round_up_nearest(49.5, 0.1, 1) == 49.5
    assert np.isclose(utilizer.round_up_nearest(7, 3, 1), 6.0, atol=0.1)  # 3*2=6.0

# -----------------------
# Tests for invalid_keyword
# -----------------------
def test_invalid_keyword_valid():
    utilizer.invalid_keyword("opt1", ["opt1", "opt2"])  # Should not raise

def test_invalid_keyword_raises():
    with pytest.raises(ValueError) as exc:
        utilizer.invalid_keyword("invalid", ["opt1", "opt2"], message=" is wrong.")
    assert "'invalid' is wrong." in str(exc.value)