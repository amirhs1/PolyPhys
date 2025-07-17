import pytest


@pytest.fixture
def make_temp_file(tmp_path):
    """
    Create an empty temporary file with a given filename and return its full
    path. Useful for testing parsers that expect a file path as input.
    """
    def _make_temp_file(filename: str):
        f = tmp_path / filename
        f.touch()
        return str(f)
    return _make_temp_file
