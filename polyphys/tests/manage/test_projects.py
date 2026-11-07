"""Tests for :mod:`polyphys.manage.projects`, the deprecated alias of
:mod:`polyphys.spec`."""
import pytest

from polyphys import spec
from polyphys.manage import projects


@pytest.mark.parametrize('name', [
    pytest.param('ProjectSpec', id='ProjectSpec'),
    pytest.param('ProbeSpec', id='ProbeSpec'),
    pytest.param('BinderSpec', id='BinderSpec'),
    pytest.param('PerWhole', id='PerWhole'),
    pytest.param('RecenterT', id='RecenterT'),
    pytest.param('CHAIN_SIZE_OBSERVABLES', id='CHAIN_SIZE_OBSERVABLES'),
    pytest.param('BINDER_SERIES', id='BINDER_SERIES'),
])
def test_projects_serves_moved_names_with_a_warning(name):
    """A name that moved to polyphys.spec is the same object from the old
    path, with a DeprecationWarning that names the new one."""
    with pytest.warns(DeprecationWarning, match=f'{name} moved to '
                      r'polyphys\.spec'):
        served = getattr(projects, name)

    assert served is getattr(spec, name)


def test_projects_rejects_an_unknown_name():
    """A name the alias does not serve raises AttributeError."""
    with pytest.raises(AttributeError, match="no attribute 'Spec'"):
        getattr(projects, 'Spec')
