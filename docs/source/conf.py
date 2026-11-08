# Configuration file for the Sphinx documentation builder.
import os
import tomllib
from pathlib import Path
from datetime import datetime
from polyphys.__version__ import __version__

# -- Project information -----------------------------------------------------
pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
with open(pyproject_path, "rb") as f:
    meta = tomllib.load(f)["project"]
project = meta["name"]
author = meta["authors"][0].get("name", "Unknown Author")
copyright = f"2022-{datetime.now().year}, {author}"
release = __version__
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.intersphinx',
    'sphinx_autodoc_typehints',
    'myst_parser',
]

# Napoleon docstring style
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True
# Render "Attributes" sections as :ivar: fields rather than standalone
# ``.. attribute::`` directives, which would collide with the same attributes
# picked up by autodoc and raise duplicate-object warnings.
napoleon_use_ivar = True
# Turn NumPy-style type text such as "int, optional", "default False", and
# "{'first', 'last'}" into emphasis and literals instead of cross-references
# that cannot resolve.
napoleon_preprocess_types = True
# ParserBase documents its class attributes and a subclass's declaration in
# sections of their own, laid out like Parameters.
napoleon_custom_sections = [
    ('Class Attributes', 'params_style'),
    ('Declaration', 'params_style'),
]
napoleon_type_aliases = {
    'callable': ':py:func:`callable`',
    'mapping': ':py:class:`~collections.abc.Mapping`',
}

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "inherited-members": False,
}
autosummary_generate = True
# List only the names in a module's __all__ in its summary tables, as autodoc
# documents them, so a table never names a private, undocumented member.
autosummary_ignore_module_all = False
typehints_fully_qualified = True
autodoc_typehints = 'description'  # Cleaner type hint formatting


def _public_scipy_sparse(module: str) -> str:
    """Name SciPy's sparse arrays by their public module, ``scipy.sparse``,
    the one SciPy's inventory lists, not the private one that defines them."""
    return 'scipy.sparse' if module.startswith('scipy.sparse._') else module


typehints_fixup_module_name = _public_scipy_sparse


# The package's runtime dependencies, NumPy, pandas, and SciPy, are
# importable wherever the docs are built. MDAnalysis and matplotlib are
# optional (the ``mda`` and ``viz`` extras) and the docs build installs the
# ``docs`` extra only, so they are mocked for the ``polyphys.probe`` and
# ``polyphys.viz`` pages.
autodoc_mock_imports: list[str] = ['MDAnalysis', 'matplotlib']

# Intersphinx mappings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None),
    'MDAnalysis': ('https://docs.mdanalysis.org/stable/', None),
    # SciPy and matplotlib are named in docstrings and type hints:
    # scipy.optimize.curve_fit, scipy.sparse.csr_array, matplotlib.axes.Axes.
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

# Cross-references that cannot resolve, each with its reason; everything
# else that nitpicky mode reports is fixed at its source.
nitpick_ignore = [
    # NumPy's inventory lists its scalar types as py:attribute, so the
    # py:class references that type hints such as NDArray[np.float64] make
    # cannot match them.
    ('py:class', 'numpy.float64'),
    ('py:class', 'numpy.int64'),
    # The dataclass signatures name NDArray by the private module that
    # defines it, numpy._typing._array_like, which no inventory lists; the
    # public name is numpy.typing.NDArray.
    ('py:class', 'numpy._typing._array_like.NDArray'),
    # The private base class of the probe analyses, which "Bases:" names.
    ('py:class', 'polyphys.probe.core._ProbeAnalysis'),
    # make_universe's return annotation names the public alias
    # MDAnalysis.Universe; MDAnalysis's inventory lists the class only as
    # MDAnalysis.core.universe.Universe.
    ('py:class', 'MDAnalysis.Universe'),
    # String forward references the typehints extension renders unqualified:
    # polyphys.spec.PostHook names 'ArrayWholes' and 'HookContext' of
    # polyphys.analyze.pipeline, which polyphys.spec cannot import, and
    # polyphys.manage.types.ParserClassT names "ParserBase".
    ('py:class', 'ArrayWholes'),
    ('py:class', 'HookContext'),
    ('py:class', 'ParserBase'),
]

# MyST configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
]

# Canonical URL: Read the Docs passes it to the build as
# READTHEDOCS_CANONICAL_URL. Local and CI builds leave it empty, so no
# canonical link is emitted.
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Paths and templates
templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# TODOs
todo_include_todos = True

# -- HTML output -------------------------------------------------------------
html_theme = 'pydata_sphinx_theme'  # Modern, well-supported theme
html_static_path = ['_static']

html_theme_options = {
    "github_url": "https://github.com/amirhs1/PolyPhys",
    "navbar_align": "content",
    "navigation_with_keys": True,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "use_edit_page_button": False,
    "show_prev_next": False,
    # "show_source": True,
}
