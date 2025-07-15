# Configuration file for the Sphinx documentation builder.
from datetime import datetime
from polyphys.__version__ import __version__
import tomli

# -- Project information -----------------------------------------------------
with open("../../pyproject.toml", "rb") as f:
    meta = tomli.load(f)["project"]
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

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "inherited-members": False,
}
autosummary_generate = True
typehints_fully_qualified = True
autodoc_typehints = 'description'  # Cleaner type hint formatting

# Mock imports for optional dependencies
autodoc_mock_imports = ["MDAnalysis", "pyarrow", "statsmodels"]

# Intersphinx mappings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None),
}

# MyST configuration
myst_enable_extensions = [
    "deflist",
    "colon_fence"
]

# Sitemap configuration
html_baseurl = "https://amirhs1.github.io/PolyPhys/"

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
