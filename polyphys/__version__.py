"""
Version information for PolyPhys, :mod:`polyphys.__version__`.

:data:`__version__` is the release of PolyPhys, a `semantic versioning`_
string, ``MAJOR.MINOR.PATCH`` with optional pre-release and build labels.
Which names the versioning promise covers, which modules are stable or
provisional, and how a name is deprecated are stated once, on the
documentation's "API stability" page (``docs/source/stability.md``).

.. _`semantic versioning`: https://semver.org/spec/v2.0.0.html

Data
----

.. autodata:: __version__

"""

# keep __version__ in separate file to avoid circular imports

#: Release of PolyPhys as a string, using `semantic versioning`_.
__version__ = "1.0.0"
