"""
Version information for PolyPhys - :mod:`PolyPhys.version`
=============================================

The version information in :mod:`PolyPhys.version` indicates the
release of PolyPhys. PolyPhys uses a `semantic versioning`.

Given a version number MAJOR.MINOR.PATCH, we increment the

1. **MAJOR** version when we make **incompatible API changes**,
2. **MINOR** version when we **add functionality** in a
   **backwards-compatible** manner, and
3. **PATCH** version when we make backwards-compatible **bug fixes**.

However, as long as the **MAJOR** number is **0** (i.e. the API has
not stabilized), even **MINOR** increases *may* introduce incompatible
API changes. As soon as we have a 1.0.0 release, the public API can
only be changed in a backward-incompatible manner with an increase in
MAJOR version.

Additional labels for pre-release and build metadata are available as
extensions to the MAJOR.MINOR.PATCH format.

.. _`semantic versioning`: http://semver.org/

Data
----

.. autodata:: __version__

"""

# keep __version__ in separate file to avoid circular imports

#: Release of PolyPhys as a string, using `semantic versioning`_.
__version__ = "0.4.0"
