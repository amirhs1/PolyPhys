API
===

Spec
----

The records that describe one study as data; no MDAnalysis needed.

.. autosummary::
   :toctree: generated

   polyphys.spec

Manage
------

.. autosummary::
   :toctree: generated

   polyphys.manage.parser
   polyphys.manage.organizer
   polyphys.manage.utils
   polyphys.manage.types

Simulate
--------

.. autosummary::
   :toctree: generated

   polyphys.simulate.core
   polyphys.simulate.log
   polyphys.simulate.dump
   polyphys.simulate.data
   polyphys.simulate.build

Analyze
-------

.. autosummary::
   :toctree: generated

   polyphys.analyze.measurer
   polyphys.analyze.stats
   polyphys.analyze.clusters
   polyphys.analyze.correlations
   polyphys.analyze.sumrule
   polyphys.analyze.density
   polyphys.analyze.pipeline

Probe
-----

Needs MDAnalysis, which the ``mda`` extra installs.

.. autosummary::
   :toctree: generated

   polyphys.probe.core
   polyphys.probe.chain_size
   polyphys.probe.binders
   polyphys.probe.histograms
   polyphys.probe.pairs

Viz
---

Needs matplotlib, which the ``viz`` extra installs; the color lists in
``colors`` do not.

.. autosummary::
   :toctree: generated

   polyphys.viz.colors
   polyphys.viz.helpers

Packages
--------

The ``simulate`` and ``probe`` packages re-export the public names of their
modules; ``viz`` imports neither of its modules on its own. A re-export is the
name documented in its module above.

.. automodule:: polyphys.simulate
   :no-members:

.. automodule:: polyphys.probe
   :no-members:

.. automodule:: polyphys.viz
   :no-members:
