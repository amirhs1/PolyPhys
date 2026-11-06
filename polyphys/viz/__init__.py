"""
The viz layer: color lists and small matplotlib helpers for the studies'
figures, which the ``viz`` extra makes usable
(``python -m pip install 'polyphys[viz]'``).

:mod:`polyphys.viz.colors` holds the color lists, which need no plotting
library, and :func:`~polyphys.viz.colors.truncate_colormap`, which imports
matplotlib when called. :mod:`polyphys.viz.helpers` needs matplotlib at
import and raises :exc:`ImportError` with that hint when it is missing. This
package imports neither module on its own; import the one you need.

.. versionadded:: 0.8.0
"""

__all__ = ["colors", "helpers"]
