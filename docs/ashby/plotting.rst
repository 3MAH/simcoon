Plotting
========

The plotting module provides the high-level ``ashby_plot()`` function for
creating Ashby material-property charts, as well as the original legacy
helper functions for drawing envelopes around point clouds.

Requires ``matplotlib`` (and ``scipy`` for convex-hull envelopes).
Install with:

.. code-block:: bash

   pip install 'simcoon[ashby]'

ashby_plot
----------

.. autofunction:: simcoon.ashby.plotting.ashby_plot

**Basic usage**

.. code-block:: python

   from simcoon.ashby import load_builtin, ashby_plot

   mats = load_builtin()
   ax = ashby_plot(mats, "density", "E")

**Pre-grouped data**

.. code-block:: python

   groups = mats.group_by("category")
   ax = ashby_plot(groups, "density", "E", envelope="convex_hull")

**Performance-index guide lines**

Guide lines represent constant values of a material index in log-log space.
For a tie-rod stiffness index :math:`E/\rho`, the slope is 1 in log-log
space:

.. code-block:: python

   ax = ashby_plot(mats, "density", "E", guidelines=[
       {"slope": 1.0, "label": "E/ρ = const", "intercepts": [-1.0, 0.0, 1.0]},
   ])

**Envelope styles**

The ``envelope`` parameter controls how material families are outlined:

- ``"ellipse"`` — covariance-based ellipse (default)
- ``"convex_hull"`` — convex hull polygon
- ``"enclose"`` — smoothed convex hull (Bézier curves)
- ``None`` — scatter points only

Category colours
----------------

The ``CATEGORY_COLORS`` dictionary maps standard material families to
colours used by ``ashby_plot()``:

.. code-block:: python

   from simcoon.ashby.plotting import CATEGORY_COLORS
   # {'Metal': '#1f77b4', 'Ceramic': '#d62728', 'Polymer': '#2ca02c', ...}

Legacy functions
----------------

The following functions are preserved verbatim from the original
``simcoon.ashby`` module for backward compatibility.  They operate directly
on ``matplotlib.pyplot`` (the current axes) and on numpy point arrays.

.. autofunction:: simcoon.ashby.plotting.poly_convexHull
.. autofunction:: simcoon.ashby.plotting.poly_enclose
.. autofunction:: simcoon.ashby.plotting.ellip_enclose
.. autofunction:: simcoon.ashby.plotting.unit_vector
.. autofunction:: simcoon.ashby.plotting.length
.. autofunction:: simcoon.ashby.plotting.uv_2
