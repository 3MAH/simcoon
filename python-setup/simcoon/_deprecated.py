"""Deprecated free projector functions, scheduled for removal in simcoon 2.0.

``Ireal``/``Ireal2``/``Ivol``/``Idev``/``Idev2`` built 4th-order projectors as raw
``arma::mat`` before the typed :class:`~simcoon.tensor.Tensor4` existed. The Tensor4
factories now carry the Voigt convention in the type tag, so these are redundant in
Python. The C++ free functions stay (the ``arma::mat`` UMAT code still uses them); only
these Python entry points are deprecated. (``Ir2``/``Ir05`` are Voigt scaling *vectors*,
not projectors, and are left untouched.)

The names exported here shadow the ones pulled in by ``from simcoon._core import *`` in
``simcoon/__init__``. Removing the deprecated surface for 2.0 is a single-file delete plus
dropping the import line.
"""

import warnings

from simcoon import _core


def _deprecate_projector(name, replacement):
    impl = getattr(_core, name)

    def shim(*args, **kwargs):
        warnings.warn(
            f"simcoon.{name}() is deprecated and will be removed in simcoon 2.0; "
            f"use {replacement} instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return impl(*args, **kwargs)

    shim.__name__ = name
    shim.__qualname__ = name
    shim.__doc__ = f"Deprecated: use {replacement}."
    return shim


Ireal = _deprecate_projector("Ireal", "Tensor4.identity('stiffness').mat")
Ireal2 = _deprecate_projector("Ireal2", "Tensor4.identity('compliance').mat")
Ivol = _deprecate_projector("Ivol", "Tensor4.volumetric('stiffness').mat")
Idev = _deprecate_projector("Idev", "Tensor4.deviatoric('stiffness').mat")
Idev2 = _deprecate_projector("Idev2", "Tensor4.deviatoric('compliance').mat")
