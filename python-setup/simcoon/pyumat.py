"""Constitutive laws written in Python, driven by the C++ solver (UMAT name ``PYEXT``).

A :class:`PythonUMAT` is a point-wise small-strain constitutive law whose
:meth:`~PythonUMAT.integrate` method is called by the C++ material-point solver
(:func:`simcoon.solver.solve`) and by the batch entry point ``sim.umat("PYEXT", ...)``
exactly like a built-in kernel. The C++ side owns the Newton loop, the step
cutting, the state rollback and the finite-strain kinematics; the Python side
only integrates one increment at one material point.

Example
-------
>>> import numpy as np, simcoon as sim
>>> class Elastic(sim.PythonUMAT):
...     nstatev = 0
...     def __init__(self, E, nu):
...         self.L = sim.L_iso([E, nu], "Enu")
...     def integrate(self, *, Etot, DEtot, sigma, statev, Wm, **kw):
...         stress = self.L @ (Etot + DEtot)
...         Wm = Wm.copy()
...         Wm[0] += 0.5 * (sigma + stress) @ DEtot
...         return stress, self.L, statev, Wm
>>> step = sim.solver.StepMeca(control=["strain"] + ["stress"] * 5,
...                            value=[0.01, 0, 0, 0, 0, 0], ninc=20)
>>> res = sim.solver.solve(step, Elastic(70000., 0.3))
"""

from __future__ import annotations

import contextlib
import threading
from abc import ABC, abstractmethod
from typing import Callable, Optional

import numpy as np

import simcoon._core as _core

#: UMAT name under which the registered Python law is dispatched by the C++ code.
UMAT_NAME = "PYEXT"

_EMPTY_PROPS = np.zeros(0)
_EMPTY_PROPS.flags.writeable = False


class StepCut(_core.StepCut):
    """Raise from :meth:`PythonUMAT.integrate` to ask the solver for a smaller increment.

    The current trial is discarded (the solver restores the start-of-increment
    state) and the increment is retried with a smaller time step.

    Parameters
    ----------
    ratio : float
        Suggested reduction factor in ``(0, 1)``, forwarded as ``tnew_dt``. The material-point
        solver applies it directly when the cut interrupts a strain-driven increment and
        replaces it by its own ``div_tnew_dt`` when a Newton loop had to be abandoned.
    msg : str
        Message attached to the exception.
    """

    def __init__(self, ratio: float = 0.5, msg: str = "step cut requested"):
        super().__init__(msg)
        self.ratio = float(ratio)


class PythonUMAT(ABC):
    """Point-wise small-strain constitutive law implemented in Python.

    Subclasses set :attr:`nstatev` (number of internal variables), optionally
    :attr:`props` (material properties forwarded to the C++ material record; may
    stay empty since the object holds its own parameters) and implement
    :meth:`integrate`.

    Contract of :meth:`integrate` (identical to the C++ kernels):

    * **stateless between calls** except through ``statev`` and ``Wm``: the solver
      re-calls the same increment several times (Newton iterations, step cuts) with
      ``statev``/``sigma``/``Wm`` reset to their start-of-increment values. A
      recurrent model must keep its hidden state inside ``statev``.
    * ``start`` is ``True`` on the first call of a block (``Time == 0``): initialise
      ``statev``/``Wm`` there. ``statev`` arrives zero-filled.
    * the solver primes the tangent at the start of **every** block with a
      zero-increment call (``DTime``, ``DT`` and ``DEtot`` all zero), and finite-element
      couplers do the same at initialisation. That probe is not a step of the loading
      history: :class:`PythonUMAT` answers it as a pure tangent query (see
      :meth:`__call__`), so ``integrate`` need not care. A bare callable registered
      instead of a :class:`PythonUMAT` must handle it itself.
    * Voigt order ``11 22 33 12 13 23``, **engineering** shear strains, quantities in
      the material (local) frame. ``Etot`` is the strain at the beginning of the
      increment, ``DEtot`` the increment, ``sigma`` the stress at the beginning of the
      increment. Under finite strain the caller feeds the logarithmic strain and
      expects the Kirchhoff stress (same convention as ELISO/EPICP).
    * ``ndi`` follows the classical convention: 3 = 3D (plane strain is 3D with a null
      out-of-plane strain), 2 = plane stress, 1 = uniaxial. ``nshr`` = number of shear
      components.
    * ``tangent_mode``: 0 = return the elastic operator as ``Lt``, 1 = continuum
      tangent, 2 = algorithmic (consistent) tangent.
    * ``Wm`` is the accumulated ``(Wm, Wm_r, Wm_ir, Wm_d)`` at the beginning of the
      increment and must be returned accumulated.
    * returns ``(sigma (6,), Lt (6,6), statev (nstatev,), Wm (4,)[, L (6,6)])`` as
      float arrays (lists / float32 / CPU torch tensors are converted). ``L`` defaults
      to ``Lt``. Raise :class:`StepCut` to request a smaller increment; any other
      exception aborts the solve and is re-raised unchanged.
    """

    #: Material properties forwarded to the C++ material record (floats). May be empty.
    #: The default is a shared read-only array: a law with properties assigns its own
    #: (``self.props = np.array([...])``) rather than mutating this one in place.
    props: np.ndarray = _EMPTY_PROPS
    #: Number of internal state variables.
    nstatev: int = 0

    @abstractmethod
    def integrate(
        self,
        *,
        Etot: np.ndarray,
        DEtot: np.ndarray,
        sigma: np.ndarray,
        DR: np.ndarray,
        props: np.ndarray,
        statev: np.ndarray,
        T: float,
        DT: float,
        Time: float,
        DTime: float,
        Wm: np.ndarray,
        ndi: int,
        nshr: int,
        start: bool,
        tangent_mode: int,
        **_,
    ):
        """Integrate one increment at one material point (see class docstring)."""

    def __call__(self, **kwargs):
        """Entry point used by the C++ bridge: :meth:`integrate`, except that the
        zero-increment tangent probe (``DTime == 0`` and ``DEtot == 0``) is a pure tangent
        query: stress, ``statev`` and ``Wm`` are returned exactly as received (a zero
        increment changes none of them), only ``Lt``/``L`` come from the evaluation."""
        probe = (float(kwargs["DTime"]) == 0.0 and float(kwargs["DT"]) == 0.0
                 and not np.any(kwargs["DEtot"]))
        if not probe:
            return self.integrate(**kwargs)
        # copies first: a law updating its inputs in place must not leak the probe
        sigma0 = np.array(kwargs["sigma"], dtype=float, copy=True)
        statev0 = np.array(kwargs["statev"], dtype=float, copy=True)
        Wm0 = np.array(kwargs["Wm"], dtype=float, copy=True)
        out = self.integrate(**kwargs)
        return (sigma0, out[1], statev0, Wm0, *out[4:])


_current: Optional[object] = None
_owner: Optional[int] = None


def current():
    """Return the currently registered Python law (or ``None``)."""
    return _current


def _as_callable(umat) -> Callable:
    if callable(umat):                  # PythonUMAT instances are callable (probe rule inside)
        return umat
    raise TypeError(
        "expected a PythonUMAT instance or a callable with the integrate() keyword "
        f"signature, got {type(umat).__name__}"
    )


def register(umat) -> None:
    """Register ``umat`` as the process-wide ``PYEXT`` law (prefer :func:`registered`).

    The slot is process-wide: registering from a second thread while another thread's
    law is in place raises instead of silently rebinding that thread's running solve.
    """
    global _current, _owner
    me = threading.get_ident()
    if _current is not None and _owner != me and _owner in {t.ident for t in threading.enumerate()}:
        raise RuntimeError("a Python UMAT is already registered by another live thread "
                           "(PYEXT serves one law per process)")
    _core.register_python_umat(_as_callable(umat))
    _current, _owner = umat, me


def unregister() -> None:
    """Remove the process-wide ``PYEXT`` law (only the registering thread may do so while
    it is alive: a solve running on that thread must not lose its law)."""
    global _current, _owner
    me = threading.get_ident()
    if _current is not None and _owner != me and _owner in {t.ident for t in threading.enumerate()}:
        raise RuntimeError("the Python UMAT was registered by another live thread")
    _core.unregister_python_umat()
    _current, _owner = None, None


@contextlib.contextmanager
def registered(umat):
    """Context manager: register ``umat`` for the block, restore the previous law on exit.

    Parameters
    ----------
    umat : PythonUMAT or callable
        The law to serve under the ``PYEXT`` name.
    """
    previous = _current
    register(umat)
    try:
        yield umat
    finally:
        if previous is None:
            unregister()
        else:
            register(previous)


__all__ = [
    "UMAT_NAME", "PythonUMAT", "StepCut", "registered", "register", "unregister",
    "current",
]
