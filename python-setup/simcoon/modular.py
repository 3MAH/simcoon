"""Modular UMAT configuration for composable constitutive models.

This module provides a Pythonic interface to the modular UMAT system.
It allows users to compose constitutive models from building blocks:
- Elasticity (isotropic, cubic, transversely isotropic, orthotropic)
- Plasticity (yield criteria + isotropic/kinematic hardening)
- Viscoelasticity (Prony series)
- Damage (linear, exponential, power-law, Weibull)

Example
-------
>>> from simcoon.modular import (
...     ModularMaterial, IsotropicElasticity,
...     Plasticity, VonMisesYield, VoceHardening
... )
>>> mat = ModularMaterial(
...     elasticity=IsotropicElasticity(C1=210000., C2=0.3, alpha=1.2e-5),
...     mechanisms=[
...         Plasticity(
...             sigma_Y=300.,
...             yield_criterion=VonMisesYield(),
...             isotropic_hardening=VoceHardening(Q=100., b=10.),
...         )
...     ]
... )
>>> props = mat.props    # numpy array for sim.umat("MODUL", ...)
>>> nstatev = mat.nstatev  # number of state variables required
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from enum import IntEnum
from typing import ClassVar, List, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from simcoon.rotation import Orientation, as_direction


__all__ = [
    # Enums
    "ElasticityType", "HyperPotential", "MuscleFibreLaw", "VolumetricPotential", "YieldType", "IsoHardType", "KinHardType",
    "DamageType", "MechanismType",
    # Elastic-constant conventions
    "IsoConvention", "CubicConvention",
    "IsotransConvention", "OrthoConvention",
    # Elasticity
    "IsotropicElasticity", "CubicElasticity",
    "TransverseIsotropicElasticity", "OrthotropicElasticity",
    "NeoHookeanElasticity", "MooneyRivlinElasticity", "YeohElasticity",
    "IsiharaElasticity", "GentThomasElasticity", "SwansonElasticity",
    "HolzapfelElasticity", "MuscleElasticity",
    # Yield criteria
    "VonMisesYield", "TrescaYield", "DruckerYield",
    "HillYield", "DFAYield", "AnisotropicYield",
    # Isotropic hardening
    "NoIsotropicHardening", "LinearIsotropicHardening",
    "PowerLawHardening", "VoceHardening", "CombinedVoceHardening",
    # Kinematic hardening
    "NoKinematicHardening", "PragerHardening",
    "ArmstrongFrederickHardening", "ChabocheHardening",
    # Mechanisms
    "Plasticity", "Viscoelasticity", "Damage",
    # Type aliases
    "Elasticity", "YieldCriterion", "IsotropicHardening",
    "KinematicHardening", "Mechanism",
    # Orchestrator
    "ModularMaterial",
    # Factory functions
    "elastic_model", "elastoplastic_model", "viscoelastic_model",
]


# ============================================================================
# Enums (mirrors C++ enums in the Modular UMAT headers)
# ============================================================================

class ElasticityType(IntEnum):
    """Type of elasticity block (the first four are linear)."""
    ISOTROPIC = 0
    CUBIC = 1
    TRANSVERSE_ISOTROPIC = 2
    ORTHOTROPIC = 3
    HYPER_INVARIANTS = 4


class IsoConvention(IntEnum):
    """Parameterization of the two isotropic elastic constants (C1, C2).

    Mirrors the convention strings of ``sim.L_iso`` (C++ ``IsoConv``):

    ========= ============== ========== ==========
    value     string         C1         C2
    ========= ============== ========== ==========
    ENU       ``Enu``        E          nu
    NUE       ``nuE``        nu         E
    KMU       ``Kmu``        K          mu (= G)
    MUK       ``muK``        mu (= G)   K
    LAMBDAMU  ``lambdamu``   lambda     mu (= G)
    MULAMBDA  ``mulambda``   mu (= G)   lambda
    ========= ============== ========== ==========
    """
    ENU = 0
    NUE = 1
    KMU = 2
    MUK = 3
    LAMBDAMU = 4
    MULAMBDA = 5


class CubicConvention(IntEnum):
    """Parameterization of the three cubic elastic constants (C1, C2, C3).

    ENUG: C1 = E, C2 = nu, C3 = G. CII: C1 = C11, C2 = C12, C3 = C44.
    """
    ENUG = 0
    CII = 1


class IsotransConvention(IntEnum):
    """Parameterization of the transversely isotropic constants.

    A single parameterization exists (EL, ET, nuTL, nuTT, GLT); the enum
    keeps the props layout uniform and future parameterizations additive.
    """
    ENUG = 0


class OrthoConvention(IntEnum):
    """Parameterization of the nine orthotropic elastic constants (C1..C9).

    ENUG: (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23).
    CII: (C11, C12, C13, C22, C23, C33, C44, C55, C66).
    """
    ENUG = 0
    CII = 1


# String spellings accepted anywhere a convention is expected — identical to
# the C++ convention strings of L_iso/L_cubic/L_ortho (constitutive.cpp),
# including their aliases ("KG" for "Kmu", etc.).
_ISO_CONV_FROM_STR = {
    "Enu": IsoConvention.ENU,
    "nuE": IsoConvention.NUE,
    "Kmu": IsoConvention.KMU, "KG": IsoConvention.KMU,
    "muK": IsoConvention.MUK, "GK": IsoConvention.MUK,
    "lambdamu": IsoConvention.LAMBDAMU, "lambdaG": IsoConvention.LAMBDAMU,
    "mulambda": IsoConvention.MULAMBDA, "Glambda": IsoConvention.MULAMBDA,
}
_CUBIC_CONV_FROM_STR = {
    "EnuG": CubicConvention.ENUG,
    "Cii": CubicConvention.CII,
}
_ISOTRANS_CONV_FROM_STR = {
    "EnuG": IsotransConvention.ENUG,
}
_ORTHO_CONV_FROM_STR = {
    "EnuG": OrthoConvention.ENUG,
    "Cii": OrthoConvention.CII,
}


def _as_convention(value, enum_cls, str_map):
    """Normalize a convention given as enum, int code, or string."""
    if isinstance(value, enum_cls):
        return value
    if isinstance(value, str):
        try:
            return str_map[value]
        except KeyError:
            raise ValueError(
                f"unknown {enum_cls.__name__} string {value!r} "
                f"(valid: {sorted(str_map)})") from None
    return enum_cls(value)  # int code; raises ValueError if out of range


class YieldType(IntEnum):
    """Type of yield criterion."""
    VON_MISES = 0
    TRESCA = 1
    DRUCKER = 2
    HILL = 3
    DFA = 4
    ANISOTROPIC = 5


class IsoHardType(IntEnum):
    """Type of isotropic hardening."""
    NONE = 0
    LINEAR = 1
    POWER_LAW = 2
    VOCE = 3
    COMBINED_VOCE = 4


class KinHardType(IntEnum):
    """Type of kinematic hardening."""
    NONE = 0
    PRAGER = 1
    ARMSTRONG_FREDERICK = 2
    CHABOCHE = 3


class DamageType(IntEnum):
    """Type of damage evolution law."""
    LINEAR = 0
    EXPONENTIAL = 1
    POWER_LAW = 2
    WEIBULL = 3


class MechanismType(IntEnum):
    """Type of strain mechanism."""
    PLASTICITY = 0
    VISCOELASTICITY = 1
    DAMAGE = 2


# ============================================================================
# Elasticity configurations
# ============================================================================

@dataclass(frozen=True)
class IsotropicElasticity:
    """Isotropic elasticity.

    The two elastic constants are ordinal slots whose meaning is set by
    ``convention`` (default ``"Enu"``: C1 = E, C2 = nu) — see
    :class:`IsoConvention`. The interpretation is done by the C++ builders
    (``L_iso``); no conversion happens in Python.

    Parameters
    ----------
    C1 : float
        First elastic constant (E, nu, K, mu or lambda per the convention).
    C2 : float
        Second elastic constant.
    alpha : float
        Coefficient of thermal expansion.
    convention : IsoConvention or str
        Parameterization of (C1, C2). Accepts the enum, its integer code,
        or the ``L_iso`` string (``"Enu"``, ``"Kmu"``, ``"lambdamu"``, ...).

    Examples
    --------
    >>> IsotropicElasticity(C1=210000., C2=0.3, alpha=1.2e-5)   # E, nu
    >>> IsotropicElasticity(C1=175000., C2=80769., convention="Kmu")  # K, mu
    """
    C1: float
    C2: float
    alpha: float = 0.0
    convention: IsoConvention = IsoConvention.ENU

    def __post_init__(self):
        object.__setattr__(self, "convention",
                           _as_convention(self.convention, IsoConvention,
                                          _ISO_CONV_FROM_STR))

    @property
    def elasticity_type(self) -> ElasticityType:
        return ElasticityType.ISOTROPIC

    def to_props(self) -> List[float]:
        """Return the props values for this elasticity."""
        return [float(self.convention), self.C1, self.C2, self.alpha]

    @property
    def nprops(self) -> int:
        return 4


@dataclass(frozen=True)
class CubicElasticity:
    """Cubic elasticity (3 independent constants).

    The constants are ordinal slots whose meaning is set by ``convention``
    (default ``"EnuG"``: C1 = E, C2 = nu, C3 = G; ``"Cii"``: C1 = C11,
    C2 = C12, C3 = C44) — see :class:`CubicConvention`.

    For cubic symmetry G is independent from E and nu (Zener ratio
    A = 2*G*(1+nu)/E != 1 in general).

    Parameters
    ----------
    C1, C2, C3 : float
        The three elastic constants, interpreted per the convention.
    alpha : float
        Coefficient of thermal expansion.
    convention : CubicConvention or str
        Parameterization of (C1, C2, C3): ``"EnuG"`` or ``"Cii"``.

    Examples
    --------
    >>> CubicElasticity(C1=185000., C2=0.28, C3=39700.)  # E, nu, G
    >>> CubicElasticity(C1=185000., C2=158000., C3=39700., convention="Cii")
    """
    C1: float
    C2: float
    C3: float
    alpha: float = 0.0
    convention: CubicConvention = CubicConvention.ENUG

    def __post_init__(self):
        object.__setattr__(self, "convention",
                           _as_convention(self.convention, CubicConvention,
                                          _CUBIC_CONV_FROM_STR))

    @property
    def elasticity_type(self) -> ElasticityType:
        return ElasticityType.CUBIC

    def to_props(self) -> List[float]:
        return [float(self.convention), self.C1, self.C2, self.C3, self.alpha]

    @property
    def nprops(self) -> int:
        return 5


@dataclass(frozen=True)
class TransverseIsotropicElasticity:
    """Transversely isotropic elasticity.

    Parameters
    ----------
    EL : float
        Longitudinal Young's modulus.
    ET : float
        Transverse Young's modulus.
    nuTL : float
        Poisson's ratio (transverse-longitudinal).
    nuTT : float
        Poisson's ratio (transverse-transverse).
    GLT : float
        Shear modulus.
    alpha_L : float
        Longitudinal CTE.
    alpha_T : float
        Transverse CTE.
    axis : int
        Axis of symmetry (1=x, 2=y, 3=z).
    convention : IsotransConvention or str
        Parameterization of the constants. A single one exists (``"EnuG"``);
        the field keeps the props layout uniform across elasticity types.
    """
    EL: float
    ET: float
    nuTL: float
    nuTT: float
    GLT: float
    alpha_L: float = 0.0
    alpha_T: float = 0.0
    axis: int = 3
    convention: IsotransConvention = IsotransConvention.ENUG

    def __post_init__(self):
        object.__setattr__(self, "convention",
                           _as_convention(self.convention, IsotransConvention,
                                          _ISOTRANS_CONV_FROM_STR))

    @property
    def elasticity_type(self) -> ElasticityType:
        return ElasticityType.TRANSVERSE_ISOTROPIC

    def to_props(self) -> List[float]:
        return [float(self.convention),
                self.EL, self.ET, self.nuTL, self.nuTT,
                self.GLT, self.alpha_L, self.alpha_T, float(self.axis)]

    @property
    def nprops(self) -> int:
        return 9


@dataclass(frozen=True)
class OrthotropicElasticity:
    """Orthotropic elasticity (9 independent elastic constants).

    The constants are ordinal slots whose meaning is set by ``convention``
    (see :class:`OrthoConvention`):

    - ``"EnuG"`` (default): C1..C9 = E1, E2, E3, nu12, nu13, nu23,
      G12, G13, G23.
    - ``"Cii"``: C1..C9 = C11, C12, C13, C22, C23, C33, C44, C55, C66.

    Parameters
    ----------
    C1, ..., C9 : float
        The nine elastic constants, interpreted per the convention.
    alpha1, alpha2, alpha3 : float
        Coefficients of thermal expansion.
    convention : OrthoConvention or str
        Parameterization of the nine constants: ``"EnuG"`` or ``"Cii"``.
    """
    C1: float
    C2: float
    C3: float
    C4: float
    C5: float
    C6: float
    C7: float
    C8: float
    C9: float
    alpha1: float = 0.0
    alpha2: float = 0.0
    alpha3: float = 0.0
    convention: OrthoConvention = OrthoConvention.ENUG

    def __post_init__(self):
        object.__setattr__(self, "convention",
                           _as_convention(self.convention, OrthoConvention,
                                          _ORTHO_CONV_FROM_STR))

    @property
    def elasticity_type(self) -> ElasticityType:
        return ElasticityType.ORTHOTROPIC

    def to_props(self) -> List[float]:
        return [float(self.convention),
                self.C1, self.C2, self.C3,
                self.C4, self.C5, self.C6,
                self.C7, self.C8, self.C9,
                self.alpha1, self.alpha2, self.alpha3]

    @property
    def nprops(self) -> int:
        return 13


# ---------------------------------------------------------------------------
# Hyperelastic elasticity blocks (isochoric invariants)
# ---------------------------------------------------------------------------
# Unlike the four linear symmetries, these are not a constant stiffness: the
# C++ block integrates the potential at the elastic strain it is handed. Under
# finite strain (NLGEOM) that is the elastic LOGARITHMIC strain, bridged to the
# potential by b_el = exp(2 eps_el) — so the composition is the
# logarithmic-strain-space form of multiplicative finite strain, exact for an
# isotropic material, and the mechanisms keep riding additively on ln V.
#
# Every potential shares its volumetric term U(J), chosen by ``volumetric``:
# "log" for kappa (J ln J - J + 1) (default) or "quadratic" for kappa/2 (J - 1)^2.
# Both have U''(1) = kappa, so ``kappa`` is the ground-state bulk modulus
# throughout. The ground-state stiffness (what the mechanisms take as their
# reference) is computed by the C++ block from the potential itself, not
# declared here.


class HyperPotential(IntEnum):
    """Isochoric-invariant potential (mirrors C++ ``HyperPotential``, hyperelastic.hpp)."""
    NEOHC = 0
    MOORI = 1
    YEOHH = 2
    ISHAH = 3
    GETHH = 4
    SWANH = 5
    HOLZA = 6
    MUSCL = 7


class MuscleFibreLaw(IntEnum):
    """Along-fibre force law of the ``MUSCL`` potential (mirrors C++ ``MuscleFibreLaw``).

    The four laws share everything but the fibre force
    :math:`f_d(\\bar{\\lambda}) = \\partial W / \\partial \\bar{\\lambda}`, so the choice
    travels as a code in the props and changes the *interpretation* of the fibre
    parameters, never their number.

    ``NONE``
        No fibre term. Activation acts only through the matrix multiplier
        :math:`s(a)`, which is Nazari et al. (2010, 2011): the contractile force is
        supplied by something else (1-D cable elements, in their ANSYS model).
    ``SIMPLE``
        :math:`f_d = a\\,\\sigma_{max}`, a constant active stress. ArtiSynth
        ``SimpleForceMuscle``.
    ``GENERIC``
        :math:`f_d = a\\,\\sigma_{max} + P_1 (e^{P_2(\\bar{\\lambda}-1)} - 1)/\\bar{\\lambda}`.
        ArtiSynth ``GenericMuscle``. **Here** :math:`P_1` **is a stress**, there is no
        optimal length, and the active term does not depend on stretch.
    ``BLEMKER``
        :math:`f_d = \\sigma_{max}(a f_a(\\hat{\\lambda}) + f_p(\\hat{\\lambda}))/\\lambda_{opt}`
        with the Hill force-length curve. Blemker et al. (2005); ArtiSynth
        ``BlemkerMuscle``; FEBio. **Here** :math:`P_1` **is dimensionless.**

    .. warning::

       ArtiSynth ships ``P1 = 0.05`` as the default for both ``GenericMuscle`` and
       ``BlemkerMuscle`` although it means different things in each. Use
       :meth:`MuscleElasticity.generic` and :meth:`MuscleElasticity.blemker`, which
       carry each law's own published values, rather than copying a number across.
    """
    NONE = 0
    SIMPLE = 1
    GENERIC = 2
    BLEMKER = 3


_FIBRE_LAW_NAMES = {"none": MuscleFibreLaw.NONE, "simple": MuscleFibreLaw.SIMPLE,
                    "generic": MuscleFibreLaw.GENERIC, "blemker": MuscleFibreLaw.BLEMKER}


class VolumetricPotential(IntEnum):
    """Volumetric term U(J) (mirrors C++ ``VolumetricPotential``, hyperelastic.hpp).

    Selected by the trailing prop of every hyperelastic law (NEOHC, MOORI, YEOHH,
    ISHAH, GETHH, SWANH, OGDEN and the modular blocks): absent or 0 for ``LOG_J``,
    1 for ``QUADRATIC``.
    """
    LOG_J = 0       # kappa (J ln J - J + 1)
    QUADRATIC = 1   # kappa / 2 (J - 1)^2


_VOLUMETRIC_NAMES = {"log": VolumetricPotential.LOG_J,
                     "quadratic": VolumetricPotential.QUADRATIC}


@dataclass(frozen=True)
class _HyperInvariantsElasticity:
    r"""Common serialization of an isochoric-invariant potential.

    Props layout: ``[potential, n_params, params..., volumetric, alpha]``, the
    volumetric selector counted in ``n_params``. The potential comes from the
    concrete subclass and the parameters are its fields in declaration order,
    which is what gives each model its own named arguments. They have no
    default; ``volumetric`` and ``alpha`` are keyword-only, so positional
    arguments always fill the potential's parameters.

    ``volumetric`` selects U(J): ``"log"`` for
    :math:`\kappa (J \ln J - J + 1)` (default) or ``"quadratic"`` for
    :math:`\frac{\kappa}{2} (J - 1)^2`.
    """
    potential: ClassVar[HyperPotential]
    volumetric: str = field(default="log", kw_only=True)
    alpha: float = field(default=0.0, kw_only=True)

    @property
    def elasticity_type(self) -> ElasticityType:
        return ElasticityType.HYPER_INVARIANTS

    def potential_params(self) -> List[float]:
        return [getattr(self, f.name) for f in fields(self)
                if f.name not in ("alpha", "volumetric")]

    def __post_init__(self):
        self.volumetric_potential   # a bad selector fails at construction, not at to_props()

    @property
    def volumetric_potential(self) -> VolumetricPotential:
        if isinstance(self.volumetric, str):
            if self.volumetric in _VOLUMETRIC_NAMES:
                return _VOLUMETRIC_NAMES[self.volumetric]
        elif self.volumetric in tuple(VolumetricPotential):   # the enum or its 0/1 code
            return VolumetricPotential(self.volumetric)
        raise ValueError(f"volumetric must be 'log' or 'quadratic' (or the VolumetricPotential "
                         f"code 0/1), got {self.volumetric!r}")

    def to_props(self) -> List[float]:
        """Return the props values for this elasticity."""
        params = list(self.potential_params()) + [float(self.volumetric_potential)]
        return [float(self.potential), float(len(params))] + params + [self.alpha]

    @property
    def nprops(self) -> int:
        return len(self.to_props())


@dataclass(frozen=True)
class NeoHookeanElasticity(_HyperInvariantsElasticity):
    r"""Compressible neo-Hookean potential (the ``NEOHC`` UMAT's).

    :math:`W = \frac{\mu}{2}(\bar{I}_1 - 3) + U(J)`

    Parameters
    ----------
    mu : float
        Ground-state shear modulus.
    kappa : float
        Ground-state bulk modulus, :math:`U''(1)`.
    volumetric : str
        ``"log"`` (:math:`U = \kappa (J \ln J - J + 1)`, default) or
        ``"quadratic"`` (:math:`U = \frac{\kappa}{2} (J - 1)^2`); keyword-only.
    alpha : float
        Coefficient of thermal expansion.
    """
    potential = HyperPotential.NEOHC
    mu: float
    kappa: float


@dataclass(frozen=True)
class MooneyRivlinElasticity(_HyperInvariantsElasticity):
    r"""Mooney-Rivlin potential (the ``MOORI`` UMAT's).

    :math:`W = C_{10}(\bar{I}_1 - 3) + C_{01}(\bar{I}_2 - 3) + U(J)`
    """
    potential = HyperPotential.MOORI
    C10: float
    C01: float
    kappa: float


@dataclass(frozen=True)
class YeohElasticity(_HyperInvariantsElasticity):
    r"""Yeoh potential (the ``YEOHH`` UMAT's).

    :math:`W = C_{10}(\bar{I}_1 - 3) + C_{20}(\bar{I}_1 - 3)^2
    + C_{30}(\bar{I}_1 - 3)^3 + U(J)`

    The ground-state shear modulus is :math:`\mu = 2 C_{10}`; the higher-order
    terms carry the upturn at large stretch that a neo-Hookean cannot fit.

    Examples
    --------
    >>> YeohElasticity(C10=0.30, C20=-0.010, C30=0.0005, kappa=1000.)
    """
    potential = HyperPotential.YEOHH
    C10: float
    C20: float
    C30: float
    kappa: float


@dataclass(frozen=True)
class IsiharaElasticity(_HyperInvariantsElasticity):
    r"""Isihara potential (the ``ISHAH`` UMAT's).

    :math:`W = C_{10}(\bar{I}_1 - 3) + C_{20}(\bar{I}_1 - 3)^2 + C_{01}(\bar{I}_2 - 3)
    + U(J)`
    """
    potential = HyperPotential.ISHAH
    C10: float
    C20: float
    C01: float
    kappa: float


@dataclass(frozen=True)
class GentThomasElasticity(_HyperInvariantsElasticity):
    r"""Gent-Thomas potential (the ``GETHH`` UMAT's).

    :math:`W = c_1(\bar{I}_1 - 3) + c_2 \ln(\bar{I}_2 / 3) + U(J)`
    """
    potential = HyperPotential.GETHH
    c1: float
    c2: float
    kappa: float


@dataclass(frozen=True)
class SwansonElasticity(_HyperInvariantsElasticity):
    """Swanson potential (the ``SWANH`` UMAT's), N terms.

    Parameters
    ----------
    terms : sequence of (A, B, alpha, beta) tuples
        One tuple per term; all four values are required.
    kappa : float
        Bulk compressibility.
    """
    potential = HyperPotential.SWANH
    terms: Tuple[Tuple[float, float, float, float], ...]
    kappa: float

    def __post_init__(self):
        super().__post_init__()
        for i, term in enumerate(self.terms):
            if len(term) != 4:
                raise TypeError(
                    f"SwansonElasticity.terms[{i}] must be (A, B, alpha, beta); "
                    f"got {len(term)} values.")

    def potential_params(self) -> List[float]:
        params: List[float] = [float(len(self.terms)), self.kappa]
        for A, B, a, b in self.terms:
            params.extend([A, B, a, b])
        return params


@dataclass(frozen=True)
class HolzapfelElasticity(_HyperInvariantsElasticity):
    r"""Gasser-Ogden-Holzapfel potential (the ``HOLZA`` UMAT's).

    An isotropic neo-Hookean ground matrix reinforced by one or two families of
    dispersed collagen fibres:

    :math:`W = C_{10}(\bar{I}_1 - 3)
    + \sum_i \frac{k_1}{2 k_2}\left[\exp\left(k_2 (\bar{I}^*_{4,i} - 1)^2\right) - 1\right]
    + U(J)`

    where :math:`\bar{I}^*_{4,i} = \kappa_d \bar{I}_1 + (1 - 3\kappa_d)\bar{I}_{4,i}`
    and :math:`\bar{I}_{4,i} = \mathbf{a}_{0,i} \cdot \bar{\mathbf{C}}\, \mathbf{a}_{0,i}`.
    The fibre term is inactive where :math:`\bar{I}^*_{4,i} < 1`. That is a condition
    on the *generalized* invariant, not on fibre compression: once
    :math:`\kappa_d > 0` a fibre with :math:`\bar{I}_{4,i} < 1` can still be active,
    through the :math:`\kappa_d \bar{I}_1` term.

    Parameters
    ----------
    C10 : float
        Neo-Hookean ground matrix; the matrix shear modulus is :math:`\mu = 2 C_{10}`.
    k1 : float
        Fibre stiffness, in stress units (MPa).
    k2 : float
        Fibre stiffening exponent, dimensionless.
    kappa_d : float
        Fibre dispersion, :math:`\kappa_d \in [0, 1/3]`. 0 gives perfectly aligned
        fibres (the Holzapfel-Gasser-Ogden 2000 model), 1/3 an isotropic
        distribution, for which the response no longer depends on ``fibres``.
    fibres : Rotation or array-like
        The fibre directions in the LOCAL material frame, as a (possibly batched)
        :class:`simcoon.Rotation` -- one entry per family -- applied to
        :math:`\mathbf{e}_1`, or directly as an ``(n, 3)`` array of components.
        Going through a ``Rotation`` keeps Euler angles, and hence gimbal lock, out
        of the path to the kernel. The solver's material orientation places the
        local frame globally, as it does for ELIST/ELORT.
    kappa : float
        Ground-state bulk modulus, :math:`U''(1)`.
    volumetric : str
        ``"log"`` (default) or ``"quadratic"``; keyword-only.
    alpha : float
        Coefficient of thermal expansion.

    Examples
    --------
    >>> import simcoon as sim
    >>> HolzapfelElasticity(
    ...     C10=0.0354, k1=0.0107, k2=7.48, kappa_d=0.0,
    ...     fibres=sim.Rotation.from_euler('zxz', [[0, 0, 40], [0, 0, -40]], degrees=True),
    ...     kappa=1000.)               # doctest: +ELLIPSIS
    HolzapfelElasticity(...)

    Notes
    -----
    Composed with :class:`Damage` -- the classic anisotropic tissue with
    softening -- this block is exact: damage subtracts no inelastic strain (it
    scales the stiffness instead), so the elastic stretch is still the total
    one, and its driving force uses the current anisotropic tangent.

    .. warning::

       Composed with :class:`Plasticity` or :class:`Viscoelasticity`, the fibre
       convection is **approximate**. The block carries the reference directions
       and pushes them forward with the *elastic* stretch, so the inelastic
       strain does not reorient the fibres. The composition is well posed and
       converges -- the return mapping gets the anisotropic tangent and its
       consistency condition holds exactly -- but it is only as good as that
       assumption, degrading as the inelastic strain reorients the fibres.
       Nothing rejects it at run time; it is a modelling choice.

       For :class:`Viscoelasticity` there is a second, separate caveat: every
       Prony branch is built as an *isotropic* ``L_iso(E_i, nu_i)``, so the
       viscous response carries none of the fibre anisotropy while the
       equilibrium response does. This is a property of the viscoelastic
       mechanism's ``(E_i, nu_i)`` parameterization, not of this block -- a
       branch cannot follow an orthotropic or transversely isotropic elasticity
       either.

    A single scalar damage variable also degrades matrix and fibres at the same
    rate; the Holzapfel damage literature uses separate variables per term.
    """
    potential = HyperPotential.HOLZA
    C10: float
    k1: float
    k2: float
    kappa_d: float
    fibres: Orientation
    kappa: float

    def __post_init__(self):
        super().__post_init__()
        if not 0.0 <= float(self.kappa_d) <= 1.0/3.0:
            raise ValueError(f"HolzapfelElasticity: kappa_d must lie in [0, 1/3], "
                             f"got {self.kappa_d!r}")
        self.directions     # a bad fibre spec fails at construction, not at to_props()

    @property
    def directions(self) -> NDArray[np.float64]:
        """The unit fibre directions as a ``(3, n_fam)`` array, one per column."""
        return as_direction(self.fibres)

    # `fibres` may hold a Rotation or an array, neither of which the dataclass-generated
    # __eq__/__hash__ can handle (an array comparison is ambiguous, and an array is
    # unhashable). Compare and hash the resolved directions instead, so this block behaves
    # like every other Elasticity. Same reason micromechanics.py carries _dataclass_eq.
    def _key(self):
        return (float(self.C10), float(self.k1), float(self.k2), float(self.kappa_d),
                float(self.kappa), self.volumetric, float(self.alpha),
                tuple(self.directions.ravel()))

    def __eq__(self, other) -> bool:
        if not isinstance(other, HolzapfelElasticity):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self) -> int:
        return hash(self._key())

    def potential_params(self) -> List[float]:
        a0 = self.directions
        return ([float(self.C10), float(self.k1), float(self.k2), float(self.kappa_d),
                 float(a0.shape[1])]
                + a0.T.ravel().tolist()
                + [float(self.kappa)])


@dataclass(frozen=True)
class MuscleElasticity(_HyperInvariantsElasticity):
    r"""Activated skeletal muscle (the ``MUSCL`` UMAT's).

    A 5-parameter Mooney-Rivlin ground matrix whose stiffness rises with the
    activation, plus an along-fibre force law chosen by ``fibre_law``:

    :math:`W = s(a) \left[ C_{10}(\bar{I}_1 - 3) + C_{01}(\bar{I}_2 - 3)
    + C_{20}(\bar{I}_1 - 3)^2 + C_{11}(\bar{I}_1 - 3)(\bar{I}_2 - 3)
    + C_{02}(\bar{I}_2 - 3)^2 \right] + \sum_i \Phi(\bar{\lambda}_i; a) + s(a) U(J)`

    with :math:`s(a) = 1 + (s_{max} - 1) a` and
    :math:`\bar{\lambda} = \sqrt{\bar{I}^*_4}` the isochoric fibre stretch. See
    :class:`MuscleFibreLaw` for :math:`\Phi`. The matrix form is Nazari et al.'s
    Eq. (1) and is a superset of neo-Hookean, Mooney-Rivlin and second-order Yeoh.

    The named constructors :meth:`nazari`, :meth:`blemker`, :meth:`generic`,
    :meth:`simple_force` and :meth:`face` carry each source's published parameters
    and are the intended entry points; the raw constructor is keyword-only, because
    seventeen positional numbers would be unreadable.

    Parameters
    ----------
    fibre_law : MuscleFibreLaw or str
        ``"none"``, ``"simple"``, ``"generic"`` or ``"blemker"``.
    C10, C01, C20, C11, C02 : float
        Ground-matrix constants (MPa). Setting only ``C10`` gives neo-Hookean,
        ``C10, C01`` Mooney-Rivlin, ``C10, C20`` the second-order Yeoh that
        Nazari uses.
    s_max : float
        Matrix stiffness multiplier at full activation, :math:`s(1) \ge 1`.
        Nazari's value is 10. Must be 1 unless ``fibre_law`` is ``NONE``.
    activation : float
        The activation :math:`a \in [0, 1]`. **This is an input that changes every
        increment**, unlike every other parameter here: rewrite it and pass the
        props again. See :attr:`activation_index`.
    sigma_max : float
        Maximum isometric fibre stress (MPa).
    lambda_opt : float
        Optimal fibre stretch. ``BLEMKER`` only; ``GENERIC`` has no optimal length.
    lambda_star : float
        Stretch at which the passive fibre curve switches from exponential to
        linear. The switch is C1 by construction.
    P1 : float
        Passive fibre coefficient. **A stress for** ``GENERIC``, **dimensionless
        for** ``BLEMKER`` -- see :class:`MuscleFibreLaw`.
    P2 : float
        Uncrimping factor (dimensionless).
    zero_below_opt : bool
        Whether the passive fibre force is zero below the optimal length. True in
        ``BlemkerMuscle``, **false** in ``GenericMuscle``, which therefore produces
        a negative fibre stress in compression. True is the physically safe choice.
    kappa_d : float
        Fibre dispersion :math:`\kappa_d \in [0, 1/3]`; 0 (perfectly aligned) for
        muscle. Shared with :class:`HolzapfelElasticity`.
    fibres : Rotation or array-like, optional
        Fibre directions in the LOCAL material frame, as for
        :class:`HolzapfelElasticity`. Required unless ``fibre_law`` is ``NONE``.
    kappa : float
        Ground-state bulk modulus, :math:`U''(1)`. It is scaled by :math:`s(a)`
        along with the matrix constants, which is what keeps Poisson's ratio fixed
        as the muscle stiffens.

    Examples
    --------
    >>> import simcoon as sim
    >>> law = sim.modular.MuscleElasticity.blemker(
    ...     fibres=sim.Rotation.from_euler('zxz', [[0, 0, 0]], degrees=True),
    ...     kappa=1000., activation=0.5)          # doctest: +ELLIPSIS
    >>> law.fibre_law
    <MuscleFibreLaw.BLEMKER: 3>

    Notes
    -----
    The activation is a *prescribed* input, not a governed internal variable, so the
    material point is thermodynamically open. At frozen activation the law is
    hyperelastic and :math:`W_{m,d} = W_{m,ir} = 0`; along a path where the
    activation varies, :math:`W_{m,r}` is the mechanical work residue rather than
    the stored energy and **can go negative**.

    ``s_max > 1`` and an active fibre law are rejected together: scaling the passive
    stiffness with activation is Nazari's surrogate for the stress-stiffening a real
    contractile fibre produces, and his own later 3-D muscle element drops it for
    exactly that reason. Composing them double-counts.

    The fibre-convection caveats of :class:`HolzapfelElasticity` apply verbatim.
    """
    potential = HyperPotential.MUSCL

    #: Index of the activation among this potential's own parameters -- that is, its
    #: index in the props of the standalone ``MUSCL`` UMAT. In a
    #: :class:`ModularMaterial` the block is preceded by the elasticity type, the
    #: potential code and the parameter count, so use
    #: :attr:`ModularMaterial.activation_index` there rather than adding 3 by hand.
    activation_index: ClassVar[int] = 7

    fibre_law: Union[MuscleFibreLaw, str] = field(default=MuscleFibreLaw.BLEMKER, kw_only=True)
    C10: float = field(default=0.0, kw_only=True)
    C01: float = field(default=0.0, kw_only=True)
    C20: float = field(default=0.0, kw_only=True)
    C11: float = field(default=0.0, kw_only=True)
    C02: float = field(default=0.0, kw_only=True)
    s_max: float = field(default=1.0, kw_only=True)
    activation: float = field(default=0.0, kw_only=True)
    sigma_max: float = field(default=0.0, kw_only=True)
    lambda_opt: float = field(default=1.0, kw_only=True)
    lambda_star: float = field(default=1.4, kw_only=True)
    P1: float = field(default=0.0, kw_only=True)
    P2: float = field(default=6.6, kw_only=True)
    zero_below_opt: bool = field(default=True, kw_only=True)
    kappa_d: float = field(default=0.0, kw_only=True)
    fibres: Optional[Orientation] = field(default=None, kw_only=True)
    kappa: float = field(default=0.0, kw_only=True)

    # ---- named constructors: one per source, with its own published parameters ----

    @classmethod
    def nazari(cls, **kwargs) -> "MuscleElasticity":
        """Nazari et al. (2010, 2011): activation raises the matrix stiffness x1 -> x10.

        No fibre term -- the contractile force comes from elsewhere. Parameters from
        Nazari et al. (2010) Table 1 (:math:`d = 0.8` MPa\\ :sup:`-1`, so
        :math:`\\kappa = 2/d`), the same numbers ArtiSynth's ``BadinFaceDemo`` carries.
        """
        return cls(**{"fibre_law": MuscleFibreLaw.NONE, "C10": 2.5e-3, "C20": 1.175e-3,
                      "s_max": 10.0, "kappa": 2.5, "volumetric": "quadratic", **kwargs})

    @classmethod
    def blemker(cls, **kwargs) -> "MuscleElasticity":
        """Blemker et al. (2005) fibre law, with ArtiSynth ``BlemkerMuscle``'s defaults.

        ``sigma_max`` is 0.3 MPa (3e5 Pa) and ``P1`` is dimensionless.
        """
        return cls(**{"fibre_law": MuscleFibreLaw.BLEMKER, "sigma_max": 0.3,
                      "lambda_opt": 1.0, "lambda_star": 1.4, "P1": 0.05, "P2": 6.6,
                      "zero_below_opt": True, **kwargs})

    @classmethod
    def generic(cls, **kwargs) -> "MuscleElasticity":
        """ArtiSynth ``GenericMuscle``: constant active stress, exponential passive fibre.

        ``sigma_max`` is 0.03 MPa (3e4 Pa) and ``P1`` is a **stress**. ``zero_below_opt``
        is False, matching the reference implementation -- which means a negative passive
        fibre stress in compression.
        """
        return cls(**{"fibre_law": MuscleFibreLaw.GENERIC, "sigma_max": 0.03,
                      "lambda_opt": 1.0, "lambda_star": 1.4, "P1": 0.05, "P2": 6.6,
                      "zero_below_opt": False, **kwargs})

    @classmethod
    def simple_force(cls, **kwargs) -> "MuscleElasticity":
        """ArtiSynth ``SimpleForceMuscle``: active stress only, no passive fibre term."""
        return cls(**{"fibre_law": MuscleFibreLaw.SIMPLE, "sigma_max": 0.03, **kwargs})

    @classmethod
    def face(cls, **kwargs) -> "MuscleElasticity":
        """The ArtiSynth face models' setting: :meth:`generic` at 0.1 MPa (100 kPa).

        ``BadinFemMuscleFaceDemo`` and ``RefFemMuscleFaceDemo`` both do
        ``setMaxStress(100000)`` on a ``GenericMuscle``, over a Mooney-Rivlin matrix
        with Nazari's constants -- which this constructor also supplies.
        """
        return cls.generic(**{"sigma_max": 0.1, "C10": 2.5e-3, "C20": 1.175e-3,
                              "kappa": 2.5, **kwargs})

    # ---- validation and serialization ----

    @property
    def law(self) -> MuscleFibreLaw:
        """``fibre_law`` resolved to the enum."""
        if isinstance(self.fibre_law, str):
            if self.fibre_law in _FIBRE_LAW_NAMES:
                return _FIBRE_LAW_NAMES[self.fibre_law]
        elif self.fibre_law in tuple(MuscleFibreLaw):
            return MuscleFibreLaw(self.fibre_law)
        raise ValueError(f"fibre_law must be one of {sorted(_FIBRE_LAW_NAMES)} (or the "
                         f"MuscleFibreLaw code 0-3), got {self.fibre_law!r}")

    def __post_init__(self):
        super().__post_init__()
        law = self.law
        if not 0.0 <= float(self.activation) <= 1.0:
            raise ValueError(f"MuscleElasticity: activation must lie in [0, 1], "
                             f"got {self.activation!r}")
        if float(self.s_max) < 1.0:
            raise ValueError(f"MuscleElasticity: s_max must be >= 1, got {self.s_max!r}")
        if float(self.s_max) > 1.0 and law is not MuscleFibreLaw.NONE:
            raise ValueError(
                "MuscleElasticity: s_max > 1 (activation-scaled matrix stiffness) and an "
                "active fibre law represent the same physics, so composing them "
                "double-counts it. Use s_max=1 with a fibre law, or fibre_law='none' "
                "with s_max > 1 (Nazari's own model).")
        if not 0.0 <= float(self.kappa_d) <= 1.0/3.0:
            raise ValueError(f"MuscleElasticity: kappa_d must lie in [0, 1/3], "
                             f"got {self.kappa_d!r}")
        if float(self.lambda_opt) <= 0.0:
            raise ValueError(f"MuscleElasticity: lambda_opt must be > 0, "
                             f"got {self.lambda_opt!r}")
        if law not in (MuscleFibreLaw.NONE, MuscleFibreLaw.SIMPLE):
            if float(self.lambda_star) <= float(self.lambda_opt):
                raise ValueError(f"MuscleElasticity: lambda_star must exceed lambda_opt, "
                                 f"got {self.lambda_star!r} <= {self.lambda_opt!r}")
            if float(self.P2) <= 0.0:
                raise ValueError(f"MuscleElasticity: P2 must be > 0, got {self.P2!r}")
        if law is not MuscleFibreLaw.NONE and self.fibres is None:
            raise ValueError(f"MuscleElasticity: fibre_law {law.name} needs a fibre "
                             f"direction; pass fibres=")
        self.directions     # a bad fibre spec fails at construction, not at to_props()

    @property
    def directions(self) -> NDArray[np.float64]:
        """The unit fibre directions as a ``(3, n_fam)`` array, one per column.

        Empty, ``(3, 0)``, when no fibres are declared -- legal only for
        ``fibre_law='none'``, whose response has no fibre term to place.
        """
        if self.fibres is None:
            return np.zeros((3, 0))
        return as_direction(self.fibres)

    # `fibres` may hold a Rotation or an array, so the dataclass-generated __eq__ and
    # __hash__ cannot be used; compare and hash the resolved directions instead, exactly
    # as HolzapfelElasticity does.
    def _key(self):
        return (int(self.law), float(self.C10), float(self.C01), float(self.C20),
                float(self.C11), float(self.C02), float(self.s_max), float(self.activation),
                float(self.sigma_max), float(self.lambda_opt), float(self.lambda_star),
                float(self.P1), float(self.P2), bool(self.zero_below_opt),
                float(self.kappa_d), float(self.kappa), self.volumetric, float(self.alpha),
                tuple(self.directions.ravel()))

    def __eq__(self, other) -> bool:
        if not isinstance(other, MuscleElasticity):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self) -> int:
        return hash(self._key())

    def potential_params(self) -> List[float]:
        a0 = self.directions
        return ([float(self.law), float(self.C10), float(self.C01), float(self.C20),
                 float(self.C11), float(self.C02), float(self.s_max), float(self.activation),
                 float(self.sigma_max), float(self.lambda_opt), float(self.lambda_star),
                 float(self.P1), float(self.P2), 1.0 if self.zero_below_opt else 0.0,
                 float(self.kappa_d), float(a0.shape[1])]
                + a0.T.ravel().tolist()
                + [float(self.kappa)])


Elasticity = Union[IsotropicElasticity, CubicElasticity,
                   TransverseIsotropicElasticity, OrthotropicElasticity,
                   NeoHookeanElasticity, MooneyRivlinElasticity, YeohElasticity,
                   IsiharaElasticity, GentThomasElasticity, SwansonElasticity,
                   HolzapfelElasticity, MuscleElasticity]

# ============================================================================
# Yield criteria
# ============================================================================

@dataclass(frozen=True)
class VonMisesYield:
    """Von Mises (J2) yield criterion. No additional parameters."""

    @property
    def yield_type(self) -> YieldType:
        return YieldType.VON_MISES

    def to_props(self) -> List[float]:
        return []

    @property
    def nprops(self) -> int:
        return 0


@dataclass(frozen=True)
class TrescaYield:
    """Tresca yield criterion. No additional parameters."""

    @property
    def yield_type(self) -> YieldType:
        return YieldType.TRESCA

    def to_props(self) -> List[float]:
        return []

    @property
    def nprops(self) -> int:
        return 0


@dataclass(frozen=True)
class DruckerYield:
    """Drucker yield criterion (J2/J3-based).

    Parameters
    ----------
    b : float
        J3 influence parameter.
    n : float
        Exponent parameter.
    """
    b: float
    n: float

    @property
    def yield_type(self) -> YieldType:
        return YieldType.DRUCKER

    def to_props(self) -> List[float]:
        return [self.b, self.n]

    @property
    def nprops(self) -> int:
        return 2


@dataclass(frozen=True)
class HillYield:
    """Hill 1948 anisotropic yield criterion.

    Parameters
    ----------
    F, G, H, L, M, N : float
        Hill anisotropy parameters.
    """
    F: float
    G: float
    H: float
    L: float
    M: float
    N: float

    @property
    def yield_type(self) -> YieldType:
        return YieldType.HILL

    def to_props(self) -> List[float]:
        return [self.F, self.G, self.H, self.L, self.M, self.N]

    @property
    def nprops(self) -> int:
        return 6


@dataclass(frozen=True)
class DFAYield:
    """Deshpande-Fleck-Ashby yield criterion.

    Parameters
    ----------
    F, G, H, L, M, N : float
        Anisotropy parameters.
    K : float
        Hydrostatic sensitivity parameter.
    """
    F: float
    G: float
    H: float
    L: float
    M: float
    N: float
    K: float

    @property
    def yield_type(self) -> YieldType:
        return YieldType.DFA

    def to_props(self) -> List[float]:
        return [self.F, self.G, self.H, self.L, self.M, self.N, self.K]

    @property
    def nprops(self) -> int:
        return 7


@dataclass(frozen=True)
class AnisotropicYield:
    """Generic anisotropic yield criterion (9 parameters).

    Parameters
    ----------
    P11, P22, P33 : float
        Normal components of the anisotropy tensor.
    P12, P13, P23 : float
        Off-diagonal components.
    P44, P55, P66 : float
        Shear components.
    """
    P11: float
    P22: float
    P33: float
    P12: float
    P13: float
    P23: float
    P44: float
    P55: float
    P66: float

    @property
    def yield_type(self) -> YieldType:
        return YieldType.ANISOTROPIC

    def to_props(self) -> List[float]:
        return [self.P11, self.P22, self.P33,
                self.P12, self.P13, self.P23,
                self.P44, self.P55, self.P66]

    @property
    def nprops(self) -> int:
        return 9


YieldCriterion = Union[VonMisesYield, TrescaYield, DruckerYield,
                       HillYield, DFAYield, AnisotropicYield]

# ============================================================================
# Isotropic hardening
# ============================================================================

@dataclass(frozen=True)
class NoIsotropicHardening:
    """No isotropic hardening."""

    @property
    def iso_hard_type(self) -> IsoHardType:
        return IsoHardType.NONE

    def to_props(self) -> List[float]:
        return []

    @property
    def nprops(self) -> int:
        return 0

    @property
    def N(self) -> int:
        return 1


@dataclass(frozen=True)
class LinearIsotropicHardening:
    """Linear isotropic hardening: R = H * p.

    Parameters
    ----------
    H : float
        Hardening modulus.
    """
    H: float

    @property
    def iso_hard_type(self) -> IsoHardType:
        return IsoHardType.LINEAR

    def to_props(self) -> List[float]:
        return [self.H]

    @property
    def nprops(self) -> int:
        return 1

    @property
    def N(self) -> int:
        return 1


@dataclass(frozen=True)
class PowerLawHardening:
    """Power-law isotropic hardening: R = k * p^m.

    Parameters
    ----------
    k : float
        Hardening coefficient.
    m : float
        Hardening exponent.
    """
    k: float
    m: float

    @property
    def iso_hard_type(self) -> IsoHardType:
        return IsoHardType.POWER_LAW

    def to_props(self) -> List[float]:
        return [self.k, self.m]

    @property
    def nprops(self) -> int:
        return 2

    @property
    def N(self) -> int:
        return 1


@dataclass(frozen=True)
class VoceHardening:
    """Voce saturation hardening: R = Q * (1 - exp(-b*p)).

    Parameters
    ----------
    Q : float
        Saturation stress.
    b : float
        Hardening rate.
    """
    Q: float
    b: float

    @property
    def iso_hard_type(self) -> IsoHardType:
        return IsoHardType.VOCE

    def to_props(self) -> List[float]:
        return [self.Q, self.b]

    @property
    def nprops(self) -> int:
        return 2

    @property
    def N(self) -> int:
        return 1


@dataclass(frozen=True)
class CombinedVoceHardening:
    """Combined Voce hardening: R = sum_i Q_i * (1 - exp(-b_i*p)).

    Parameters
    ----------
    terms : tuple of (Q, b) tuples
        Each tuple is (saturation stress, hardening rate) for one Voce term.
    """
    terms: Tuple[Tuple[float, float], ...] = ()

    @property
    def iso_hard_type(self) -> IsoHardType:
        return IsoHardType.COMBINED_VOCE

    def to_props(self) -> List[float]:
        props = []
        for Q, b in self.terms:
            props.extend([Q, b])
        return props

    @property
    def nprops(self) -> int:
        return 2 * len(self.terms)

    @property
    def N(self) -> int:
        return len(self.terms)


IsotropicHardening = Union[NoIsotropicHardening, LinearIsotropicHardening,
                           PowerLawHardening, VoceHardening, CombinedVoceHardening]

# ============================================================================
# Kinematic hardening
# ============================================================================

@dataclass(frozen=True)
class NoKinematicHardening:
    """No kinematic hardening."""

    @property
    def kin_hard_type(self) -> KinHardType:
        return KinHardType.NONE

    def to_props(self) -> List[float]:
        return []

    @property
    def nprops(self) -> int:
        return 0

    @property
    def N(self) -> int:
        return 1

    @property
    def num_backstresses(self) -> int:
        return 0


@dataclass(frozen=True)
class PragerHardening:
    """Linear Prager kinematic hardening: X = (2/3)*C*alpha.

    Parameters
    ----------
    C : float
        Kinematic hardening modulus.
    """
    C: float

    def __post_init__(self):
        if not isinstance(self.C, (int, float)):
            raise TypeError(
                "PragerHardening.C must be a scalar. "
                "For multiple Prager-like terms, use ChabocheHardening with D=0 in each term."
            )

    @property
    def kin_hard_type(self) -> KinHardType:
        return KinHardType.PRAGER

    def to_props(self) -> List[float]:
        return [self.C]

    @property
    def nprops(self) -> int:
        return 1

    @property
    def N(self) -> int:
        return 1

    @property
    def num_backstresses(self) -> int:
        return 1


@dataclass(frozen=True)
class ArmstrongFrederickHardening:
    """Armstrong-Frederick kinematic hardening: dX = (2/3)*C*dep - D*X*dp.

    Parameters
    ----------
    C : float
        Hardening parameter.
    D : float
        Dynamic recovery parameter.
    """
    C: float
    D: float

    def __post_init__(self):
        if not (isinstance(self.C, (int, float)) and isinstance(self.D, (int, float))):
            raise TypeError(
                "ArmstrongFrederickHardening(C, D) takes scalars. "
                "For multi-term AF (Chaboche), use "
                "ChabocheHardening(terms=[(C1, D1), (C2, D2), ...])."
            )

    @property
    def kin_hard_type(self) -> KinHardType:
        return KinHardType.ARMSTRONG_FREDERICK

    def to_props(self) -> List[float]:
        return [self.C, self.D]

    @property
    def nprops(self) -> int:
        return 2

    @property
    def N(self) -> int:
        return 1

    @property
    def num_backstresses(self) -> int:
        return 1


@dataclass(frozen=True)
class ChabocheHardening:
    """Chaboche kinematic hardening with multiple backstress terms.

    Each term i follows: dX_i = (2/3)*C_i*dep - D_i*X_i*dp.

    Parameters
    ----------
    terms : tuple of (C, D) tuples
        Each tuple is (hardening parameter, dynamic recovery parameter)
        for one Armstrong-Frederick backstress term.
    """
    terms: Tuple[Tuple[float, float], ...] = ()

    @property
    def kin_hard_type(self) -> KinHardType:
        return KinHardType.CHABOCHE

    def to_props(self) -> List[float]:
        props = []
        for C, D in self.terms:
            props.extend([C, D])
        return props

    @property
    def nprops(self) -> int:
        return 2 * len(self.terms)

    @property
    def N(self) -> int:
        return len(self.terms)

    @property
    def num_backstresses(self) -> int:
        return len(self.terms)


KinematicHardening = Union[NoKinematicHardening, PragerHardening,
                           ArmstrongFrederickHardening, ChabocheHardening]

# ============================================================================
# Strain mechanisms
# ============================================================================

@dataclass(frozen=True)
class Plasticity:
    """Plasticity mechanism combining yield + isotropic + kinematic hardening.

    Parameters
    ----------
    sigma_Y : float
        Initial yield stress.
    yield_criterion : YieldCriterion
        Yield criterion (default: VonMisesYield).
    isotropic_hardening : IsotropicHardening
        Isotropic hardening law (default: NoIsotropicHardening).
    kinematic_hardening : KinematicHardening
        Kinematic hardening law (default: NoKinematicHardening).
    """
    sigma_Y: float
    yield_criterion: YieldCriterion = field(default_factory=VonMisesYield)
    isotropic_hardening: IsotropicHardening = field(default_factory=NoIsotropicHardening)
    kinematic_hardening: KinematicHardening = field(default_factory=NoKinematicHardening)

    @property
    def mechanism_type(self) -> MechanismType:
        return MechanismType.PLASTICITY

    def to_props(self) -> List[float]:
        """Return props values for this mechanism (excluding mechanism_type header)."""
        yc = self.yield_criterion
        ih = self.isotropic_hardening
        kh = self.kinematic_hardening
        # Header: yield_type, iso_type, kin_type, N_iso, N_kin
        header = [
            float(yc.yield_type),
            float(ih.iso_hard_type),
            float(kh.kin_hard_type),
            float(ih.N),
            float(kh.N),
        ]
        # sigma_Y + yield params + iso params + kin params
        params = [self.sigma_Y] + yc.to_props() + ih.to_props() + kh.to_props()
        return header + params

    @property
    def nstatev(self) -> int:
        """Number of state variables for this mechanism: p(1) + EP(6) + backstresses(6*N)."""
        return 7 + 6 * self.kinematic_hardening.num_backstresses


@dataclass(frozen=True)
class Viscoelasticity:
    """Generalized-Maxwell viscoelastic mechanism (Prony branches).

    Each branch is a Maxwell element with its own branch stiffness L_i(E, nu)
    and viscosity tensor H_i(etaB, etaS). The flow rate in branch i is

        dEV_i/dt = invH_i . L_i . (eps - EV_i)

    and the mechanism contribution to the total inelastic strain is

        eps^{in,visco} = sum_i (M_0 . L_i) . EV_i

    where M_0 is the compliance at the reference stiffness, i.e. the stiffness of
    the elasticity block of the material.

    That block is the **instantaneous** (glassy) stiffness, not the long-term one:
    at t = 0 every EV_i is zero and the response is L_0, while at t -> infinity each
    EV_i saturates at the total strain and the response relaxes to

        L_infinity = L_0 - sum_i L_i

    which is the usual Prony series E(t) = E_inf + sum_i E_i exp(-t/tau_i) written
    with E(0) = E_inf + sum_i E_i. Two consequences:

    * the branch moduli must satisfy ``sum_i E_i < E_0``, otherwise the long-term
      stiffness is negative and the stress crosses zero during a hold (measured:
      E_0 = 1, branches 1.0 and 0.5, a hold at 1 % strain relaxes from 0.01 to
      -0.005 MPa). Nothing validates this today;
    * over a **hyperelastic** elasticity block, the potential plays the
      instantaneous role and the relaxed part is subtracted through the
      ground-state compliance M_0, which is linear. The long-term response is then
      not itself a hyperelastic potential, so the classical rubber model (an
      equilibrium hyperelastic spring carrying Maxwell branches) is out of reach
      of this mechanism.

    Parameters
    ----------
    terms : sequence of (E, nu, etaB, etaS) tuples
        For each Prony branch: Young's modulus, Poisson ratio, bulk viscosity,
        shear viscosity. All four are required per branch.
    """
    terms: Tuple[Tuple[float, float, float, float], ...] = ()

    def __post_init__(self):
        for i, term in enumerate(self.terms):
            if len(term) != 4:
                raise TypeError(
                    f"Viscoelasticity.terms[{i}] must be (E, nu, etaB, etaS); "
                    f"got {len(term)} values. Previous (g, tau) layout is no "
                    "longer supported — port to the Prony_Nfast form."
                )

    @property
    def mechanism_type(self) -> MechanismType:
        return MechanismType.VISCOELASTICITY

    def to_props(self) -> List[float]:
        """Return props values for this mechanism (excluding mechanism_type header)."""
        header = [float(len(self.terms))]  # N_prony
        params: List[float] = []
        for E, nu, etaB, etaS in self.terms:
            params.extend([E, nu, etaB, etaS])
        return header + params

    @property
    def nstatev(self) -> int:
        """Per Prony branch: 1 scalar lead variable v_i + 6-Voigt EV_i = 7."""
        return 7 * len(self.terms)


@dataclass(frozen=True)
class Damage:
    """Scalar damage mechanism.

    Parameters
    ----------
    Y_0 : float
        Damage threshold (no damage below this energy).
    Y_c : float
        Critical damage driving force.
    damage_type : DamageType
        Type of damage evolution law (default: LINEAR).
    A : float, optional
        Scale parameter (for EXPONENTIAL and WEIBULL).
    n : float, optional
        Exponent parameter (for POWER_LAW and WEIBULL).
    """
    Y_0: float
    Y_c: float
    damage_type: DamageType = DamageType.LINEAR
    A: float = 0.0
    n: float = 1.0

    @property
    def mechanism_type(self) -> MechanismType:
        return MechanismType.DAMAGE

    def to_props(self) -> List[float]:
        """Return props values for this mechanism (excluding mechanism_type header)."""
        # Header: damage_type
        header = [float(self.damage_type)]
        # Common params: Y_0, Y_c
        params = [self.Y_0, self.Y_c]
        # Type-specific params
        if self.damage_type == DamageType.EXPONENTIAL:
            params.append(self.A)
        elif self.damage_type == DamageType.POWER_LAW:
            params.append(self.n)
        elif self.damage_type == DamageType.WEIBULL:
            params.extend([self.A, self.n])
        return header + params

    @property
    def nstatev(self) -> int:
        """Number of state variables: D(1) + Y_max(1)."""
        return 2


Mechanism = Union[Plasticity, Viscoelasticity, Damage]

# ============================================================================
# Modular material orchestrator
# ============================================================================

class ModularMaterial:
    """Composable constitutive model using the modular UMAT system.

    This class assembles an elasticity module and one or more strain mechanisms
    into a complete constitutive model. It builds the flat props array expected
    by the C++ ``umat_modular`` function, which can be called via
    ``simcoon.umat("MODUL", ...)``.

    Parameters
    ----------
    elasticity : Elasticity
        Elasticity configuration.
    mechanisms : list of Mechanism, optional
        List of strain mechanisms (plasticity, viscoelasticity, damage).

    Examples
    --------
    Pure elasticity:

    >>> mat = ModularMaterial(elasticity=IsotropicElasticity(C1=210000., C2=0.3))
    >>> mat.props  # array([0., 210000., 0.3, 0., 0.])

    Elastoplastic with Voce hardening:

    >>> mat = ModularMaterial(
    ...     elasticity=IsotropicElasticity(C1=210000., C2=0.3, alpha=1.2e-5),
    ...     mechanisms=[
    ...         Plasticity(
    ...             sigma_Y=300.,
    ...             yield_criterion=VonMisesYield(),
    ...             isotropic_hardening=VoceHardening(Q=100., b=10.),
    ...         )
    ...     ]
    ... )

    Coupled plasticity + viscoelasticity:

    >>> mat = ModularMaterial(
    ...     elasticity=IsotropicElasticity(C1=70000., C2=0.33, alpha=2.3e-5),
    ...     mechanisms=[
    ...         Plasticity(sigma_Y=200., isotropic_hardening=PowerLawHardening(k=500., m=0.3)),
    ...         Viscoelasticity(terms=(
    ...             (70000., 0.33, 7e5, 3e5),    # branch 1: (E, nu, etaB, etaS)
    ...             (35000., 0.33, 3.5e5, 1.5e5),
    ...         )),
    ...     ]
    ... )

    Notes
    -----
    Under finite strain (NLGEOM control types 2-6), the composition becomes a
    Hencky hyperelastic law: the solver kinematics feed the model the
    logarithmic strain, so the elasticity block acts as a stored-energy
    function of ln V and the mechanisms ride additively on that measure.
    This holds only when the accumulated corotational strain is exactly the
    logarithmic strain, so the solver requires ``corate_type=3`` (log_R) for
    ``"MODUL"`` under NLGEOM and raises a ``RuntimeError`` for any other
    corate (which would degrade the model to a non-integrable hypoelastic
    rate with spurious dissipation in closed cycles).
    """

    def __init__(
        self,
        elasticity: Elasticity,
        mechanisms: Sequence[Mechanism] = (),
    ):
        self._elasticity = elasticity
        self._mechanisms = list(mechanisms)
        self._props: NDArray[np.float64] | None = None

    @property
    def elasticity(self) -> Elasticity:
        """The elasticity configuration."""
        return self._elasticity

    @property
    def mechanisms(self) -> List[Mechanism]:
        """The list of strain mechanisms."""
        return list(self._mechanisms)

    @property
    def umat_name(self) -> str:
        """UMAT name string for use with ``simcoon.umat()``."""
        return "MODUL"

    @property
    def activation_index(self) -> int:
        """Index of the muscle activation in :attr:`props`.

        The activation of a :class:`MuscleElasticity` block is the one parameter that
        is *driven*: a caller rewrites it every increment and passes the props again,
        which is how a time-varying, per-point activation reaches the kernel without
        any change to the UMAT interface. Read the index from here rather than
        hard-coding it::

            props = np.tile(mat.props[:, None], (1, n_points))
            props[mat.activation_index, :] = activation_field
            sim.umat("MODUL", ..., props=props, ...)

        Raises
        ------
        TypeError
            if the elasticity block is not a :class:`MuscleElasticity`.
        """
        if not isinstance(self._elasticity, MuscleElasticity):
            raise TypeError(f"activation_index: the elasticity block is a "
                            f"{type(self._elasticity).__name__}, which has no activation")
        # props[0] is the elasticity type, props[1] the potential code and props[2] the
        # parameter count; the potential's own parameters start at 3.
        return 3 + MuscleElasticity.activation_index

    @property
    def props(self) -> NDArray[np.float64]:
        """Build and return the flat props array for the C++ modular UMAT.

        The props format follows the ``configure_from_props()`` convention:

        - ``props[0]``: elasticity type
        - ``props[1..N_el]``: elasticity block — for the linear types the
          elastic-constant convention code (see :class:`IsoConvention` and
          friends) then the constants; for the hyperelastic blocks
          ``[potential, n_params, params..., alpha]``
        - ``props[N_el+1]``: number of mechanisms
        - For each mechanism: type code + mechanism-specific parameters
        """
        if self._props is not None:
            return self._props

        values: List[float] = []

        # Elasticity type + parameters
        values.append(float(self._elasticity.elasticity_type))
        values.extend(self._elasticity.to_props())

        # Number of mechanisms
        values.append(float(len(self._mechanisms)))

        # Each mechanism
        for mech in self._mechanisms:
            values.append(float(mech.mechanism_type))
            values.extend(mech.to_props())

        self._props = np.array(values, dtype=np.float64)
        self._props.flags.writeable = False
        return self._props

    @property
    def nstatev(self) -> int:
        """Compute the total number of state variables required.

        Accounts for:
        - T_init (1 scalar)
        - Each mechanism's internal variables
        """
        count = 1  # T_init
        for mech in self._mechanisms:
            count += mech.nstatev
        return count

    @property
    def nprops(self) -> int:
        """Total number of material properties."""
        return len(self.props)

    def __repr__(self) -> str:
        mechs = ", ".join(m.__class__.__name__ for m in self._mechanisms)
        return (f"ModularMaterial(elasticity={self._elasticity.__class__.__name__}, "
                f"mechanisms=[{mechs}], nprops={self.nprops}, nstatev={self.nstatev})")

    def summary(self) -> str:
        """Return a human-readable summary of the material configuration."""
        lines = ["ModularMaterial:"]
        el = self._elasticity
        lines.append(f"  Elasticity: {el.__class__.__name__}")
        if isinstance(el, IsotropicElasticity):
            lines.append(f"    C1={el.C1}, C2={el.C2}, alpha={el.alpha} "
                         f"[{el.convention.name}]")
        elif isinstance(el, CubicElasticity):
            lines.append(f"    C1={el.C1}, C2={el.C2}, C3={el.C3}, alpha={el.alpha} "
                         f"[{el.convention.name}]")
        elif isinstance(el, TransverseIsotropicElasticity):
            lines.append(f"    EL={el.EL}, ET={el.ET}, nuTL={el.nuTL}, nuTT={el.nuTT}")
            lines.append(f"    GLT={el.GLT}, alpha_L={el.alpha_L}, alpha_T={el.alpha_T}, axis={el.axis}")
        elif isinstance(el, OrthotropicElasticity):
            lines.append(f"    C1={el.C1}, C2={el.C2}, C3={el.C3}")
            lines.append(f"    C4={el.C4}, C5={el.C5}, C6={el.C6}")
            lines.append(f"    C7={el.C7}, C8={el.C8}, C9={el.C9} [{el.convention.name}]")
            lines.append(f"    alpha1={el.alpha1}, alpha2={el.alpha2}, alpha3={el.alpha3}")
        elif isinstance(el, SwansonElasticity):
            lines.append(f"    {len(el.terms)} Swanson terms (A, B, alpha, beta), "
                         f"kappa={el.kappa}, volumetric={el.volumetric}, alpha={el.alpha}")
            for k, term in enumerate(el.terms):
                A, B, a, b = term
                lines.append(f"      [{k}] A={A}, B={B}, alpha={a}, beta={b}")
        elif isinstance(el, _HyperInvariantsElasticity):
            named = [f.name for f in fields(el) if f.name != "alpha"]
            params = ", ".join(f"{n}={getattr(el, n)}" for n in named)
            lines.append(f"    {params}, alpha={el.alpha}")

        if not self._mechanisms:
            lines.append("  Mechanisms: (none - pure elastic)")
        else:
            lines.append(f"  Mechanisms ({len(self._mechanisms)}):")
            for i, mech in enumerate(self._mechanisms):
                lines.append(f"    [{i}] {mech.__class__.__name__}")
                if isinstance(mech, Plasticity):
                    lines.append(f"        sigma_Y={mech.sigma_Y}")
                    lines.append(f"        yield: {mech.yield_criterion.__class__.__name__}")
                    lines.append(f"        iso_hard: {mech.isotropic_hardening.__class__.__name__}")
                    lines.append(f"        kin_hard: {mech.kinematic_hardening.__class__.__name__}")
                elif isinstance(mech, Viscoelasticity):
                    lines.append(f"        {len(mech.terms)} Prony terms "
                                 "(E_i, nu_i, etaB_i, etaS_i)")
                    for k, term in enumerate(mech.terms):
                        E_i, nu_i, etaB_i, etaS_i = term
                        lines.append(
                            f"          [{k}] E={E_i}, nu={nu_i}, "
                            f"etaB={etaB_i}, etaS={etaS_i}")
                elif isinstance(mech, Damage):
                    lines.append(f"        type: {mech.damage_type.name}")
                    lines.append(f"        Y_0={mech.Y_0}, Y_c={mech.Y_c}")

        lines.append(f"  nprops={self.nprops}, nstatev={self.nstatev}")
        return "\n".join(lines)


# ============================================================================
# Factory functions for common material models
# ============================================================================

def elastic_model(E: float, nu: float, alpha: float = 0.0) -> ModularMaterial:
    """Create an isotropic elastic material.

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.
    alpha : float
        Coefficient of thermal expansion (default: 0).

    Returns
    -------
    ModularMaterial
    """
    return ModularMaterial(
        elasticity=IsotropicElasticity(C1=E, C2=nu, alpha=alpha)
    )


def elastoplastic_model(
    E: float,
    nu: float,
    sigma_Y: float,
    k: float = 0.0,
    m: float = 1.0,
    alpha: float = 0.0,
) -> ModularMaterial:
    """Create an isotropic elastoplastic material with power-law hardening.

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.
    sigma_Y : float
        Initial yield stress.
    k : float
        Power-law hardening coefficient.
    m : float
        Power-law hardening exponent.
    alpha : float
        Coefficient of thermal expansion.

    Returns
    -------
    ModularMaterial
    """
    return ModularMaterial(
        elasticity=IsotropicElasticity(C1=E, C2=nu, alpha=alpha),
        mechanisms=[
            Plasticity(
                sigma_Y=sigma_Y,
                isotropic_hardening=PowerLawHardening(k=k, m=m),
            )
        ],
    )


def viscoelastic_model(
    E: float,
    nu: float,
    prony_terms: Sequence[Tuple[float, float, float, float]],
    alpha: float = 0.0,
) -> ModularMaterial:
    """Create an isotropic generalized-Maxwell viscoelastic material.

    Parameters
    ----------
    E : float
        **Instantaneous** (glassy) Young's modulus, i.e. E(0) of the Prony series:
        the long-term modulus is ``E - sum_i E_i`` and must stay positive. Pass
        ``E_inf + sum_i E_i`` when calibration gives the long-term modulus.
    nu : float
        Poisson's ratio of that same reference stiffness.
    prony_terms : sequence of (E_i, nu_i, etaB_i, etaS_i) tuples
        Per-branch parameters: branch modulus, branch Poisson ratio, bulk
        viscosity, shear viscosity.
    alpha : float
        Coefficient of thermal expansion.

    Returns
    -------
    ModularMaterial
    """
    return ModularMaterial(
        elasticity=IsotropicElasticity(C1=E, C2=nu, alpha=alpha),
        mechanisms=[
            Viscoelasticity(terms=tuple(prony_terms)),
        ],
    )
