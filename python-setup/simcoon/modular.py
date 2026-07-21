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

from dataclasses import dataclass, field
from enum import IntEnum
from typing import List, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray


__all__ = [
    # Enums
    "ElasticityType", "YieldType", "IsoHardType", "KinHardType",
    "DamageType", "MechanismType",
    # Elastic-constant conventions
    "IsoConvention", "CubicConvention",
    "IsotransConvention", "OrthoConvention",
    # Elasticity
    "IsotropicElasticity", "CubicElasticity",
    "TransverseIsotropicElasticity", "OrthotropicElasticity",
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
    """Type of linear elasticity."""
    ISOTROPIC = 0
    CUBIC = 1
    TRANSVERSE_ISOTROPIC = 2
    ORTHOTROPIC = 3


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


Elasticity = Union[IsotropicElasticity, CubicElasticity,
                   TransverseIsotropicElasticity, OrthotropicElasticity]

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

    where M_0 is the compliance at the reference (undamaged, long-term)
    stiffness.

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
    def props(self) -> NDArray[np.float64]:
        """Build and return the flat props array for the C++ modular UMAT.

        The props format follows the ``configure_from_props()`` convention:

        - ``props[0]``: elasticity type
        - ``props[1]``: elastic-constant convention code (see
          :class:`IsoConvention` and friends)
        - ``props[2..N_el]``: elasticity parameters
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
        Reference (long-term) Young's modulus.
    nu : float
        Reference Poisson's ratio.
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
