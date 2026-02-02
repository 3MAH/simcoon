"""
Python 0D Solver for Simcoon.

This module provides a Python-based material point solver that mirrors the C++ solver
architecture while enabling flexible material simulations with Python control flow.
It leverages the existing C++ UMAT implementations via pybind11.

The solver uses Newton-Raphson iteration to solve mixed strain/stress boundary
conditions at a material point.

Note: This implementation minimizes array copies following the carma copy=false pattern.
Arrays are passed directly to bindings where possible, and np.copyto() is used for
in-place copying instead of creating new arrays.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import List, Optional, Literal

import numpy as np
from numpy.linalg import norm
import simcoon._core as scc


# =============================================================================
# Control Type and Corate Type Mappings
# =============================================================================

CONTROL_TYPES = {
    'small_strain': 1,
    'green_lagrange': 2,
    'logarithmic': 3,
    'biot': 4,
    'F': 5,
    'gradU': 6,
}


# =============================================================================
# Lightweight History Point (optimized for minimal memory allocation)
# =============================================================================

@dataclass(slots=True)
class HistoryPoint:
    """
    Lightweight history point storing only essential state variables.

    This class is optimized for history storage, containing only the fields
    typically accessed after simulation: strain, stress, work, and state variables.
    Using this instead of full StateVariablesM copies reduces memory allocation
    by ~8x per history point.

    Uses __slots__ for faster attribute access and lower memory footprint.

    Attributes
    ----------
    Etot : np.ndarray
        Green-Lagrange strain tensor in Voigt notation (6,)
    sigma : np.ndarray
        Cauchy stress tensor in Voigt notation (6,)
    Wm : np.ndarray
        Mechanical work components [Wm, Wm_r, Wm_ir, Wm_d] (4,)
    statev : np.ndarray
        Internal state variables vector (nstatev,)
    R : np.ndarray
        Rotation tensor (3,3) - for objective rate analysis
    T : float
        Current temperature
    """
    Etot: np.ndarray
    sigma: np.ndarray
    Wm: np.ndarray
    statev: np.ndarray
    R: np.ndarray
    T: float

    @classmethod
    def from_state(cls, sv: 'StateVariablesM') -> 'HistoryPoint':
        """Create a HistoryPoint from a StateVariablesM (copies only essential fields)."""
        return cls(
            Etot=sv.Etot.copy(),
            sigma=sv.sigma.copy(),
            Wm=sv.Wm.copy(),
            statev=sv.statev.copy(),
            R=sv.R.copy(),
            T=sv.T,
        )

CORATE_TYPES = {
    'jaumann': 0,
    'green_naghdi': 1,
    'logarithmic': 2,
    'logarithmic_R': 3,
    'truesdell': 4,
    'logarithmic_F': 5,
}


# =============================================================================
# State Variable Classes
# =============================================================================

@dataclass
class StateVariables:
    """
    Base state variables class mirroring C++ state_variables.

    Stores all mechanical state variables (strains, stresses, deformation gradients)
    in both current and reference configurations. Supports finite strain formulations
    with various stress measures (Cauchy, Kirchhoff, 2nd Piola-Kirchhoff) and strain
    measures (Green-Lagrange, logarithmic).

    Note: This class uses in-place operations (np.copyto) to minimize memory allocation.
    Arrays are owned by the instance and modified in-place where possible.

    Attributes
    ----------
    Etot : np.ndarray
        Green-Lagrange strain tensor in Voigt notation (6,)
    DEtot : np.ndarray
        Increment of Green-Lagrange strain (6,)
    etot : np.ndarray
        Logarithmic (Hencky) strain tensor in Voigt notation (6,)
    Detot : np.ndarray
        Increment of logarithmic strain (6,)
    PKII : np.ndarray
        2nd Piola-Kirchhoff stress tensor in Voigt notation (6,)
    PKII_start : np.ndarray
        2nd Piola-Kirchhoff stress at start of increment (6,)
    tau : np.ndarray
        Kirchhoff stress tensor in Voigt notation (6,)
    tau_start : np.ndarray
        Kirchhoff stress at start of increment (6,)
    sigma : np.ndarray
        Cauchy stress tensor in Voigt notation (6,)
    sigma_start : np.ndarray
        Cauchy stress at start of increment (6,)
    F0 : np.ndarray
        Deformation gradient at start of increment (3,3)
    F1 : np.ndarray
        Deformation gradient at end of increment (3,3)
    U0 : np.ndarray
        Right stretch tensor at start of increment (3,3)
    U1 : np.ndarray
        Right stretch tensor at end of increment (3,3)
    R : np.ndarray
        Rotation tensor (3,3)
    DR : np.ndarray
        Increment of rotation tensor (3,3)
    T : float
        Current temperature
    DT : float
        Temperature increment
    nstatev : int
        Number of internal state variables
    statev : np.ndarray
        Internal state variables vector (nstatev,)
    statev_start : np.ndarray
        Internal state variables at start of increment (nstatev,)
    """

    # Strain measures
    Etot: np.ndarray = None
    DEtot: np.ndarray = None
    etot: np.ndarray = None
    Detot: np.ndarray = None

    # Stress measures
    PKII: np.ndarray = None
    PKII_start: np.ndarray = None
    tau: np.ndarray = None
    tau_start: np.ndarray = None
    sigma: np.ndarray = None
    sigma_start: np.ndarray = None

    # Deformation
    F0: np.ndarray = None
    F1: np.ndarray = None
    U0: np.ndarray = None
    U1: np.ndarray = None
    R: np.ndarray = None
    DR: np.ndarray = None

    # Temperature
    T: float = 293.15
    DT: float = 0.0

    # Internal state variables
    nstatev: int = 0
    statev: np.ndarray = None
    statev_start: np.ndarray = None

    def __post_init__(self):
        """Initialize arrays to default values if None."""
        if self.Etot is None:
            self.Etot = np.zeros(6)
        if self.DEtot is None:
            self.DEtot = np.zeros(6)
        if self.etot is None:
            self.etot = np.zeros(6)
        if self.Detot is None:
            self.Detot = np.zeros(6)
        if self.PKII is None:
            self.PKII = np.zeros(6)
        if self.PKII_start is None:
            self.PKII_start = np.zeros(6)
        if self.tau is None:
            self.tau = np.zeros(6)
        if self.tau_start is None:
            self.tau_start = np.zeros(6)
        if self.sigma is None:
            self.sigma = np.zeros(6)
        if self.sigma_start is None:
            self.sigma_start = np.zeros(6)
        if self.F0 is None:
            self.F0 = np.eye(3, order='F')
        if self.F1 is None:
            self.F1 = np.eye(3, order='F')
        if self.U0 is None:
            self.U0 = np.eye(3, order='F')
        if self.U1 is None:
            self.U1 = np.eye(3, order='F')
        if self.R is None:
            self.R = np.eye(3, order='F')
        if self.DR is None:
            self.DR = np.eye(3, order='F')
        if self.statev is None:
            self.statev = np.zeros(max(1, self.nstatev))
        if self.statev_start is None:
            self.statev_start = np.zeros(max(1, self.nstatev))

    def copy(self) -> 'StateVariables':
        """Create a copy of this StateVariables object (optimized, avoids deepcopy)."""
        return StateVariables(
            Etot=self.Etot.copy(),
            DEtot=self.DEtot.copy(),
            etot=self.etot.copy(),
            Detot=self.Detot.copy(),
            PKII=self.PKII.copy(),
            PKII_start=self.PKII_start.copy(),
            tau=self.tau.copy(),
            tau_start=self.tau_start.copy(),
            sigma=self.sigma.copy(),
            sigma_start=self.sigma_start.copy(),
            F0=self.F0.copy(),
            F1=self.F1.copy(),
            U0=self.U0.copy(),
            U1=self.U1.copy(),
            R=self.R.copy(),
            DR=self.DR.copy(),
            T=self.T,
            DT=self.DT,
            nstatev=self.nstatev,
            statev=self.statev.copy(),
            statev_start=self.statev_start.copy(),
        )

    def copy_to(self, other: 'StateVariables'):
        """
        Copy all values from this object to another (in-place).

        Parameters
        ----------
        other : StateVariables
            Target object to copy values into
        """
        np.copyto(other.Etot, self.Etot)
        np.copyto(other.DEtot, self.DEtot)
        np.copyto(other.etot, self.etot)
        np.copyto(other.Detot, self.Detot)
        np.copyto(other.PKII, self.PKII)
        np.copyto(other.PKII_start, self.PKII_start)
        np.copyto(other.tau, self.tau)
        np.copyto(other.tau_start, self.tau_start)
        np.copyto(other.sigma, self.sigma)
        np.copyto(other.sigma_start, self.sigma_start)
        np.copyto(other.F0, self.F0)
        np.copyto(other.F1, self.F1)
        np.copyto(other.U0, self.U0)
        np.copyto(other.U1, self.U1)
        np.copyto(other.R, self.R)
        np.copyto(other.DR, self.DR)
        other.T = self.T
        other.DT = self.DT
        other.nstatev = self.nstatev
        np.copyto(other.statev, self.statev)
        np.copyto(other.statev_start, self.statev_start)

    def to_start(self):
        """
        Reset current values TO start-of-increment values (for rollback).

        Matches C++ state_variables::to_start() - used to reset trial solution
        when NR iteration fails or when restarting an increment attempt.
        """
        np.copyto(self.PKII, self.PKII_start)
        np.copyto(self.tau, self.tau_start)
        np.copyto(self.sigma, self.sigma_start)
        np.copyto(self.statev, self.statev_start)

    def set_start(self, corate_type: int = 0):
        """
        SET _start values from current converged values and advance state.

        Matches C++ state_variables::set_start() - called after a converged
        increment to update _start values and advance strain/rotation.

        Parameters
        ----------
        corate_type : int
            Corotational rate type (0=jaumann, 1=green_naghdi, etc.)
        """
        # For small strain (corate_type not used), simple copy
        np.copyto(self.PKII_start, self.PKII)
        np.copyto(self.tau_start, self.tau)
        np.copyto(self.sigma_start, self.sigma)
        np.copyto(self.statev_start, self.statev)

        # Advance strain
        self.Etot += self.DEtot
        self.etot += self.Detot
        self.T += self.DT

        # Update deformation tensors
        np.copyto(self.F0, self.F1)
        np.copyto(self.U0, self.U1)


@dataclass
class StateVariablesM(StateVariables):
    """
    Mechanical state variables class mirroring C++ state_variables_M.

    Extends StateVariables with mechanical-specific fields including
    internal stress, mechanical work components, and tangent stiffness matrices.

    Attributes
    ----------
    sigma_in : np.ndarray
        Internal stress vector (6,)
    sigma_in_start : np.ndarray
        Internal stress at start of increment (6,)
    Wm : np.ndarray
        Mechanical work components [Wm, Wm_r, Wm_ir, Wm_d] (4,)
    Wm_start : np.ndarray
        Mechanical work at start of increment (4,)
    L : np.ndarray
        Elastic stiffness matrix (6,6)
    Lt : np.ndarray
        Tangent modulus matrix (6,6)
    """

    sigma_in: np.ndarray = None
    sigma_in_start: np.ndarray = None
    Wm: np.ndarray = None
    Wm_start: np.ndarray = None
    L: np.ndarray = None
    Lt: np.ndarray = None

    def __post_init__(self):
        """Initialize arrays to default values if None."""
        super().__post_init__()
        if self.sigma_in is None:
            self.sigma_in = np.zeros(6)
        if self.sigma_in_start is None:
            self.sigma_in_start = np.zeros(6)
        if self.Wm is None:
            self.Wm = np.zeros(4)
        if self.Wm_start is None:
            self.Wm_start = np.zeros(4)
        if self.L is None:
            self.L = np.zeros((6, 6), order='F')
        if self.Lt is None:
            self.Lt = np.zeros((6, 6), order='F')

    def copy(self) -> 'StateVariablesM':
        """Create a copy of this StateVariablesM object (optimized, avoids deepcopy)."""
        return StateVariablesM(
            Etot=self.Etot.copy(),
            DEtot=self.DEtot.copy(),
            etot=self.etot.copy(),
            Detot=self.Detot.copy(),
            PKII=self.PKII.copy(),
            PKII_start=self.PKII_start.copy(),
            tau=self.tau.copy(),
            tau_start=self.tau_start.copy(),
            sigma=self.sigma.copy(),
            sigma_start=self.sigma_start.copy(),
            F0=self.F0.copy(),
            F1=self.F1.copy(),
            U0=self.U0.copy(),
            U1=self.U1.copy(),
            R=self.R.copy(),
            DR=self.DR.copy(),
            T=self.T,
            DT=self.DT,
            nstatev=self.nstatev,
            statev=self.statev.copy(),
            statev_start=self.statev_start.copy(),
            sigma_in=self.sigma_in.copy(),
            sigma_in_start=self.sigma_in_start.copy(),
            Wm=self.Wm.copy(),
            Wm_start=self.Wm_start.copy(),
            L=self.L.copy(),
            Lt=self.Lt.copy(),
        )

    def copy_to(self, other: 'StateVariablesM'):
        """Copy all values from this object to another (in-place)."""
        super().copy_to(other)
        np.copyto(other.sigma_in, self.sigma_in)
        np.copyto(other.sigma_in_start, self.sigma_in_start)
        np.copyto(other.Wm, self.Wm)
        np.copyto(other.Wm_start, self.Wm_start)
        np.copyto(other.L, self.L)
        np.copyto(other.Lt, self.Lt)

    def to_start(self):
        """
        Reset current values TO start-of-increment values (for rollback).

        Matches C++ state_variables_M::to_start().
        """
        super().to_start()
        np.copyto(self.sigma_in, self.sigma_in_start)
        np.copyto(self.Wm, self.Wm_start)

    def set_start(self, corate_type: int = 0):
        """
        SET _start values from current converged values and advance state.

        Matches C++ state_variables_M::set_start().
        """
        super().set_start(corate_type)
        np.copyto(self.sigma_in_start, self.sigma_in)
        np.copyto(self.Wm_start, self.Wm)


@dataclass
class StateVariablesT(StateVariables):
    """
    Thermomechanical state variables class mirroring C++ state_variables_T.

    Extends StateVariables with thermal-specific fields including
    heat quantities and coupled thermomechanical tangent matrices.

    Attributes
    ----------
    sigma_in : np.ndarray
        Internal stress vector (6,)
    sigma_in_start : np.ndarray
        Internal stress at start of increment (6,)
    Wm : np.ndarray
        Mechanical work components (4,)
    Wt : np.ndarray
        Thermal work components (4,)
    Wm_start : np.ndarray
        Mechanical work at start of increment (4,)
    Wt_start : np.ndarray
        Thermal work at start of increment (4,)
    dSdE : np.ndarray
        Mechanical tangent dStress/dStrain (6,6)
    dSdEt : np.ndarray
        Coupling tangent (6,6)
    dSdT : np.ndarray
        dStress/dTemperature (6,1)
    Q : float
        Heat flux
    r : float
        Heat source
    r_in : float
        Internal heat source
    drdE : np.ndarray
        dr/dStrain (1,6)
    drdT : np.ndarray
        dr/dTemperature (1,1)
    """

    sigma_in: np.ndarray = None
    sigma_in_start: np.ndarray = None
    Wm: np.ndarray = None
    Wt: np.ndarray = None
    Wm_start: np.ndarray = None
    Wt_start: np.ndarray = None

    # Thermomechanical tangents
    dSdE: np.ndarray = None
    dSdEt: np.ndarray = None
    dSdT: np.ndarray = None

    # Heat quantities
    Q: float = 0.0
    r: float = 0.0
    r_in: float = 0.0
    drdE: np.ndarray = None
    drdT: np.ndarray = None

    def __post_init__(self):
        """Initialize arrays to default values if None."""
        super().__post_init__()
        if self.sigma_in is None:
            self.sigma_in = np.zeros(6)
        if self.sigma_in_start is None:
            self.sigma_in_start = np.zeros(6)
        if self.Wm is None:
            self.Wm = np.zeros(4)
        if self.Wt is None:
            self.Wt = np.zeros(4)
        if self.Wm_start is None:
            self.Wm_start = np.zeros(4)
        if self.Wt_start is None:
            self.Wt_start = np.zeros(4)
        if self.dSdE is None:
            self.dSdE = np.zeros((6, 6))
        if self.dSdEt is None:
            self.dSdEt = np.zeros((6, 6))
        if self.dSdT is None:
            self.dSdT = np.zeros((6, 1))
        if self.drdE is None:
            self.drdE = np.zeros((1, 6))
        if self.drdT is None:
            self.drdT = np.zeros((1, 1))

    def copy_to(self, other: 'StateVariablesT'):
        """Copy all values from this object to another (in-place)."""
        super().copy_to(other)
        np.copyto(other.sigma_in, self.sigma_in)
        np.copyto(other.sigma_in_start, self.sigma_in_start)
        np.copyto(other.Wm, self.Wm)
        np.copyto(other.Wt, self.Wt)
        np.copyto(other.Wm_start, self.Wm_start)
        np.copyto(other.Wt_start, self.Wt_start)
        np.copyto(other.dSdE, self.dSdE)
        np.copyto(other.dSdEt, self.dSdEt)
        np.copyto(other.dSdT, self.dSdT)
        other.Q = self.Q
        other.r = self.r
        other.r_in = self.r_in
        np.copyto(other.drdE, self.drdE)
        np.copyto(other.drdT, self.drdT)

    def to_start(self):
        """
        Reset current values TO start-of-increment values (for rollback).

        Matches C++ state_variables_T::to_start().
        """
        super().to_start()
        np.copyto(self.sigma_in, self.sigma_in_start)
        np.copyto(self.Wm, self.Wm_start)
        np.copyto(self.Wt, self.Wt_start)

    def set_start(self, corate_type: int = 0):
        """
        SET _start values from current converged values and advance state.

        Matches C++ state_variables_T::set_start().
        """
        super().set_start(corate_type)
        np.copyto(self.sigma_in_start, self.sigma_in)
        np.copyto(self.Wm_start, self.Wm)
        np.copyto(self.Wt_start, self.Wt)


# =============================================================================
# Step Classes
# =============================================================================

@dataclass
class Step:
    """
    Base class for a loading step.

    A step defines a loading increment with targets for strain/stress components
    and control mode for each component.

    Attributes
    ----------
    Dn_init : int
        Initial number of sub-increments
    Dn_mini : int
        Minimum number of sub-increments
    Dn_inc : int
        Maximum number of sub-increments
    control : List[str]
        Control mode per component ('strain' or 'stress'), length 6
    time : float
        Time duration for this step
    """

    Dn_init: int = 1
    Dn_mini: int = 1
    Dn_inc: int = 100
    control: List[str] = None
    time: float = 1.0

    def __post_init__(self):
        if self.control is None:
            self.control = ['strain'] * 6

    def get_cBC_meca(self) -> np.ndarray:
        """
        Get mechanical boundary condition control array.

        Returns
        -------
        np.ndarray
            Array of shape (6,) where 1 = stress controlled, 0 = strain controlled
        """
        return np.array([1 if c == 'stress' else 0 for c in self.control], dtype=int)


@dataclass
class StepMeca(Step):
    """
    Mechanical loading step.

    Defines a mechanical loading increment with strain and/or stress targets.

    Attributes
    ----------
    DEtot_end : np.ndarray
        Target strain increment (6,), for strain-controlled components
    Dsigma_end : np.ndarray
        Target stress increment (6,), for stress-controlled components
    """

    DEtot_end: np.ndarray = None
    Dsigma_end: np.ndarray = None

    def __post_init__(self):
        super().__post_init__()
        if self.DEtot_end is None:
            self.DEtot_end = np.zeros(6)
        if self.Dsigma_end is None:
            self.Dsigma_end = np.zeros(6)


@dataclass
class StepThermomeca(StepMeca):
    """
    Thermomechanical loading step.

    Extends StepMeca with temperature control.

    Attributes
    ----------
    DT_end : float
        Target temperature increment
    Q_end : float
        Target heat flux
    thermal_control : str
        Thermal control mode: 'temperature', 'heat_flux', or 'convection'
    """

    DT_end: float = 0.0
    Q_end: float = 0.0
    thermal_control: str = 'temperature'

    def get_cBC_T(self) -> int:
        """
        Get thermal boundary condition control flag.

        Returns
        -------
        int
            0 for temperature controlled, 1 for heat flux controlled
        """
        return 0 if self.thermal_control == 'temperature' else 1


# =============================================================================
# Block Class
# =============================================================================

@dataclass
class Block:
    """
    A block containing multiple steps with shared settings.

    A block groups steps that share the same UMAT, control type, and corotational
    formulation settings.

    Attributes
    ----------
    steps : List[Step]
        List of steps in this block
    nstatev : int
        Number of internal state variables for the UMAT
    umat_name : str
        Name of the UMAT to use (e.g., 'ELISO', 'EPICP')
    umat_type : str
        Type of UMAT: 'mechanical' or 'thermomechanical'
    props : np.ndarray
        Material properties array for the UMAT
    control_type : str
        Control type: 'small_strain', 'green_lagrange', 'logarithmic', 'biot', 'F', 'gradU'
    corate_type : str
        Corotational rate type: 'jaumann', 'green_naghdi', 'logarithmic', 'truesdell'
    ncycle : int
        Number of cycles to repeat the steps
    """

    steps: List[Step] = None
    nstatev: int = 0
    umat_name: str = "ELISO"
    umat_type: str = "mechanical"
    props: np.ndarray = None
    control_type: str = 'small_strain'
    corate_type: str = 'jaumann'
    ncycle: int = 1

    def __post_init__(self):
        if self.steps is None:
            self.steps = []
        if self.props is None:
            self.props = np.array([])

    def add_step(self, step: Step):
        """Add a step to this block."""
        self.steps.append(step)

    def get_control_type_int(self) -> int:
        """Get the integer control type code."""
        return CONTROL_TYPES.get(self.control_type, 1)

    def get_corate_type_int(self) -> int:
        """Get the integer corotational type code."""
        return CORATE_TYPES.get(self.corate_type, 0)


# =============================================================================
# Jacobian Helper Functions
# =============================================================================

def Lt_2_K(Lt: np.ndarray, cBC_meca: np.ndarray, lambda_solver: float,
           K: np.ndarray = None) -> np.ndarray:
    """
    Build 6x6 Jacobian for mechanical solver.

    Constructs the Jacobian matrix for mixed strain/stress boundary conditions.

    Parameters
    ----------
    Lt : np.ndarray
        Tangent modulus matrix (6,6)
    cBC_meca : np.ndarray
        Boundary condition control array (6,) where 1 = stress controlled
    lambda_solver : float
        Effective stiffness for strain-controlled components
    K : np.ndarray, optional
        Pre-allocated output array (6,6). If None, a new array is created.

    Returns
    -------
    np.ndarray
        Jacobian matrix K (6,6)
    """
    if K is None:
        K = np.zeros((6, 6))
    else:
        K.fill(0.0)

    for i in range(6):
        if cBC_meca[i]:
            K[i, :] = Lt[i, :]
        else:
            K[i, i] = lambda_solver
    return K


def Lth_2_K(dSdE: np.ndarray, dSdT: np.ndarray, dQdE: np.ndarray, dQdT: float,
            cBC_meca: np.ndarray, cBC_T: int, lambda_solver: float,
            K: np.ndarray = None) -> np.ndarray:
    """
    Build 7x7 Jacobian for thermomechanical solver.

    Constructs the coupled thermomechanical Jacobian matrix.

    Parameters
    ----------
    dSdE : np.ndarray
        Mechanical tangent (6,6)
    dSdT : np.ndarray
        Stress-temperature coupling (6,1)
    dQdE : np.ndarray
        Heat-strain coupling (1,6)
    dQdT : float
        Thermal tangent
    cBC_meca : np.ndarray
        Mechanical boundary condition control (6,)
    cBC_T : int
        Thermal boundary condition control (0=temp, 1=heat flux)
    lambda_solver : float
        Effective stiffness for controlled components
    K : np.ndarray, optional
        Pre-allocated output array (7,7). If None, a new array is created.

    Returns
    -------
    np.ndarray
        Jacobian matrix K (7,7)
    """
    if K is None:
        K = np.zeros((7, 7))
    else:
        K.fill(0.0)

    K[0:6, 0:6] = dSdE
    K[0:6, 6:7] = dSdT.reshape(6, 1) if dSdT.ndim == 1 else dSdT
    K[6:7, 0:6] = dQdE.reshape(1, 6) if dQdE.ndim == 1 else dQdE
    K[6, 6] = dQdT

    for i in range(6):
        if cBC_meca[i] == 0:
            K[i, :] = 0.0
            K[i, i] = lambda_solver

    if cBC_T == 0:
        K[6, :] = 0.0
        K[6, 6] = lambda_solver

    return K


# =============================================================================
# Main Solver Class
# =============================================================================

class Solver:
    """
    0D material point solver with Newton-Raphson iterations.

    This solver handles mechanical and thermomechanical simulations at a material
    point, with support for mixed strain/stress boundary conditions and various
    finite strain formulations.

    Note: This implementation minimizes array copies. Arrays are passed directly
    to C++ bindings where possible (carma copy=false pattern), and in-place
    operations are used throughout.

    Parameters
    ----------
    blocks : List[Block]
        List of loading blocks to simulate
    max_iter : int
        Maximum Newton-Raphson iterations per increment
    tol : float
        Convergence tolerance for Newton-Raphson
    lambda_solver : float
        Effective stiffness for strain-controlled components

    Attributes
    ----------
    history : List[HistoryPoint]
        History of essential state variables at each converged increment

    Examples
    --------
    >>> import numpy as np
    >>> from simcoon.solver import Solver, Block, StepMeca, StateVariablesM
    >>>
    >>> # Material properties for ELISO (E, nu)
    >>> props = np.array([210000.0, 0.3])
    >>>
    >>> # Uniaxial tension step
    >>> step = StepMeca(
    ...     DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    ...     control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    ...     Dn_init=10
    ... )
    >>>
    >>> block = Block(
    ...     steps=[step],
    ...     umat_name="ELISO",
    ...     props=props,
    ...     nstatev=1
    ... )
    >>>
    >>> solver = Solver(blocks=[block])
    >>> history = solver.solve()
    """

    def __init__(self, blocks: List[Block] = None,
                 max_iter: int = 10, tol: float = 1e-9,
                 lambda_solver: float = 10000.0,
                 div_tnew_dt: float = 0.5, mul_tnew_dt: float = 2.0):
        self.blocks = blocks or []
        self.max_iter = max_iter
        self.tol = tol
        self.lambda_solver = lambda_solver
        self.div_tnew_dt = div_tnew_dt
        self.mul_tnew_dt = mul_tnew_dt
        self.history = []

        # Pre-allocate work arrays for Newton-Raphson
        self._K = np.zeros((6, 6))
        self._residual = np.zeros(6)
        self._Delta = np.zeros(6)

        # Pre-allocate UMAT batch arrays (Fortran-contiguous for C++ binding)
        # These are modified in-place by umat_inplace for zero-copy performance
        self._etot_batch = np.zeros((6, 1), order='F')
        self._Detot_batch = np.zeros((6, 1), order='F')
        self._sigma_batch = np.zeros((6, 1), order='F')
        self._F0_batch = np.zeros((3, 3, 1), order='F')
        self._F1_batch = np.zeros((3, 3, 1), order='F')
        self._DR_batch = np.zeros((3, 3, 1), order='F')
        self._Wm_batch = np.zeros((4, 1), order='F')
        self._Lt_batch = np.zeros((6, 6, 1), order='F')  # Tangent modulus
        self._props_batch = None  # Allocated per-block (variable size)
        self._statev_batch = None  # Allocated per-block (variable size)

        # Cache function references for hot path (avoids repeated lookups)
        self._umat_inplace = scc.umat_inplace
        self._np_copyto = np.copyto
        self._np_fill_diagonal = np.fill_diagonal
        self._norm = norm

    def solve(self, sv_init: StateVariables = None) -> List[HistoryPoint]:
        """
        Run the full simulation.

        Parameters
        ----------
        sv_init : StateVariables, optional
            Initial state variables. If None, creates default StateVariablesM.

        Returns
        -------
        List[HistoryPoint]
            History of essential state variables at each converged increment.
            Each HistoryPoint contains: Etot, sigma, Wm, statev, R, T.
        """
        self.history = []

        # Initialize state if not provided
        if sv_init is None:
            # Determine nstatev from first block
            nstatev = self.blocks[0].nstatev if self.blocks else 1
            sv = StateVariablesM(nstatev=nstatev)
        else:
            # Use the provided state directly (no copy)
            sv = sv_init

        # Store initial state (lightweight copy for history)
        self.history.append(HistoryPoint.from_state(sv))

        Time = 0.0
        start = True

        # Process each block
        for block in self.blocks:
            control_type_int = block.get_control_type_int()
            corate_type_int = block.get_corate_type_int()

            # Initialize with zero increment to get initial tangent
            if start:
                self._initialize_umat(block, sv, Time)
                start = False

            # Process cycles
            for _ in range(block.ncycle):
                # Process steps
                for step in block.steps:
                    Time = self._solve_step(block, step, sv, Time,
                                            control_type_int, corate_type_int)

        return self.history

    def _initialize_umat(self, block: Block, sv: StateVariables, Time: float):
        """
        Initialize the UMAT by calling it with zero increment.

        This gets the initial tangent stiffness matrix.
        Modifies sv in-place.
        """
        DTime = 0.0
        sv.DEtot.fill(0.0)
        sv.Detot.fill(0.0)
        sv.DT = 0.0
        sv.DR.fill(0.0)
        np.fill_diagonal(sv.DR, 1.0)

        # Call UMAT (modifies sv in-place)
        self._call_umat(block, sv, Time, DTime)
        # Set _start values from current (C++ set_start pattern)
        # With DEtot=0, this just saves initial state without advancing
        sv.set_start(0)

    def _solve_step(self, block: Block, step: Step, sv: StateVariables,
                    Time: float, control_type_int: int,
                    corate_type_int: int) -> float:
        """
        Solve a single step with adaptive sub-incrementation.

        Modifies sv in-place and returns updated time.
        """
        cBC_meca = step.get_cBC_meca()
        nK = np.sum(cBC_meca)  # Number of stress-controlled components

        # Get targets (views, no copy)
        if isinstance(step, StepMeca):
            DEtot_target = step.DEtot_end
            Dsigma_target = step.Dsigma_end
        else:
            DEtot_target = np.zeros(6)
            Dsigma_target = np.zeros(6)

        # Thermomechanical targets
        DT_target = 0.0
        cBC_T = 0
        if isinstance(step, StepThermomeca):
            DT_target = step.DT_end
            cBC_T = step.get_cBC_T()

        # Sub-incrementation
        ninc = step.Dn_init
        tinc = 0.0  # Fraction of step completed

        while tinc < 1.0:
            # Calculate increment fraction
            Dtinc = min(1.0 / ninc, 1.0 - tinc)
            DTime = Dtinc * step.time

            # _start values are already set from previous set_start() or initialization
            # No explicit save needed here (C++ pattern)

            # Try to solve this increment
            converged = self._solve_increment(
                block, sv, Time, DTime, Dtinc,
                DEtot_target, Dsigma_target, DT_target,
                cBC_meca, cBC_T, nK, control_type_int, corate_type_int
            )

            if converged:
                # Accept increment
                tinc += Dtinc
                Time += DTime

                # Advance state: set _start from current + update strain/rotation
                # (C++ set_start pattern - must be called before recording history)
                sv.set_start(corate_type_int)

                # Store converged state (lightweight copy for history)
                self.history.append(HistoryPoint.from_state(sv))

                # Try to increase step size
                if ninc > step.Dn_mini:
                    ninc = max(step.Dn_mini, int(ninc * self.div_tnew_dt))
            else:
                # Reject increment, reset current TO _start values (C++ to_start pattern)
                sv.to_start()
                ninc = min(step.Dn_inc, int(ninc * self.mul_tnew_dt))

                if ninc >= step.Dn_inc:
                    raise RuntimeError(
                        f"Step failed to converge after reaching maximum "
                        f"sub-increments ({step.Dn_inc})"
                    )

        return Time

    def _solve_increment(self, block: Block, sv: StateVariables,
                         Time: float, DTime: float, Dtinc: float,
                         DEtot_target: np.ndarray, Dsigma_target: np.ndarray,
                         DT_target: float, cBC_meca: np.ndarray, cBC_T: int,
                         nK: int, control_type_int: int,
                         corate_type_int: int) -> bool:
        """
        Solve a single increment using Newton-Raphson iteration.

        Modifies sv in-place and returns convergence status.
        """
        # If fully strain controlled (nK == 0), single UMAT call suffices
        if nK == 0:
            self._apply_strain_increment(
                sv, Dtinc, DEtot_target, DT_target,
                control_type_int, corate_type_int, DTime
            )
            self._call_umat(block, sv, Time, DTime)
            # Strain advancement is done by set_start() in _solve_step after convergence
            return True

        # Mixed control: Newton-Raphson iteration
        # Initialize strain increment (in-place)
        sv.DEtot.fill(0.0)
        sv.Detot.fill(0.0)
        sv.DT = Dtinc * DT_target

        # Compute initial residual (reuse pre-allocated array)
        self._compute_residual(
            sv, Dtinc, DEtot_target, Dsigma_target, cBC_meca, control_type_int
        )
        error = norm(self._residual)

        compteur = 0
        while error > self.tol and compteur < self.max_iter:
            # Build Jacobian (reuse pre-allocated array)
            self._build_jacobian(sv, cBC_meca, control_type_int, corate_type_int)

            # Solve for correction
            np.copyto(self._Delta, np.linalg.solve(self._K, -self._residual))

            # Update strain (in-place)
            if control_type_int == 1:  # small_strain
                sv.DEtot += self._Delta
            elif control_type_int == 3:  # logarithmic
                sv.Detot += self._Delta
            else:
                sv.DEtot += self._Delta

            # Update kinematics for finite strain (modifies sv in-place)
            self._update_kinematics(sv, control_type_int, corate_type_int, DTime)

            # Reset state to start-of-increment values before UMAT call
            # This is critical for NR convergence: each UMAT call should start from
            # the same initial state (stress, statev, Wm) and only DEtot changes
            np.copyto(sv.sigma, sv.sigma_start)
            np.copyto(sv.statev, sv.statev_start)
            if isinstance(sv, (StateVariablesM, StateVariablesT)):
                np.copyto(sv.Wm, sv.Wm_start)

            # Call UMAT (modifies sv in-place)
            self._call_umat(block, sv, Time, DTime)

            # Compute new residual
            self._compute_residual(
                sv, Dtinc, DEtot_target, Dsigma_target, cBC_meca, control_type_int
            )
            error = norm(self._residual)
            compteur += 1

        # Strain advancement is done by set_start() in _solve_step after convergence
        return error <= self.tol

    def _compute_residual(self, sv: StateVariables, Dtinc: float,
                          DEtot_target: np.ndarray, Dsigma_target: np.ndarray,
                          cBC_meca: np.ndarray, control_type_int: int):
        """Compute the residual vector for Newton-Raphson (stores in self._residual)."""
        for k in range(6):
            if cBC_meca[k]:  # Stress controlled
                if control_type_int == 1:  # small_strain - Cauchy stress
                    self._residual[k] = sv.sigma[k] - sv.sigma_start[k] - Dtinc * Dsigma_target[k]
                elif control_type_int == 2:  # green_lagrange - PKII stress
                    self._residual[k] = sv.PKII[k] - sv.PKII_start[k] - Dtinc * Dsigma_target[k]
                elif control_type_int == 3:  # logarithmic - Cauchy stress
                    self._residual[k] = sv.sigma[k] - sv.sigma_start[k] - Dtinc * Dsigma_target[k]
                else:
                    self._residual[k] = sv.sigma[k] - sv.sigma_start[k] - Dtinc * Dsigma_target[k]
            else:  # Strain controlled
                if control_type_int == 3:  # logarithmic
                    self._residual[k] = self.lambda_solver * (sv.Detot[k] - Dtinc * DEtot_target[k])
                else:
                    self._residual[k] = self.lambda_solver * (sv.DEtot[k] - Dtinc * DEtot_target[k])

    def _build_jacobian(self, sv: StateVariables, cBC_meca: np.ndarray,
                        control_type_int: int, corate_type_int: int):
        """Build the Jacobian matrix (stores in self._K)."""
        if isinstance(sv, StateVariablesM):
            Lt = sv.Lt
        else:
            Lt = np.zeros((6, 6))

        # For small strain, directly use tangent
        # For finite strain, tangent transformations would be needed
        Lt_2_K(Lt, cBC_meca, self.lambda_solver, self._K)

    def _apply_strain_increment(self, sv: StateVariables, Dtinc: float,
                                DEtot_target: np.ndarray, DT_target: float,
                                control_type_int: int, corate_type_int: int,
                                DTime: float):
        """Apply strain increment for fully strain-controlled case (modifies sv in-place)."""
        sv.DT = Dtinc * DT_target
        sv.DR.fill(0.0)
        np.fill_diagonal(sv.DR, 1.0)

        if control_type_int == 1:  # small_strain
            np.copyto(sv.DEtot, DEtot_target)
            sv.DEtot *= Dtinc
        elif control_type_int == 3:  # logarithmic
            np.copyto(sv.Detot, DEtot_target)
            sv.Detot *= Dtinc
            self._update_kinematics(sv, control_type_int, corate_type_int, DTime)
        else:
            np.copyto(sv.DEtot, DEtot_target)
            sv.DEtot *= Dtinc

    def _update_kinematics(self, sv: StateVariables, control_type_int: int,
                           corate_type_int: int, DTime: float):
        """Update kinematic quantities for finite strain formulations (modifies sv in-place)."""
        if control_type_int == 1:  # small_strain
            # No kinematic update needed
            sv.DR.fill(0.0)
            np.fill_diagonal(sv.DR, 1.0)
            return

        if control_type_int == 2:  # green_lagrange
            # F from E and R (results written directly, carma copy=false used internally)
            sv.F0 = scc.ER_to_F(scc.v2t_strain(sv.Etot), sv.R, copy=False)
            sv.F1 = scc.ER_to_F(scc.v2t_strain(sv.Etot + sv.DEtot), sv.R, copy=False)
        elif control_type_int == 3:  # logarithmic
            # F from logarithmic strain and R
            sv.F0 = scc.eR_to_F(scc.v2t_strain(sv.etot), sv.R, copy=False)
            sv.F1 = scc.eR_to_F(scc.v2t_strain(sv.etot + sv.Detot), sv.R, copy=False)
            # Update Green-Lagrange from F
            GL = scc.Green_Lagrange(sv.F1, copy=False)
            GL_vec = scc.t2v_strain(GL, copy=False)
            np.copyto(sv.DEtot, GL_vec.ravel())  # Flatten in case of 2D return
            sv.DEtot -= sv.Etot

        # Compute objective rate quantities
        if DTime > 1e-12 and control_type_int > 1:
            # objective_rate returns (D, DR, Omega)
            D, DR, Omega = scc.objective_rate(
                self._get_corate_name(corate_type_int),
                sv.F0, sv.F1, DTime
            )
            np.copyto(sv.DR, DR)

    def _get_corate_name(self, corate_type_int: int) -> str:
        """Get corotational rate name from integer code."""
        corate_names = {0: 'jaumann', 1: 'green_naghdi', 2: 'logarithmic', 4: 'truesdell'}
        return corate_names.get(corate_type_int, 'jaumann')

    def _call_umat(self, block: Block, sv: StateVariables,
                   Time: float, DTime: float):
        """
        Call the UMAT via pybind11 binding (zero-copy version).

        Updates sv in-place using reshaped views - no array copies needed.
        The reshape operation creates views that share memory with the original arrays.
        """
        control_type_int = block.get_control_type_int()

        # Create reshaped views (no copy - shares memory with sv arrays)
        if control_type_int == 1:  # small_strain - use Green-Lagrange
            etot_view = sv.Etot.reshape(6, 1, order='F')
            Detot_view = sv.DEtot.reshape(6, 1, order='F')
        else:  # finite strain - use logarithmic strain
            etot_view = sv.etot.reshape(6, 1, order='F')
            Detot_view = sv.Detot.reshape(6, 1, order='F')

        sigma_view = sv.sigma.reshape(6, 1, order='F')
        F0_view = sv.F0.reshape(3, 3, 1, order='F')
        F1_view = sv.F1.reshape(3, 3, 1, order='F')
        DR_view = sv.DR.reshape(3, 3, 1, order='F')
        statev_view = sv.statev.reshape(-1, 1, order='F')

        # Props need to be copied (block.props may not be contiguous)
        nprops = len(block.props)
        if self._props_batch is None or self._props_batch.shape[0] != nprops:
            self._props_batch = np.zeros((nprops, 1), order='F')
        self._props_batch[:, 0] = block.props

        # Wm and Lt views (for StateVariablesM/T)
        if isinstance(sv, (StateVariablesM, StateVariablesT)):
            Wm_view = sv.Wm.reshape(4, 1, order='F')
            Lt_view = sv.Lt.reshape(6, 6, 1, order='F')
        else:
            Wm_view = self._Wm_batch
            Lt_view = self._Lt_batch

        # Temperature array for UMAT (single value reshaped for batch interface)
        temp_arr = np.array([sv.T], dtype=np.float64)

        # Call UMAT in-place - modifies sigma, statev, Wm, Lt through views
        self._umat_inplace(
            block.umat_name,
            etot_view, Detot_view,
            F0_view, F1_view,
            sigma_view, DR_view,
            self._props_batch, statev_view,
            Time, DTime,
            Wm_view, Lt_view,
            temp_arr,  # temp - pass actual temperature
            3,     # ndi
            1      # n_threads
        )
        # No copy needed - sv.sigma, sv.statev, sv.Wm, sv.Lt already modified!
        # NOTE: Strain totals are NOT updated here - they are updated after NR convergence
        # in _solve_increment to avoid accumulating strain during NR iterations

        # Update deformation for finite strain
        if block.get_control_type_int() > 1:
            np.copyto(sv.F0, sv.F1)
