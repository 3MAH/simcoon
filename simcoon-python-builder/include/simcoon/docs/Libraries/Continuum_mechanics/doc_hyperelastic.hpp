#pragma once

namespace simcoon_docs {

constexpr auto isochoric_invariants = R"pbdoc(
    Provides the isochoric strain invariants from the left Cauchy-Green deformation tensor \( \mathbf{b} \).

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 3x3 matrix representing the left Cauchy-Green deformation tensor \( \mathbf{b} \).
    J : double, optional
        The determinant of the transformation gradient \( \mathbf{F} \). Default is 0.
    copy : bool, optional
        Whether to copy the input data. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A column vector of dimension 3 containing the three isochoric invariants.

    Notes
    -----
    The isochoric invariants are defined as:

    .. math::

        \bar{I}_1 = \textrm{tr} \bar{\mathbf{b}} \\
        \bar{I}_2 = \frac{1}{2} \left( \left(\textrm{tr} \bar{\mathbf{b}} \right)^2 - \textrm{tr} \bar{\mathbf{b}}^2 \right) \\
        \bar{I}_3 = \textrm{det} \bar{\mathbf{b}} = 1

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        b = np.random.rand(3, 3)
        J = np.linalg.det(b)
        invariants = sim.isochoric_invariants(b, J)
        print(invariants)
)pbdoc";

constexpr auto isochoric_pstretch = R"pbdoc(
    Provides the isochoric principal stretches from the left Cauchy-Green tensor \( \mathbf{b} \) or the Eulerian stretch tensor \( \mathbf{v} \).

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 3x3 matrix representing the left Cauchy-Green tensor \( \mathbf{b} \) or the Eulerian stretch tensor \( \mathbf{v} \).
    input_tensor : str, optional
        A string indicating the type of input tensor ("b" for \( \mathbf{b} \) or "V" for \( \mathbf{v} \)). Default is "V".
    J : double, optional
        The determinant of the transformation gradient \( \mathbf{F} \). Default is 0.
    copy : bool, optional
        Whether to copy the input data. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A column vector of dimension 3 containing the three isochoric principal stretches.

    Notes
    -----
    The isochoric principal stretches are defined as:

    .. math::

        \bar{\lambda}_i = J^{-1/3} \lambda_i

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        b = np.random.rand(3, 3)
        J = np.linalg.det(b)
        pstretch = sim.isochoric_pstretch(b, "b", J)
        print(pstretch)
)pbdoc";

constexpr auto tau_iso_hyper_invariants = R"pbdoc(
    Computes the isochoric part of the Kirchhoff stress tensor using the invariants-based formulation.

    For an incompressible hyperelastic material with strain energy density :math:`W(\bar{I}_1, \bar{I}_2)`,
    the isochoric Kirchhoff stress is:

    .. math::

        \boldsymbol{\tau}_{iso} = 2 \frac{\partial W}{\partial \bar{I}_1} \text{dev}(\bar{\mathbf{b}}) 
        + 2 \frac{\partial W}{\partial \bar{I}_2} \left( \text{tr}(\bar{\mathbf{b}}) \text{dev}(\bar{\mathbf{b}}) - \text{dev}(\bar{\mathbf{b}}^2) \right)

    Parameters
    ----------
    dWdI_1_bar : float
        Derivative of the strain energy with respect to the first isochoric invariant.
    dWdI_2_bar : float
        Derivative of the strain energy with respect to the second isochoric invariant.
    input : numpy.ndarray
        Left Cauchy-Green tensor b (3x3 matrix).
    J : float, optional
        Jacobian determinant. If 0 or not provided, it is computed from the input tensor.

    Returns
    -------
    numpy.ndarray
        Isochoric Kirchhoff stress tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Neo-Hookean material: W = (mu/2)(I1_bar - 3)
        # dW/dI1_bar = mu/2, dW/dI2_bar = 0
        mu = 1000.0  # Shear modulus in Pa
        
        F = np.array([[1.2, 0.0, 0.0],
                      [0.0, 1.1, 0.0],
                      [0.0, 0.0, 0.9]])
        b = F @ F.T
        
        tau_iso = sim.tau_iso_hyper_invariants(mu/2, 0.0, b)
        print("Isochoric Kirchhoff stress:")
        print(tau_iso)
)pbdoc";

constexpr auto sigma_iso_hyper_invariants = R"pbdoc(
    Computes the isochoric part of the Cauchy stress tensor using the invariants-based formulation.

    The isochoric Cauchy stress is related to the isochoric Kirchhoff stress by:

    .. math::

        \boldsymbol{\sigma}_{iso} = \frac{1}{J} \boldsymbol{\tau}_{iso}

    Parameters
    ----------
    dWdI_1_bar : float
        Derivative of the strain energy with respect to the first isochoric invariant.
    dWdI_2_bar : float
        Derivative of the strain energy with respect to the second isochoric invariant.
    input : numpy.ndarray
        Left Cauchy-Green tensor b (3x3 matrix).
    J : float, optional
        Jacobian determinant. If 0 or not provided, it is computed from the input tensor.

    Returns
    -------
    numpy.ndarray
        Isochoric Cauchy stress tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Neo-Hookean material: W = (mu/2)(I1_bar - 3)
        mu = 1000.0  # Shear modulus in Pa
        
        F = np.array([[1.2, 0.0, 0.0],
                      [0.0, 1.1, 0.0],
                      [0.0, 0.0, 0.9]])
        b = F @ F.T
        
        sigma_iso = sim.sigma_iso_hyper_invariants(mu/2, 0.0, b)
        print("Isochoric Cauchy stress:")
        print(sigma_iso)
)pbdoc";

constexpr auto tau_vol_hyper = R"pbdoc(
    Computes the volumetric part of the Kirchhoff stress tensor.

    For a hyperelastic material with volumetric strain energy :math:`U(J)`,
    the volumetric Kirchhoff stress is:

    .. math::

        \boldsymbol{\tau}_{vol} = J \frac{dU}{dJ} \mathbf{I}

    Parameters
    ----------
    dUdJ : float
        Derivative of the volumetric strain energy with respect to J.
    input : numpy.ndarray
        Left Cauchy-Green tensor b (3x3 matrix), used to compute J if not provided.
    J : float, optional
        Jacobian determinant. If 0 or not provided, it is computed from the input tensor.

    Returns
    -------
    numpy.ndarray
        Volumetric Kirchhoff stress tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Volumetric penalty: U(J) = (kappa/2)(J-1)^2
        # dU/dJ = kappa(J-1)
        kappa = 10000.0  # Bulk modulus in Pa
        
        F = np.array([[1.2, 0.0, 0.0],
                      [0.0, 1.1, 0.0],
                      [0.0, 0.0, 0.9]])
        J = np.linalg.det(F)
        b = F @ F.T
        
        dUdJ = kappa * (J - 1)
        tau_vol = sim.tau_vol_hyper(dUdJ, b, J)
        print("Volumetric Kirchhoff stress:")
        print(tau_vol)
)pbdoc";

constexpr auto sigma_vol_hyper = R"pbdoc(
    Computes the volumetric part of the Cauchy stress tensor.

    The volumetric Cauchy stress is:

    .. math::

        \boldsymbol{\sigma}_{vol} = \frac{dU}{dJ} \mathbf{I}

    Parameters
    ----------
    dUdJ : float
        Derivative of the volumetric strain energy with respect to J.
    input : numpy.ndarray
        Left Cauchy-Green tensor b (3x3 matrix), used to compute J if not provided.
    J : float, optional
        Jacobian determinant. If 0 or not provided, it is computed from the input tensor.

    Returns
    -------
    numpy.ndarray
        Volumetric Cauchy stress tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Volumetric penalty: U(J) = (kappa/2)(J-1)^2
        # dU/dJ = kappa(J-1)
        kappa = 10000.0  # Bulk modulus in Pa
        
        F = np.array([[1.2, 0.0, 0.0],
                      [0.0, 1.1, 0.0],
                      [0.0, 0.0, 0.9]])
        J = np.linalg.det(F)
        b = F @ F.T
        
        dUdJ = kappa * (J - 1)
        sigma_vol = sim.sigma_vol_hyper(dUdJ, b, J)
        print("Volumetric Cauchy stress:")
        print(sigma_vol)
)pbdoc";

} // namespace simcoon_docs