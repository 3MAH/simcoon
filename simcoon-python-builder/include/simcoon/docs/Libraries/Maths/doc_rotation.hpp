#pragma once

namespace simcoon_docs {

constexpr auto CppRotation_class = R"pbdoc(
    Internal C++ rotation backend using unit quaternions (scalar-last).

    End users should use ``simcoon.Rotation`` instead, which inherits from
    ``scipy.spatial.transform.Rotation`` and delegates mechanics operations
    to this class.
)pbdoc";

constexpr auto CppRotation_from_quat = R"pbdoc(
    Create a rotation from a quaternion in scalar-last convention.

    Parameters
    ----------
    quat : numpy.ndarray
        A 1D array of 4 elements [qx, qy, qz, qw] (scalar-last).

    Returns
    -------
    _CppRotation
        A rotation object.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as smc

        q = np.array([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])  # 90deg around z
        r = smc.Rotation.from_quat(q)
)pbdoc";

constexpr auto as_voigt_stress_rotation = R"pbdoc(
    Get the 6x6 rotation matrix for stress tensors in Voigt notation.

    The stress rotation matrix :math:`Q_\sigma` satisfies
    :math:`\boldsymbol{\sigma}' = Q_\sigma \boldsymbol{\sigma}`.

    Parameters
    ----------
    active : bool, optional
        If True (default), returns the active (alibi) rotation matrix.
        If False, returns the passive (alias) rotation matrix.

    Returns
    -------
    numpy.ndarray
        A (6, 6) rotation matrix for stress tensors.
        For batch rotations, returns (N, 6, 6).

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/4, 3)
        QS = r.as_voigt_stress_rotation()  # (6, 6)

        # Batch: N rotations
        rots = smc.Rotation.random(100)
        QS_batch = rots.as_voigt_stress_rotation()  # (100, 6, 6)
)pbdoc";

constexpr auto as_voigt_strain_rotation = R"pbdoc(
    Get the 6x6 rotation matrix for strain tensors in Voigt notation.

    The strain rotation matrix :math:`Q_\varepsilon` satisfies
    :math:`\boldsymbol{\varepsilon}' = Q_\varepsilon \boldsymbol{\varepsilon}`.

    Parameters
    ----------
    active : bool, optional
        If True (default), returns the active (alibi) rotation matrix.
        If False, returns the passive (alias) rotation matrix.

    Returns
    -------
    numpy.ndarray
        A (6, 6) rotation matrix for strain tensors.
        For batch rotations, returns (N, 6, 6).

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/4, 3)
        QE = r.as_voigt_strain_rotation()  # (6, 6)
)pbdoc";

constexpr auto apply_tensor = R"pbdoc(
    Apply the rotation to a 3x3 tensor: :math:`R \cdot T \cdot R^T`.

    Parameters
    ----------
    m : numpy.ndarray
        A (3, 3) tensor. For batch rotations, shape (3, 3, N).
    inverse : bool, optional
        If True, applies the inverse rotation. Default is False.

    Returns
    -------
    numpy.ndarray
        The rotated tensor, same shape as input.

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/2, 3)
        T = np.diag([1.0, 2.0, 3.0])
        T_rot = r.apply_tensor(T)
)pbdoc";

constexpr auto apply_stress = R"pbdoc(
    Apply the rotation to a stress vector in Voigt notation.

    Voigt convention: :math:`[\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{12}, \sigma_{13}, \sigma_{23}]`.

    Parameters
    ----------
    sigma : numpy.ndarray
        A 1D array of 6 elements (single rotation) or shape (6, N) for batch.
    active : bool, optional
        If True (default), active rotation. If False, passive rotation.

    Returns
    -------
    numpy.ndarray
        The rotated stress vector, same shape as input.

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/4, 3)
        sigma = np.array([100.0, 50.0, 25.0, 10.0, 5.0, 2.0])
        sigma_rot = r.apply_stress(sigma)
)pbdoc";

constexpr auto apply_strain = R"pbdoc(
    Apply the rotation to a strain vector in Voigt notation.

    Voigt convention: :math:`[\varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33}, 2\varepsilon_{12}, 2\varepsilon_{13}, 2\varepsilon_{23}]`.

    Parameters
    ----------
    epsilon : numpy.ndarray
        A 1D array of 6 elements (single rotation) or shape (6, N) for batch.
    active : bool, optional
        If True (default), active rotation. If False, passive rotation.

    Returns
    -------
    numpy.ndarray
        The rotated strain vector, same shape as input.

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/4, 3)
        epsilon = np.array([0.01, -0.005, -0.005, 0.002, 0.001, 0.0])
        epsilon_rot = r.apply_strain(epsilon)
)pbdoc";

constexpr auto apply_stiffness = R"pbdoc(
    Apply the rotation to a 6x6 stiffness matrix.

    Computes :math:`L' = Q_\sigma \cdot L \cdot Q_\sigma^T`.

    Parameters
    ----------
    L : numpy.ndarray
        A (6, 6) stiffness matrix (single) or (6, 6, N) for batch.
    active : bool, optional
        If True (default), active rotation. If False, passive rotation.

    Returns
    -------
    numpy.ndarray
        The rotated stiffness matrix, same shape as input.

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/4, 3)
        L = smc.L_iso(np.array([210000, 0.3]), "Enu")
        L_rot = r.apply_stiffness(L)
)pbdoc";

constexpr auto apply_compliance = R"pbdoc(
    Apply the rotation to a 6x6 compliance matrix.

    Computes :math:`M' = Q_\varepsilon \cdot M \cdot Q_\varepsilon^T`.

    Parameters
    ----------
    M : numpy.ndarray
        A (6, 6) compliance matrix (single) or (6, 6, N) for batch.
    active : bool, optional
        If True (default), active rotation. If False, passive rotation.

    Returns
    -------
    numpy.ndarray
        The rotated compliance matrix, same shape as input.

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_axis_angle(np.pi/4, 3)
        M = smc.M_iso(np.array([210000, 0.3]), "Enu")
        M_rot = r.apply_compliance(M)
)pbdoc";

constexpr auto apply_strain_concentration = R"pbdoc(
    Apply the rotation to a 6x6 strain concentration tensor.

    Computes :math:`A' = Q_\varepsilon \cdot A \cdot Q_\sigma^T`.

    Parameters
    ----------
    A : numpy.ndarray
        A (6, 6) strain concentration tensor (single) or (6, 6, N) for batch.
    active : bool, optional
        If True (default), active rotation. If False, passive rotation.

    Returns
    -------
    numpy.ndarray
        The rotated strain concentration tensor, same shape as input.
)pbdoc";

constexpr auto apply_stress_concentration = R"pbdoc(
    Apply the rotation to a 6x6 stress concentration tensor.

    Computes :math:`B' = Q_\sigma \cdot B \cdot Q_\varepsilon^T`.

    Parameters
    ----------
    B : numpy.ndarray
        A (6, 6) stress concentration tensor (single) or (6, 6, N) for batch.
    active : bool, optional
        If True (default), active rotation. If False, passive rotation.

    Returns
    -------
    numpy.ndarray
        The rotated stress concentration tensor, same shape as input.
)pbdoc";

constexpr auto dR_drotvec = R"pbdoc(
    Compute the derivatives of the rotation matrix w.r.t. rotation vector components.

    Uses the exact differentiation of the Rodrigues formula
    (Gallego & Yezzi, J. Math. Imaging Vis., 2015).

    Returns
    -------
    numpy.ndarray
        A (3, 3, 3) array where ``result[k]`` is :math:`\partial R / \partial \omega_k`.

    Examples
    --------
    .. code-block:: python

        import simcoon as smc
        import numpy as np

        r = smc.Rotation.from_rotvec([0.1, 0.2, 0.3])
        dR = r.dR_drotvec()  # (3, 3, 3)
        # dR[0] = dR/d(omega_0), dR[1] = dR/d(omega_1), dR[2] = dR/d(omega_2)
)pbdoc";

constexpr auto dR_drotvec_free = R"pbdoc(
    Compute the derivatives of R(omega) w.r.t. rotation vector components.

    Free function version — takes a rotation vector directly.

    Parameters
    ----------
    rotvec : numpy.ndarray
        A 1D array of 3 elements: rotation vector (axis * angle).

    Returns
    -------
    numpy.ndarray
        A (3, 3, 3) array where ``result[k]`` is :math:`\partial R / \partial \omega_k`.

    References
    ----------
    Gallego & Yezzi, "A Compact Formula for the Derivative of a 3-D
    Rotation in Exponential Coordinates", J. Math. Imaging Vis., 2015.
)pbdoc";

} // namespace simcoon_docs
