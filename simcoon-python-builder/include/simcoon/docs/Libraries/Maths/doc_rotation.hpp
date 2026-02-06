#pragma once

namespace simcoon_docs {

constexpr auto fillR_angle = R"pbdoc(
    Generate a 3x3 rotation matrix for rotation around a single axis.

    Parameters
    ----------
    angle : float
        Rotation angle in radians.
    axis : int
        Axis of rotation: 1=x, 2=y, 3=z.
    active : bool, optional
        If True (default), active (alibi) rotation; if False, passive (alias) rotation.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        3x3 rotation matrix.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # 45 degree rotation around z-axis
        R = sim.fillR_angle(np.pi/4, 3)
        print(R)
)pbdoc";

constexpr auto fillR_euler = R"pbdoc(
    Generate a 3x3 rotation matrix from Euler angles.

    Parameters
    ----------
    psi : float
        First Euler angle in radians.
    theta : float
        Second Euler angle in radians.
    phi : float
        Third Euler angle in radians.
    active : bool
        If True, active (alibi) rotation; if False, passive (alias) rotation.
    conv : str
        Euler angle convention. Supported conventions:
        - Proper Euler angles: "zxz", "xyx", "yzy", "zyz", "xzx", "yxy"
        - Tait-Bryan angles: "xyz", "yzx", "zxy", "xzy", "zyx", "yxz"
        - Custom: "user" (uses axis_psi, axis_theta, axis_phi from parameter.hpp)
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        3x3 rotation matrix.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Rotation using ZXZ Euler angles
        psi, theta, phi = np.pi/6, np.pi/4, np.pi/3
        R = sim.fillR_euler(psi, theta, phi, True, "zxz")
        print(R)

        # Rotation using XYZ Tait-Bryan angles
        R_xyz = sim.fillR_euler(0.1, 0.2, 0.3, True, "xyz")
        print(R_xyz)
)pbdoc";

constexpr auto fillQS_angle = R"pbdoc(
    Generate a 6x6 rotation matrix for stress tensors in Voigt notation.

    Parameters
    ----------
    angle : float
        Rotation angle in radians.
    axis : int
        Axis of rotation: 1=x, 2=y, 3=z.
    active : bool, optional
        If True (default), active rotation; if False, passive rotation.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        6x6 stress rotation matrix (QS).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # 45 degree rotation matrix for stress
        QS = sim.fillQS_angle(np.pi/4, 3)
        print(QS)
)pbdoc";

constexpr auto fillQS_R = R"pbdoc(
    Generate a 6x6 rotation matrix for stress tensors from a 3x3 rotation matrix.

    Parameters
    ----------
    R : ndarray
        3x3 rotation matrix.
    active : bool, optional
        If True (default), active rotation; if False, passive rotation.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        6x6 stress rotation matrix (QS).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        R = sim.fillR_angle(np.pi/4, 3)
        QS = sim.fillQS_R(R)
        print(QS)
)pbdoc";

constexpr auto fillQE_angle = R"pbdoc(
    Generate a 6x6 rotation matrix for strain tensors in Voigt notation.

    Parameters
    ----------
    angle : float
        Rotation angle in radians.
    axis : int
        Axis of rotation: 1=x, 2=y, 3=z.
    active : bool, optional
        If True (default), active rotation; if False, passive rotation.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        6x6 strain rotation matrix (QE).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # 45 degree rotation matrix for strain
        QE = sim.fillQE_angle(np.pi/4, 3)
        print(QE)
)pbdoc";

constexpr auto fillQE_R = R"pbdoc(
    Generate a 6x6 rotation matrix for strain tensors from a 3x3 rotation matrix.

    Parameters
    ----------
    R : ndarray
        3x3 rotation matrix.
    active : bool, optional
        If True (default), active rotation; if False, passive rotation.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        6x6 strain rotation matrix (QE).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        R = sim.fillR_angle(np.pi/4, 3)
        QE = sim.fillQE_R(R)
        print(QE)
)pbdoc";

constexpr auto rotate_vec_R = R"pbdoc(
    Rotate a 3D vector using a rotation matrix.

    Parameters
    ----------
    input : ndarray
        3D vector to rotate.
    R : ndarray
        3x3 rotation matrix.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        Rotated 3D vector.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        v = np.array([1.0, 0.0, 0.0])
        R = sim.fillR_angle(np.pi/4, 3)
        v_rot = sim.rotate_vec_R(v, R)
        print(v_rot)
)pbdoc";

constexpr auto rotate_vec_angle = R"pbdoc(
    Rotate a 3D vector by an angle around a given axis.

    Parameters
    ----------
    input : ndarray
        3D vector to rotate.
    angle : float
        Rotation angle in radians.
    axis : int
        Axis of rotation: 1=x, 2=y, 3=z.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        Rotated 3D vector.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        v = np.array([1.0, 0.0, 0.0])
        v_rot = sim.rotate_vec_angle(v, np.pi/4, 3)
        print(v_rot)
)pbdoc";

constexpr auto rotate_mat_R = R"pbdoc(
    Rotate a 3x3 matrix using a rotation matrix.

    Parameters
    ----------
    input : ndarray
        3x3 matrix to rotate.
    R : ndarray
        3x3 rotation matrix.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        Rotated 3x3 matrix: R * input * R^T.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        m = np.eye(3)
        R = sim.fillR_angle(np.pi/4, 3)
        m_rot = sim.rotate_mat_R(m, R)
        print(m_rot)
)pbdoc";

constexpr auto rotate_mat_angle = R"pbdoc(
    Rotate a 3x3 matrix by an angle around a given axis.

    Parameters
    ----------
    input : ndarray
        3x3 matrix to rotate.
    angle : float
        Rotation angle in radians.
    axis : int
        Axis of rotation: 1=x, 2=y, 3=z.
    copy : bool, optional
        If True (default), return a copy of the result.

    Returns
    -------
    ndarray
        Rotated 3x3 matrix.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        m = np.eye(3)
        m_rot = sim.rotate_mat_angle(m, np.pi/4, 3)
        print(m_rot)
)pbdoc";

} // namespace simcoon_docs
