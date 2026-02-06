#pragma once

namespace simcoon_docs {

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
        R = sim.Rotation.from_axis_angle(np.pi/4, 3).as_matrix()
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
        R = sim.Rotation.from_axis_angle(np.pi/4, 3).as_matrix()
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
