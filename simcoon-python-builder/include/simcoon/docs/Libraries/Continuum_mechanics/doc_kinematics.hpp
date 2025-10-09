#pragma once

namespace simcoon_docs {

constexpr auto ER_to_F = R"pbdoc(
    Provides the transformation gradient from the Green-Lagrange strain and the rotation.

    Parameters
    ----------
    E : pybind11::array_t<double>
        Green-Lagrange strain tensor (3x3 matrix).
    R : pybind11::array_t<double>
        Rotation tensor (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        E = np.random.rand(3, 3)
        R = np.eye(3)
        F = sim.ER_to_F(E, R)
        print(F)
)pbdoc";

constexpr auto eR_to_F = R"pbdoc(
    Provides the transformation gradient from the logarithmic strain and the rotation.

    Parameters
    ----------
    e : pybind11::array_t<double>
        Logarithmic strain tensor (3x3 matrix).
    R : pybind11::array_t<double>
        Rotation tensor (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        e = np.random.rand(3, 3)
        R = np.eye(3)
        F = sim.eR_to_F(e, R)
        print(F)
)pbdoc";

constexpr auto G_UdX = R"pbdoc(
    Computes the gradient of displacement (Lagrangian) from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Gradient of displacement (Lagrangian) (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        GradU = sim.G_UdX(F)
        print(GradU)
)pbdoc";

constexpr auto G_Udx = R"pbdoc(
    Computes the gradient of displacement (Eulerian) from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Gradient of displacement (Eulerian) (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        gradU = sim.G_Udx(F)
        print(gradU)
)pbdoc";

constexpr auto R_Cauchy_Green = R"pbdoc(
    Computes the Right Cauchy-Green tensor from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Right Cauchy-Green tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        C = sim.R_Cauchy_Green(F)
        print(C)
)pbdoc";

constexpr auto L_Cauchy_Green = R"pbdoc(
    Computes the Left Cauchy-Green tensor from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Left Cauchy-Green tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        B = sim.L_Cauchy_Green(F)
        print(B)
)pbdoc";

constexpr auto RU_decomposition = R"pbdoc(
    Computes the RU decomposition of the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    tuple
        Rotation matrix (3x3) and right stretch tensor (3x3).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        R, U = sim.RU_decomposition(F)
        print(R, U)
)pbdoc";

constexpr auto VR_decomposition = R"pbdoc(
    Computes the VR decomposition of the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    tuple
        Left stretch tensor (3x3) and rotation matrix (3x3).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        V, R = sim.VR_decomposition(F)
        print(V, R)
)pbdoc";

constexpr auto Inv_X = R"pbdoc(
    Computes the invariants of a symmetric tensor.

    Parameters
    ----------
    X : pybind11::array_t<double>
        Symmetric tensor (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Vector containing the three invariants of the tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        X = np.random.rand(3, 3)
        invariants = sim.Inv_X(X)
        print(invariants)
)pbdoc";

constexpr auto Cauchy = R"pbdoc(
    Computes the Cauchy tensor from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Cauchy tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        b = sim.Cauchy(F)
        print(b)
)pbdoc";

constexpr auto Green_Lagrange = R"pbdoc(
    Computes the Green-Lagrange tensor from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Green-Lagrange tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        E = sim.Green_Lagrange(F)
        print(E)
)pbdoc";

constexpr auto Euler_Almansi = R"pbdoc(
    Computes the Euler-Almansi tensor from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Euler-Almansi tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        A = sim.Euler_Almansi(F)
        print(A)
)pbdoc";

constexpr auto Log_strain = R"pbdoc(
    Computes the logarithmic strain tensor from the transformation gradient.

    Parameters
    ----------
    F : pybind11::array_t<double>
        Transformation gradient (3x3 matrix).

    Returns
    -------
    pybind11::array_t<double>
        Logarithmic strain tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F = np.random.rand(3, 3)
        e = sim.Log_strain(F)
        print(e)
)pbdoc";

constexpr auto finite_L = R"pbdoc(
    Computes the Eulerian velocity tensor from the transformation gradient at two different times.

    Parameters
    ----------
    F0 : pybind11::array_t<double>
        Transformation gradient at time t0 (3x3 matrix).
    F1 : pybind11::array_t<double>
        Transformation gradient at time t1 (3x3 matrix).
    DTime : double
        Time difference (t1 - t0).

    Returns
    -------
    pybind11::array_t<double>
        Eulerian velocity tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F0 = np.random.rand(3, 3)
        F1 = np.random.rand(3, 3)
        DTime = 0.1
        L = sim.finite_L(F0, F1, DTime)
        print(L)
)pbdoc";

constexpr auto finite_D = R"pbdoc(
    Computes the Eulerian symmetric rate tensor from the transformation gradient at two different times.

    Parameters
    ----------
    F0 : pybind11::array_t<double>
        Transformation gradient at time t0 (3x3 matrix).
    F1 : pybind11::array_t<double>
        Transformation gradient at time t1 (3x3 matrix).
    DTime : double
        Time difference (t1 - t0).

    Returns
    -------
    pybind11::array_t<double>
        Eulerian symmetric rate tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F0 = np.random.rand(3, 3)
        F1 = np.random.rand(3, 3)
        DTime = 0.1
        D = sim.finite_D(F0, F1, DTime)
        print(D)
)pbdoc";

constexpr auto finite_W = R"pbdoc(
    Computes the Eulerian antisymmetric spin tensor from the transformation gradient at two different times.

    Parameters
    ----------
    F0 : pybind11::array_t<double>
        Transformation gradient at time t0 (3x3 matrix).
    F1 : pybind11::array_t<double>
        Transformation gradient at time t1 (3x3 matrix).
    DTime : double
        Time difference (t1 - t0).

    Returns
    -------
    pybind11::array_t<double>
        Eulerian antisymmetric spin tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F0 = np.random.rand(3, 3)
        F1 = np.random.rand(3, 3)
        DTime = 0.1
        W = sim.finite_W(F0, F1, DTime)
        print(W)
)pbdoc";

constexpr auto finite_Omega = R"pbdoc(
    Computes the rigid-body rotation spin tensor from the transformation gradient at two different times.

    Parameters
    ----------
    F0 : pybind11::array_t<double>
        Transformation gradient at time t0 (3x3 matrix).
    F1 : pybind11::array_t<double>
        Transformation gradient at time t1 (3x3 matrix).
    DTime : double
        Time difference (t1 - t0).

    Returns
    -------
    pybind11::array_t<double>
        Rigid-body rotation spin tensor (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        F0 = np.random.rand(3, 3)
        F1 = np.random.rand(3, 3)
        DTime = 0.1
        Omega = sim.finite_Omega(F0, F1, DTime)
        print(Omega)
)pbdoc";

constexpr auto finite_DQ = R"pbdoc(
    Computes the Hughes-Winget approximation of the increment of rotation or transformation.

    Parameters
    ----------
    Omega0 : pybind11::array_t<double>
        Spin/velocity at time t0 (3x3 matrix).
    Omega1 : pybind11::array_t<double>
        Spin/velocity at time t1 (3x3 matrix).
    DTime : double
        Time difference (t1 - t0).

    Returns
    -------
    pybind11::array_t<double>
        Increment of rotation or transformation (3x3 matrix).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        Omega0 = np.random.rand(3, 3)
        Omega1 = np.random.rand(3, 3)
        DTime = 0.1
        DQ = sim.finite_DQ(Omega0, Omega1, DTime)
        print(DQ)
)pbdoc";

} // namespace simcoon_docs