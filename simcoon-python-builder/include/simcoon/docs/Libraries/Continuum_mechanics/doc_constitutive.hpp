#pragma once

namespace simcoon_docs {

constexpr auto Ireal = R"pbdoc(
    Returns the fourth order identity tensor :math:`I_{real}` written in Voigt notation.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The fourth order identity tensor :math:`I_{real}` in Voigt notation.

    Notes
    -----
    The tensor is defined as:

    .. math::

        I_{real} = \begin{bmatrix}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for real
        Ireal = sim.Ireal()
        print(Ireal)
)pbdoc";

constexpr auto Ivol = R"pbdoc(
    Returns the volumetric part of the identity tensor :math:`I_{vol}` written in Voigt notation.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The volumetric part of the identity tensor :math:`I_{vol}` in Voigt notation.

    Notes
    -----
    The tensor is defined as:

    .. math::

        I_{vol} = \begin{bmatrix}
        \frac{1}{3} & \frac{1}{3} & \frac{1}{3} & 0 & 0 & 0 \\
        \frac{1}{3} & \frac{1}{3} & \frac{1}{3} & 0 & 0 & 0 \\
        \frac{1}{3} & \frac{1}{3} & \frac{1}{3} & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Ivol
        Ivol = sim.Ivol()
        print(Ivol)
)pbdoc";

constexpr auto Idev = R"pbdoc(
    Returns the deviatoric part of the identity tensor :math:`I_{dev}` written in Voigt notation.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The deviatoric part of the identity tensor :math:`I_{dev}` in Voigt notation.

    Notes
    -----
    The tensor is defined as:

    .. math::

        I_{dev} = \begin{bmatrix}
        \frac{2}{3} & -\frac{1}{3} & -\frac{1}{3} & 0 & 0 & 0 \\
        -\frac{1}{3} & \frac{2}{3} & -\frac{1}{3} & 0 & 0 & 0 \\
        -\frac{1}{3} & -\frac{1}{3} & \frac{2}{3} & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Idev
        Idev = sim.Idev()
        print(Idev)
)pbdoc";

constexpr auto Ireal2 = R"pbdoc(
    Returns the fourth order identity tensor :math:`\widehat{I}` written in Voigt notation.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The fourth order identity tensor :math:`\widehat{I}` in Voigt notation.

    Notes
    -----
    The tensor is defined as:

    .. math::

        \widehat{I} = \begin{bmatrix}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Ireal2
        Ireal2 = sim.Ireal2()
        print(Ireal2)
)pbdoc";

constexpr auto Idev2 = R"pbdoc(
    Returns the deviatoric part of the identity tensor :math:`\widehat{I}_{dev}` written in Voigt notation.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The deviatoric part of the identity tensor :math:`\widehat{I}_{dev}` in Voigt notation.

    Notes
    -----
    The tensor is defined as:

    .. math::

        \widehat{I}_{dev} = \begin{bmatrix}
        \frac{2}{3} & -\frac{1}{3} & -\frac{1}{3} & 0 & 0 & 0 \\
        -\frac{1}{3} & \frac{2}{3} & -\frac{1}{3} & 0 & 0 & 0 \\
        -\frac{1}{3} & -\frac{1}{3} & \frac{2}{3} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Idev2
        Idev2 = sim.Idev2()
        print(Idev2)
)pbdoc";

constexpr auto Ith = R"pbdoc(
    Returns the expansion vector :math:`I_{th}`.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the vector is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The expansion vector :math:`I_{th}`.

    Notes
    -----
    The vector is defined as:

    .. math::

        I_{th} = \begin{bmatrix}
        1 \\
        1 \\
        1 \\
        0 \\
        0 \\
        0
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Ith
        Ith = sim.Ith()
        print(Ith)
)pbdoc";

constexpr auto Ir2 = R"pbdoc(
    Returns the operator from stress to strain in Voigt notation :math:`I_{r2}`.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the operator is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The operator :math:`I_{r2}` in Voigt notation.

    Notes
    -----
    The operator is defined as:

    .. math::

        I_{r2} = \begin{bmatrix}
        1 \\
        1 \\
        1 \\
        2 \\
        2 \\
        2
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Ir2
        Ir2 = sim.Ir2()
        print(Ir2)
)pbdoc";

constexpr auto Ir05 = R"pbdoc(
    Returns the operator from strain to stress in Voigt notation :math:`I_{r05}`.

    Parameters
    ----------
    copy : bool, optional
        If true, a copy of the operator is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The operator :math:`I_{r05}` in Voigt notation.

    Notes
    -----
    The operator is defined as:

    .. math::

        I_{r05} = \begin{bmatrix}
        1 \\
        1 \\
        1 \\
        0.5 \\
        0.5 \\
        0.5
        \end{bmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Ir05
        Ir05 = sim.Ir05()
        print(Ir05)
)pbdoc";

constexpr auto L_iso = R"pbdoc(
    Provides the elastic stiffness tensor for an isotropic material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the material properties (e.g., Lamé coefficients).
    conv : str
        A string specifying the convention used for the material properties. Possible values include:
        'Enu', 'nuE', 'Kmu', 'muK', 'KG', 'GK', 'lambdamu', 'mulambda', 'lambdaG', 'Glambda'.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 stiffness matrix in Voigt notation.

    Notes
    -----
    The stiffness tensor is computed based on the provided material properties and convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for L_iso
        props = np.array([210000, 0.3])
        L_iso = sim.L_iso(props, "Enu")
        print(L_iso)
)pbdoc";

constexpr auto M_iso = R"pbdoc(
    Provides the elastic compliance tensor for an isotropic material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the material properties (e.g., Lamé coefficients).
    conv : str
        A string specifying the convention used for the material properties. Possible values include:
        'Enu', 'nuE', 'Kmu', 'muK', 'KG', 'GK', 'lambdamu', 'mulambda', 'lambdaG', 'Glambda'.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 compliance matrix in Voigt notation.

    Notes
    -----
    The compliance tensor is computed based on the provided material properties and convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for M_iso
        props = np.array([210000, 0.3])    
        M_iso = sim.M_iso(props, "Enu")
        print(M_iso)
)pbdoc";

constexpr auto L_cubic = R"pbdoc(
    Provides the elastic stiffness tensor for a cubic material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the stiffness coefficients (C11, C12, C44).
    conv : str
        A string specifying the convention used for the material properties. Possible values include:
        'Cii'.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 stiffness matrix in Voigt notation.

    Notes
    -----
    The stiffness tensor is computed based on the provided stiffness coefficients and convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for L_cubic
        props_cubic = np.array([100000, 50000, 30000])
        L_cubic = sim.L_cubic(props_cubic, "Cii")
        print(L_cubic)
)pbdoc";

constexpr auto M_cubic = R"pbdoc(
    Provides the elastic compliance tensor for a cubic material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the stiffness coefficients (C11, C12, C44).
    conv : str
        A string specifying the convention used for the material properties. Possible values include:
        'Cii'.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 compliance matrix in Voigt notation.

    Notes
    -----
    The compliance tensor is computed based on the provided stiffness coefficients and convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for M_cubic
        props_cubic = np.array([100000, 50000, 30000])    
        M_cubic = sim.M_cubic(props_cubic, "Cii")
        print(M_cubic)
)pbdoc";

constexpr auto L_ortho = R"pbdoc(
    Provides the elastic stiffness tensor for an orthotropic material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the stiffness coefficients or material parameters.
    conv : str
        A string specifying the convention used for the material properties. Possible values include:
        'Cii', 'EnuG'.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 stiffness matrix in Voigt notation.

    Notes
    -----
    The stiffness tensor is computed based on the provided coefficients or material parameters and convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for L_ortho
        props_ortho = np.array([100000, 50000, 30000, 0.3, 0.3, 0.3, 40000, 40000, 40000])
        L_ortho = sim.L_ortho(props_ortho, "EnuG")
        print(L_ortho)
)pbdoc";

constexpr auto M_ortho = R"pbdoc(
    Provides the elastic compliance tensor for an orthotropic material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the stiffness coefficients or material parameters.
    conv : str
        A string specifying the convention used for the material properties. Possible values include:
        'Cii', 'EnuG'.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 compliance matrix in Voigt notation.

    Notes
    -----
    The compliance tensor is computed based on the provided coefficients or material parameters and convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for M_ortho
        props_ortho = np.array([100000, 50000, 30000, 0.3, 0.3, 0.3, 40000, 40000, 40000])    
        M_ortho = sim.M_ortho(props_ortho, "EnuG")
        print(M_ortho)
)pbdoc";

constexpr auto L_isotrans = R"pbdoc(
    Provides the elastic stiffness tensor for an isotropic transverse material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the material properties (e.g., EL, ET, nuTL, nuTT, GLT).
    axis : int
        The axis of symmetry.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 stiffness matrix in Voigt notation.

    Notes
    -----
    The stiffness tensor is computed based on the provided material properties and axis of symmetry.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for L_isotrans
        props_isotrans = np.array([210000, 70000, 0.3, 0.3, 50000])
        axis = 1
        L_isotrans = sim.L_isotrans(props_isotrans, axis)
        print(L_isotrans)
)pbdoc";

constexpr auto M_isotrans = R"pbdoc(
    Provides the elastic compliance tensor for an isotropic transverse material.

    Parameters
    ----------
    props : pybind11::array_t<double>
        A 1D array containing the material properties (e.g., EL, ET, nuTL, nuTT, GLT).
    axis : int
        The axis of symmetry.
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 compliance matrix in Voigt notation.

    Notes
    -----
    The compliance tensor is computed based on the provided material properties and axis of symmetry.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for M_isotrans
        props_isotrans = np.array([210000, 70000, 0.3, 0.3, 50000])
        axis = 1    
        M_isotrans = sim.M_isotrans(props_isotrans, axis)
        print(M_isotrans)
)pbdoc";

constexpr auto H_iso = R"pbdoc(
    Provides the viscous tensor for an isotropic material.

    Parameters
    ----------
    props_py : pybind11::array_t<double>
        A 1D array containing the viscous coefficients (bulk and shear).
    copy : bool, optional
        If true, a copy of the tensor is returned. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        The 6x6 viscous matrix in Voigt notation.

    Notes
    -----
    The viscous tensor is computed based on the provided bulk and shear coefficients.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for H_iso
        props_viscous = np.array([0.1, 0.05])
        H_iso = sim.H_iso(props_viscous)
        print(H_iso)
)pbdoc";

constexpr auto umat = R"pbdoc(
    Integrate a constitutive law over one increment at a batch of material points.

    This is the point-level entry used by finite-element couplers (fedoo).

    Parameters
    ----------
    umat_name : str
        5-letter model name (see the UMAT catalog).
    etot, Detot : numpy.ndarray, shape (6, N)
        Strain at the start of the increment and its increment (Voigt, engineering
        shear). Under finite strain: the corotational logarithmic strain.
    F0, F1 : numpy.ndarray, shape (3, 3, N)
        Deformation gradient at the start and the end of the increment. May be
        empty for a small-strain call.
    sigma : numpy.ndarray, shape (6, N)
        CAUCHY stress at the start of the increment.
    DR : numpy.ndarray, shape (3, 3, N)
        Rotation increment of the objective rate.
    props, statev, Wm : numpy.ndarray
        Material parameters, state variables and the accumulated work terms
        (Wm, Wm_r, Wm_ir, Wm_d).
    time, dtime : float
        Time at the start of the increment and its increment.
    temp : numpy.ndarray, optional
        Temperature per point.
    ndi : int
        Number of direct stress components (3, 2 or 1).
    n_threads : int
        Point-loop threads: 1 runs serially on every platform; otherwise it
        caps the worker threads on Windows (0 = all hardware threads), while
        macOS (GCD) and Linux (OpenMP) size the pool themselves. Batches of
        100 points or fewer always run serially.
    tangent_mode : int
        0 none, 1 continuum, 2 algorithmic (default).
    corate : int
        Objective rate the returned tangent is expressed in: 0 Jaumann,
        1 Green-Naghdi, 2 XBM (logarithmic), 3 log_R (default), 4 Truesdell,
        5 log_F. Choosing the coupler's own rate here spares it a tangent
        conversion. Kernels fed the corotated strain (the small-strain and
        log-strain boxes, MODUL, HYPOO) are in-rate already and ignore it.

    Returns
    -------
    tuple
        (sigma, statev, Wm, Lt): the CAUCHY stress at the end of the increment,
        the updated state variables and work terms, and the tangent Lt, shape
        (6, 6, N).

    Notes
    -----
    Stress measures. Inside simcoon the finite-strain route carries the KIRCHHOFF
    stress tau: every native kernel takes and returns tau, and Wm is accumulated
    per reference volume on tau. The Cauchy stress is formed only at the
    boundaries, and this function is one of them: for a Kirchhoff kernel with F0
    and F1 given, sigma is multiplied by det(F0) on the way in and divided by
    det(F1) on the way out, so the caller always exchanges Cauchy. Without F0/F1
    (or with a degenerate F) no conversion is applied, which is exact at small
    strain. Only the plugin adapters (UMEXT, UMABA) are Cauchy-native and pass
    through unconverted.

    Lt is NOT rescaled to Cauchy: it is the Kirchhoff box tangent
    d(tau_hat)/d(De) in the requested corate, with no J. Rescaling it by 1/J
    would break Lt_convert, which consumes exactly this object.

    A kernel may request a step cut instead of integrating a too-large increment
    (the modular engine on a non-finite or runaway return mapping, the SMR* SMA
    and LLDM0 damage laws on a failed local iteration, a Python law raising
    simcoon.StepCut). The simcoon solver retries automatically; this batch entry
    cannot subdivide, so it raises simcoon.StepCut (a RuntimeError) carrying the
    smallest requested ``ratio``. The input arrays are untouched: discard the
    call and retry with a smaller increment. umat_T behaves the same way.
)pbdoc";

} // namespace simcoon_docs