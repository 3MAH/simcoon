#pragma once

namespace simpy_docs {

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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Ireal
    >>> Ireal = sim.Ireal()
    >>> print(Ireal)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Ivol
    >>> Ivol = sim.Ivol()
    >>> print(Ivol)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Idev
    >>> Idev = sim.Idev()
    >>> print(Idev)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Ireal2
    >>> Ireal2 = sim.Ireal2()
    >>> print(Ireal2)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Idev2
    >>> Idev2 = sim.Idev2()
    >>> print(Idev2)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Ith
    >>> Ith = sim.Ith()
    >>> print(Ith)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Ir2
    >>> Ir2 = sim.Ir2()
    >>> print(Ir2)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for Ir05
    >>> Ir05 = sim.Ir05()
    >>> print(Ir05)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for L_iso
    >>> props = np.array([210000, 0.3])
    >>> L_iso = sim.L_iso(props, "Enu")
    >>> print(L_iso)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for M_iso
    >>> props = np.array([210000, 0.3])    
    >>> M_iso = sim.M_iso(props, "Enu")
    >>> print(M_iso)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for L_cubic
    >>> props_cubic = np.array([100000, 50000, 30000])
    >>> L_cubic = sim.L_cubic(props_cubic, "Cii")
    >>> print(L_cubic)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for M_cubic
    >>> props_cubic = np.array([100000, 50000, 30000])    
    >>> M_cubic = sim.M_cubic(props_cubic, "Cii")
    >>> print(M_cubic)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for L_ortho
    >>> props_ortho = np.array([100000, 50000, 30000, 0.3, 0.3, 0.3, 40000, 40000, 40000])
    >>> L_ortho = sim.L_ortho(props_ortho, "EnuG")
    >>> print(L_ortho)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for M_ortho
    >>> props_ortho = np.array([100000, 50000, 30000, 0.3, 0.3, 0.3, 40000, 40000, 40000])    
    >>> M_ortho = sim.M_ortho(props_ortho, "EnuG")
    >>> print(M_ortho)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for L_isotrans
    >>> props_isotrans = np.array([210000, 70000, 0.3, 0.3, 50000])
    >>> axis = 1
    >>> L_isotrans = sim.L_isotrans(props_isotrans, axis)
    >>> print(L_isotrans)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for M_isotrans
    >>> props_isotrans = np.array([210000, 70000, 0.3, 0.3, 50000])
    >>> axis = 1    
    >>> M_isotrans = sim.M_isotrans(props_isotrans, axis)
    >>> print(M_isotrans)
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
    >>> import numpy as np
    >>> import simcoon as sim

    # Example for H_iso
    >>> props_viscous = np.array([0.1, 0.05])
    >>> H_iso = sim.H_iso(props_viscous)
    >>> print(H_iso)
)pbdoc";

} // namespace simpy_docs