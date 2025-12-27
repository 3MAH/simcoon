#pragma once

namespace simcoon_docs {

constexpr auto Eshelby_sphere = R"pbdoc(
    Provides the Eshelby tensor of a spherical inclusion for isotropic linear elasticity.

    Parameters
    ----------
    nu : float
        the Poisson ratio

    Returns
    -------
    pybind11::array_t<double>
        Provides the Eshelby tensor of a spherical inclusion for isotropic linear elasticity as a numpy array, 
        according to the conventions of a localisation tensor, as a function of the Poisson ratio :math:`\nu`

    Notes
    -----
    The Eshelby tensor for spherical inclusions is defined as:

    .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{matrix}\right)

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        nu = 0.3
        S = sim.Eshelby_sphere(nu)
        print(S)
)pbdoc";

constexpr auto Eshelby_cylinder = R"pbdoc(
    Provides the Eshelby tensor of a cylindrical inclusion for isotropic linear elasticity

    Parameters
    ----------
    nu : float
        the Poisson ratio

    Returns
    -------
    pybind11::array_t<double>
        Provides the Eshelby tensor of a cylindrical inclusion for isotropic linear elasticity as a numpy array
        according to the conventions of a localisation tensor, as a function of the Poisson ratio :math:`\nu`

    Notes
    -----
    The Eshelby tensor for spherical inclusions is defined as:

    .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{matrix}\right)

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        nu = 0.3
        S = sim.Eshelby_sphere(nu)
        print(S)
)pbdoc";

constexpr auto Eshelby_prolate = R"pbdoc(
    Provides the Eshelby tensor of a prolate spheroid inclusion for isotropic linear elasticity

    Parameters
    ----------
    nu : float
        the Poisson ratio
    a_r : float
        Aspect ratio defined as :math:`a_r = \frac{a1}{a2} = \frac{a1}{a3}`

    Returns
    -------
    pybind11::array_t<double>
        Provides the Eshelby tensor of a prolate spheroid inclusion for isotropic linear elasticity as a numpy array
        according to the conventions of a localisation tensor, as a function of the Poisson ratio :math:`\nu`
        and the aspect ratio :math:`a_r = \frac{a1}{a2} = \frac{a1}{a3}`

    Notes
    -----
    The Eshelby tensor for prolate spheroidal inclusions is defined as:

    .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        S_{11} & S_{12} & S_{12} & 0 & 0 & 0 \\ 
        S_{21} & S_{22} & S_{23} & 0 & 0 & 0 \\
        S_{21} & S_{23} & S_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & S_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & S_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & S_{66} \end{matrix}\right)

    with

    .. math::    

        \begin{matrix}
        S_{11} & = \frac{1}{2(1-\nu)}\left(1-2\nu+\frac{3a_r^2-1}{a_r^2-1}-g\left(1-2\nu+\frac{3a_r^2}{a_r^2-1}\right)\right) \\
        S_{12} & = \frac{-1}{2(1-\nu)}\left(1-2\nu+\frac{1}{a_r^2-1}+g\left(1-2\nu+\frac{3}{a_r^2-1}\right)\right) \\
        S_{21} & = \frac{-a_r^2}{2(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(\frac{3a_r^2}{a_r^2-1}-\left(1-2\nu\right)\right) \\
        S_{22} & = \frac{3a_r^2}{8(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(1-2\nu-\frac{9}{4\left(a_r^2-1\right)}\right) \\
        S_{23} & = \frac{1}{4(1-\nu)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}-g\left(1-2\nu+\frac{3}{4\left(a_r^2-1\right)}\right)\right) \\
        S_{44} & = \frac{2}{4\left(1-\nu\right)}\left(1-2\nu-\frac{a_r^2+1}{a_r^2-1}-\frac{g}{2}\left(1-2\nu-\frac{3a_r^2+1}{a_r^2-1}\right)\right) \\
        S_{66} & = \frac{2}{4\left(1-\nu\right)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}+g\left(1-2\nu-\frac{3}{4\left(a_r^2-1\right)}\right)\right)
        \end{matrix}

    considering :math:`g = a_r\frac{a_r\sqrt{a_r^2-1}}{\left(a_r^2-1\right)^{\frac{3}{2}}} - acos(a_r)`

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        nu = 0.3
        a_r = 2.0
        S = sim.Eshelby_prolate(nu,a_r)
        print(S)
)pbdoc";

constexpr auto Eshelby_oblate = R"pbdoc(
    Provides the Eshelby tensor of a oblate spheroid inclusion for isotropic linear elasticity

    Parameters
    ----------
    nu : float
        the Poisson ratio
    a_r : float
        Aspect ratio defined as :math:`a_r = \frac{a1}{a2} = \frac{a1}{a3}`

    Returns
    -------
    pybind11::array_t<double>
        Provides the Eshelby tensor of a oblate spheroid inclusion for isotropic linear elasticity as a numpy array
        according to the conventions of a localisation tensor, as a function of the Poisson ratio :math:`\nu`
        and the aspect ratio :math:`a_r = \frac{a1}{a2} = \frac{a1}{a3}`

    Notes
    -----
    The Eshelby tensor for oblate spheroidal inclusions is defined as:

    .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        S_{11} & S_{12} & S_{12} & 0 & 0 & 0 \\ 
        S_{21} & S_{22} & S_{23} & 0 & 0 & 0 \\
        S_{21} & S_{23} & S_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & S_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & S_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & S_{66} \end{matrix}\right)

    with

    .. math::    

        \begin{matrix}
        S_{11} &= \frac{1}{2(1-\nu)}\left(1-2\nu+\frac{3a_r^2-1}{a_r^2-1}-g\left(1-2\nu+\frac{3a_r^2}{a_r^2-1}\right)\right) \\
        S_{12} &= \frac{-1}{2(1-\nu)}\left(1-2\nu+\frac{1}{a_r^2-1}+g\left(1-2\nu+\frac{3}{a_r^2-1}\right)\right) \\
        S_{21} &= \frac{-a_r^2}{2(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(\frac{3a_r^2}{a_r^2-1}-\left(1-2\nu\right)\right) \\
        S_{22} &= \frac{3a_r^2}{8(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(1-2\nu-\frac{9}{4\left(a_r^2-1\right)}\right) \\
        S_{23} &= \frac{1}{4(1-\nu)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}-g\left(1-2\nu+\frac{3}{4\left(a_r^2-1\right)}\right)\right) \\
        S_{44} &= \frac{2}{4\left(1-\nu\right)}\left(1-2\nu-\frac{a_r^2+1}{a_r^2-1}-\frac{g}{2}\left(1-2\nu-\frac{3a_r^2+1}{a_r^2-1}\right)\right) \\
        S_{66} &= \frac{2}{4\left(1-\nu\right)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}+g\left(1-2\nu-\frac{3}{4\left(a_r^2-1\right)}\right)\right)
        \end{matrix}

    considering :math:`g = a_r\frac{-a_r\sqrt{1-a_r^2}}{\left(1-a_r^2\right)^{\frac{3}{2}}} - acos(a_r)`

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        nu = 0.3
        a_r = 2.0
        S = sim.Eshelby_oblate(nu,a_r)
        print(S)
)pbdoc";

constexpr auto Eshelby = R"pbdoc(
    Provides a numerical estimation of the Eshelby tensor of an ellipsoisal, possibly anisotropic inclusion, 

    Parameters
    ----------
    L : numpy array
        the stiffness tensor of the considered material
    a1 : float
        the aspect ratio in the direction :math:`1`
    a2 : float
        the aspect ratio in the direction :math:`2`
    a3 : float
        the aspect ratio in the direction :math:`3`
    mp : int
        the number of integration points in the direction :math:`1`
    np : int
        the number of integration points in the direction :math:`2`

    Returns
    -------
    pybind11::array_t<double>
        Provides the numerical estimation of the Eshelby tensor of an ellispoid in the general case of anisotropic media,
        as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
        :math:`1`, :math:`2` and :math:`3`, respectively. It also requires the list of integration points and their respective weight
        for the numerical integration, as well as the number of integration points in the :math:`1` and :math:`2` directions.
        The points and weights are calculated using the *point* function that require to be called previously.

    Notes
    -----
    The Eshelby tensor is defined as:

    .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        S_{11} & S_{12} & S_{12} & 0 & 0 & 0 \\ 
        S_{21} & S_{22} & S_{23} & 0 & 0 & 0 \\
        S_{21} & S_{23} & S_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & S_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & S_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & S_{66} \end{matrix}\right)

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        mat L = sim.L_iso(210000., 0.3, "Enu");
        double a1 = 2.0;
        double a2 = 1.0;
        double a3 = 1.0;
        int mp = 4;
        int np = 4;
        mat S = Eshelby(L, a1, a2, a3, mp, np);
        print(S)
)pbdoc";

constexpr auto T_II = R"pbdoc(
    Provides a numerical estimation of the Hill interaction tensor of an ellipsoisal, possibly anisotropic inclusion, 

    Parameters
    ----------
    L : numpy array
        the stiffness tensor of the considered material
    a1 : float
        the aspect ratio in the direction :math:`1`
    a2 : float
        the aspect ratio in the direction :math:`2`
    a3 : float
        the aspect ratio in the direction :math:`3`
    mp : int
        the number of integration points in the direction :math:`1`
    np : int
        the number of integration points in the direction :math:`2`

    Returns
    -------
    pybind11::array_t<double>
        Provides the numerical estimation of the Hill interaction tensor of an ellispoid in the general case of anisotropic media,
        as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
        :math:`1`, :math:`2` and :math:`3`, respectively. It also requires the list of integration points and their respective weight
        for the numerical integration, as well as the number of integration points in the :math:`1` and :math:`2` directions.
        The points and weights are calculated using the *point* function that require to be called previously.

    Notes
    -----
    The Hill interaction tensor is defined as:

    .. math::

        \boldsymbol{T}=\left(\begin{matrix}
        T_{11} & T_{12} & T_{12} & 0 & 0 & 0 \\ 
        T_{21} & T_{22} & T_{23} & 0 & 0 & 0 \\
        T_{21} & T_{23} & T_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & T_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & T_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & T_{66} \end{matrix}\right)

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        mat L = sim.L_iso(210000., 0.3, "Enu");
        double a1 = 2.0;
        double a2 = 1.0;
        double a3 = 1.0;
        int mp = 4;
        int np = 4;
        mat S = T_II(L, a1, a2, a3, mp, np);
        print(S)
)pbdoc";

} // namespace simcoon_docs