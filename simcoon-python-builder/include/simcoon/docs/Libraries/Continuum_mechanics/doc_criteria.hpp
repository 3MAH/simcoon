#pragma once

namespace simpy_docs {

constexpr auto Prager_stress = R"pbdoc(
    Provides the Prager equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress vector.
    b : double
        Parameter that defines the equivalent stress.
    n : double
        Parameter that defines the equivalent stress.

    Returns
    -------
    double
        The Prager equivalent stress.

    Notes
    -----
    Returns the Prager equivalent stress \( \sigma^{P} \), considering:

    .. math::

        \sigma^{P} = \sigma^{VM} \left(\frac{1 + b \cdot J_3}{\left(J_2\right)^{3/2}} \right)^{m}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        b = 1.2
        n = 0.5
        sigma_Prager = sim.Prager_stress(sigma, b, n)
        print(sigma_Prager)
)pbdoc";

constexpr auto Tresca_stress = R"pbdoc(
    Provides the Tresca equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress vector.

    Returns
    -------
    double
        The Tresca equivalent stress.

    Notes
    -----
    Returns the Tresca equivalent stress \( \sigma^{T} \), considering:

    .. math::

        \sigma^{T} = \sigma_{I} - \sigma_{III}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        sigma_Tresca = sim.Tresca_stress(sigma)
        print(sigma_Tresca)
)pbdoc";

constexpr auto dPrager_stress = R"pbdoc(
    Provides the derivative of the Prager equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress vector.
    b : double
        Parameter that defines the equivalent stress.
    n : double
        Parameter that defines the equivalent stress.

    Returns
    -------
    pybind11::array_t<double>
        The derivative of the Prager equivalent stress.

    Notes
    -----
    Returns the derivative of the Prager equivalent stress with respect to stress. It is mainly used to define evolution equations for strain based on an associated rule of a convex yield surface.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        b = 1.2
        n = 0.5
        dsigma_Prager = sim.dPrager_stress(sigma, b, n)
        print(dsigma_Prager)
)pbdoc";

constexpr auto dTresca_stress = R"pbdoc(
    Provides the derivative of the Tresca equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress vector.

    Returns
    -------
    pybind11::array_t<double>
        The derivative of the Tresca equivalent stress.

    Notes
    -----
    Returns the derivative of the Tresca equivalent stress with respect to stress. This is primarily used to define evolution equations for strain based on the Tresca yield surface.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        dsigma_Tresca = sim.dTresca_stress(sigma)
        print(dsigma_Tresca)
)pbdoc";

constexpr auto P_Ani = R"pbdoc(
    Returns an anisotropic configurational tensor \( P \) in the Voigt format (6x6 matrix).

    Parameters
    ----------
    P_params : pybind11::array_t<double>
        A vector of parameters defining the tensor.

    Returns
    -------
    pybind11::array_t<double>
        The anisotropic configurational tensor.

    Notes
    -----
    The vector of parameters must contain 9 values, corresponding to specific tensor components.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        P_params = np.array([1., 1.2, 1.3, -0.2, -0.2, -0.33, 1., 1., 1.4])
        P = sim.P_Ani(P_params)
        print(P)
)pbdoc";

constexpr auto P_Hill = R"pbdoc(
    Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion.

    Parameters
    ----------
    P_params : pybind11::array_t<double>
        A vector of parameters defining the Hill criterion.

    Returns
    -------
    pybind11::array_t<double>
        The Hill configurational tensor.

    Notes
    -----
    The tensor is computed based on the provided parameters \( F, G, H, L, M, N \).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        P_params = np.array([0.5, 0.6, 0.7, 3.0, 3.0, 3.2])
        P = sim.P_Hill(P_params)
        print(P)
)pbdoc";

constexpr auto Hill_stress = R"pbdoc(
    Provides the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t<double>
        The stress vector.
    P_params : pybind11::array_t<double>
        A vector of parameters defining the Hill48 criterion.

    Returns
    -------
    double
        The Hill48 anisotropic equivalent stress.

    Notes
    -----
    The equivalent stress is computed based on the Hill48 criterion.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        P_params = np.array([0.5, 0.6, 0.7, 1.5, 1.5, 1.6])
        sigma_Hill = sim.Hill_stress(sigma, P_params)
        print(sigma_Hill)
)pbdoc";

constexpr auto Eq_stress = R"pbdoc(
    Provides an equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress, and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t<double>
        The stress vector.
    type_eq : str
        The type of equivalent stress ("Mises", "Tresca", "Prager", "Hill", "Ani").
    P_params : pybind11::array_t<double>, optional
        A vector of parameters defining the equivalent stress. Default is an empty vector.

    Returns
    -------
    double
        The equivalent stress.

    Notes
    -----
    The equivalent stress is computed based on the specified type and parameters.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        type_eq = "Mises"
        P_params = np.array([1.0, 1.2])
        sigma_eq = sim.Eq_stress(sigma, type_eq, P_params)
        print(sigma_eq)
)pbdoc";

constexpr auto dEq_stress = R"pbdoc(
    Provides the derivative of the equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress, and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t<double>
        The stress vector.
    type_eq : str
        The type of equivalent stress ("Mises", "Tresca", "Prager", "Hill", "Ani").
    P_params : pybind11::array_t<double>, optional
        A vector of parameters defining the equivalent stress. Default is an empty vector.

    Returns
    -------
    pybind11::array_t<double>
        The derivative of the equivalent stress.

    Notes
    -----
    The derivative is computed based on the specified type and parameters.

    Examples
    --------
        .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        type_eq = "Mises"
        P_params = np.array([1.0, 1.2])
        dsigma_eq = sim.dEq_stress(sigma, type_eq, P_params)
        print(dsigma_eq)
)pbdoc";

// Add more documentation strings for other functions in criteria.hpp as needed.

} // namespace simpy_docs