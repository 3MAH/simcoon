#pragma once

namespace simcoon_docs
{

    // Drucker stress
    constexpr auto Drucker_stress = R"pbdoc(
    Provides the Drucker equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t[double]
        Input stress vector.
    b : float
        Parameter that defines the equivalent stress.
    n : float
        Parameter that defines the equivalent stress.

    Returns
    -------
    float
        The Drucker equivalent stress.

    Notes
    -----
    Returns the Drucker equivalent stress :math:`\sigma^{P}`, considering the input stress :math:`\mathbf{\sigma}`, :math:`\sigma^{VM}` is the Von Mises computed equivalent stress, and :math:`b` and :math:`n` are parameters that define the equivalent stress:

    .. math::

        \sigma^{P} = \sigma^{VM} \left(\frac{1 + b \cdot J_3(\mathbf{\sigma})}{\left(J_2(\mathbf{\sigma})\right)^{3/2}} \right)^{n}

    where :math:`J_2` and :math:`J_3` are the second and third invariants of the deviatoric stress, and :math:`\sigma^{VM}` is the Von Mises equivalent stress:

    .. math::

        J_2 = \frac{1}{2} \mathrm{tr}(\mathbf{s}^2), \qquad J_3 = \frac{1}{3} \mathrm{tr}(\mathbf{s}^3)

    Note that if :math:`n > 10`, the Drucker criterion is sufficiently close to the Mises criterion that the Mises norm is used to avoid numerical instabilities from high-power computations.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        b = 1.2
        n = 0.5
        sigma_Drucker = sim.Drucker_stress(sigma, b, n)
        print(sigma_Drucker)
)pbdoc";

    // Derivative of Drucker stress
    constexpr auto dDrucker_stress = R"pbdoc(
    Provides the derivative of the Drucker equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t[double]
        Input stress vector.
    b : float
        Parameter that defines the equivalent stress.
    n : float
        Parameter that defines the equivalent stress.

    Returns
    -------
    pybind11::array_t[double]
        The derivative of the Drucker equivalent stress.

    Notes
    -----
    Returns the derivative of the Drucker equivalent stress with respect to stress. Its main use is to define evolution equations for strain based on an associated rule of a convex yield surface.

    .. math::

        \frac{\partial \sigma^{P}}{\partial \mathbf{\sigma}} =
        n \left(\frac{1 + b J_3}{J_2^{3/2}}\right)^{n-1}
        \left[
            \frac{b \, \sigma^{VM}}{J_2^{3/2}} \frac{\partial J_3}{\partial \mathbf{\sigma}}
            - \frac{3}{2} \frac{1 + b J_3}{J_2^{5/2}} \frac{\partial J_2}{\partial \mathbf{\sigma}}
        \right]
        + \left(\frac{1 + b J_3}{J_2^{3/2}}\right)^{n} \frac{\partial \sigma^{VM}}{\partial \mathbf{\sigma}}

    where :math:`J_2` and :math:`J_3` are the second and third invariants of the deviatoric stress, and :math:`\sigma^{VM}` is the Von Mises equivalent stress.

    Considering the input stress :math:`\mathbf{\sigma}`, :math:`\sigma^{VM}` is the Von Mises computed equivalent stress, and :math:`b` and :math:`n` are parameters that define the equivalent stress.

    Note that if :math:`n > 10`, the Drucker criterion is sufficiently close to the Mises criterion that the derivative of the Mises norm is used to avoid numerical instabilities from high-power computations.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        b = 1.2
        n = 0.5
        dsigma_Drucker = sim.dDrucker_stress(sigma, b, n)
        print(dsigma_Drucker)
)pbdoc";

    // Tresca stress
    constexpr auto Tresca_stress = R"pbdoc(
    Provides the Tresca equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t[double]
        Input stress vector.

    Returns
    -------
    float
        The Tresca equivalent stress.

    Notes
    -----
    Returns the Tresca equivalent stress :math:`\sigma^{T}`, considering the input stress :math:`\mathbf{\sigma}`:

    .. math::

        \sigma^{T} = \sigma_{I} - \sigma_{III}

    where :math:`\sigma_{I} \geq \sigma_{II} \geq \sigma_{III}` are the principal stresses.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        sigma_Tresca = sim.Tresca_stress(sigma)
        print(sigma_Tresca)
)pbdoc";

    // Derivative of Tresca stress
    constexpr auto dTresca_stress = R"pbdoc(
    Provides the derivative of the Tresca equivalent stress, given its vector representation.

    Parameters
    ----------
    v : pybind11::array_t[double]
        Input stress vector.

    Returns
    -------
    pybind11::array_t[double]
        The derivative of the Tresca equivalent stress.

    Notes
    -----
    Returns the derivative of the Tresca equivalent stress with respect to stress. Its main use is to define evolution equations for strain based on an associated rule of a Tresca convex but not smooth yield surface.

    .. warning::
        Note that so far the correct derivative is not implemented! Only the stress flow

    .. math::

        \eta_{stress} = \frac{3}{2} \frac{\sigma_{dev}}{\sigma_{Mises}}

    is returned.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        dsigma_Tresca = sim.dTresca_stress(sigma)
        print(dsigma_Tresca)
)pbdoc";

    // Anisotropic tensor P (general)
    constexpr auto P_Ani = R"pbdoc(
    Returns an anisotropic configurational tensor :math:`P` in the Voigt format (6x6 matrix), given its vector representation.

    Parameters
    ----------
    P_params : pybind11::array_t[double]
        A vector of 9 parameters: :math:`P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44},P_{55},P_{66}`.

    Returns
    -------
    pybind11::array_t[double]
        The anisotropic configurational tensor.

    Notes
    -----
    The tensor is constructed as:

    .. math::

        P_{ani} = \begin{pmatrix}
            P_{11} & P_{12} & P_{13} & 0 & 0 & 0 \\
            P_{12} & P_{22} & P_{23} & 0 & 0 & 0 \\
            P_{13} & P_{23} & P_{33} & 0 & 0 & 0 \\
            0 & 0 & 0 & 2P_{44} & 0 & 0 \\
            0 & 0 & 0 & 0 & 2P_{55} & 0 \\
            0 & 0 & 0 & 0 & 0 & 2P_{66}
        \end{pmatrix}

    The equivalent stress is :math:`\sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} }`.

    It reduces to:

    .. math::

        \sigma^{ani} = \left( P_{11}\,\sigma_{11}^2 + P_{22}\,\sigma_{22}^2 + P_{33}\,\sigma_{33}^2 + 2\,P_{12}\,\sigma_{11}\,\sigma_{22} + 2\,P_{13}\,\sigma_{11}\,\sigma_{33} + 2\,P_{23}\,\sigma_{22} \sigma_{33} + 2\,P_{44}\,\sigma_{12}^2 + 2\,P_{55}\,\sigma_{13}^2 + 2\,P_{66}\,\sigma_{23}^2 \right)^{1/2}

    For the Mises equivalent stress, :math:`P` reduces to:

    .. math::

        P_{ani} = \begin{pmatrix}
            1 & -1/2 & -1/2 & 0 & 0 & 0 \\
            -1/2 & 1 & -1/2 & 0 & 0 & 0 \\
            -1/2 & -1/2 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 3 & 0 & 0 \\
            0 & 0 & 0 & 0 & 3 & 0 \\
            0 & 0 & 0 & 0 & 0 & 3
        \end{pmatrix}

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        P_params = np.array([1., 1.2, 1.3, -0.2, -0.2, -0.33, 1., 1., 1.4])
        P = sim.P_Ani(P_params)
        print(P)
)pbdoc";

    // Hill tensor
    constexpr auto P_Hill = R"pbdoc(
    Provides an anisotropic configurational tensor for the quadratic Hill yield criterion in the Voigt format (6x6 matrix), given its vector representation.

    Parameters
    ----------
    P_params : pybind11::array_t[double]
        A vector of 6 parameters: :math:`F, G, H, L, M, N`.

    Returns
    -------
    pybind11::array_t[double]
        The Hill configurational tensor.

    Notes
    -----
    The tensor is constructed as:

    .. math::

        P_{Hill48} = \begin{pmatrix}
            G+H & -H & -G & 0 & 0 & 0 \\
            -H & F+H & -F & 0 & 0 & 0 \\
            -G & -F & F+G & 0 & 0 & 0 \\
            0 & 0 & 0 & 2N & 0 & 0 \\
            0 & 0 & 0 & 0 & 2M & 0 \\
            0 & 0 & 0 & 0 & 0 & 2L
        \end{pmatrix}

    The equivalent stress is :math:`\sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} }`.

    It reduces to:

    .. math::

        \sigma^{H48} = \left( H\, (\sigma_{11} - \sigma_{22})^2 + G\, (\sigma_{11} - \sigma_{33})^2 + F\, (\sigma_{22} - \sigma_{33})^2 + 2N\,\sigma_{12}^2 + 2M\,\sigma_{13}^2 + 2L\,\sigma_{23}^2 \right)^{1/2}

    For the Mises equivalent stress, :math:`P` reduces to:

    .. math::

        P_{Hill48} = \begin{pmatrix}
            1 & -1/2 & -1/2 & 0 & 0 & 0 \\
            -1/2 & 1 & -1/2 & 0 & 0 & 0 \\
            -1/2 & -1/2 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 3 & 0 & 0 \\
            0 & 0 & 0 & 0 & 3 & 0 \\
            0 & 0 & 0 & 0 & 0 & 3
        \end{pmatrix}

    So that :math:`F = H = G = 1/2`, :math:`L = M = N = 3/2`.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        P_params = np.array([0.5, 0.6, 0.7, 3.0, 3.0, 3.2])
        P = sim.P_Hill(P_params)
        print(P)
)pbdoc";

    // DFA tensor
    constexpr auto P_DFA = R"pbdoc(
    Provides an anisotropic configurational tensor for the quadratic Deshpande–Fleck–Ashby (DFA) yield criterion in the Voigt format (6x6 matrix), given its vector representation.

    Parameters
    ----------
    P_params : pybind11::array_t[double]
        A vector of 7 parameters: :math:`F, G, H, L, M, N, K`.

    Returns
    -------
    pybind11::array_t[double]
        The DFA configurational tensor.

    Notes
    -----
    The tensor is constructed as:

    .. math::

        P_{DFA} = \begin{pmatrix}
            G+H+K/9 & -H+K/9 & -G+K/9 & 0 & 0 & 0 \\
            -H+K/9 & F+H+K/9 & -F+K/9 & 0 & 0 & 0 \\
            -G+K/9 & -F+K/9 & F+G+K/9 & 0 & 0 & 0 \\
            0 & 0 & 0 & 2N & 0 & 0 \\
            0 & 0 & 0 & 0 & 2M & 0 \\
            0 & 0 & 0 & 0 & 0 & 2L
        \end{pmatrix}

    The equivalent stress is :math:`\sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} }`.

    It reduces to:

    .. math::

        \sigma^{DFA} = \left( H\, (\sigma_{11} - \sigma_{22})^2 + G\, (\sigma_{11} - \sigma_{33})^2 + F\, (\sigma_{22} - \sigma_{33})^2 + 2L\,\sigma_{12}^2 + 2M\,\sigma_{13}^2 + 2N\,\sigma_{23}^2 + K\,p^2 \right)^{1/2}


    with : :math:`p = \frac{1}{3} \textrm{tr} (\sigma)`

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        P_params = np.array([0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.0])
        P = sim.P_DFA(P_params)
        print(P)
)pbdoc";

    // Hill48 equivalent stress
    constexpr auto Hill_stress = R"pbdoc(
    Provides the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    P_params : pybind11::array_t[double]
        A vector of parameters :math:`F, G, H, L, M, N`.

    Returns
    -------
    float
        The Hill48 anisotropic equivalent stress.

    Notes
    -----
    Returns the Hill48 anisotropic equivalent stress :math:`\sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} }`.

    See the function :func:`P_Hill` for more details to obtain the tensor :math:`P` from the set of parameters :math:`(F,G,H,L,M,N)`.

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

    // Derivative of Hill48 equivalent stress
    constexpr auto dHill_stress = R"pbdoc(
    Provides the derivative of the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    P_params : pybind11::array_t[double]
        A vector of parameters :math:`F, G, H, L, M, N`.

    Returns
    -------
    pybind11::array_t[double]
        The derivative of the Hill48 anisotropic equivalent stress.

    Notes
    -----
    Returns the derivative of the Hill48 anisotropic equivalent stress

    .. math::

        \frac{\partial \sigma^{eq}_{ani}}{\partial \mathbf{\sigma}} = \frac{P \cdot \mathbf{\sigma}}{\sqrt{\mathbf{\sigma} : P : \mathbf{\sigma}}}

    where :math:`P` is constructed from the parameters as in :func:`P_Hill`.

    See the function :func:`P_Hill` for more details to obtain the tensor :math:`P` from the set of parameters :math:`(F,G,H,L,M,N)`.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        P_params = np.array([0.5, 0.6, 0.7, 1.5, 1.5, 1.6])
        dsigma_Hill = sim.dHill_stress(sigma, P_params)
        print(dsigma_Hill)
)pbdoc";

    // DFA equivalent stress
    constexpr auto DFA_stress = R"pbdoc(
    Provides the DFA anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    P_params : pybind11::array_t[double]
        A vector of parameters :math:`F, G, H, L, M, N, K`.

    Returns
    -------
    float
        The DFA anisotropic equivalent stress.

    Notes
    -----
    Returns the DFA anisotropic equivalent stress :math:`\sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} }`.

    See the function :func:`P_DFA` for more details to obtain the tensor :math:`P` from the set of parameters :math:`(F,G,H,L,M,N,K)`.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        P_params = np.array([0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2])
        sigma_DFA = sim.DFA_stress(sigma, P_params)
        print(sigma_DFA)
)pbdoc";

    // Derivative of DFA equivalent stress
    constexpr auto dDFA_stress = R"pbdoc(
    Provides the derivative of the DFA anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    P_params : pybind11::array_t[double]
        A vector of parameters :math:`F, G, H, L, M, N, K`.

    Returns
    -------
    pybind11::array_t[double]
        The derivative of the DFA anisotropic equivalent stress.

    Notes
    -----
    Returns the derivative of the DFA anisotropic equivalent stress

    .. math::

        \frac{\partial \sigma^{eq}_{ani}}{\partial \mathbf{\sigma}} = \frac{P \cdot \mathbf{\sigma}}{\sqrt{\mathbf{\sigma} : P : \mathbf{\sigma}}}

    where :math:`P` is constructed from the parameters as in :func:`P_DFA`.

    See the function :func:`P_DFA` for more details to obtain the tensor :math:`P` from the set of parameters :math:`(F,G,H,L,M,N,K)`.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        P_params = np.array([0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2])
        dsigma_DFA = sim.dDFA_stress(sigma, P_params)
        print(dsigma_DFA)
)pbdoc";

    // General anisotropic equivalent stress
    constexpr auto Ani_stress = R"pbdoc(
    Provides the anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    P_params : pybind11::array_t[double]
        A vector of parameters :math:`P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44},P_{55},P_{66}`.

    Returns
    -------
    float
        The anisotropic equivalent stress.

    Notes
    -----
    Returns the anisotropic equivalent stress

    .. math::

        \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} }

    where :math:`P` is constructed from the parameters as in :func:`P_Ani`.

    See the function :func:`P_Ani` for more details to obtain the tensor :math:`P` from the set of parameters.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        P_params = np.array([1.0, 1.2, 1.3, -0.2, -0.2, -0.33, 1.0, 1.0, 1.4])
        sigma_Ani = sim.Ani_stress(sigma, P_params)
        print(sigma_Ani)
)pbdoc";

    // Derivative of general anisotropic equivalent stress
    constexpr auto dAni_stress = R"pbdoc(
    Provides the derivative of the anisotropic equivalent stress, given the stress in a vector format and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    P_params : pybind11::array_t[double]
        A vector of parameters :math:`P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44},P_{55},P_{66}`.

    Returns
    -------
    pybind11::array_t[double]
        The derivative of the anisotropic equivalent stress.

    Notes
    -----
    Returns the derivative of the anisotropic equivalent stress

    .. math::

        \frac{\partial \sigma^{eq}_{ani}}{\partial \mathbf{\sigma}} = \frac{P \cdot \mathbf{\sigma}}{\sqrt{\mathbf{\sigma} : P : \mathbf{\sigma}}}

    where :math:`P` is constructed from the parameters as in :func:`P_Ani`.

    See the function :func:`P_Ani` for more details to obtain the tensor :math:`P` from the set of parameters.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.random.rand(6)
        P_params = np.array([1.0, 1.2, 1.3, -0.2, -0.2, -0.33, 1.0, 1.0, 1.4])
        dsigma_Ani = sim.dAni_stress(sigma, P_params)
        print(dsigma_Ani)
)pbdoc";

    // General equivalent stress selector
    constexpr auto Eq_stress = R"pbdoc(
    Provides an equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress, and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    type_eq : str
        The type of equivalent stress ("Mises", "Tresca", "Drucker", "Hill", "Ani").
    P_params : pybind11::array_t[double], optional
        A vector of parameters defining the equivalent stress.

    Returns
    -------
    float
        The equivalent stress.

    Notes
    -----
    Returns the equivalent stress, depending on the equivalent type provided: "Mises", "Tresca", "Drucker", "Hill", "Ani".

    See the detailed function for each equivalent stress to determine the number and type of parameters to provide, in the form of a numpy array (for the Drucker criterion, P_params = [b, n] for instance).

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

    // Derivative of general equivalent stress selector
    constexpr auto dEq_stress = R"pbdoc(
    Provides the derivative of the equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress, and a vector of parameters.

    Parameters
    ----------
    sigma : pybind11::array_t[double]
        The stress vector.
    type_eq : str
        The type of equivalent stress ("Mises", "Tresca", "Drucker", "Hill", "Ani").
    P_params : pybind11::array_t[double], optional
        A vector of parameters defining the equivalent stress.

    Returns
    -------
    pybind11::array_t[double]
        The derivative of the equivalent stress.

    Notes
    -----
    Returns the derivative of the equivalent stress, depending on the equivalent type provided: "Mises", "Tresca", "Drucker", "Hill", "Ani".

    See the detailed function for each equivalent stress to determine the number and type of parameters to provide, in the form of a numpy array (for the Drucker criterion, P_params = [b, n] for instance).

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

} // namespace simcoon_docs