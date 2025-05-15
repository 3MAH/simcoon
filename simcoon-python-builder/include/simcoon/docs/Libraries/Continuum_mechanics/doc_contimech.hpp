#pragma once

namespace simpy_docs {

constexpr auto dev = R"pbdoc(
    Returns the deviatoric part of a 3x3 matrix.

    Parameters
    ----------
    m : pybind11::array_t<double>
        Input 3x3 matrix.

    Returns
    -------
    pybind11::array_t<double>
        Deviatoric part of the input matrix.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for dev
        m = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dev_m = sim.dev(m)
        print(dev_m)
)pbdoc";

constexpr auto sph = R"pbdoc(
    Returns the spherical part of a 3x3 matrix.

    Parameters
    ----------
    m : pybind11::array_t<double>
        Input 3x3 matrix.

    Returns
    -------
    pybind11::array_t<double>
        Spherical part of the input matrix.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for sph
        m = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sph_m = sim.sph(m)
        print(sph_m)
)pbdoc";

constexpr auto tr = R"pbdoc(
    Returns the trace of a tensor expressed in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Trace of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for tr
        v = np.array([1, 1, 1, 0, 0, 0])
        trace = sim.tr(v)
        print(trace)
)pbdoc";

constexpr auto dev_voigt = R"pbdoc(
    Returns the deviatoric part of a tensor expressed in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Deviatoric part of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for dev_voigt
        v = np.array([1, 1, 1, 0, 0, 0])
        dev_v = sim.dev_voigt(v)
        print(dev_v)
)pbdoc";

constexpr auto Mises_stress = R"pbdoc(
    Provides the Von Mises stress of a second-order stress tensor written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Von Mises stress of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Mises_stress
        v = np.array([1, 1, 1, 0, 0, 0])
        stress = sim.Mises_stress(v)
        print(stress)
)pbdoc";

constexpr auto eta_stress = R"pbdoc(
    Provides the stress flow of a second-order stress tensor written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Stress flow of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for eta_stress
        v = np.array([1, 1, 1, 0, 0, 0])
        flow = sim.eta_stress(v)
        print(flow)
)pbdoc";

constexpr auto eta_norm_stress = R"pbdoc(
    Provides the strain flow (direction) from a stress tensor (Euclidian norm), according to the Voigt convention for strains.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Strain flow of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for eta_norm_stress
        v = np.array([1, 1, 1, 0, 0, 0])
        flow = sim.eta_norm_stress(v)
        print(flow)
)pbdoc";

constexpr auto eta_norm_strain = R"pbdoc(
    Provides the strain flow (direction) from a strain tensor (Euclidian norm), according to the Voigt convention for strains.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input strain tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Strain flow of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for eta_norm_strain
        v = np.array([1, 1, 1, 0, 0, 0])
        flow = sim.eta_norm_strain(v)
        print(flow)
)pbdoc";

constexpr auto norm_stress = R"pbdoc(
    Provides the Euclidian norm of a stress tensor.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Euclidian norm of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for norm_stress
        v = np.array([1, 1, 1, 0, 0, 0])
        norm = sim.norm_stress(v)
        print(norm)
)pbdoc";

constexpr auto norm_strain = R"pbdoc(
    Provides the Euclidian norm of a strain tensor.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input strain tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Euclidian norm of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for norm_strain
        v = np.array([1, 1, 1, 0, 0, 0])
        norm = sim.norm_strain(v)
        print(norm)
)pbdoc";

constexpr auto Mises_strain = R"pbdoc(
    Provides the Von Mises strain of a second-order strain tensor.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input strain tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Von Mises strain of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Mises_strain
        v = np.array([1, 1, 1, 0, 0, 0])
        strain = sim.Mises_strain(v)
        print(strain)
)pbdoc";

constexpr auto eta_strain = R"pbdoc(
    Provides the strain flow of a second-order strain tensor written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input strain tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Strain flow of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for eta_strain
        v = np.array([1, 1, 1, 0, 0, 0])
        flow = sim.eta_strain(v)
        print(flow)
)pbdoc";

constexpr auto J2_stress = R"pbdoc(
    Provides the second invariant of the deviatoric part of a second-order stress tensor, written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Second invariant of the deviatoric part of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for J2_stress
        v = np.array([1, 1, 1, 0, 0, 0])
        invariant = sim.J2_stress(v)
        print(invariant)
)pbdoc";

constexpr auto J2_strain = R"pbdoc(
    Provides the second invariant of the deviatoric part of a second-order strain tensor, written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input strain tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Second invariant of the deviatoric part of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for J2_strain
        v = np.array([1, 1, 1, 0, 0, 0])
        invariant = sim.J2_strain(v)
        print(invariant)
)pbdoc";

constexpr auto J3_stress = R"pbdoc(
    Provides the third invariant of the deviatoric part of a second-order stress tensor, written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input stress tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Third invariant of the deviatoric part of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for J3_stress
        v = np.array([1, 1, 1, 0, 0, 0])
        invariant = sim.J3_stress(v)
        print(invariant)
)pbdoc";

constexpr auto J3_strain = R"pbdoc(
    Provides the third invariant of the deviatoric part of a second-order strain tensor, written as a vector in Voigt notation.

    Parameters
    ----------
    v : pybind11::array_t<double>
        Input strain tensor in Voigt notation.

    Returns
    -------
    pybind11::array_t<double>
        Third invariant of the deviatoric part of the input tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for J3_strain
        v = np.array([1, 1, 1, 0, 0, 0])
        invariant = sim.J3_strain(v)
        print(invariant)
)pbdoc";

constexpr auto Macaulay_p = R"pbdoc(
    Provides the results of the MacCaulay brackets operator <>+.

    Parameters
    ----------
    d : pybind11::array_t<double>
        Input tensor.

    Returns
    -------
    pybind11::array_t<double>
        Result of the MacCaulay brackets operator <>+.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Macaulay_p
        d = np.array([1, -1, 0])
        result = sim.Macaulay_p(d)
        print(result)
)pbdoc";

constexpr auto Macaulay_n = R"pbdoc(
    Provides the results of the MacCaulay brackets operator <>-.

    Parameters
    ----------
    d : pybind11::array_t<double>
        Input tensor.

    Returns
    -------
    pybind11::array_t<double>
        Result of the MacCaulay brackets operator <>-.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for Macaulay_n
        d = np.array([1, -1, 0])
        result = sim.Macaulay_n(d)
        print(result)
)pbdoc";

constexpr auto sign = R"pbdoc(
    Provides the results of the sign operator.

    Parameters
    ----------
    d : pybind11::array_t<double>
        Input tensor.

    Returns
    -------
    pybind11::array_t<double>
        Result of the sign operator.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for sign
        d = np.array([1, -1, 0])
        result = sim.sign(d)
        print(result)
)pbdoc";

constexpr auto normal_ellipsoid = R"pbdoc(
    Provides the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3.

    Parameters
    ----------
    u : pybind11::array_t<double>
        Input vector u.
    v : pybind11::array_t<double>
        Input vector v.
    a1 : pybind11::array_t<double>
        Length of semi-principal axis a1.
    a2 : pybind11::array_t<double>
        Length of semi-principal axis a2.
    a3 : pybind11::array_t<double>
        Length of semi-principal axis a3.

    Returns
    -------
    pybind11::array_t<double>
        Normalized vector normal to the ellipsoid.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for normal_ellipsoid
        u = np.array([1, 0, 0])
        v = np.array([0, 1, 0])
        a1 = np.array([1])
        a2 = np.array([1])
        a3 = np.array([1])
        normal = sim.normal_ellipsoid(u, v, a1, a2, a3)
        print(normal)
)pbdoc";

constexpr auto curvature_ellipsoid = R"pbdoc(
    Provides the curvature of an ellipsoid with semi-principal axes of length a1, a2, a3 at the angle u,v.

    Parameters
    ----------
    u : pybind11::array_t<double>
        Input vector u.
    v : pybind11::array_t<double>
        Input vector v.
    a1 : pybind11::array_t<double>
        Length of semi-principal axis a1.
    a2 : pybind11::array_t<double>
        Length of semi-principal axis a2.
    a3 : pybind11::array_t<double>
        Length of semi-principal axis a3.

    Returns
    -------
    pybind11::array_t<double>
        Curvature of the ellipsoid.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for curvature_ellipsoid
        u = np.array([1, 0, 0])
        v = np.array([0, 1, 0])
        a1 = np.array([1])
        a2 = np.array([1])
        a3 = np.array([1])
        curvature = sim.curvature_ellipsoid(u, v, a1, a2, a3)
        print(curvature)
)pbdoc";

constexpr auto sigma_int = R"pbdoc(
    Provides the normal and tangent components of the traction vector in the normal direction n to an ellipsoid with axes a1, a2, a3 from an input inside stress.

    Parameters
    ----------
    sigma_in : pybind11::array_t<double>
        Input stress tensor.
    u : pybind11::array_t<double>
        Input vector u.
    v : pybind11::array_t<double>
        Input vector v.
    a1 : pybind11::array_t<double>
        Length of semi-principal axis a1.
    a2 : pybind11::array_t<double>
        Length of semi-principal axis a2.
    a3 : pybind11::array_t<double>
        Length of semi-principal axis a3.

    Returns
    -------
    pybind11::array_t<double>
        Normal and tangent components of the traction vector.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for sigma_int
        sigma_in = np.array([1, 1, 1, 0, 0, 0])
        u = np.array([1, 0, 0])
        v = np.array([0, 1, 0])
        a1 = np.array([1])
        a2 = np.array([1])
        a3 = np.array([1])
        result = sim.sigma_int(sigma_in, u, v, a1, a2, a3)
        print(result)
)pbdoc";

constexpr auto p_ikjl = R"pbdoc(
    Computes the Hill interfacial operator according to a normal a.

    Parameters
    ----------
    a : pybind11::array_t<double>
        Input normal vector.

    Returns
    -------
    pybind11::array_t<double>
        Hill interfacial operator.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for p_ikjl
        a = np.array([1, 0, 0])
        result = sim.p_ikjl(a)
        print(result)
)pbdoc";

constexpr auto auto_sym_dyadic = R"pbdoc(
    Provides the dyadic product of a symmetric tensor with itself.

    Parameters
    ----------
    a : pybind11::array_t<double>
        Input symmetric tensor.

    Returns
    -------
    pybind11::array_t<double>
        Dyadic product of the input tensor with itself.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for auto_sym_dyadic
        a = np.array([1, 1, 1, 0, 0, 0])
        result = sim.auto_sym_dyadic(a)
        print(result)
)pbdoc";

constexpr auto sym_dyadic = R"pbdoc(
    Provides the dyadic product of two symmetric tensors.

    Parameters
    ----------
    a : pybind11::array_t<double>
        Input symmetric tensor a.
    b : pybind11::array_t<double>
        Input symmetric tensor b.

    Returns
    -------
    pybind11::array_t<double>
        Dyadic product of the input tensors.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for sym_dyadic
        a = np.array([1, 1, 1, 0, 0, 0])
        b = np.array([1, 1, 1, 0, 0, 0])
        result = sim.sym_dyadic(a, b)
        print(result)
)pbdoc";

constexpr auto auto_dyadic = R"pbdoc(
    Provides the dyadic product of a tensor with itself.

    Parameters
    ----------
    a : pybind11::array_t<double>
        Input tensor.

    Returns
    -------
    pybind11::array_t<double>
        Dyadic product of the input tensor with itself.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for auto_dyadic
        a = np.array([1, 1, 1, 0, 0, 0])
        result = sim.auto_dyadic(a)
        print(result)
)pbdoc";

constexpr auto dyadic_4vectors_sym = R"pbdoc(
    Provides the dyadic product of four vectors to provide a symmetric 4th order tensor.

    Parameters
    ----------
    n_a : pybind11::array_t<double>
        Input vector n_a.
    n_b : pybind11::array_t<double>
        Input vector n_b.
    conv : std::string
        Convention for the dyadic product.

    Returns
    -------
    pybind11::array_t<double>
        Symmetric 4th order tensor.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim        

        # Example for dyadic_4vectors_sym
        n_a = np.array([1, 0, 0])
        n_b = np.array([0, 1, 0])
        conv = "some_convention"
        result = sim.dyadic_4vectors_sym(n_a, n_b, conv)
        print(result)
)pbdoc";

constexpr auto auto_sym_dyadic_operator = R"pbdoc(
    Provides the symmetric 4th-order dyadic product of a symmetric tensor with itself.

    Parameters
    ----------
    a : pybind11::array_t<double>
        Input symmetric tensor.

    Returns
    -------
    pybind11::array_t<double>
        Symmetric 4th-order dyadic product of the input tensor with itself.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for auto_sym_dyadic_operator
        a = np.array([1, 1, 1, 0, 0, 0])
        result = sim.auto_sym_dyadic_operator(a)
        print(result)
)pbdoc";

constexpr auto sym_dyadic_operator = R"pbdoc(
    Provides the symmetric 4th-order dyadic product of two symmetric tensors.

    Parameters
    ----------
    a : pybind11::array_t<double>
        Input symmetric tensor a.
    b : pybind11::array_t<double>
        Input symmetric tensor b.

    Returns
    -------
    pybind11::array_t<double>
        Symmetric 4th-order dyadic product of the input tensors.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for sym_dyadic_operator
        a = np.array([1, 1, 1, 0, 0, 0])
        b = np.array([1, 1, 1, 0, 0, 0])
        result = sim.sym_dyadic_operator(a, b)
        print(result)
)pbdoc";

constexpr auto B_klmn = R"pbdoc(
    Provides the symmetric 4th-order tensor B_klmn.

    Parameters
    ----------
    b_i : pybind11::array_t<double>
        Input tensor b_i.
    b_j : pybind11::array_t<double>
        Input tensor b_j.

    Returns
    -------
    pybind11::array_t<double>
        Symmetric 4th-order tensor B_klmn.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        # Example for B_klmn
        b_i = np.array([1, 1, 1, 0, 0, 0])
        b_j = np.array([1, 1, 1, 0, 0, 0])
        result = sim.B_klmn(b_i, b_j)
        print(result)
)pbdoc";

} // namespace simpy_docs
