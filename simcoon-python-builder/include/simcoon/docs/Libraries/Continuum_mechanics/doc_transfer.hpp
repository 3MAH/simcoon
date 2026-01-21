#pragma once

namespace simcoon_docs {

constexpr auto v2t_strain = R"pbdoc(
    Transform a strain Voigt vector into a 3x3 strain matrix.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6-element Voigt vector representing strain in the order
        [e11, e22, e33, e23, e13, e12].
    copy : bool, optional
        If true, a copy of the returned array is made. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A 3x3 matrix with the full small-strain tensor corresponding to the
        input Voigt vector. Shear components are placed using engineering
        shear convention (e.g., e12 stored in position (0,1) and (1,0)).

    Notes
    -----
    This wrapper calls the underlying C++ helper which converts a 6-component
    Voigt column vector into a symmetric 3x3 matrix.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        v = np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.0])
        M = sim.v2t_strain(v)
        print(M)
)pbdoc";

constexpr auto t2v_strain = R"pbdoc(
    Transform a 3x3 strain matrix into a 6-element strain Voigt vector.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 3x3 symmetric strain matrix.
    copy : bool, optional
        If true, a copy of the returned array is made. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A 6-element Voigt column vector in the order
        [e11, e22, e33, e23, e13, e12].

    Notes
    -----
    The conversion uses the engineering shear components convention. The
    resulting vector is suitable for functions that expect Voigt-form strains.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        M = np.array([[0.1, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        v = sim.t2v_strain(M)
        print(v)
)pbdoc";

constexpr auto v2t_stress = R"pbdoc(
    Transform a stress Voigt vector into a 3x3 stress matrix.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6-element Voigt vector representing stress in the order
        [s11, s22, s33, s23, s13, s12].
    copy : bool, optional
        If true, a copy of the returned array is made. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A 3x3 symmetric stress matrix.

    Notes
    -----
    Shear components follow the engineering convention. This wrapper simply
    maps the Voigt components into the 3x3 matrix used by the C++ core.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        s = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        S = sim.v2t_stress(s)
        print(S)
)pbdoc";

constexpr auto t2v_stress = R"pbdoc(
    Transform a 3x3 stress matrix into a 6-element stress Voigt vector.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 3x3 symmetric stress matrix.
    copy : bool, optional
        If true, a copy of the returned array is made. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A 6-element Voigt column vector in the order
        [s11, s22, s33, s23, s13, s12].

    Notes
    -----
    This conversion is the inverse of `v2t_stress` and uses the same
    engineering shear convention.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        S = np.zeros((3,3))
        S[0,0] = 100.0
        v = sim.t2v_stress(S)
        print(v)
)pbdoc";

} // namespace simcoon_docs
