#pragma once

namespace simcoon_docs {

constexpr auto isochoric_invariants = R"pbdoc(
    Provides the isochoric strain invariants from the left Cauchy-Green deformation tensor \( \mathbf{b} \).

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 3x3 matrix representing the left Cauchy-Green deformation tensor \( \mathbf{b} \).
    J : double, optional
        The determinant of the transformation gradient \( \mathbf{F} \). Default is 0.
    copy : bool, optional
        Whether to copy the input data. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A column vector of dimension 3 containing the three isochoric invariants.

    Notes
    -----
    The isochoric invariants are defined as:

    .. math::

        \bar{I}_1 = \textrm{tr} \bar{\mathbf{b}} \\
        \bar{I}_2 = \frac{1}{2} \left( \left(\textrm{tr} \bar{\mathbf{b}} \right)^2 - \textrm{tr} \bar{\mathbf{b}}^2 \right) \\
        \bar{I}_3 = \textrm{det} \bar{\mathbf{b}} = 1

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        b = np.random.rand(3, 3)
        J = np.linalg.det(b)
        invariants = sim.isochoric_invariants(b, J)
        print(invariants)
)pbdoc";

constexpr auto isochoric_pstretch = R"pbdoc(
    Provides the isochoric principal stretches from the left Cauchy-Green tensor \( \mathbf{b} \) or the Eulerian stretch tensor \( \mathbf{v} \).

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 3x3 matrix representing the left Cauchy-Green tensor \( \mathbf{b} \) or the Eulerian stretch tensor \( \mathbf{v} \).
    input_tensor : str, optional
        A string indicating the type of input tensor ("b" for \( \mathbf{b} \) or "V" for \( \mathbf{v} \)). Default is "V".
    J : double, optional
        The determinant of the transformation gradient \( \mathbf{F} \). Default is 0.
    copy : bool, optional
        Whether to copy the input data. Default is true.

    Returns
    -------
    pybind11::array_t<double>
        A column vector of dimension 3 containing the three isochoric principal stretches.

    Notes
    -----
    The isochoric principal stretches are defined as:

    .. math::

        \bar{\lambda}_i = J^{-1/3} \lambda_i

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        b = np.random.rand(3, 3)
        J = np.linalg.det(b)
        pstretch = sim.isochoric_pstretch(b, "b", J)
        print(pstretch)
)pbdoc";

} // namespace simcoon_docs