#pragma once

namespace simcoon_docs {

constexpr auto check_symetries = R"pbdoc(
    Check material symmetries and recover elastic properties from a stiffness tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 stiffness matrix in Voigt notation.
    tol : float, optional
        Tolerance used when testing symmetry conditions. Defaults to 0.0.

    Returns
    -------
    dict
        A python dictionary with keys:
        - "umat_type": str, short code for symmetry (e.g. 'ELISO','ELIST','ELCUB','ELORTO')
        - "axis": int, principal axis index (1,2 or 3) when relevant, 0 otherwise
        - "maj_sym": int, 1 if major symmetry (L = L^T) holds within tolerance, 0 otherwise
        - "props": numpy.ndarray, column vector of recovered material properties (size depends on symmetry)

    Notes
    -----
    The function uses rotation/reflection checks and value-equality checks to classify the
    stiffness tensor symmetry (triclinic, monoclinic, orthotropic, cubic, transversely isotropic, isotropic).

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        L = sim.L_iso(np.array([210000, 0.3]), "Enu")
        info = sim.check_symetries(L, tol=1e-8)
        print(info["umat_type"], info["props"])
)pbdoc";

constexpr auto L_iso_props = R"pbdoc(
    Recover isotropic elastic properties (E, nu) from a stiffness tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 stiffness matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-2 array [:math:`E`, :math:`\nu`] where :math:`E` is Young's modulus and :math:`\nu` is Poisson's ratio.

    Notes
    -----
    Uses averaged shear and Lame-like combinations of the stiffness components to compute :math:`E` and :math:`\nu`.

    Examples
    --------
    .. code-block:: python

        import simcoon as sim
        props = sim.L_iso_props(L)
        E, nu = props
)pbdoc";

constexpr auto M_iso_props = R"pbdoc(
    Recover isotropic elastic properties (E, nu) from a compliance tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 compliance matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-2 array [:math:`E`, :math:`\nu`] where :math:`E` is Young's modulus and :math:`\nu` is Poisson's ratio.

    Examples
    --------
    .. code-block:: python

        import simcoon as sim
        props = sim.M_iso_props(M)
)pbdoc";

constexpr auto L_isotrans_props = R"pbdoc(
    Recover transversely isotropic elastic properties from a stiffness tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 stiffness matrix in Voigt notation.
    axis : int
        Axis of transverse isotropy (1, 2 or 3).

    Returns
    -------
    numpy.ndarray
        A length-5 array [:math:`E_{L}`, :math:`E_{T}`, :math:`\nu_{TL}`, :math:`\nu_{TT}`, :math:`G_{LT}`].

    Notes
    -----
    The returned properties follow the ordering used by the C++ implementation.

    Examples
    --------
    .. code-block:: python

        props = sim.L_isotrans_props(L, axis=3)
)pbdoc";

constexpr auto M_isotrans_props = R"pbdoc(
    Recover transversely isotropic elastic properties from a compliance tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 compliance matrix in Voigt notation.
    axis : int
        Axis of transverse isotropy (1, 2 or 3).

    Returns
    -------
    numpy.ndarray
        A length-5 array [:math:`E_{L}`, :math:`E_{T}`, :math:`\nu_{TL}`, :math:`\nu_{TT}`, :math:`G_{LT}`].
)pbdoc";

constexpr auto L_cubic_props = R"pbdoc(
    Recover cubic elastic properties (E, nu, G) from a stiffness tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 stiffness matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-3 array [:math:`E`, :math:`\nu`, :math:`G`].
)pbdoc";

constexpr auto M_cubic_props = R"pbdoc(
    Recover cubic elastic properties (E, nu, G) from a compliance tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 compliance matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-3 array [:math:`E`, :math:`\nu`, :math:`G`].
)pbdoc";

constexpr auto L_ortho_props = R"pbdoc(
    Recover orthotropic elastic properties from a stiffness tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 stiffness matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-9 array [:math:`E_{1}`, :math:`E_{2}`, :math:`E_{3}`, :math:`\nu_{12}`, :math:`\nu_{13}`, :math:`\nu_{23}`, :math:`G_{12}`, :math:`G_{13}`, :math:`G_{23}`].
)pbdoc";

constexpr auto M_ortho_props = R"pbdoc(
    Recover orthotropic elastic properties from a compliance tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 compliance matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-9 array [:math:`E_{1}`, :math:`E_{2}`, :math:`E_{3}`, :math:`\nu_{12}`, :math:`\nu_{13}`, :math:`\nu_{23}`, :math:`G_{12}`, :math:`G_{13}`, :math:`G_{23}`].
)pbdoc";

constexpr auto M_aniso_props = R"pbdoc(
    Recover a full set of anisotropic elastic properties from a compliance tensor.

    Parameters
    ----------
    input : pybind11::array_t<double>
        A 6x6 compliance matrix in Voigt notation.

    Returns
    -------
    numpy.ndarray
        A length-21 array [E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, :math:`\eta_{14}`, :math:`\eta_{15}`, :math:`\eta_{16}`,
        :math:`\eta_{24}`, :math:`\eta_{25}`, :math:`\eta_{26}`, :math:`\eta_{34}`, :math:`\eta_{35}`, :math:`\eta_{36}`, :math:`\eta_{45}`, :math:`\eta_{46}`, :math:`\eta_{56}`].

    Notes
    -----
    The :math:`\eta_{ij}` parameters are coupling coefficients that appear for a general anisotropic
    compliance description and are returned in the order used by the C++ implementation.
)pbdoc";

} // namespace simcoon_docs
