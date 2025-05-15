#pragma once

namespace simpy_docs {

constexpr auto damage_weibull = R"pbdoc(
    Provides the damage evolution \( \delta D \) considering a Weibull damage law.

    Parameters
    ----------
    stress : pybind11::array_t<double>
        The stress vector \( \sigma \).
    damage : double
        The old damage \( D_{old} \).
    alpha : double
        The shape parameter \( \alpha \).
    beta : double
        The scale parameter \( \beta \).
    DTime : double
        The time increment \( \Delta T \).
    criterion : str, optional
        The criterion (default is "vonmises").

    Returns
    -------
    double
        The damage evolution \( \delta D \).

    Notes
    -----
    The damage evolution is given by:

    .. math::

        \Delta D = (1-D_{old}) \cdot \Big(1-\exp\big(-1(\frac{crit}{\beta})^{\alpha}\big)\Big)

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        stress = np.random.rand(6)
        damage = 0.1
        alpha = 2.0
        beta = 1.5
        DTime = 0.01
        criterion = "vonmises"
        delta_D = sim.damage_weibull(stress, damage, alpha, beta, DTime, criterion)
        print(delta_D)
)pbdoc";

constexpr auto damage_kachanov = R"pbdoc(
    Provides the damage evolution \( \delta D \) considering a Kachanov damage law.

    Parameters
    ----------
    stress : pybind11::array_t<double>
        The stress vector \( \sigma \).
    strain : pybind11::array_t<double>
        The strain vector \( \epsilon \).
    damage : double
        The old damage \( D_{old} \).
    A0 : double
        The material properties characteristic of creep damage \( A_0 \).
    r : double
        The parameter \( r \).
    criterion : str
        The criterion.

    Returns
    -------
    double
        The damage evolution \( \delta D \).

    Notes
    -----
    The damage evolution is given by:

    .. math::

        \delta D = \Big(\frac{crit}{A_0(1-D_{old})}\Big)^r

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        stress = np.random.rand(6)
        strain = np.random.rand(6)
        damage = 0.1
        A0 = 1.0
        r = 2.0
        criterion = "vonmises"
        delta_D = sim.damage_kachanov(stress, strain, damage, A0, r, criterion)
        print(delta_D)
)pbdoc";

constexpr auto damage_miner = R"pbdoc(
    Provides the constant damage evolution \( \Delta D \) considering a Woehler-Miner's damage law.

    Parameters
    ----------
    S_max : double
        The maximum stress value \( \sigma_{Max} \).
    S_mean : double
        The mean stress value \( \sigma_{Mean} \).
    S_ult : double
        The "ultimate" stress value \( \sigma_{ult} \).
    b : double
        The parameter \( b \).
    B0 : double
        The parameter \( B_0 \).
    beta : double
        The parameter \( \beta \).
    Sl_0 : double, optional
        The parameter \( Sl_0 \). Default is 0.0.

    Returns
    -------
    double
        The damage evolution \( \Delta D \).

    Notes
    -----
    The damage evolution is given by:

    .. math::

        \Delta D = \left(\frac{S_{Max}-S_{Mean}+Sl_0\cdot(1-b\cdot S_{Mean})}{S_{ult}-S_{Max}}\right)\cdot\left(\frac{S_{Max}-S_{Mean}}{B_0\cdot(1-b\cdot S_{Mean})}\right)^\beta

    Examples
    --------
    .. code-block:: python

        import simcoon as sim

        S_max = 500.0
        S_mean = 300.0
        S_ult = 1000.0
        b = 0.1
        B0 = 200.0
        beta = 2.0
        Sl_0 = 0.0
        delta_D = sim.damage_miner(S_max, S_mean, S_ult, b, B0, beta, Sl_0)
        print(delta_D)
)pbdoc";

constexpr auto damage_manson = R"pbdoc(
    Provides the constant damage evolution \( \Delta D \) considering a Coffin-Manson's damage law.

    Parameters
    ----------
    S_amp : double
        The stress amplitude \( \sigma_{Amp} \).
    C2 : double
        The parameter \( C_2 \).
    gamma2 : double
        The parameter \( \gamma_2 \).

    Returns
    -------
    double
        The damage evolution \( \Delta D \).

    Notes
    -----
    The damage evolution is given by:

    .. math::

        \Delta D = \left(\frac{\sigma_{Amp}}{C_{2}}\right)^{\gamma_2}

    Examples
    --------
    .. code-block:: python

        import simcoon as sim

        S_amp = 150.0
        C2 = 100.0
        gamma2 = 1.5
        delta_D = sim.damage_manson(S_amp, C2, gamma2)
        print(delta_D)
)pbdoc";

// Add more documentation strings for other functions in damage.hpp as needed.

} // namespace simpy_docs