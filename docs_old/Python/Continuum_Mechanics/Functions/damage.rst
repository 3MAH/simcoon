The Damage Library
==================

.. default-domain:: cpp

.. function:: double damage_weibull(const vec &stress, const double &damage, const double &alpha, const double &beta, const double &DTime, const string &criterion)

    Provides the damage evolution :math:`\delta D` considering a Weibull damage law.
    It is given by : :math:`\delta D = (1-D_{old})*\Big(1-exp\big(-1(\frac{crit}{\beta})^{\alpha}\big)\Big)`
    Parameters of this function are: the stress vector :math:`\sigma`, the old damage :math:`D_{old}`, the shape parameter :math:`\alpha`, the scale parameter :math:`\beta`, the time increment :math:`\Delta T` and the criterion (which is a string).

    The criterion possibilities are :
    “vonmises” : :math:`crit = \sigma_{Mises}`
    “hydro” : :math:`crit = tr(\sigma)`
    “J3” : :math:`crit = J3(\sigma)`
    Default value of the criterion is “vonmises”.

    .. code-block:: cpp

        double varD = damage_weibull(stress, damage, alpha, beta, DTime, criterion);

.. function:: double damage_kachanov(const vec &stress, const vec &strain, const double &damage, const double &A0, const double &r, const string &criterion)

    Provides the damage evolution :math:`\delta D` considering a Kachanov’s creep damage law.
    It is given by : :math:`\delta D = \Big(\frac{crit}{A_0(1-D_{old})}\Big)^r`
    Parameters of this function are: the stress vector :math:`\sigma`, the strain vector :math:`\epsilon`, the old damage :math:`D_{old}`, the material properties characteristic of creep damage :math:`(A_0,r)` and the criterion (which is a string).

    The criterion possibilities are :
    “vonmises” : :math:`crit = (\sigma*(1+\varepsilon))_{Mises}`
    “hydro” : :math:`crit = tr(\sigma*(1+\varepsilon))`
    “J3” : :math:`crit = J3(\sigma*(1+\varepsilon))`
    Here, the criterion has no default value.

    .. code-block:: cpp

        double varD = damage_kachanov(stress, strain, damage, A0, r, criterion);

.. function:: double damage_miner(const double &S_max, const double &S_mean, const double &S_ult, const double &b, const double &B0, const double &beta, const double &Sl_0)

    Provides the constant damage evolution :math:`\Delta D` considering a Woehler- Miner’s damage law.
    It is given by : :math:`\Delta D = \big(\frac{S_{Max}-S_{Mean}+Sl_0*(1-b*S_{Mean})}{S_{ult}-S_{Max}}\big)*\big(\frac{S_{Max}-S_{Mean}}{B_0*(1-b*S_{Mean})}\big)^\beta`
    Parameters of this function are: the max stress value :math:`\sigma_{Max}`, the mean stress value :math:`\sigma_{Mean}`, the “ult” stress value :math:`\sigma_{ult}`, the :math:`b`, the :math:`B_0`, the :math:`\beta` and the :math:`Sl_0`.

    Default value of :math:`Sl_0` is 0.0.

    .. code-block:: cpp

        double varD = damage_minerl(S_max, S_mean, S_ult, b, B0, beta, Sl_0);

.. function:: double damage_manson(const double &S_amp, const double &C2, const double &gamma2)

    Provides the constant damage evolution :math:`\Delta D` considering a Coffin-Manson’s damage law.
    It is given by : :math:`\Delta D = \big(\frac{\sigma_{Amp}}{C_{2}}\big)^{\gamma_2}`
    Parameters of this function are: the “amp” stress value :math:`\sigma_{Amp}`, the :math:`C_2` and the :math:`\gamma_2`.

    .. code-block:: cpp

        double varD = damage_manson(S_amp, C2, gamma2);
