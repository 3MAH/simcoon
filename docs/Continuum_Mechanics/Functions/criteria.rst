The Criteria Library
========================

.. default-domain:: cpp

.. function:: double Prager_stress(const vec &v, const double &b, const double &n)

    Returns the Prager equivalent stress :math:`\boldsymbol{\sigma}^{P}`, considering

    .. math::

        \boldsymbol{\sigma}^{P} = \boldsymbol{\sigma}^{VM} \left(\frac{1 + b \cdot J_3 \left(\boldsymbol{\sigma} \right)}{\left(J_2 \left(\boldsymbol{\sigma} \right) \right)^{3/2} } \right)^{m}

    considering the input stress :math:`\boldsymbol{\sigma}`, :math:`\boldsymbol{\sigma}^{VM}` is the Von Mises computed equivalent stress, and :math:`b` and :math:`m` are parameter that define the equivalent stress.

    .. code-block:: cpp

        vec sigma(6);
        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        double sigma_Prager = Prager_stress(sigma, b, n);

.. function:: vec dPrager_stress(const vec &v, const double &b, const double &n)

    Returns the derivative of the Prager equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a convex yield surface

   .. code-block:: cpp

        vec sigma(6);
        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        vec dsigma_Pragerdsigma = dPrager_stress(sigma, b, n);
