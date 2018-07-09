The Criteria Library
========================

.. default-domain:: cpp

.. function:: double Prager_stress(const vec &v, const double &b, const double &n)

    Returns the Prager equivalent stress :math:`\boldsymbol{\sigma}^{P}`, considering

    .. math::

        \sigma^{P} = \sigma^{VM} \left(\frac{1 + b \cdot J_3 \left(\boldsymbol{\sigma} \right)}{\left(J_2 \left(\boldsymbol{\sigma} \right) \right)^{3/2} } \right)^{m}

    considering the input stress :math:`\boldsymbol{\sigma}`, :math:`\boldsymbol{\sigma}^{VM}` is the Von Mises computed equivalent stress, and :math:`b` and :math:`m` are parameter that define the equivalent stress.

    .. code-block:: cpp

        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        double sigma_Prager = Prager_stress(sigma, b, n);

.. function:: vec dPrager_stress(const vec &v, const double &b, const double &n)

    Returns the derivative of the Prager equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a convex yield surface

   .. code-block:: cpp

        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        vec dsigma_Pragerdsigma = dPrager_stress(sigma, b, n);

.. function:: double Tresca(const vec &v)

    Returns the Tresca equivalent stress :math:`\boldsymbol{\sigma}^{T}`, considering

    .. math::

        \sigma^{T} = \sigma_{I} - \sigma_{III},

    where \sigma_{I} and \sigma_{III} are the highest and lowest principal stress values, respectively.

    .. code-block:: cpp

        vec sigma = randu(6);
        double sigma_Prager = Tresca_stress(sigma);

.. function:: vec dTresca_stress(const vec &v)

    Returns the derivative of the Tresca equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a convex yield surface.

    .. warning:: Note that so far that the correct derivative it is not implemented! Only stress flow :math:`\eta_{stress}=\frac{3/2\sigma_{dev}}{\sigma_{Mises}}` is returned

    .. code-block:: cpp

        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        vec dsigma_Pragerdsigma = dPrager_stress(sigma, b, n);

.. function:: mat P_Ani(const vec &params);

    Returns an anisotropic configurational tensor in the Voigt format (6x6 matrix)

    The vector of parameters must be constituted of 9 values, respectively:
    :math:`P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44}=P_{1212},P_{55}=P_{1313},P_66=P_{2323}`

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
        mat P = P_Ani(P_params);

.. function:: mat P_Hill(const vec &params);

    Returns an anisotropic configurational tensor considering the quadratic Hill yield criterion [Hill48].

    The vector of parameters must be constituted of 5 values, respectively:
    :math:`F^*,G^*,H^*,L,M,N`

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,0.95,0.8,1.2};
        mat P = P_Hill(P_params);

    Note that the values of :math:`F^*,G^*,H^*` have been scaled up so that math:`F^*=\frac{1}{3}F,G^*=\frac{1}{3}G,H^*=\frac{1}{3}H`.
    The reason is that if :math:`F^*=G^*=H^*=L=M=N=1`, the Mises equivalent stress is retrieved when defining an equivalent stress based on the obtained configurational tensor (see below).

.. function:: double Ani_stress(const vec &v, const mat &H)

    Returns an anisotropic equivalent stress, providing a configurational tensor

    .. math::

        \sigma^{Ani} = \sqrt{\frac{3.}{2.} \boldsymbol{\sigma} \cdot \boldsymbol{H} \cdot \boldsymbol{\sigma}}

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,0.95,0.8,1.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        double sigma_ani = Ani_stress(sigma,P_Hill);

.. function:: double dAni_stress(const vec &v, const mat &H)

    Returns an the derivative (with respect to stress) of an anisotropic equivalent stress, providing a configurational tensor

    .. warning:: Might be not stable for pure deviatoric criteria

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,0.95,0.8,1.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        vec dsigma_anidsigma = dAni_stress(sigma,P_Hill);
}

.. function:: double Hill_stress(const vec &v, const vec &params)

    Returns an the Hill equivalent stress, providing a set of Parameters

    .. seealso:: The definition of the *P_Hill* function: :func:`P_Hill`.

    .. math::

        \sigma^{Ani} = \sqrt{\frac{3.}{2.} \boldsymbol{\sigma} \cdot \boldsymbol{H} \cdot \boldsymbol{\sigma}}

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,0.95,0.8,1.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        double sigma_ani = Ani_stress(sigma,P_Hill);

.. function:: double Hill_stress(const vec &v, const vec &params)

    Returns an the Hill equivalent stress, providing a configurational tensor

    .. math::

        \sigma^{Ani} = \sqrt{\frac{3.}{2.} \boldsymbol{\sigma} \cdot \boldsymbol{H} \cdot \boldsymbol{\sigma}}

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,0.95,0.8,1.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        double sigma_ani = Ani_stress(sigma,P_Hill);

.. function:: double Hill_stress(const vec &v, const vec &params)

    Returns an the Hill equivalent stress, providing a configurational tensor

    .. math::

        \sigma^{Ani} = \sqrt{\frac{3.}{2.} \boldsymbol{\sigma} \cdot \boldsymbol{H} \cdot \boldsymbol{\sigma}}

    .. code-block:: cpp

        vec P_params = {1.,1.2,1.3,0.95,0.8,1.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        double sigma_ani = Ani_stress(sigma,P_Hill);


double Hill_stress(const vec &v, const vec &params) {

.. rubric:: References

[Hill48] Hill R. A theory of the yielding and plastic fow of anisotropic materials. Proc R Soc. 1947;(193):281–97.

