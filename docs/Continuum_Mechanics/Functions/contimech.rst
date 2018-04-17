The Continuum Mechanics Library
===============================

.. default-domain:: cpp

.. function:: double tr(const vec &v)

    Provides the trace of a second order tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double trace = tr(v);

.. function:: vec dev(const vec &v)

    Provides the deviatoric part of a second order tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        vec deviatoric = dev(v);

.. function:: double Mises_stress(const vec &v)

    Provides the Von Mises stress :math:`\sigma^{Mises}` of a second order stress tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double Mises_sig = Mises_stress(v);

.. function:: vec eta_stress(const vec &v)

    Provides the stress flow :math:`\eta_{stress}=\frac{3/2\sigma_{dev}}{\sigma_{Mises}}` from a second order stress tensor written as a vector v in the SMART+ formalism (i.e. the shear terms are multiplied by 2, providing shear angles).

    .. code-block:: cpp

        vec v = randu(6);
        vec sigma_f = eta_stress(v);

.. function:: double Mises_strain(const vec &v)

    Provides the Von Mises strain :math:`\varepsilon^{Mises}` of a second order stress tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double Mises_eps = Mises_strain(v);

.. function:: vec eta_strain(const vec &v)

    Provides the strain flow :math:`\eta_{strain}=\frac{2/3\varepsilon_{dev}}{\varepsilon_{Mises}}` from a second order strain tensor written as a vector v in the SMART+ formalism (i.e. the shear terms are multiplied by 2, providing shear angles).

    .. code-block:: cpp

        vec v = randu(6);
        vec eps_f = eta_strain(v);

.. function:: mat v2t_strain(const vec &v)

    Converts a second order strain tensor written as a vector v in the SMART+ formalism into a second order strain tensor written as a matrix m.

    .. code-block:: cpp

        vec v = randu(6);
        mat m = v2t_strain(v);

.. function:: vec t2v_strain (const mat &strain)

    Converts a second order strain tensor written as a matrix m in the SMART+ formalism into a second order strain tensor written as a vector v.

    .. code-block:: cpp

        mat m = randu(6,6);
        vec v = t2v_strain(m);

.. function:: mat v2t_stress(const vec &v)

    Converts a second order stress tensor written as a vector v in the SMART+ formalism into a second order stress tensor written as a matrix m.

    .. code-block:: cpp

        vec v = randu(6);
        mat m = v2t_stress(v);

.. function:: vec t2v_stress (const mat &stress)

    Converts a second order stress tensor written as a matrix m in the SMART+ formalism into a second order stress tensor written as a vector v.

    .. code-block:: cpp

        mat m = randu(6,6);
        vec v = t2v_stress(m);

.. function:: double J2_stress(const vec &v)

    Provides the second invariant of a second order stress tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double J2 = J2_stress(v);

.. function:: double J2_strain(const vec &v)

    Provides the second invariant of a second order strain tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double J2 = J2_strain(v);

.. function:: double J3_stress(const vec &v)

    Provides the third invariant of a second order stress tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double J3 = J3_stress(v);

.. function:: double J3_strain(const vec &v)

    Provides the third invariant of a second order strain tensor written as a vector v in the SMART+ formalism.

    .. code-block:: cpp

        vec v = randu(6);
        double J3 = J3_strain(v);

.. function:: double Macaulay_p(const double &d)

   This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)

.. function:: double Macaulay_n(const double &d)

   This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)

.. function:: double sign(const double &d)

   This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)

.. function:: vec normal_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3)

    Provides the normalized vector to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u and v. These 2 angles correspond to the rotations in the plan defined by the center of the ellipsoid, a1 and a2 directions for u, a1 and a3 ones for v. u = 0 corresponds to a1 direction and v = 0 correspond to a3 one. So the normal vector is set at the parametrized position :

    .. math::

        \begin{align}
        x & = a_{1} cos(u) sin(v) \\
        y & = a_{2} sin(u) sin(v) \\
        z & = a_{3} cos(v)
        \end{align}

    .. code-block:: cpp

        const double Pi = 3.14159265358979323846

        double u = (double)rand()/(double)(RAND_MAX) % 2*Pi - 2*Pi;
        double v = (double)rand()/(double)(RAND_MAX) % Pi - Pi;
        double a1 = (double)rand();
        double a2 = (double)rand();
        double a3 = (double)rand();
        vec v = normal_ellipsoid(u, v, a1, a2, a3);

.. function:: vec sigma_int(const vec &sigma_in, const double &a1, const double &a2, const double &a3, const double &u, const double &v)

    Provides the normal and tangent components of a stress vector σin in accordance with the normal direction n to an ellipsoid with axes a1, a2, a3. The normal vector is set at the parametrized position :

    .. math::

        \begin{align}
        x & = a_{1} cos(u) sin(v) \\
        y & = a_{2} sin(u) sin(v) \\
        z & = a_{3} cos(v)
        \end{align}

    .. code-block:: cpp

        vec sigma_in = randu(6);
        double u = (double)rand()/(double)(RAND_MAX) % Pi - Pi/2;
        double v = (double)rand()/(double)(RAND_MAX) % 2*Pi - Pi;
        double a1 = (double)rand();
        double a2 = (double)rand();
        double a3 = (double)rand();
        vec sigma_i = sigma_int(sigma_in, a1, a2, a3, u, v));

.. function:: mat p_ikjl(const vec &a)

    Provides the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer Ph.D. dissertation).

    .. code-block:: cpp

        vec v = randu(6);
        mat H = p_ikjl(v);
