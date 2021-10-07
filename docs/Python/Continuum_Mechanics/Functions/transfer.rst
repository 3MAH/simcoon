The Transfer Library
========================

.. default-domain:: cpp

.. function:: mat v2t_strain(const vec &v)

    Converts a second order strain tensor written as a vector v in the 'simcoon' formalism into a second order strain tensor written as a matrix m.

    .. code-block:: cpp

        vec v = randu(6);
        mat m = v2t_strain(v);

.. function:: vec t2v_strain (const mat &strain)

    Converts a second order strain tensor written as a matrix m in the 'simcoon' formalism into a second order strain tensor written as a vector v.

    .. code-block:: cpp

        mat m = randu(6,6);
        vec v = t2v_strain(m);

.. function:: mat v2t_stress(const vec &v)

    Converts a second order stress tensor written as a vector v in the 'simcoon' formalism into a second order stress tensor written as a matrix m.

    .. code-block:: cpp

        vec v = randu(6);
        mat m = v2t_stress(v);

.. function:: vec t2v_stress (const mat &stress)

    Converts a second order stress tensor written as a matrix m in the 'simcoon' formalism into a second order stress tensor written as a vector v.

    .. code-block:: cpp

        mat m = randu(6,6);
        vec v = t2v_stress(m);
