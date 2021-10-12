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

.. function:: mat v2t_stress(const vec &v)

    Converts a second order stress tensor written as a vector v in the 'simcoon' formalism into a second order stress tensor written as a matrix m.

.. code-block:: cpp

    vec v = randu(6);
    mat m = v2t_stress(v);

.. function:: vec t2v_sym(const mat &m)

    Converts a 3x3 symmetric matrix into a 6 component vector {11,22,33,12,13,23}

.. code-block:: cpp

    mat m = randu(3,3);
    vec v = t2v_sym(m);

.. function:: mat v2t_sym(const vec &v)

    Converts a 6 component vector {11,22,33,12,13,23} into a 3x3 symmetric matrix

.. code-block:: cpp

    vec v = randu(6);
    mat m = t2v_sym(m);

.. function:: mat v2t_skewsym(const vec &v)

    Converts a 6 component vector {11,22,33,12,13,23} into a 3x3 antisymmetric matrix, while keeping the diagonal components

.. math::

    m = \left( \begin{array}{ccc}
    v_1 & v_4 & v_5 \\
    -v_4 & v_2 & v_6 \\
    v_5 & -v_6 & v_3 \end{array} \right)

.. code-block:: cpp

    vec v = randu(6);
    mat m = t2v_sym(m);
    
.. function:: mat v2t(const vec &v)

    Converts a 9 component vector {11,12,13,21,22,23,31,32,33} into a 3x3 symmetric matrix

.. code-block:: cpp

    vec v = randu(9);
    mat m = t2v(m);

