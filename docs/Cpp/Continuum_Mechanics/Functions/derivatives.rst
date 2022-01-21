The Kinematics Library
========================

.. default-domain:: cpp

.. function:: dI1DS(const mat &S)

    Provides the derivative of the first invariant (trace) of a 2nd order tensor :math:`\mathbf{S}`. Such derivative returns the identity matrix :math:`\mathbf{I}`:
    
    .. math::

        \frac{\partial I_1}{\partial \mathbf{S}} = \mathbf{I}

    .. code-block:: cpp

        mat S = randu(3,3);
        mat dI1 = dI1DS(S);

.. function:: mat dI2DS(const mat &S) {

    Provides the derivative of the second invariant :math:`I_2 = \frac{1}{2} S_{ij} S_{ij}` of a 2nd order tensor :math:`\mathbf{S}`. Such derivative returns the tensor :math:`\mathbf{S}`:

    .. math::

        \frac{\partial I_2}{\partial \mathbf{S}} = \mathbf{S}

    .. code-block:: cpp

        mat S = randu(3,3);
        mat dI2 = dI2DS(S);

.. function:: dI3DS(const mat &S)

    Provides the derivative of the third invariant :math:`I_3 = \frac{1}{3} S_{ij} S_{jk} S_{ki}` of a 2nd order tensor :math:`\mathbf{S}`. Such derivative returns the tensor :math:`\left(\mathbf{S} \cdot \mathbf{S}\right)^T`

    .. math::

        \frac{\partial I_3}{\partial \mathbf{S}} = \left(\mathbf{S} \cdot \mathbf{S}\right)^T

    .. code-block:: cpp

        mat S = randu(3,3);
        mat dI3 = dI3DS(S);

.. function:: dtrSdS(const mat &S)

    Provides the derivative of the trace of a 2nd order tensor :math:`\mathbf{S}`. Such derivative returns the identity matrix :
    
    .. math::

        \frac{\partial tr(\mathbf{S})}{\partial \mathbf{S}} = \mathbf{I}

    .. code-block:: cpp

        mat S = randu(3,3);
        mat dtrS = dtrSdS(S);
                
.. function:: mat ddetSdS(const mat &S)

    Provides the derivative of the determinant of a 2nd order tensor :math:`\mathbf{S}`:

    .. math::

        \mathbf{C} = \textrm{det} (\mathbf{S}) \cdot \mathbf{S}^{-T}

    .. code-block:: cpp

        mat S = randu(3,3);
        mat ddetS = ddetSdS(S);

.. function:: mat dinvSdS(const mat &S)

    Provides the derivative of the inverse of a 2nd order tensor :math:`\mathbf{S}`:
    
    .. code-block:: cpp

        mat S = randu(3,3);
        mat dinvS = dinvSdS(S);
