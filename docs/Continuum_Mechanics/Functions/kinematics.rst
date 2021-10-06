The Kinematics Library
========================

.. default-domain:: cpp

.. function:: mat ER_to_F(const mat &E, const mat &R) {

    Provides the transformation gradient :math:`\mathbf{F}`, from the Green-Lagrange strain :math:`\mathbf{E}` and the rotation :math:`\mathbf{R}`:

    .. math::

        \mathbf{F} = \mathbf{R} \cdot \mathbf{U} \quad \mathbf{E} = \frac{1}{2} \left( \sqrt{\mathbf{U}^2} - \mathbf{I} \right)

    .. code-block:: cpp

        mat E = randu(3,3);
        mat R = eye(3,3);
        mat F = ER_to_F(E, R);

.. function:: mat eR_to_F(const mat &e, const mat &R) {

    Provides the transformation gradient :math:`\mathbf{F}`, from the logarithmic strain :math:`\mathbf{e}` and the rotation :math:`\mathbf{R}`:

    .. math::

        \mathbf{F} = \mathbf{V} \cdot \mathbf{R} \quad \mathbf{e} = \textrm{ln} \mathbf{V}

    .. code-block:: cpp

        mat e = randu(3,3);
        mat R = eye(3,3);
        mat F = eR_to_F(e, R);
        
.. function:: G_UdX(const mat &F) {

    Provides the gradient of the displacement (Lagrangian) from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \nabla_X \mathbf{U} = \mathbf{F} - \mathbf{I}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat GradU = G_UdX(F);

.. function:: G_Udx(const mat &F) {

    Provides the gradient of the displacement (Eulerian) from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \nabla_x \mathbf{U} =  \mathbf{I} - \mathbf{F}^{-1}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat GradU = G_UdX(F);