The Kinematics Library
========================

.. default-domain:: cpp

.. function:: mat ER_to_F(const mat &E, const mat &R)

    Provides the transformation gradient :math:`\mathbf{F}`, from the Green-Lagrange strain :math:`\mathbf{E}` and the rotation :math:`\mathbf{R}`:

    .. math::

        \mathbf{F} = \mathbf{R} \cdot \mathbf{U} \quad \mathbf{E} = \frac{1}{2} \left( \sqrt{\mathbf{U}^2} - \mathbf{I} \right)

    .. code-block:: cpp

        mat E = randu(3,3);
        mat R = eye(3,3);
        mat F = ER_to_F(E, R);

.. function:: mat eR_to_F(const mat &e, const mat &R)

    Provides the transformation gradient :math:`\mathbf{F}`, from the logarithmic strain :math:`\mathbf{e}` and the rotation :math:`\mathbf{R}`:

    .. math::

        \mathbf{F} = \mathbf{V} \cdot \mathbf{R} \quad \mathbf{e} = \textrm{ln} \mathbf{V}

    .. code-block:: cpp

        mat e = randu(3,3);
        mat R = eye(3,3);
        mat F = eR_to_F(e, R);
        
.. function:: mat G_UdX(const mat &F)

    Provides the gradient of the displacement (Lagrangian) from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \nabla_X \mathbf{U} = \mathbf{F} - \mathbf{I}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat GradU = G_UdX(F);

.. function:: mat G_Udx(const mat &F)

    Provides the gradient of the displacement (Eulerian) from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \nabla_x \mathbf{U} =  \mathbf{I} - \mathbf{F}^{-1}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat gradU = G_UdX(F);

.. function:: mat R_Cauchy_Green(const mat &F)

    Provides the Right Cauchy-Green tensor :math:`\mathbf{C}`: from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \mathbf{C} =  \mathbf{F}^T \cdot \mathbf{F}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat C = R_Cauchy_Green(F);

.. function:: mat L_Cauchy_Green(const mat &F)

    Provides the Left Cauchy-Green tensor :math:`\mathbf{B}`: from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \mathbf{B} =  \mathbf{F} \cdot \mathbf{F}^T

    .. code-block:: cpp

        mat F = randu(3,3);
        mat B = L_Cauchy_Green(F);

.. function:: RU_decomposition(mat &R, mat &U, const mat &F)

    Provides the RU decomposition of the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \mathbf{F} = \mathbf{R} \cdot \mathbf{U} \quad \mathbf{U} = \sqrt{\mathbf{F}^T \cdot \mathbf{F}} \quad \mathbf{R} = \mathbf{F} \cdot \mathbf{U}^{-1}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat R = zeros(3,3);
        mat U = zeros(3,3);
        RU_decomposition(R, U, F);

.. function:: VR_decomposition(mat &R, mat &V, const mat &F)

Provides the VR decomposition of the transformation gradient :math:`\mathbf{F}`:

.. math::

    \mathbf{F} = \mathbf{V} \cdot \mathbf{R} \quad \mathbf{V} = \sqrt{\mathbf{F} \cdot \mathbf{F}^T} \quad \mathbf{R} = \mathbf{V}^{-1} \cdot \mathbf{F}

.. code-block:: cpp

    mat F = randu(3,3);
    mat R = zeros(3,3);
    mat V = zeros(3,3);
    VR_decomposition(R, V, F);

.. function:: vec Inv_X(const mat &X)

Provides the invariants of a symmetric tensor :math:`\mathbf{X}`:

    .. math::

        \mathbf{I}_1 = \textrm{trace} \left( X \right) \quad \mathbf{I}_2 = \frac{1}{2} \left( \textrm{trace} \left( X \right)^2 - \textrm{trace} \left( X^2 \right) \right) \quad \mathbf{I}_3 = \textrm{det} \left( X \right)

    .. code-block:: cpp

        mat F = randu(3,3);
        mat C = R_Cauchy_Green(F);
        vec I = Inv_X(F);
        

.. function:: mat Cauchy(const mat &F)

Provides the Cauchy tensor :math:`\mathbf{b}`: from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \mathbf{b} = \left( \mathbf{F} \cdot \mathbf{F}^T \right)^{-1}

    .. code-block:: cpp

        mat F = randu(3,3);
        mat b = Cauchy(F);

.. function:: mat Green_Lagrange(const mat &F)

Provides the Green-Lagrange tensor :math:`\mathbf{E}`: from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \mathbf{E} = \frac{1}{2} \left( \mathbf{F}^T \cdot \mathbf{F} - \mathbf{I} \right)

    .. code-block:: cpp

        mat F = randu(3,3);
        mat E = Green_Lagrange(F);

.. function:: mat Euler_Almansi(const mat &F)

Provides the Euler_Almansi tensor :math:`\mathbf{e}`: from the transformation gradient :math:`\mathbf{F}`:

    .. math::

        \mathbf{e} = \frac{1}{2} \left( \mathbf{I} - \left( \mathbf{F} \cdot \mathbf{F}^T \right)^T \right)

    .. code-block:: cpp

        mat F = randu(3,3);
        mat e = Euler_Almansi(F);

.. function:: mat Log_strain(const mat &F)

Provides the logarithmic strain tensor :math:`\mathbf{h}`: from the transformation gradient :math:`\mathbf{F}`:

.. math::

    \mathbf{h} = \textrm{ln} \left( V \right) = \frac{1}{2} \textrm{ln} \left( V^2 \right) = \frac{1}{2} \textrm{ln} \left( \mathbf{F} \cdot \mathbf{F}^T \right)

.. code-block:: cpp

    mat F = randu(3,3);
    mat h = Log_strain(F);

.. function:: mat finite_L(const mat &F0, const mat &F1, const double &DTime)

Provides the approximation of the Eulerian velocity tensor :math:`\mathbf{L}`: from the transformation gradient :math:`\mathbf{F}_0`: at time :math:`t_0`:, :math:`\mathbf{F}_1`: at time :math:`t_1`: and the time difference :math:`\Delta t = t_1 - t_0`:

.. math::

    \mathbf{L} = \frac{1}{\Delta t} \left( \mathbf{F}_1 - \mathbf{F}_0 \right) \cdot \mathbf{F}_1^{-1}

.. code-block:: cpp

    mat F0 = randu(3,3);
    mat F1 = randu(3,3);
    mat DTime = 0.1;
    mat L = finite_L(F0, F1, DTime);

.. function:: mat finite_D(const mat &F0, const mat &F1, const double &DTime)

Provides the approximation of the Eulerian symmetric rate tensor :math:`\mathbf{D}`: from the transformation gradient :math:`\mathbf{F}_0`: at time :math:`t_0`:, :math:`\mathbf{F}_1`: at time :math:`t_1`: and the time difference :math:`\Delta t = t_1 - t_0`: This is commonly referred as the rate of deformation (this necessitates although a specific discussion)

.. math::

    \mathbf{W} = \frac{1}{2} \left( \mathbf{L} - \mathbf{L}^T \right)

.. code-block:: cpp

    mat F0 = randu(3,3);
    mat F1 = randu(3,3);
    mat DTime = 0.1;
    mat D = finite_D(F0, F1, DTime);

.. function:: mat finite_W(const mat &F0, const mat &F1, const double &DTime)

Provides the approximation of the Eulerian antisymmetric spin tensor :math:`\mathbf{W}`: from the transformation gradient :math:`\mathbf{F}_0`: at time :math:`t_0`:, :math:`\mathbf{F}_1`: at time :math:`t_1`: and the time difference :math:`\Delta t = t_1 - t_0`: . This correspond to the Jaumann corotationnal rate:

.. math::

    \mathbf{W} = \frac{1}{2} \left( \mathbf{L} - \mathbf{L}^T \right)

.. code-block:: cpp

    mat F0 = randu(3,3);
    mat F1 = randu(3,3);
    mat DTime = 0.1;
    mat W = finite_W(F0, F1, DTime);

.. function:: mat finite_Omega(const mat &F0, const mat &F1, const double &DTime)

Provides the approximation of the Eulerian rigid-body rotation spin tensor :math:`\mathbf{\Omega}`: from the transformation gradient :math:`\mathbf{F}_0`: at time :math:`t_0`:, :math:`\mathbf{F}_1`: at time :math:`t_1`: and the time difference :math:`\Delta t = t_1 - t_0`: . This correspond to the Green-Naghdi corotationnal rate:

.. math::

    \mathbf{\Omega} = \frac{1}{\Delta t} \left( \mathbf{R}_1 - \mathbf{R}_0 \right) \cdot \mathbf{R}_1^{T}

Note that the rotation matrix is obtained from a RU decomposition of the transformation gradient :math:`\mathbf{F}`:

.. code-block:: cpp

    mat F0 = randu(3,3);
    mat F1 = randu(3,3);
    mat DTime = 0.1;
    mat Omega = finite_Omega(F0, F1, DTime);

.. function:: mat finite_DQ(const mat &Omega0, const mat &Omega1, const double &DTime)

Provides the Hughes-Winget approximation of a increment of rotation or transformation from the spin/velocity :math:`\mathbf{\Omega}_0`: at time :math:`t_0`:, :math:`\mathbf{\Omega}_1`: at time :math:`t_1`: and the time difference :math:`\Delta t = t_1 - t_0`:

.. math::

    \mathbf{\Delta Q} = \left( \mathbf{I} + \frac{1}{2} \Delta t \, \mathbf{\Omega}_0 \right) \cdot \left( \mathbf{I} - \frac{1}{2} \Delta t \, \mathbf{\Omega}_1 \right)^{-1}

.. code-block:: cpp

    mat Omega0 = randu(3,3);
    mat Omega1 = randu(3,3);
    mat DTime = 0.1;
    mat DQ = finite_DQ(Omega0, Omega1, DTime);
