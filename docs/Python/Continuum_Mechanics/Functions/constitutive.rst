The Constitutive Library
========================

.. default-domain:: python

.. function:: Import of the simmit module

    To import the simmit module as *sim*

    .. code-block:: python
    
        from simcoon import simmit as sim

.. function:: np.ndarray Ireal()

    Provides the fourth order identity tensor written in Voigt notation :math:`I_{real}`, where :

    .. math::

        I_{real} = \left( \begin{array}{cccccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)

    .. code-block:: python

        Ir = sim.Ireal()

.. function:: mat Ivol()

    Provides the volumic of the identity tensor :math:`I_{vol}` written in the Simcoon formalism. So :

    .. math::

        I_{vol} = \left( \begin{array}{ccc}
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \end{array} \right)

   .. code-block:: python

       Iv = sim.Ivol()

.. function:: mat Idev()

    Provides the deviatoric of the identity tensor :math:`I_{dev}` written in the Simcoon formalism. So :
    
     .. math:: 
     
     	I_{dev} = I_{real} - I_{vol} = \left( \begin{array}{ccc}
        2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)

    .. code-block:: python

        Id = sim.Idev()

.. function:: mat Ireal2()

    Provides the fourth order identity tensor :math:`\widehat{I}` written in the form. So :

    .. math::

        \widehat{I} = \left( \begin{array}{ccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)

   For example, this tensor allows to obtain : :math:`L*\widehat{M}=I` or :math:`\widehat{L}*M=I`, where a matrix :math:`\widehat{A}` is set by :math:`\widehat{A}=\widehat{I}A\widehat{I}`

   .. code-block:: python

        Ir2 = sim.Ireal2()

.. function:: mat Idev2()

    Provides the deviatoric of the identity tensor :math:`\widehat{I}` written in the Simcoon formalism. So :

    .. math::

        I_{dev2} = \left( \begin{array}{ccc}
        2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)

    .. code-block:: python

        Id2 = sim.Idev2()

.. function:: vec Ith()

    Provide the vector :math:`I_{th} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    0 \\
    0 \\
    0 \end{array} \right)`

    .. code-block:: python

        It = sim.Ith()

.. function:: vec Ir2()

    Provide the vector :math:`I_{r2} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    2 \\
    2 \\
    2 \end{array} \right)`

    .. code-block:: python

        I2 = sim.Ir2()

.. function:: vec Ir05()

    Provide the vector :math:`I_{r05} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    0.5 \\
    0.5 \\
    0.5 \end{array} \right)`

    .. code-block:: python

        I05 = sim.Ir05()

.. function:: mat L_iso(const double &C1, const double &C2, const std::string &conv)

    Provides the elastic stiffness tensor for an isotropic material.
    The two first arguments are a couple of elastic properties. The third argument specifies which couple has been provided and the nature and order of coefficients.
    Exhaustive list of possible third argument :
    ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.

    .. code-block:: python

        E = 210000.0
        nu = 0.3;
        Liso = sim.L_iso(E, nu, "Enu")

.. function:: mat M_iso(const double &C1, const double &C2, const string &conv)

    Provides the elastic compliance tensor for an isotropic material.
    The two first arguments are a couple of elastic properties. The third argument specify which couple has been provided and the nature and order of coefficients.
    Exhaustive list of possible third argument :
    ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.

    .. code-block:: python

        E = 210000.0
        nu = 0.3
        Miso = sim.M_iso(E, nu, "Enu")

.. function:: mat L_cubic(const double &C1, const double &C2, const double &C4, const string &conv)

    Provides the elastic stiffness tensor for a cubic material.
    The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.

    .. code-block:: python

        E = 70000.0
        nu = 0.3
        G = 23000.0
        Lcubic = sim.L_cubic(E, nu, G, "EnuG")

        import numpy as np
        C11 = np.random.uniform(10000., 100000.)
        C12 = np.random.uniform(10000., 100000.)
        C44 = np.random.uniform(10000., 100000.)
        Lcubic = sim.L_cubic(C11, C12, C44, "Cii")

.. function:: mat M_cubic(const double &C1, const double &C2, const double &C4, const string &conv)

    Provides the elastic compliance tensor for a cubic material.
    The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.

    .. code-block:: python

        E = 70000.0
        nu = 0.3
        G = 23000.0
        Lcubic = sim.L_cubic(E, nu, G, "EnuG")

        C11 = np.random.uniform(10000., 100000.)
        C12 = np.random.uniform(10000., 100000.)
        C44 = np.random.uniform(10000., 100000.)
        Mcubic = M_cubic(C11, C12, C44, "Cii")

.. function:: mat L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv)

    Provides the elastic stiffness tensor for an orthotropic material.
    Arguments could be all the stiffness coefficients or the material parameter. For an orthotropic material the material parameters should be : Ex,Ey,Ez,nuxy,nuyz,nxz,Gxy,Gyz,Gxz.

    The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.

    .. code-block:: python

        C11 = np.random.uniform(10000., 100000.)
        C12 = np.random.uniform(10000., 100000.)
        C13 = np.random.uniform(10000., 100000.)
        C22 = np.random.uniform(10000., 100000.)
        C23 = np.random.uniform(10000., 100000.)
        C33 = np.random.uniform(10000., 100000.)
        C44 = np.random.uniform(10000., 100000.)
        C55 = np.random.uniform(10000., 100000.)
        C66 = np.random.uniform(10000., 100000.)
        Lortho = sim.L_ortho(C11, C12, C13, C22, C23, C33, C44, C55, C66, "Cii")

.. function:: mat M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv)


    Provides the elastic compliance tensor for an orthotropic material.
    Arguments could be all the stiffness coefficients or the material parameter. For an orthotropic material the material parameters should be : Ex,Ey,Ez,nuxy,nuyz,nxz,Gxy,Gyz,Gxz.

    The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.

   .. code-block:: python

       C11 = np.random.uniform(10000., 100000.)
       C12 = np.random.uniform(10000., 100000.)
       C13 = np.random.uniform(10000., 100000.)
       C22 = np.random.uniform(10000., 100000.)
       C23 = np.random.uniform(10000., 100000.)
       C33 = np.random.uniform(10000., 100000.)
       C44 = np.random.uniform(10000., 100000.)
       C55 = np.random.uniform(10000., 100000.)
       C66 = np.random.uniform(10000., 100000.)
       Mortho = sim.M_ortho(C11, C12, C13, C22, C23, C33, C44, C55, C66, "Cii")

.. function:: mat L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis)

    Provides the elastic stiffness tensor for an isotropic transverse material.
    Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.

    .. code-block:: python

        EL = np.random.uniform(10000., 100000.)
        ET = np.random.uniform(10000., 100000.)
        nuTL = np.random.uniform(0., 0.5)
        nuTT = np.random.uniform(0., 0.5)
        GLT = np.random.uniform(10000., 100000.)
        axis = 1
        Lisotrans = sim.L_isotrans(EL, ET, nuTL, nuTT, GLT, axis)

.. function:: mat M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis)

    Provides the elastic compliance tensor for an isotropic transverse material.
    Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.

    .. code-block:: python

        EL = np.random.uniform(10000., 100000.)
        ET = np.random.uniform(10000., 100000.)
        nuTL = np.random.uniform(0., 0.5)
        nuTT = np.random.uniform(0., 0.5)
        GLT = np.random.uniform(10000., 100000.)
        axis = 1
        Misotrans = sim.M_isotrans(EL, ET, nuTL, nuTT, GLT, axis)

.. function:: mat H_iso(const double &etaB, const double &etaS)

    Provides the viscoelastic tensor H, providing Bulk viscosity etaB and shear viscosity etaS. 
    It actually returns :
    
    .. math::

        H_iso = \left( \begin{array}{ccc}
        \eta_B & \eta_B & \eta_B & 0 & 0 & 0 \\
        \eta_B & \eta_B & \eta_B & 0 & 0 & 0 \\
        \eta_B & \eta_B & \eta_B & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)
    
    
    .. code-block:: python

        etaB = np.random.uniform(0., 1.)
        etaS = np.random.uniform(0., 1.)
        Hiso = sim.H_iso(etaB, etaS)
