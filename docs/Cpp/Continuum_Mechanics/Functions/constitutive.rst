The Constitutive Library
========================

.. default-domain:: cpp

.. cpp:function:: mat Ireal()

    :parameter: None

    :Description:
        Provides the fourth order identity tensor written in Voigt notation :math:`I_{real}`, where :

    .. math::

        I_{real} = \left( \begin{array}{cccccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)

    :return: The above 6x6 mat (arma::mat)

    :example:

    .. code-block:: cpp

        mat Ir = Ireal();

.. cpp:function:: mat Ivol()

    :parameter: None

    :Description:

        Provides the volumic of the identity tensor :math:`I_{vol}` written in the Simcoon formalism. So :

    .. math::

        I_{vol} = \left( \begin{array}{ccc}
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \end{array} \right)

    :return: The above 6x6 mat (arma::mat)

    :example:

   .. code-block:: cpp

       mat Iv = Ivol();

.. cpp:function:: mat Idev()

    :parameter: None

    :Description:
    
        Provides the deviatoric of the identity tensor :math:`I_{dev}` written in the Simcoon formalism. So :
    
    .. math:: 
     
     	I_{dev} = I_{real} - I_{vol} = \left( \begin{array}{ccc}
        2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)

    :return: The above 6x6 mat (arma::mat)

    :example:

    .. code-block:: cpp

        mat Id = Idev();

.. cpp:function:: mat Ireal2()

    :parameter: None

    :Description:

        Provides the fourth order identity tensor :math:`\widehat{I}` written in the form. So :

    .. math::

        \widehat{I} = \left( \begin{array}{ccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)

   For example, this tensor allows to obtain : :math: `L*\widehat{M}=I` or :math:`\widehat{L}*M=I`, where a matrix :math:`\widehat{A}` is set by :math:`\widehat{A}=\widehat{I}\,A\,\widehat{I}`

    :return: The above 6x6 mat (arma::mat)

    :example: 

   .. code-block:: cpp

        mat Ir2 = Ireal2();

.. cpp:function:: mat Idev2()

    Provides the deviatoric of the identity tensor :math: `\widehat{I}` written in the Simcoon formalism. So :

    .. math::

        I_{dev2} = \left( \begin{array}{ccc}
        2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)

    .. code-block:: cpp

        mat Id2 = Idev2();

.. function:: vec Ith()

    Provide the vector :math:`I_{th} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    0 \\
    0 \\
    0 \end{array} \right)`

    .. code-block:: cpp

        vec It = Ith();

.. function:: vec Ir2()

    Provide the vector :math:`I_{r2} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    2 \\
    2 \\
    2 \end{array} \right)`

    .. code-block:: cpp

        vec I2 = Ir2();

.. function:: vec Ir05()

    Provide the vector :math:`I_{r05} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    0.5 \\
    0.5 \\
    0.5 \end{array} \right)`

    .. code-block:: cpp

        vec I05 = Ir05();

.. function:: mat L_iso(const double &C1, const double &C2, const std::string &conv)

    Provides the elastic stiffness tensor for an isotropic material.
    The two first arguments are a couple of elastic properties. The third argument specifies which couple has been provided and the nature and order of coefficients.
    Exhaustive list of possible third argument :
    ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.

    .. code-block:: cpp

        double E = 210000;
        double nu = 0.3;
        mat Liso = L_iso(E, nu, "Enu");

.. function:: mat M_iso(const double &C1, const double &C2, const string &conv)

    Provides the elastic compliance tensor for an isotropic material.
    The two first arguments are a couple of elastic properties. The third argument specify which couple has been provided and the nature and order of coefficients.
    Exhaustive list of possible third argument :
    ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.

    .. code-block:: cpp

        double E = 210000;
        double nu = 0.3;
        mat Miso = M_iso(E, nu, "Enu");

.. function:: mat L_cubic(const double &C1, const double &C2, const double &C4, const string &conv)

    Provides the elastic stiffness tensor for a cubic material.
    Arguments are the stiffness coefficients C11, C12 and C44.

    .. code-block:: cpp

        double C11 = alead(10000., 100000.);
        double C12 = alead(10000., 100000.);
        double C44 = alead(10000., 100000.);
        mat Lcubic = L_cubic(C11, C12, C44, "Cii");

.. function:: mat M_cubic(const double &C1, const double &C2, const double &C4, const string &conv)

    Provides the elastic compliance tensor for a cubic material.
    Arguments are the stiffness coefficients C11, C12 and C44.

    .. code-block:: cpp

        double C11 = alead(10000., 100000.);
        double C12 = alead(10000., 100000.);
        double C44 = alead(10000., 100000.);
        mat Mcubic = M_cubic(C11,C12,C44);

.. function:: mat L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv)

    Provides the elastic stiffness tensor for an orthotropic material.
    Arguments could be all the stiffness coefficients or the material parameter. For an orthotropic material the material parameters should be : Ex,Ey,Ez,nuxy,nuyz,nxz,Gxy,Gyz,Gxz.

    The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.

    .. code-block:: cpp

        double C11 = alead(10000., 100000.);
        double C12 = alead(10000., 100000.);
        double C13 = alead(10000., 100000.);
        double C22 = alead(10000., 100000.);
        double C23 = alead(10000., 100000.);
        double C33 = alead(10000., 100000.);
        double C44 = alead(10000., 100000.);
        double C55 = alead(10000., 100000.);
        double C66 = alead(10000., 100000.);
        mat Lortho = L_ortho(C11, C12, C13, C22, C23, C33, C44, C55, C66,"Cii");

.. function:: mat M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv)


    Provides the elastic compliance tensor for an orthotropic material.
    Arguments could be all the stiffness coefficients or the material parameter. For an orthotropic material the material parameters should be : Ex,Ey,Ez,nuxy,nuyz,nxz,Gxy,Gyz,Gxz.

    The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.

   .. code-block:: cpp

       double C11 = alead(10000., 100000.);
       double C12 = alead(10000., 100000.);
       double C13 = alead(10000., 100000.);
       double C22 = alead(10000., 100000.);
       double C23 = alead(10000., 100000.);
       double C33 = alead(10000., 100000.);
       double C44 = alead(10000., 100000.);
       double C55 = alead(10000., 100000.);
       double C66 = alead(10000., 100000.);
       mat Mortho = M_ortho(C11, C12, C13, C22, C23, C33, C44, C55, C66,"Cii");

.. function:: mat L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis)

    Provides the elastic stiffness tensor for an isotropic transverse material.
    Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.

    .. code-block:: cpp

        double EL = alead(10000., 100000.);
        double ET = alead(10000., 100000.);
        double nuTL = alead(0., 0.5);
        double nuTT = alead(0.5, 0.5);
        double GLT = alead(10000., 100000.);
        double axis = 1;
        mat Lisotrans = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);

.. function:: mat M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis)

    Provides the elastic compliance tensor for an isotropic transverse material.
    Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.

    .. code-block:: cpp

        double EL = alead(10000., 100000.);
        double ET = alead(10000., 100000.);
        double nuTL = alead(0., 0.5);
        double nuTT = alead(0., 0.5);
        double GLT = alead(10000., 100000.);
        double axis = 1;
        mat Misotrans = M_isotrans(EL, ET, nuTL, nuTT, GLT, axis);

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
    
    
    .. code-block:: cpp

        double etaB = alead(0., 0.1);
        double etaS = alead(0., 0.1);
        mat Hiso = H_iso(etaB, etaS);

.. function:: void el_pred(see below)

    Provides the stress tensor from an elastic prediction
    There are two possible ways:

    1. From the elastic stiffness tensor and the trial elastic strain:
    parameters : L : Stiffness matrix; Eel ; elastic strain vector, ndi (optional, default = 3): number of dimensions
    (const mat &L, const vec &E_el, const int &ndi)

    .. code-block:: cpp
        
        mat L = L_iso(70000, 0.3,"Enu");
        vec Eel;
        Eel.randu(6);
        int ndi = 3;
        vec sigma =  el_pred(L, Eel, ndi);

    2. From the previous stress increment, providing the elastic stiffness tensor and the trial elastic strain increment:
    parameters : sigma_start: The previous stress, L : Stiffness matrix; Eel : elastic strain vector, ndi (optional, default = 3): number of dimensions
    (const vec &sigma_start, const mat &L, const vec &DE_el, const int &ndi)

    .. code-block:: cpp
        
        vec sigma_start = zeros(6);
        sigma_start.randu(6);
        mat L = L_iso(70000, 0.3,"Enu");
        vec Eel;
        Eel.randu(6);
        int ndi = 3;
        vec sigma =  el_pred(sigma_start,L, Eel, ndi);

.. function:: mat Isotropize(const mat &Lt)

    Provides an isotropized version of an anisotropic stiffness tensor. Such isotropic tensor is called consistent since for any given strain it return the same stress as the anisotropic version.

    .. code-block:: cpp

        double EL = (double)rand();
        double ET = (double)rand();
        double nuTL = (double)rand();
        double nuTT = (double)rand();
        double GLT = (double)rand();
        double axis = 1;
        mat L_isotrans = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);
        mat L_iso = Isotropize(Lisotrans);
