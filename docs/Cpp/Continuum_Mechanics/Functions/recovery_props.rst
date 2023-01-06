The Recovery Props Library
========================

The recovery props library provides a set of function to check and evaluate the properties of stiffness and compliance tensors.

.. default-domain:: cpp

.. function:: void check_symetries(const mat &L, std::string &umat_type, int &axis, vec &props, int &maj_sym)

    Check the symmetries of a stiffness matrix L, and fill the vector of material properties. 
    Depending on the symmetry found, the string umat_type, the axis of symmetry (if applicable) the vector of material properties, and the  major symmetry maj_sym (L_ij = L_ji ?).
    If the major symmetry condition is not fulfilled, the check of symmetries if performed on the symmetric part of L. For fully anisotropic and monoclinic symmetries, the vector f parameters is not returned (the full stiffness tensor is generally directly utilized)
    
	.. csv-table:: Material Symmetries considered
   		:header: "Symmetry", "umat_type", "axis", "size of props"
   		:widths: 60, 20, 30,20

		"Fully anisotropic", "ELANI", "0", "N/A"
   		"Monoclinic", "ELMON", "1,2 or 3", "N/A"   
   		"Orthotropic", "ELORT", "0", "9"      
   		"Cubic", "ELCUB", "0", "3"  
   		"Transversely isotropic", "ELITR", "1,2 or 3", "5"  
   		"Isotropic", "ELISO", "0", "2"        
    
    .. code-block:: cpp

        check_symetries(L, umat_name, axis, props, maj_sym);
        
.. function:: vec L_iso_props(const mat &L)

Returns a vector containing the Young's modulus and the Poisson ratio :math:`\left(E, \nu \right)` of a linear elastic isotropic material, providing the stiffness matrix :math:`\mathbf{L}`. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to an isotropic material)

.. code-block:: cpp

    mat L = L_iso(70000., 0.3, 'Enu');
    vec eliso_props = L_iso_props(L);

.. function:: vec M_iso_props(const mat &M)

Returns a vector containing the Young's modulus and the Poisson ratio :math:`\left(E, \nu \right)` of a linear elastic isotropic material, providing the compliance matrix :math:`\mathbf{M}`. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to an isotropic material)

.. code-block:: cpp

    mat M = M_iso(70000., 0.3, 'Enu');
    vec eliso_props = M_iso_props(M);

.. function:: vec L_isotrans_props(const mat &L)

Returns a vector containing the material properties :math:`\left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right)` of a linear elastic transversely isotropic material, providing the stiffness matrix :math:`\mathbf{L}` and the axis of symmetry. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    int axis = 1;
    double E_L = 4500;
    double E_T = 2300;
    double nu_TL = 0.05;
    double nu_TT = 0.3;
    double G_LT = 2700;
    mat L = L_isotrans(E_L, E_T, nu_TL, nu_TT, G_LT., axis);
    vec isotrans_props = L_isotrans_props(L, axis);
    
.. function:: vec M_isotrans_props(const mat &M)

Returns a vector containing the material properties :math:`\left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right)` of a linear elastic transversely isotropic material, providing the compliance matrix :math:`\mathbf{M}` and the axis of symmetry. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    int axis = 1;
    double E_L = 4500;
    double E_T = 2300;
    double nu_TL = 0.05;
    double nu_TT = 0.3;
    double G_LT = 2700;
    mat M = M_isotrans(E_L, E_T, nu_TL, nu_TT, G_LT., axis);
    vec isotrans_props = M_isotrans_props(M, axis);
    
.. function:: vec L_cubic_props(const mat &L)

Returns a vector containing the material properties :math:`\left(E, \nu, G \right)` of a linear elastic, providing the stiffness matrix :math:`\mathbf{L}`. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    mat L = L_cubic(185000., 158000., 39700., 'Cii') //C11, C12, C44
    vec cubic_props = L_cubic_props(L);

.. function:: vec M_cubic_props(const mat &M)

Returns a vector containing the material properties :math:`\left(E, \nu, G \right)` of a linear elastic, providing the compliance matrix :math:`\mathbf{M}`. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    mat M = M_cubic(185000., 158000., 39700., 'Cii') //C11, C12, C44
    vec cubic_props = M_cubic_props(M);

.. function:: vec L_ortho_props(const mat &L)

Returns a vector containing the material properties :math:`\left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23} \right)` of a linear elastic orthotropic material, providing the stiffness matrix :math:`\mathbf{L}` and the axis of symmetry. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    double E_1 = 4500;
    double E_2 = 2300;
    double E_3 = 2700;
    double nu_12 = 0.06;
    double nu_13 = 0.08;
    double nu_23 = 0.3;
    double G_12 = 2200;
    double G_13 = 2100;
    double G_23 = 2400;
    mat L = L_ortho(E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23);
    vec ortho_props = L_ortho_props(L);
    
.. function:: vec M_ortho_props(const mat &M)

Returns a vector containing the material properties :math:`\left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23} \right)` of a linear elastic orthotropic material, providing the stiffness matrix :math:`\mathbf{L}` and the axis of symmetry. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    double E_1 = 4500;
    double E_2 = 2300;
    double E_3 = 2700;
    double nu_12 = 0.06;
    double nu_13 = 0.08;
    double nu_23 = 0.3;
    double G_12 = 2200;
    double G_13 = 2100;
    double G_23 = 2400;
    mat L = L_ortho(E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23);
    vec ortho_props = L_ortho_props(L);
    
.. function:: vec M_aniso_props(const mat &M)

Returns a vector containing the material properties :math:`\left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23}, \eta_{14}, \eta_{15}, \eta_{16}, \eta_{24}, \eta_{25}, \eta_{26}, \eta_{34}, \eta_{35}, \eta_{36}, \eta_{45}, \eta_{46}, \eta_{56} \right)` of a linear elastic fully anisotropic material, providing the stiffness matrix :math:`\mathbf{L}` and the axis of symmetry. Note that an averaging over the component is operated (usefull when the provided matrix do not exactly correspond to a transversely isotropic material)

.. code-block:: cpp

    string umat_name;
    string path_data = "data";
    string materialfile = "material.dat";
    
    unsigned int nprops = 0;
    unsigned int nstatev = 0;
    vec props;
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    double T_init = 273.15;

    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    phase_characteristics rve;
    
    rve.construct(0,1);
    natural_basis nb;
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(3,3), zeros(3,3), eye(3,3), eye(3,3),T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
    
    auto sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);

    //Second we call a recursive method that find all the elastic moduli iof the phases
    get_L_elastic(rve);
    mat M = arma::inv(sv_M->Lt);
    vec aniso_props = M_aniso_props(M);
