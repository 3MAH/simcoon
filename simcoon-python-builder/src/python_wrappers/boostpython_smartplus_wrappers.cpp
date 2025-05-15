
#include <armadillo>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/constitutive.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/contimech.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/stress.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/criteria.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/damage.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/hyperelastic.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/recovery_props.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Leff.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/objective_rates.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/umat.hpp>
//#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/RunUmat.hpp>

#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/lagrange.hpp>
#include <simcoon/python_wrappers/Libraries/Material/ODF.hpp>
#include <simcoon/python_wrappers/Libraries/Homogenization/eshelby.hpp>

#include <simcoon/python_wrappers/Libraries/Solver/read.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver.hpp>
//#include <simcoon/python_wrappers/Libraries/Solver/step_meca.hpp>
//#include <simcoon/python_wrappers/Libraries/Solver/step_thermomeca.hpp>

#include <simcoon/python_wrappers/Libraries/Identification/identification.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/constants.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/parameters.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/optimize.hpp>

#include <simcoon/docs/Libraries/Continuum_mechanics/doc_constitutive.hpp>

using namespace std;
using namespace arma;
using namespace simpy;

using namespace pybind11::literals;

PYBIND11_MODULE(simmit, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    // Register the from-python converters for constitutive.hpp
    m.def("Ireal", &Ireal, "copy"_a=true, simcoon_docs::Ireal);
    m.def("Ivol", &Ivol, "copy"_a=true, simcoon_docs::Ivol);
    m.def("Idev", &Idev, "copy"_a=true, simcoon_docs::Idev);
    m.def("Ireal2", &Ireal2, "copy"_a=true, simcoon_docs::Ireal2);
    m.def("Idev2", &Idev2, "copy"_a=true, simcoon_docs::Idev2);
    m.def("Ith", &Ith, "copy"_a=true, simcoon_docs::Ith);
    m.def("Ir2", &Ir2, "copy"_a=true, simcoon_docs::Ir2);
    m.def("Ir05", &Ir05, "copy"_a=true, simcoon_docs::Ir05);
    m.def("L_iso", &L_iso, "props"_a, "conv"_a, "copy"_a=true, simcoon_docs::L_iso);
    m.def("M_iso", &M_iso, "props"_a, "conv"_a, "copy"_a=true, simcoon_docs::M_iso);
    m.def("L_cubic", &L_cubic, "props"_a, "conv"_a, "copy"_a=true, simcoon_docs::L_cubic);
    m.def("M_cubic", &M_cubic, "props"_a, "conv"_a, "copy"_a=true, simcoon_docs::M_cubic);
    m.def("L_ortho", &L_ortho, "props"_a, "conv"_a, "copy"_a=true, simcoon_docs::L_ortho);
    m.def("M_ortho", &M_ortho, "props"_a, "conv"_a, "copy"_a=true, simcoon_docs::M_ortho);
    m.def("L_isotrans", &L_isotrans, "props"_a, "axis"_a, "copy"_a=true, simcoon_docs::L_isotrans);
    m.def("M_isotrans", &M_isotrans, "props"_a, "axis"_a, "copy"_a=true, simcoon_docs::M_isotrans);
    m.def("H_iso", &H_iso, "props"_a, "copy"_a=true, simcoon_docs::H_iso);

    // Register the from-python converters for contimech
    m.def("tr", &tr, "input"_a, "This function returns the trace of a tensor");
    m.def("dev", &dev, "input"_a, "copy"_a=true, "This function returns the deviatoric part of a tensor");
    m.def("Mises_stress", &Mises_stress, "input"_a, "This function determines the Mises equivalent of a stress tensor");
    m.def("eta_stress", &eta_stress, "input"_a, "copy"_a=true, "This function determines the strain flow (direction) from a stress tensor");
    m.def("Mises_strain", &Mises_strain, "input"_a, "This function determines the Mises equivalent of a strain tensor");
    m.def("eta_strain", &eta_strain, "input"_a, "copy"_a=true, "This function determines the strain flow (direction) from a strain tensor");
    m.def("J2_stress", &J2_stress, "input"_a, "Returns the second invariant of the deviatoric part of a second order stress tensor");
    m.def("J2_strain", &J2_strain, "input"_a, "Returns the second invariant of the deviatoric part of a second order strain tensor");
    m.def("J3_stress", &J3_stress, "input"_a, "Returns the third invariant of the deviatoric part of a second order stress tensor");
    m.def("J3_strain", &J3_strain, "input"_a, "Returns the third invariant of the deviatoric part of a second order stress tensor");
    m.def("Macaulay_p", &Macaulay_p, "value"_a, "This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)");
    m.def("Macaulay_n", &Macaulay_n, "value"_a, "This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)");
    m.def("sign", &simpy::sign, "value"_a, "This function returns the sign of a double");
    m.def("normal_ellipsoid", &normal_ellipsoid, "u"_a, "v"_a, "a1"_a, "a2"_a, "a3"_a, "copy"_a=true, "Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3");
    m.def("curvature_ellipsoid", &curvature_ellipsoid, "u"_a, "v"_a, "a1"_a, "a2"_a, "a3"_a, "Provides the curvature of an ellipsoid with semi-principal axes of length a1, a2, a3 at the angle u,v.");        
    m.def("sigma_int", &sigma_int, "input"_a, "u"_a, "v"_a, "a1"_a, "a2"_a, "a3"_a, "copy"_a=true, "Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3");
    m.def("p_ikjl", &p_ikjl, "normal"_a, "copy"_a=true, "This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)");
    m.def("auto_sym_dyadic", &auto_sym_dyadic, "input"_a, "copy"_a=true, "Provides the dyadic product of a symmetric tensor with itself (auto dyadic product)");
    m.def("sym_dyadic", &sym_dyadic, "a"_a, "b"_a, "copy"_a=true, "Provides the dyadic product of two symmetric tensors"); 
    m.def("auto_dyadic", &auto_dyadic, "input"_a, "copy"_a=true, "Provides the dyadic product of a tensor with itself (auto dyadic product)");
    m.def("dyadic", &dyadic, "a"_a, "b"_a, "copy"_a=true, "Provides the dyadic product of of two symmetric tensors");

    // Register the from-python converters for criteria
    m.def("Prager_stress", &Prager_stress, "input"_a, "props"_a, "This function returns the Prager equivalent stress");
    m.def("dPrager_stress", &dPrager_stress, "input"_a, "props"_a, "copy"_a=true, "This function returns the derivative of the Prager equivalent stress");
    m.def("Tresca_stress", &Tresca_stress, "input"_a, "This function returns the Tresca equivalent stress");
    m.def("dTresca_stress", &dTresca_stress, "input"_a, "copy"_a=true, "This function returns the derivative of the Tresca equivalent stress");
    m.def("P_Ani", &P_Ani, "props"_a, "copy"_a=true, "Returns an anisotropic configurational tensor P in the Voigt format (6x6 numpy array), given its vector representation");
    m.def("P_Hill", &P_Hill, "props"_a, "copy"_a=true, "Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion in the Voigt format (6x6 numpy array), given its vector representation");
    m.def("Hill_stress", &Hill_stress, "input"_a, "props"_a, "This function returns the Hill equivalent stress");
    m.def("dHill_stress", &dHill_stress, "input"_a, "props"_a, "copy"_a=true, "This function returns the derivative of the Hill equivalent stress");
    m.def("Ani_stress", &Ani_stress, "input"_a, "props"_a, "This function returns the Ani equivalent stress");
    m.def("dAni_stress", &dAni_stress, "input"_a, "props"_a, "copy"_a=true, "This function returns the derivative of the Ani equivalent stress");
    m.def("DFA_stress", &DFA_stress, "input"_a, "props"_a, "This function returns the DFA equivalent stress");
    m.def("dDFA_stress", &dDFA_stress, "input"_a, "props"_a, "copy"_a=true, "This function returns the derivative of the DFA equivalent stress");
    m.def("Eq_stress", &Eq_stress, "input"_a, "criteria"_a, "props"_a, "This function computes the selected equivalent stress function");
    m.def("dEq_stress", &dEq_stress, "input"_a, "criteria"_a, "props"_a, "copy"_a=true, "This function computes the deriavtive of the selected equivalent stress function");

    // Register the from-python converters for recovery_props
    m.def("check_symetries", &check_symetries, "input"_a, "tol"_a, "Check the material symetries and the type of elastic response for a given stiffness tensor");
    m.def("L_iso_props", &L_iso_props, "input"_a, "Return a list of elastic properties for the isotropic case (E,nu) from a stiffness tensor");
    m.def("M_iso_props", &M_iso_props, "input"_a, "Return a list of elastic properties for the isotropic case (E,nu) from a compliance tensor");
    m.def("L_isotrans_props", &L_isotrans_props, "input"_a, "axis"_a, "Return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a stiffness tensor");
    m.def("M_isotrans_props", &M_isotrans_props, "input"_a, "axis"_a, "Return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a compliance tensor");
    m.def("L_cubic_props", &L_cubic_props, "input"_a, "Return a list of elastic properties for the cubic case (E,nu,G) from a stiffness tensor");
    m.def("M_cubic_props", &M_cubic_props, "input"_a, "Return a list of elastic properties for the cubic case (E,nu,G) from a compliance tensor");
    m.def("L_ortho_props", &L_ortho_props, "input"_a, "Return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a stiffness tensor");
    m.def("M_ortho_props", &M_ortho_props, "input"_a, "Return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a compliance tensor");
    m.def("M_aniso_props", &M_aniso_props, "input"_a, "Return a list of elastic properties for the anisotropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,deviations) from a compliance tensor");

    // Register the L_eff for composites
    m.def("L_eff", &L_eff, "umat_name"_a, "props"_a, "nstatev"_a, "psi_rve"_a=0., "theta_rve"_a=0., "phi_rve"_a=0., "Return the elastic stiffness tensor of a composite material");

    //Register the from-python converters for kinematics
    m.def("ER_to_F", &ER_to_F, "E"_a, "R"_a, "copy"_a=true, "Provides the transformation gradient, from the Green-Lagrange strain and the rotation");
    m.def("eR_to_F", &eR_to_F, "e"_a, "R"_a, "copy"_a=true, "Provides the transformation gradient, from the logarithmic strain and the rotation");
    m.def("G_UdX", &G_UdX, "F"_a, "copy"_a=true, "This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor");
    m.def("G_Udx", &G_Udx, "F"_a, "copy"_a=true, "This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor");
    m.def("R_Cauchy_Green", &R_Cauchy_Green, "F"_a, "copy"_a=true, "This function computes the Right Cauchy-Green C");
    m.def("L_Cauchy_Green", &L_Cauchy_Green, "F"_a, "copy"_a=true, "This function computes the Left Cauchy-Green B");
    m.def("RU_decomposition", &RU_decomposition, "F"_a, "copy"_a=true, "Provides the RU decomposition of the transformation gradient F");
    m.def("VR_decomposition", &VR_decomposition, "F"_a, "copy"_a=true, "Provides the VR decomposition of the transformation gradient F");    
    m.def("Inv_X", &Inv_X, "input"_a, "copy"_a=true, "This function computes the common Right (or Left) Cauchy-Green invariants");
    m.def("Cauchy", &Cauchy, "F"_a, "copy"_a=true, "This function computes the Cauchy deformation tensor c from the transformation gradient F");
    m.def("Green_Lagrange", &Green_Lagrange, "F"_a, "copy"_a=true, "This function computes the Green-Lagrange finite strain tensor E");
    m.def("Euler_Almansi", &Euler_Almansi, "F"_a, "copy"_a=true, "This function computes the Euler-Almansi finite strain tensor A");
    m.def("Log_strain", &Log_strain, "F"_a, "voigt_form"_a=false, "copy"_a=true, "This function computes the logarithmic strain ln[V] = 1/2 ln[B] (B is the left Cauchy-Green Tensor)");
    m.def("finite_L", &finite_L, "F0"_a, "F1"_a, "DTime"_a,  "copy"_a=true, "This function computes the velocity difference (F0,F1,DTime)");
    m.def("finite_D", &finite_D, "F0"_a, "F1"_a, "DTime"_a,  "copy"_a=true, "This function computes the deformation rate D (F0,F1,DTime)");
    m.def("finite_W", &finite_W, "F0"_a, "F1"_a, "DTime"_a,  "copy"_a=true, "This function computes the spin tensor W (correspond to Jaumann rate) (F0,F1,DTime)");
    m.def("finite_Omega", &finite_Omega, "F0"_a, "F1"_a, "DTime"_a,  "copy"_a=true, "This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)");
    m.def("finite_DQ", &finite_DQ, "Omega0"_a, "Omega0"_a, "DTime"_a,  "copy"_a=true, "This function computes the increment of finite rotation (Omega0, Omega1, DTime)");

    //register the objective rates library
    m.def("logarithmic", &logarithmic, "F0"_a, "F1"_a, "DTime"_a, "copy"_a=true, "This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment");
    m.def("logarithmic_R", &logarithmic_R, "F0"_a, "F1"_a, "DTime"_a, "copy"_a=true, "This function computes the logarithmic strain velocity and the Green-Naghdi spin, along with the correct rotation increment");    
    m.def("Delta_log_strain", &Delta_log_strain, "D"_a, "Omega"_a, "DTime"_a, "copy"_a=true, "This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor");
    m.def("objective_rate", &objective_rate, "corate_name"_a,"F0"_a, "F1"_a, "dtime"_a, "return_de"_a=false, "n_threads"_a = 4, "This function computes the strain velocity and the spin, along with the correct rotation increment for the specified objective erivative");
    m.def("Lt_convert", &Lt_convert, "Lt"_a, "F"_a, "stress"_a, "converter_key"_a);

    //register the damage library
    m.def("damage_weibull", &damage_weibull, "stress"_a, "damage"_a, "alpha"_a, "beta"_a, "DTime"_a, "criterion"_a = "vonmises", "This function returns damage evolution (/dt) considering a Weibull damage law");
    m.def("damage_kachanov", &damage_kachanov, "stress"_a, "strain"_a, "damage"_a, "A0"_a, "r"_a, "criterion"_a, "This function returns damage evolution (/dt) considering Kachanov's creep damage law");
    m.def("damage_miner", &damage_miner, "S_max"_a, "S_mean"_a, "S_ult"_a, "b"_a, "B0"_a, "beta"_a, "Sl_0"_a = 0., "This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law");
    m.def("damage_manson", &damage_manson, "S_amp"_a, "C2"_a, "gamma2"_a, "This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law");

    //register the hyperelastic library
    m.def("isochoric_invariants", &isochoric_invariants, "input"_a, "J"_a=0., "copy"_a=true, "This function computes the isochoric invariants");
    m.def("isochoric_pstretch", &isochoric_pstretch, "input"_a, "input_tensor"_a="V", "J"_a=0., "copy"_a=true, "This function computes the isochoric invariants");    

    //register the transfer library
    m.def("v2t_strain", &v2t_strain, "input"_a, "copy"_a=true, "This function transforms the strain Voigt vector into a 3*3 strain matrix");
    m.def("t2v_strain", &t2v_strain, "input"_a, "copy"_a=true, "This function transforms a 3*3 strain matrix into a strain Voigt vector");    
    m.def("v2t_stress", &v2t_stress, "input"_a, "copy"_a=true, "This function transforms the stress Voigt vector into a 3*3 stress matrix");
    m.def("v2t_stress", &v2t_stress, "input"_a, "copy"_a=true, "This function transforms a 3*3 stress matrix into a stress Voigt vector");

    // Register the from-python converters for eshelby
    m.def("Eshelby_sphere", &Eshelby_sphere, "nu"_a, "copy"_a=true, "Eshelby tensor for a sphere");
    m.def("Eshelby_cylinder", &Eshelby_cylinder, "nu"_a, "copy"_a=true, "Eshelby tensor for a cylinder. The cylinder is oriented in such a way that the axis direction is the 1 direction. a2=a3 here");
    m.def("Eshelby_prolate", &Eshelby_prolate, "nu"_a, "aspect_ratio"_a, "copy"_a=true, "Eshelby tensor for a prolate ellipsoid. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here");
    m.def("Eshelby_oblate", &Eshelby_oblate, "nu"_a, "aspect_ratio"_a, "copy"_a=true, "Eshelby tensor for an oblate ellipsoid. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here");
    m.def("Eshelby", &Eshelby, "L"_a, "a1"_a=1., "a2"_a=1., "a3"_a=1., "mp"_a=50, "np"_a=50, "copy"_a=true, "Numerical Eshelby tensor determination");
    m.def("T_II", &T_II, "L"_a, "a1"_a=1., "a2"_a=1., "a3"_a=1., "mp"_a=50, "np"_a=50, "copy"_a=true, "Numerical Hill Interaction tensor determination");

    //Register the rotation library
    m.def("rotate_vec_R", &rotate_vec_R, "input"_a, "R"_a, "copy"_a=true, "This function returns a rotated vector (3) according to a rotation matrix");
    m.def("rotate_vec_angle", &rotate_vec_angle, "input"_a, "angle"_a, "axis"_a, "copy"_a=true, "This function returns a rotated vector (3) according to an angle and an axis");
    m.def("rotate_mat_R", &rotate_mat_R, "input"_a, "R"_a, "copy"_a=true, "This function returns a rotated matrix (3x3) according to a rotation matrix");        
    m.def("rotate_mat_angle", &rotate_mat_angle, "input"_a, "angle"_a, "axis"_a, "copy"_a=true, "This function returns a rotated matrix (3x3) according to an angle and an axis");
    m.def("fillR_angle", &fillR_angle, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "This function returns the 3*3 rotation matrix according to an angle, an axis and depending if it is active or passive rotation");
    m.def("fillR_euler", &fillR_euler, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "conv"_a="zxz", "copy"_a=true, "This function returns the 3*3 rotation matrix according to the three Euler angles, depending if it is active or passive rotation and the Euler convention (ex :zxz)");
    m.def("fillQS_angle", &fillQS_angle, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress' from an angle and an axis");
    m.def("fillQS_R", &fillQS_R, "R"_a, "active"_a=true, "copy"_a=true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress' from a rotation matrix");
    m.def("fillQE_angle", &fillQE_angle, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain' from an angle and an axis");
    m.def("fillQE_R", &fillQE_R, "R"_a, "active"_a=true, "copy"_a=true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain' from a rotation matrix");
    m.def("rotateL_angle", &rotateL_angle, "input"_a, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 stiffness matrix according to an angle and an axis");
    m.def("rotateL_R", &rotateL_R, "input"_a, "R"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 stiffness matrix according to a rotation matrix");
    m.def("rotate_l2g_L", &rotate_l2g_L, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 stiffness matrix from local to global frame");
    m.def("rotate_g2l_L", &rotate_g2l_L, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 stiffness matrix from global to local frame");
    m.def("rotateM_angle", &rotateM_angle, "input"_a, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 compliance matrix according to an angle and an axis");
    m.def("rotateM_R", &rotateM_R, "input"_a, "R"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 compliance matrix according to a rotation matrix");
    m.def("rotate_l2g_M", &rotate_l2g_M, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 compliance matrix from local to global frame");
    m.def("rotate_g2l_M", &rotate_g2l_M, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 compliance matrix from global to local frame");
    m.def("rotateA_angle", &rotateA_angle, "input"_a, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 strain concentration matrix according to an angle and an axis");
    m.def("rotateA_R", &rotateA_R, "input"_a, "R"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 strain concentration matrix according to a rotation matrix");
    m.def("rotate_l2g_A", &rotate_l2g_A, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 strain concentration matrix from local to global frame");
    m.def("rotate_g2l_A", &rotate_g2l_A, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 strain concentration matrix from global to local frame");
    m.def("rotateB_angle", &rotateB_angle, "input"_a, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 stress concentration matrix according to an angle and an axis");
    m.def("rotateB_R", &rotateB_R, "input"_a, "R"_a, "active"_a=true, "copy"_a=true, "Return the rotated 6*6 stress concentration matrix according to a rotation matrix");
    m.def("rotate_l2g_B", &rotate_l2g_B, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 stress concentration matrix from local to global frame");
    m.def("rotate_g2l_B", &rotate_g2l_B, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a=true, "Return the rotated 6*6 stress concentration matrix from global to local frame");
    m.def("rotate_strain_angle", &rotate_strain_angle, "input"_a, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=false, "Return the rotated strain matrix using voigt notations according to an angle and an axis");
    m.def("rotate_strain_R", &rotate_strain_R, "input"_a, "R"_a, "active"_a=true, "copy"_a=false, "Return the rotated strain matrix using voigt notations according to a rotation matrix");
    m.def("rotate_l2g_strain", &rotate_l2g_strain, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a=true, "Return the rotated strain matrix using voigt notations from local to global frame");
    m.def("rotate_g2l_strain", &rotate_g2l_strain, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a=true, "Return the rotated strain matrix using voigt notations from global to local frame");
    m.def("rotate_stress_angle", &rotate_stress_angle, "input"_a, "angle"_a, "axis"_a, "active"_a=true, "copy"_a=false, "Return the rotated stress matrix using voigt notations according to an angle and an axis");
    m.def("rotate_stress_R", &rotate_stress_R, "input"_a, "R"_a, "active"_a=true, "copy"_a=false, "Return the rotated stress matrix using voigt notations according to a rotation matrix");
    m.def("rotate_l2g_stress", &rotate_l2g_stress, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a=true, "Return the rotated stress matrix using voigt notations from local to global frame");
    m.def("rotate_g2l_stress", &rotate_g2l_stress, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a=true, "Return the rotated stress matrix using voigt notations from global to local frame");

    //Register the from-python converters for lagrange
    m.def("lagrange_exp", &lagrange_exp, "This function is used to determine an exponential Lagrange Multiplier (like contact in Abaqus)");
    m.def("dlagrange_exp", &dlagrange_exp, "This function is used to determine the first derivative of an exponential Lagrange Multiplier");
    m.def("lagrange_pow_0", &lagrange_pow_0, "This function is used to determine a power-law Lagrange Multiplier for problem such x >= 0");
    m.def("dlagrange_pow_0", &dlagrange_pow_0, "This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x >= 0");
    m.def("lagrange_pow_1", &lagrange_pow_1, "This function is used to determine a power-law Lagrange Multiplier for problem such x <= 1");
    m.def("dlagrange_pow_1", &dlagrange_pow_1, "This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x <= 1");
    m.def("d2lagrange_pow_1", &d2lagrange_pow_1, "This function is used to determine the SECOND derivative of a power-law Lagrange Multiplier for problem such x <= 1");

    //Register the from-python converters for stress
    m.def("stress_convert", &stress_convert, "sigma"_a, "F"_a, "converter_key"_a, "J"_a=0., "copy"_a=true, "Provides the first Piola Kirchoff stress tensor from the Cauchy stress tensor");

    //umat
    m.def("umat", &launch_umat, "umat_name"_a, "etot"_a, "Detot"_a, "F0"_a, "F1"_a, "sigma"_a, "DR"_a, "props"_a, "statev"_a, "time"_a, "dtime"_a, "Wm"_a, "temp"_a = pybind11::none(), "ndi"_a = 3, "n_threads"_a = 4);

    // Register the from-python converters for read and solver
    m.def("read_matprops", &read_matprops);
    m.def("read_path", &read_path);
    m.def("solver", &solver);        

    // Register the from-python converters for ODF functions
    m.def("get_densities_ODF", &get_densities_ODF);
    m.def("ODF_discretization", &ODF_discretization);

    // Register the from-python converters for identification
    m.def("identification", &identification);
    m.def("calc_cost", &calc_cost, "nfiles"_a, "data_num_name"_a); 

}
