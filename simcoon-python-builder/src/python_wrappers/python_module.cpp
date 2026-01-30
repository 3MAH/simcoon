
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

#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/lagrange.hpp>
#include <simcoon/python_wrappers/Libraries/Material/ODF.hpp>
#include <simcoon/python_wrappers/Libraries/Homogenization/eshelby.hpp>

#include <simcoon/docs/Libraries/Continuum_mechanics/doc_constitutive.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_contimech.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_criteria.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_damage.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_kinematics.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_recovery_props.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_transfer.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_stress.hpp>
#include <simcoon/docs/Libraries/Continuum_mechanics/doc_hyperelastic.hpp>

#include <simcoon/docs/Libraries/Homogenization/doc_eshelby.hpp>

using namespace std;
using namespace arma;
using namespace simpy;

using namespace pybind11::literals;

PYBIND11_MODULE(_core, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    // Register the from-python converters for constitutive.hpp
    m.def("Ireal", &Ireal, "copy"_a = true, simcoon_docs::Ireal);
    m.def("Ivol", &Ivol, "copy"_a = true, simcoon_docs::Ivol);
    m.def("Idev", &Idev, "copy"_a = true, simcoon_docs::Idev);
    m.def("Ireal2", &Ireal2, "copy"_a = true, simcoon_docs::Ireal2);
    m.def("Idev2", &Idev2, "copy"_a = true, simcoon_docs::Idev2);
    m.def("Ith", &Ith, "copy"_a = true, simcoon_docs::Ith);
    m.def("Ir2", &Ir2, "copy"_a = true, simcoon_docs::Ir2);
    m.def("Ir05", &Ir05, "copy"_a = true, simcoon_docs::Ir05);
    m.def("L_iso", &L_iso, "props"_a, "conv"_a, "copy"_a = true, simcoon_docs::L_iso);
    m.def("M_iso", &M_iso, "props"_a, "conv"_a, "copy"_a = true, simcoon_docs::M_iso);
    m.def("L_cubic", &L_cubic, "props"_a, "conv"_a, "copy"_a = true, simcoon_docs::L_cubic);
    m.def("M_cubic", &M_cubic, "props"_a, "conv"_a, "copy"_a = true, simcoon_docs::M_cubic);
    m.def("L_ortho", &L_ortho, "props"_a, "conv"_a, "copy"_a = true, simcoon_docs::L_ortho);
    m.def("M_ortho", &M_ortho, "props"_a, "conv"_a, "copy"_a = true, simcoon_docs::M_ortho);
    m.def("L_isotrans", &L_isotrans, "props"_a, "axis"_a, "copy"_a = true, simcoon_docs::L_isotrans);
    m.def("M_isotrans", &M_isotrans, "props"_a, "axis"_a, "copy"_a = true, simcoon_docs::M_isotrans);
    m.def("H_iso", &H_iso, "props"_a, "copy"_a = true, simcoon_docs::H_iso);

    // Register the from-python converters for contimech
    m.def("sph", &sph, "input"_a, "copy"_a = true, simcoon_docs::sph);
    m.def("dev", &dev, "input"_a, "copy"_a = true, simcoon_docs::dev);    
    m.def("tr", &tr, "input"_a, simcoon_docs::tr);
    m.def("Mises_stress", &Mises_stress, "input"_a, simcoon_docs::Mises_stress);
    m.def("eta_stress", &eta_stress, "input"_a, "copy"_a = true, simcoon_docs::eta_stress);
    m.def("eta_norm_stress", &eta_norm_stress, "input"_a, "copy"_a = true, simcoon_docs::eta_norm_stress);
    m.def("eta_norm_strain", &eta_norm_strain, "input"_a, "copy"_a = true, simcoon_docs::eta_norm_strain);
    m.def("norm_stress", &norm_stress, "input"_a, simcoon_docs::norm_stress);
    m.def("norm_strain", &norm_strain, "input"_a, simcoon_docs::norm_strain);
    m.def("Mises_strain", &Mises_strain, "input"_a, simcoon_docs::Mises_strain);
    m.def("eta_strain", &eta_strain, "input"_a, "copy"_a = true, simcoon_docs::eta_strain);
    m.def("J2_stress", &J2_stress, "input"_a, simcoon_docs::J2_stress);
    m.def("J2_strain", &J2_strain, "input"_a, simcoon_docs::J2_strain);
    m.def("J3_stress", &J3_stress, "input"_a, simcoon_docs::J3_stress);
    m.def("J3_strain", &J3_strain, "input"_a, simcoon_docs::J3_strain);
    m.def("Macaulay_p", &Macaulay_p, "value"_a, simcoon_docs::Macaulay_p);
    m.def("Macaulay_n", &Macaulay_n, "value"_a, simcoon_docs::Macaulay_n);
    m.def("sign", &simpy::sign, "value"_a, simcoon_docs::sign);
    m.def("normal_ellipsoid", &normal_ellipsoid, "u"_a, "v"_a, "a1"_a, "a2"_a, "a3"_a, "copy"_a = true, simcoon_docs::normal_ellipsoid);
    m.def("curvature_ellipsoid", &curvature_ellipsoid, "u"_a, "v"_a, "a1"_a, "a2"_a, "a3"_a, simcoon_docs::curvature_ellipsoid);
    m.def("sigma_int", &sigma_int, "input"_a, "u"_a, "v"_a, "a1"_a, "a2"_a, "a3"_a, "copy"_a = true, simcoon_docs::sigma_int);
    m.def("p_ikjl", &p_ikjl, "normal"_a, "copy"_a = true, simcoon_docs::p_ikjl);
    m.def("auto_sym_dyadic", &auto_sym_dyadic, "input"_a, "copy"_a = true, simcoon_docs::auto_sym_dyadic);
    m.def("sym_dyadic", &sym_dyadic, "a"_a, "b"_a, "copy"_a = true, simcoon_docs::sym_dyadic);
    m.def("auto_dyadic", &auto_dyadic, "input"_a, "copy"_a = true, simcoon_docs::auto_dyadic);
    m.def("dyadic_4vectors_sym", &dyadic_4vectors_sym, "n_a"_a, "n_b"_a, "conv"_a, "copy"_a = true, simcoon_docs::dyadic_4vectors_sym);
    m.def("dyadic", &dyadic, "a"_a, "b"_a, "copy"_a = true, simcoon_docs::dyadic);
    m.def("auto_sym_dyadic_operator", &auto_sym_dyadic_operator, "input"_a, "copy"_a = true, simcoon_docs::auto_sym_dyadic_operator);
    m.def("sym_dyadic_operator", &sym_dyadic_operator, "a"_a, "b"_a, "copy"_a = true, simcoon_docs::sym_dyadic_operator);
    m.def("dyadic", &dyadic, "a"_a, "b"_a, "copy"_a = true, simcoon_docs::dyadic);

    // Register the from-python converters for criteria
    m.def("Drucker_stress", &Drucker_stress, "input"_a, "props"_a, simcoon_docs::Drucker_stress);
    m.def("dDrucker_stress", &dDrucker_stress, "input"_a, "props"_a, "copy"_a = true, simcoon_docs::dDrucker_stress);
    m.def("Tresca_stress", &Tresca_stress, "input"_a, simcoon_docs::Tresca_stress);
    m.def("dTresca_stress", &dTresca_stress, "input"_a, "copy"_a = true, simcoon_docs::dTresca_stress);
    m.def("P_Ani", &P_Ani, "props"_a, "copy"_a = true, simcoon_docs::P_Ani);
    m.def("P_Hill", &P_Hill, "props"_a, "copy"_a = true, simcoon_docs::P_Hill);
    m.def("P_DFA", &P_DFA, "props"_a, "copy"_a = true, simcoon_docs::P_DFA);
    m.def("Hill_stress", &Hill_stress, "input"_a, "props"_a, simcoon_docs::Hill_stress);
    m.def("dHill_stress", &dHill_stress, "input"_a, "props"_a, "copy"_a = true, simcoon_docs::dHill_stress);
    m.def("Ani_stress", &Ani_stress, "input"_a, "props"_a, simcoon_docs::Ani_stress);
    m.def("dAni_stress", &dAni_stress, "input"_a, "props"_a, "copy"_a = true, simcoon_docs::dAni_stress);
    m.def("DFA_stress", &DFA_stress, "input"_a, "props"_a, simcoon_docs::DFA_stress);
    m.def("dDFA_stress", &dDFA_stress, "input"_a, "props"_a, "copy"_a = true, simcoon_docs::dDFA_stress);
    m.def("Eq_stress", &Eq_stress, "input"_a, "criteria"_a, "props"_a, simcoon_docs::Eq_stress);
    m.def("dEq_stress", &dEq_stress, "input"_a, "criteria"_a, "props"_a, "copy"_a = true, simcoon_docs::dEq_stress);

    // register the damage library
    m.def("damage_weibull", &damage_weibull, "stress"_a, "damage"_a, "alpha"_a, "beta"_a, "DTime"_a, "criterion"_a = "vonmises", simcoon_docs::damage_weibull);
    m.def("damage_kachanov", &damage_kachanov, "stress"_a, "strain"_a, "damage"_a, "A0"_a, "r"_a, "criterion"_a, simcoon_docs::damage_kachanov);
    m.def("damage_miner", &damage_miner, "S_max"_a, "S_mean"_a, "S_ult"_a, "b"_a, "B0"_a, "beta"_a, "Sl_0"_a = 0., simcoon_docs::damage_miner);
    m.def("damage_manson", &damage_manson, "S_amp"_a, "C2"_a, "gamma2"_a, simcoon_docs::damage_manson);

    // Register the from-python converters for recovery_props
    m.def("check_symetries", &check_symetries, "input"_a, "tol"_a, simcoon_docs::check_symetries);
    m.def("L_iso_props", &L_iso_props, "input"_a, simcoon_docs::L_iso_props);
    m.def("M_iso_props", &M_iso_props, "input"_a, simcoon_docs::M_iso_props);
    m.def("L_isotrans_props", &L_isotrans_props, "input"_a, "axis"_a, simcoon_docs::L_isotrans_props);
    m.def("M_isotrans_props", &M_isotrans_props, "input"_a, "axis"_a, simcoon_docs::M_isotrans_props);
    m.def("L_cubic_props", &L_cubic_props, "input"_a, simcoon_docs::L_cubic_props);
    m.def("M_cubic_props", &M_cubic_props, "input"_a, simcoon_docs::M_cubic_props);
    m.def("L_ortho_props", &L_ortho_props, "input"_a, simcoon_docs::L_ortho_props);
    m.def("M_ortho_props", &M_ortho_props, "input"_a, simcoon_docs::M_ortho_props);
    m.def("M_aniso_props", &M_aniso_props, "input"_a, simcoon_docs::M_aniso_props);

    // Register the L_eff for composites
    m.def("L_eff", &L_eff, "umat_name"_a, "props"_a, "nstatev"_a, "psi_rve"_a = 0., "theta_rve"_a = 0., "phi_rve"_a = 0., "Return the elastic stiffness tensor of a composite material");

    // Register the from-python converters for kinematics
    m.def("ER_to_F", &ER_to_F, "E"_a, "R"_a, "copy"_a = true, simcoon_docs::ER_to_F);
    m.def("eR_to_F", &eR_to_F, "e"_a, "R"_a, "copy"_a = true, simcoon_docs::eR_to_F);
    m.def("G_UdX", &G_UdX, "F"_a, "copy"_a = true, simcoon_docs::G_UdX);
    m.def("G_Udx", &G_Udx, "F"_a, "copy"_a = true, simcoon_docs::G_Udx);
    m.def("R_Cauchy_Green", &R_Cauchy_Green, "F"_a, "copy"_a = true, simcoon_docs::R_Cauchy_Green);
    m.def("L_Cauchy_Green", &L_Cauchy_Green, "F"_a, "copy"_a = true, simcoon_docs::L_Cauchy_Green);
    m.def("RU_decomposition", &RU_decomposition, "F"_a, "copy"_a = true, simcoon_docs::RU_decomposition);
    m.def("VR_decomposition", &VR_decomposition, "F"_a, "copy"_a = true, simcoon_docs::VR_decomposition);
    m.def("Inv_X", &Inv_X, "input"_a, "copy"_a = true, simcoon_docs::Inv_X);
    m.def("Cauchy", &Cauchy, "F"_a, "copy"_a = true, simcoon_docs::Cauchy);
    m.def("Green_Lagrange", &Green_Lagrange, "F"_a, "copy"_a = true, simcoon_docs::Green_Lagrange);
    m.def("Euler_Almansi", &Euler_Almansi, "F"_a, "copy"_a = true, simcoon_docs::Euler_Almansi);
    m.def("Log_strain", &Log_strain, "F"_a, "voigt_form"_a = false, "copy"_a = true, simcoon_docs::Log_strain);
    m.def("finite_L", &finite_L, "F0"_a, "F1"_a, "DTime"_a, "copy"_a = true, simcoon_docs::finite_L);
    m.def("finite_D", &finite_D, "F0"_a, "F1"_a, "DTime"_a, "copy"_a = true, simcoon_docs::finite_D);
    m.def("finite_W", &finite_W, "F0"_a, "F1"_a, "DTime"_a, "copy"_a = true, simcoon_docs::finite_W);
    m.def("finite_Omega", &finite_Omega, "F0"_a, "F1"_a, "DTime"_a, "copy"_a = true, simcoon_docs::finite_Omega);
    m.def("finite_DQ", &finite_DQ, "Omega0"_a, "Omega1"_a, "DTime"_a, "copy"_a = true, simcoon_docs::finite_DQ);

    // register the objective rates library
    m.def("logarithmic", &logarithmic, "F0"_a, "F1"_a, "DTime"_a, "copy"_a = true, "This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment");
    m.def("logarithmic_R", &logarithmic_R, "F0"_a, "F1"_a, "DTime"_a, "copy"_a = true, "This function computes the logarithmic strain velocity and the Green-Naghdi spin, along with the correct rotation increment");
    m.def("Delta_log_strain", &Delta_log_strain, "D"_a, "Omega"_a, "DTime"_a, "copy"_a = true, "This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor");
    m.def("objective_rate", &objective_rate, "corate_name"_a, "F0"_a, "F1"_a, "dtime"_a, "return_de"_a = false, "n_threads"_a = 4, "This function computes the strain velocity and the spin, along with the correct rotation increment for the specified objective erivative");
    m.def("Lt_convert", &Lt_convert, "Lt"_a, "F"_a, "stress"_a, "converter_key"_a);

    // register the hyperelastic library
    m.def("isochoric_invariants", &isochoric_invariants, "input"_a, "J"_a = 0., "copy"_a = true, simcoon_docs::isochoric_invariants);
    m.def("isochoric_pstretch", &isochoric_pstretch, "input"_a, "input_tensor"_a = "V", "J"_a = 0., "copy"_a = true, simcoon_docs::isochoric_pstretch);
    m.def("tau_iso_hyper_invariants", &tau_iso_hyper_invariants, "dWdI_1_bar"_a, "dWdI_2_bar"_a, "input"_a, "J"_a = 0., "copy"_a = true, simcoon_docs::tau_iso_hyper_invariants);
    m.def("sigma_iso_hyper_invariants", &sigma_iso_hyper_invariants, "dWdI_1_bar"_a, "dWdI_2_bar"_a, "input"_a, "J"_a = 0., "copy"_a = true, simcoon_docs::sigma_iso_hyper_invariants);
    m.def("tau_vol_hyper", &tau_vol_hyper, "dUdJ"_a, "input"_a, "J"_a = 0., "copy"_a = true, simcoon_docs::tau_vol_hyper);
    m.def("sigma_vol_hyper", &sigma_vol_hyper, "dUdJ"_a, "input"_a, "J"_a = 0., "copy"_a = true, simcoon_docs::sigma_vol_hyper);

    // register the transfer library
    m.def("v2t_strain", &v2t_strain, "input"_a, "copy"_a = true, simcoon_docs::v2t_strain);
    m.def("t2v_strain", &t2v_strain, "input"_a, "copy"_a = true, simcoon_docs::t2v_strain);
    m.def("v2t_stress", &v2t_stress, "input"_a, "copy"_a = true, simcoon_docs::v2t_stress);
    m.def("t2v_stress", &t2v_stress, "input"_a, "copy"_a = true, simcoon_docs::t2v_stress); 

    // Register the from-python converters for eshelby
    m.def("Eshelby_sphere", &Eshelby_sphere, "nu"_a, "copy"_a = true, simcoon_docs::Eshelby_sphere);
    m.def("Eshelby_cylinder", &Eshelby_cylinder, "nu"_a, "copy"_a = true, simcoon_docs::Eshelby_cylinder);
    m.def("Eshelby_prolate", &Eshelby_prolate, "nu"_a, "aspect_ratio"_a, "copy"_a = true, simcoon_docs::Eshelby_prolate);
    m.def("Eshelby_oblate", &Eshelby_oblate, "nu"_a, "aspect_ratio"_a, "copy"_a = true, simcoon_docs::Eshelby_oblate);
    m.def("Eshelby_penny", &Eshelby_penny, "nu"_a, "copy"_a = true, simcoon_docs::Eshelby_penny);
    m.def("Eshelby", &Eshelby, "L"_a, "a1"_a = 1., "a2"_a = 1., "a3"_a = 1., "mp"_a = 50, "np"_a = 50, "copy"_a = true, simcoon_docs::Eshelby);
    m.def("T_II", &T_II, "L"_a, "a1"_a = 1., "a2"_a = 1., "a3"_a = 1., "mp"_a = 50, "np"_a = 50, "copy"_a = true, simcoon_docs::T_II);

    // Register the rotation library
    m.def("rotate_vec_R", &rotate_vec_R, "input"_a, "R"_a, "copy"_a = true, "This function returns a rotated vector (3) according to a rotation matrix");
    m.def("rotate_vec_angle", &rotate_vec_angle, "input"_a, "angle"_a, "axis"_a, "copy"_a = true, "This function returns a rotated vector (3) according to an angle and an axis");
    m.def("rotate_mat_R", &rotate_mat_R, "input"_a, "R"_a, "copy"_a = true, "This function returns a rotated matrix (3x3) according to a rotation matrix");
    m.def("rotate_mat_angle", &rotate_mat_angle, "input"_a, "angle"_a, "axis"_a, "copy"_a = true, "This function returns a rotated matrix (3x3) according to an angle and an axis");
    m.def("fillR_angle", &fillR_angle, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "This function returns the 3*3 rotation matrix according to an angle, an axis and depending if it is active or passive rotation");
    m.def("fillR_euler", &fillR_euler, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "conv"_a = "zxz", "copy"_a = true, "This function returns the 3*3 rotation matrix according to the three Euler angles, depending if it is active or passive rotation and the Euler convention (ex :zxz)");
    m.def("fillQS_angle", &fillQS_angle, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress' from an angle and an axis");
    m.def("fillQS_R", &fillQS_R, "R"_a, "active"_a = true, "copy"_a = true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress' from a rotation matrix");
    m.def("fillQE_angle", &fillQE_angle, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain' from an angle and an axis");
    m.def("fillQE_R", &fillQE_R, "R"_a, "active"_a = true, "copy"_a = true, "This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain' from a rotation matrix");
    m.def("rotateL_angle", &rotateL_angle, "input"_a, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 stiffness matrix according to an angle and an axis");
    m.def("rotateL_R", &rotateL_R, "input"_a, "R"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 stiffness matrix according to a rotation matrix");
    m.def("rotate_l2g_L", &rotate_l2g_L, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 stiffness matrix from local to global frame");
    m.def("rotate_g2l_L", &rotate_g2l_L, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 stiffness matrix from global to local frame");
    m.def("rotateM_angle", &rotateM_angle, "input"_a, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 compliance matrix according to an angle and an axis");
    m.def("rotateM_R", &rotateM_R, "input"_a, "R"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 compliance matrix according to a rotation matrix");
    m.def("rotate_l2g_M", &rotate_l2g_M, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 compliance matrix from local to global frame");
    m.def("rotate_g2l_M", &rotate_g2l_M, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 compliance matrix from global to local frame");
    m.def("rotateA_angle", &rotateA_angle, "input"_a, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 strain concentration matrix according to an angle and an axis");
    m.def("rotateA_R", &rotateA_R, "input"_a, "R"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 strain concentration matrix according to a rotation matrix");
    m.def("rotate_l2g_A", &rotate_l2g_A, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 strain concentration matrix from local to global frame");
    m.def("rotate_g2l_A", &rotate_g2l_A, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 strain concentration matrix from global to local frame");
    m.def("rotateB_angle", &rotateB_angle, "input"_a, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 stress concentration matrix according to an angle and an axis");
    m.def("rotateB_R", &rotateB_R, "input"_a, "R"_a, "active"_a = true, "copy"_a = true, "Return the rotated 6*6 stress concentration matrix according to a rotation matrix");
    m.def("rotate_l2g_B", &rotate_l2g_B, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 stress concentration matrix from local to global frame");
    m.def("rotate_g2l_B", &rotate_g2l_B, "input"_a, "psi"_a, "theta"_a, "phi"_a, "active"_a = true, "Return the rotated 6*6 stress concentration matrix from global to local frame");
    m.def("rotate_strain_angle", &rotate_strain_angle, "input"_a, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = false, "Return the rotated strain matrix using voigt notations according to an angle and an axis");
    m.def("rotate_strain_R", &rotate_strain_R, "input"_a, "R"_a, "active"_a = true, "copy"_a = false, "Return the rotated strain matrix using voigt notations according to a rotation matrix");
    m.def("rotate_l2g_strain", &rotate_l2g_strain, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a = true, "Return the rotated strain matrix using voigt notations from local to global frame");
    m.def("rotate_g2l_strain", &rotate_g2l_strain, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a = true, "Return the rotated strain matrix using voigt notations from global to local frame");
    m.def("rotate_stress_angle", &rotate_stress_angle, "input"_a, "angle"_a, "axis"_a, "active"_a = true, "copy"_a = false, "Return the rotated stress matrix using voigt notations according to an angle and an axis");
    m.def("rotate_stress_R", &rotate_stress_R, "input"_a, "R"_a, "active"_a = true, "copy"_a = false, "Return the rotated stress matrix using voigt notations according to a rotation matrix");
    m.def("rotate_l2g_stress", &rotate_l2g_stress, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a = true, "Return the rotated stress matrix using voigt notations from local to global frame");
    m.def("rotate_g2l_stress", &rotate_g2l_stress, "input"_a, "psi"_a, "theta"_a, "phi"_a, "copy"_a = true, "Return the rotated stress matrix using voigt notations from global to local frame");

    // Register the from-python converters for lagrange
    m.def("lagrange_exp", &lagrange_exp, "This function is used to determine an exponential Lagrange Multiplier (like contact in Abaqus)");
    m.def("dlagrange_exp", &dlagrange_exp, "This function is used to determine the first derivative of an exponential Lagrange Multiplier");
    m.def("lagrange_pow_0", &lagrange_pow_0, "This function is used to determine a power-law Lagrange Multiplier for problem such x >= 0");
    m.def("dlagrange_pow_0", &dlagrange_pow_0, "This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x >= 0");
    m.def("lagrange_pow_1", &lagrange_pow_1, "This function is used to determine a power-law Lagrange Multiplier for problem such x <= 1");
    m.def("dlagrange_pow_1", &dlagrange_pow_1, "This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x <= 1");
    m.def("d2lagrange_pow_1", &d2lagrange_pow_1, "This function is used to determine the SECOND derivative of a power-law Lagrange Multiplier for problem such x <= 1");

    // Register the from-python converters for stress
    m.def("stress_convert", &stress_convert, "sigma"_a, "F"_a, "converter_key"_a, "J"_a = 0., "copy"_a = true, simcoon_docs::stress_convert);

    // umat
    m.def("umat", &launch_umat, "umat_name"_a, "etot"_a, "Detot"_a, "F0"_a, "F1"_a, "sigma"_a, "DR"_a, "props"_a, "statev"_a, "time"_a, "dtime"_a, "Wm"_a, "temp"_a = pybind11::none(), "ndi"_a = 3, "n_threads"_a = 4);

    // ODF functions
    m.def("get_densities_ODF", &get_densities_ODF);
    m.def("ODF_discretization", &ODF_discretization);
}
