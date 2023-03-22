
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <pybind11/pybind11.h>

#include <CGAL/Simple_cartesian.h>

#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/numpy_cgal.hpp>

#include <simcoon/Simulation/Phase/state_variables.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/constitutive.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/contimech.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/criteria.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/recovery_props.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Leff.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/objective_rates.hpp>
//#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/RunUmat.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Umat_fedoo.hpp>

#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/lagrange.hpp>
#include <simcoon/python_wrappers/Libraries/Material/ODF.hpp>
#include <simcoon/python_wrappers/Libraries/Homogenization/eshelby.hpp>

#include <simcoon/python_wrappers/Libraries/Solver/read.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/step_meca.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/step_thermomeca.hpp>

#include <simcoon/python_wrappers/Libraries/Identification/identification.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/constants.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/parameters.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/optimize.hpp>

#include <simcoon/python_wrappers/Libraries/Unit_cell/unit_cell.hpp>

#include <simcoon/python_wrappers/Libraries/Phase/state_variables.hpp>
#include <simcoon/python_wrappers/Libraries/Phase/state_variables_M.hpp>
#include <simcoon/python_wrappers/Libraries/Phase/state_variables_T.hpp>

//#include <simcoon/python_wrappers/Libraries/Abaqus/write.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace simpy;

using namespace pybind11::literals;

PYBIND11_MODULE(simmitpybind, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    // Register the from-python converters for constitutive.hpp
    m.def("Ireal", &Ireal, "copy"_a=True, "Returns the fourth order identity tensor written in Voigt notation Ireal");
    m.def("Ivol", &Ivol, "copy"_a=True, "Returns the volumic of the identity tensor Ireal written in Voigt notation");
    m.def("Idev", &Idev, "copy"_a=True, "Returns the deviatoric of the identity tensor Ireal written in Voigt notation");        
    m.def("Ireal2", &Ireal2, "copy"_a=True, "Returns the fourth order identity tensor Iˆ written in Voigt notation");
    m.def("Idev2", &Idev2, "copy"_a=True, "Returns the deviatoric of the identity tensor Iˆ written in Voigt notation");    
    m.def("Ith", &Ith, "copy"_a=True, "Returns the expansion vector");        
    m.def("Ir2", &Ir2, "copy"_a=True, "Returns the stress 2 strain operator");
    m.def("Ir05", &Ir05, "copy"_a=True, "Returns the strain 2 stress operator");        
    m.def("L_iso", &L_iso, "copy"_a=True, "Provides the elastic stiffness tensor for an isotropic material");
    m.def("M_iso", &M_iso, "copy"_a=True, "Provides the elastic compliance tensor for an isotropic material");
    m.def("L_cubic", &L_cubic, "copy"_a=True, "Returns the elastic stiffness tensor for a cubic material");
    m.def("M_cubic", &M_cubic, "copy"_a=True, "Returns the elastic compliance tensor for an isotropic material");
    m.def("L_ortho", &L_ortho, "copy"_a=True, "Returns the elastic stiffness tensor for an orthotropic material");
    m.def("M_ortho", &M_ortho, "copy"_a=True, "Returns the elastic compliance tensor for an orthotropic material");
    m.def("L_isotrans", &L_isotrans, "copy"_a=True, "Returns the elastic stiffness tensor for an isotropic transverse material");
    m.def("M_isotrans", &M_isotrans, "copy"_a=True, "Returns the elastic compliance tensor for an isotropic transverse material");
    m.def("H_iso", &H_iso, "copy"_a=True, "Provides the viscous tensor H an isotropic material");

    // Register the from-python converters for contimech
    m.def(("tr", &tr, "This function returns the trace of a tensor");
    m.def("dev", &dev, "copy"_a=True, "This function returns the deviatoric part of a tensor");
    m.def("Mises_stress", &Mises_stress, "This function determines the Mises equivalent of a stress tensor");
    m.def("eta_stress", &eta_stress, "copy"_a=True, "This function determines the strain flow (direction) from a stress tensor");
    m.def("Mises_strain", &Mises_strain, "This function determines the Mises equivalent of a strain tensor");
    m.def("eta_strain", &eta_strain, "copy"_a=True, "This function determines the strain flow (direction) from a strain tensor");
    m.def("J2_stress", &J2_stress, "Returns the second invariant of the deviatoric part of a second order stress tensor");
    m.def("J2_strain", &J2_strain, "Returns the second invariant of the deviatoric part of a second order strain tensor");
    m.def("J3_stress", &J3_stress, "Returns the third invariant of the deviatoric part of a second order stress tensor");
    m.def("J3_strain", &J3_strain, "Returns the third invariant of the deviatoric part of a second order stress tensor");
    m.def("Macaulay_p", &Macaulay_p, "This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)");
    m.def("Macaulay_n", &Macaulay_n, "This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)");
    m.def("sign", &simpy::sign, "This function returns the sign of a double");
    m.def("normal_ellipsoid", "copy"_a=True, &normal_ellipsoid, "Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3.");
    m.def("sigma_int", &sigma_int, "copy"_a=True, "Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3.");
    m.def("p_ikjl", &p_ikjl, "copy"_a=True, "This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)");

    // Register the from-python converters for criteria
    m.def("Prager_stress", &Prager_stress, "This function returns the Prager equivalent stress");
    m.def("dPrager_stress", &dPrager_stress, "copy"_a=True, "This function returns the derivative of the Prager equivalent stress");
    m.def("Tresca_stress", &Tresca_stress, "This function returns the Tresca equivalent stress");
    m.def("dTresca_stress", &dTresca_stress, "copy"_a=True, "This function returns the derivative of the Tresca equivalent stress");
    m.def("P_ani", &P_ani, "copy"_a=True, "Returns an anisotropic configurational tensor P in the Voigt format (6x6 numpy array), given its vector representation");
    m.def("P_Hill", &P_Hill, "copy"_a=True, "Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion in the Voigt format (6x6 numpy array), given its vector representation");
    m.def("Hill_stress", &Hill_stress, "This function returns the Hill equivalent stress");
    m.def("dHill_stress", &dHill_stress, "copy"_a=True, "This function returns the derivative of the Hill equivalent stress");
    m.def("Ani_stress", &Ani_stress, "This function returns the Ani equivalent stress");
    m.def("dAni_stress", &dAni_stress, "copy"_a=True, "This function returns the derivative of the Ani equivalent stress");
    m.def("Eq_stress", &Eq_stress, "This function computes the selected equivalent stress function");
    m.def("dEq_stress", &dEq_stress, "copy"_a=True, "This function computes the deriavtive of the selected equivalent stress function");

    // Register the from-python converters for recovery_props
    m.def("check_symetries", &check_symetries, "Check the material symetries and the type of elastic response for a given stiffness tensor");
    m.def("L_iso_props", &L_iso_props, "Return a list of elastic properties for the isotropic case (E,nu) from a stiffness tensor");
    m.def("M_iso_props", &M_iso_props, "Return a list of elastic properties for the isotropic case (E,nu) from a compliance tensor");
    m.def("L_isotrans_props", &L_isotrans_props, "Return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a stiffness tensor");
    m.def("M_isotrans_props", &M_isotrans_props, "Return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a compliance tensor");
    m.def("L_cubic_props", &L_cubic_props, "Return a list of elastic properties for the cubic case (E,nu,G) from a stiffness tensor");
    m.def("M_cubic_props", &M_cubic_props, "Return a list of elastic properties for the cubic case (E,nu,G) from a compliance tensor");
    m.def("L_ortho_props", &L_ortho_props, "Return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a stiffness tensor");
    m.def("M_ortho_props", &M_ortho_props, "Return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a compliance tensor");
    m.def("M_aniso_props", &M_aniso_props, "//Return a list of elastic properties for the anisotropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,deviations) from a compliance tensor");

}


BOOST_PYTHON_MODULE(simmit) {

    Py_Initialize();
    bn::initialize();
    
/*    // Register the from-python converters for constitutive
    bp::def("Ireal", Ireal);
    bp::def("Ivol", Ivol);
    bp::def("Idev", Idev);
    bp::def("Ireal2", Ireal2);
    bp::def("Idev2", Idev2);
    bp::def("Ith", Ith);
    bp::def("Ir2", Ir2);
    bp::def("Ir05", Ir05);
    bp::def("L_iso", L_iso);
    bp::def("M_iso", M_iso);
    bp::def("L_cubic", L_cubic);
    bp::def("M_cubic", M_cubic);
    bp::def("L_ortho", L_ortho);
    bp::def("M_ortho", M_ortho);
    bp::def("L_isotrans", L_isotrans);
    bp::def("M_isotrans", M_isotrans);
    bp::def("H_iso", H_iso);
*/
    
/*    // Register the from-python converters for contimech
    bp::def("tr", tr);
    bp::def("dev", dev);
    bp::def("Mises_stress", Mises_stress);
    bp::def("eta_stress", eta_stress);
    bp::def("Mises_strain", Mises_strain);
    bp::def("eta_strain", eta_strain);
    bp::def("v2t_strain", v2t_strain);
    bp::def("t2v_strain", t2v_strain);
    bp::def("v2t_stress", v2t_stress);
    bp::def("t2v_stress", t2v_stress);
    bp::def("J2_stress", J2_stress);
    bp::def("J2_strain", J2_strain);
    bp::def("J3_stress", J3_stress);
    bp::def("J3_strain", J3_strain);
    bp::def("Macaulay_p", Macaulay_p);
    bp::def("Macaulay_n", Macaulay_n);
    bp::def("sign", simpy::sign);
    bp::def("normal_ellipsoid", normal_ellipsoid);
    bp::def("sigma_int", sigma_int);
    bp::def("p_ikjl", p_ikjl);
*/
    
/*    // Register the from-python converters for criteria
    bp::def("Prager_stress", Prager_stress);
    bp::def("dPrager_stress", dPrager_stress);
    bp::def("Tresca_stress", Tresca_stress);
    bp::def("dTresca_stress", dTresca_stress);
    bp::def("Eq_stress", Eq_stress);
    bp::def("dEq_stress", dEq_stress);
*/
    
    // Register the from-python converters for recovery_props
/*    bp::def("check_symetries", check_symetries);
    bp::def("L_iso_props", L_iso_props);
    bp::def("M_iso_props", M_iso_props);
    bp::def("L_isotrans_props", L_isotrans_props);
    bp::def("M_isotrans_props", M_isotrans_props);
    bp::def("L_cubic_props", L_cubic_props);
    bp::def("M_cubic_props", M_cubic_props);
    bp::def("L_ortho_props", L_ortho_props);
    bp::def("M_ortho_props", M_ortho_props);
    bp::def("M_aniso_props", M_aniso_props);
    */
    
    // Register the L_eff for composites
    bp::def("L_eff", L_eff);
    
    //register the kinematics library
    bp::def("G_UdX", G_UdX);
    bp::def("G_Udx", G_Udx);
    bp::def("R_Cauchy_Green", R_Cauchy_Green);
    bp::def("L_Cauchy_Green", L_Cauchy_Green);
    bp::def("Inv_X", Inv_X);
    bp::def("Cauchy", Cauchy);
    bp::def("Green_Lagrange", Green_Lagrange);
    bp::def("Euler_Almansi", Euler_Almansi);
    bp::def("Log_strain", Log_strain);
    bp::def("finite_L", finite_L);
    bp::def("finite_D", finite_D);
    bp::def("finite_W", finite_W);
    bp::def("finite_Omega", finite_Omega);
    bp::def("finite_DQ", finite_DQ);
    
    //register the objective rates library
    bp::def("logarithmic", logarithmic);
    bp::def("Delta_log_strain", Delta_log_strain);
    
    // Register the from-python converters for eshelby
    bp::def("Eshelby_sphere", Eshelby_sphere);
    bp::def("Eshelby_cylinder", Eshelby_cylinder);
    bp::def("Eshelby_prolate", Eshelby_prolate);
    bp::def("Eshelby_oblate", Eshelby_oblate);
    bp::def("Eshelby", Eshelby);
    bp::def("T_II", T_II);
    
    // Register the from-python converters for read and solver
    bp::def("read_matprops", read_matprops);
    bp::def("read_path", read_path);
    bp::def("solver", solver);
	
    //Wrapper fedoo
 /*   bp::class_<Umat_fedoo>("Umat_fedoo", bp::init < std::string, bn::ndarray, int, int, int> ())
        .def("compute_Detot", &Umat_fedoo::compute_Detot)
        .def("Run", &Umat_fedoo::Run)
        .def("Initialize", &Umat_fedoo::Initialize)
        .def("to_start", &Umat_fedoo::to_start)
        .def("set_start", &Umat_fedoo::set_start)
        .def("set_T", &Umat_fedoo::set_T)
        .def_readwrite("corate", &Umat_fedoo::corate)
        .def_readonly("Time", &Umat_fedoo::Time)
        .def_readonly("DTime", &Umat_fedoo::DTime)
        .def_readonly("nb_points", &Umat_fedoo::nb_points)
        .def_readonly("nlgeom", &Umat_fedoo::nlgeom)
        .add_property("T", &Umat_fedoo::Get_T)
        .add_property("props", &Umat_fedoo::Get_props)
        .add_property("Kirchhoff", &Umat_fedoo::Get_Kirchhoff)
        .add_property("Cauchy", &Umat_fedoo::Get_Cauchy)
        .add_property("PKII", &Umat_fedoo::Get_PKII)
        .add_property("etot", &Umat_fedoo::Get_etot)
        .add_property("Detot", &Umat_fedoo::Get_Detot)
        .add_property("statev", &Umat_fedoo::Get_statev)
        .add_property("L", &Umat_fedoo::Get_L)
        .add_property("Lt", &Umat_fedoo::Get_Lt)
        .add_property("DR", &Umat_fedoo::Get_DR)
        .add_property("Wm", &Umat_fedoo::Get_Wm)
        .add_property("F0", &Umat_fedoo::Get_F0)
        .add_property("F1", &Umat_fedoo::Get_F1)
        ;
        
    bp::class_<state_variables_py>("state_variables", bp::init <>())
        .def(bp::init <const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const double&, const double&> ())
        .def_readwrite("T", &state_variables_py::T)
        .def_readwrite("DT", &state_variables_py::DT)
        .def_readwrite("nstatev", &state_variables_py::nstatev)
        .add_property("F0", &state_variables_py::Get_F0, &state_variables_py::Set_F0)
        .add_property("F1", &state_variables_py::Get_F1, &state_variables_py::Set_F1)
        .add_property("etot", &state_variables_py::Get_etot)
        .add_property("Detot", &state_variables_py::Get_Detot)
        .add_property("Etot", &state_variables_py::Get_Etot)
        .add_property("DEtot", &state_variables_py::Get_DEtot)
        .add_property("statev", &state_variables_py::Get_statev)
        .add_property("R", &state_variables_py::Get_R)
        .add_property("DR", &state_variables_py::Get_DR)
        .def("to_start", &state_variables_py::to_start)
        .def("set_start", &state_variables_py::set_start)
        .def("rotate_l2g", &state_variables_py::rotate_l2g_py)
        .def("rotate_g2l", &state_variables_py::rotate_g2l_py)
        ;
    
    bp::class_<state_variables_M_py>("state_variables_M", bp::init <>())
        .def(bp::init <const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const double&, const double&> ())
        .def_readwrite("T", &state_variables_M_py::T)
        .def_readwrite("DT", &state_variables_M_py::DT)
        .def_readwrite("nstatev", &state_variables_M_py::nstatev)
        .add_property("F0", &state_variables_M_py::Get_F0, &state_variables_M_py::Set_F0)
        .add_property("F1", &state_variables_M_py::Get_F1, &state_variables_M_py::Set_F1)
        .add_property("etot", &state_variables_M_py::Get_etot)
        .add_property("Detot", &state_variables_M_py::Get_Detot)
        .add_property("Etot", &state_variables_M_py::Get_Etot)
        .add_property("DEtot", &state_variables_M_py::Get_DEtot)
        .add_property("statev", &state_variables_M_py::Get_statev)
        .add_property("R", &state_variables_M_py::Get_R)
        .add_property("DR", &state_variables_M_py::Get_DR)
        .add_property("Wm", &state_variables_M_py::Get_Wm)
        .add_property("L", &state_variables_M_py::Get_L)
        .add_property("Lt", &state_variables_M_py::Get_Lt)
        .def("to_start", &state_variables_M_py::to_start)
        .def("set_start", &state_variables_M_py::set_start)
        .def("rotate_l2g", &state_variables_M_py::rotate_l2g_py)
        .def("rotate_g2l", &state_variables_M_py::rotate_l2g_py)
        ;
    
    bp::class_<state_variables_T_py>("state_variables_T", bp::init <>())
        .def(bp::init <const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const double&, const double&> ())
        .def_readwrite("T", &state_variables_T_py::T)
        .def_readwrite("DT", &state_variables_T_py::DT)
        .def_readwrite("nstatev", &state_variables_T_py::nstatev)
        .add_property("F0", &state_variables_T_py::Get_F0, &state_variables_T_py::Set_F0)
        .add_property("F1", &state_variables_T_py::Get_F1, &state_variables_T_py::Set_F1)
        .add_property("etot", &state_variables_T_py::Get_etot)
        .add_property("Detot", &state_variables_T_py::Get_Detot)
        .add_property("Etot", &state_variables_T_py::Get_Etot)
        .add_property("DEtot", &state_variables_T_py::Get_DEtot)
        .add_property("statev", &state_variables_T_py::Get_statev)
        .add_property("R", &state_variables_T_py::Get_R)
        .add_property("DR", &state_variables_T_py::Get_DR)
        .add_property("Wm", &state_variables_T_py::Get_Wm)
        .add_property("Wt", &state_variables_T_py::Get_Wt)
        .add_property("dSdE", &state_variables_T_py::Get_dSdE)
        .add_property("dSdEt", &state_variables_T_py::Get_dSdEt)
        .add_property("dSdT", &state_variables_T_py::Get_dSdT)
        .add_property("drdE", &state_variables_T_py::Get_drdE)
        .add_property("drdT", &state_variables_T_py::Get_drdT)
        .def("to_start", &state_variables_T_py::to_start)
        .def("set_start", &state_variables_T_py::set_start)
        .def("rotate_l2g", &state_variables_T_py::rotate_l2g_py)
        .def("rotate_g2l", &state_variables_T_py::rotate_l2g_py)
        ;
    
    bp::class_<step_meca_py>("step_meca", bp::init <>())
        .def(bp::init <const int &, const double &, const double &, const double &, const int &, const unsigned int &, const bn::ndarray&, const bn::ndarray&, const double&, const int&, const bn::ndarray&, const bn::ndarray&> ())
        .def_readwrite("number", &step_meca_py::number)
        .def_readwrite("Dn_init",&step_meca_py::Dn_init)
        .def_readwrite("Dn_mini", &step_meca_py::Dn_mini)
        .def_readwrite("Dn_inc", &step_meca_py::Dn_inc)
        .def_readwrite("mode", &step_meca_py::mode)
        .def_readwrite("control_type", &step_meca_py::control_type)
        .add_property("times", &step_meca_py::Get_times)
        .def_readwrite("BC_Time", &step_meca_py::BC_Time)
        .def_readwrite("file", &step_meca_py::file)
        .add_property("cBC_meca", &step_meca_py::Get_cBC_meca, &step_meca_py::Set_cBC_meca)
        .add_property("BC_meca", &step_meca_py::Get_BC_meca, &step_meca_py::Set_BC_meca)
        .add_property("mecas", &step_meca_py::Get_mecas)
        .add_property("BC_mecas", &step_meca_py::Get_BC_mecas)
        .add_property("BC_w", &step_meca_py::Get_BC_w, &step_meca_py::Set_BC_w)
        .add_property("BC_R", &step_meca_py::Get_BC_R, &step_meca_py::Set_BC_R)
        .add_property("Ts", &step_meca_py::Get_Ts)
        .add_property("BC_Ts", &step_meca_py::Get_BC_Ts)
        .def("generate", &step_meca_py::generate)
        .def("generate_kin", &step_meca_py::generate_kin)
        ;
    
    bp::class_<step_thermomeca_py>("step_thermomeca", bp::init <>())
        .def(bp::init <const int &, const double &, const double &, const double &, const int &, const unsigned int &, const bn::ndarray&, const bn::ndarray&, const double&, const int&, const bn::ndarray&, const bn::ndarray&> ())
        .def_readwrite("number", &step_thermomeca_py::number)
        .def_readwrite("Dn_init",&step_thermomeca_py::Dn_init)
        .def_readwrite("Dn_mini", &step_thermomeca_py::Dn_mini)
        .def_readwrite("Dn_inc", &step_thermomeca_py::Dn_inc)
        .def_readwrite("mode", &step_thermomeca_py::mode)
        .def_readwrite("control_type", &step_thermomeca_py::control_type)
        .add_property("times", &step_thermomeca_py::Get_times)
        .def_readwrite("BC_Time", &step_thermomeca_py::BC_Time)
        .def_readwrite("file", &step_thermomeca_py::file)
        .add_property("cBC_meca", &step_thermomeca_py::Get_cBC_meca, &step_thermomeca_py::Set_cBC_meca)
        .add_property("BC_meca", &step_thermomeca_py::Get_BC_meca, &step_thermomeca_py::Set_BC_meca)
        .add_property("mecas", &step_thermomeca_py::Get_mecas)
        .add_property("BC_mecas", &step_thermomeca_py::Get_BC_mecas)
        .add_property("BC_w", &step_thermomeca_py::Get_BC_w, &step_thermomeca_py::Set_BC_w)
        .add_property("BC_R", &step_thermomeca_py::Get_BC_R, &step_thermomeca_py::Set_BC_R)
        .add_property("Ts", &step_thermomeca_py::Get_Ts)
        .add_property("BC_Ts", &step_meca_py::Get_BC_Ts)
        .def("generate", &step_thermomeca_py::generate)
        .def("generate_kin", &step_thermomeca_py::generate_kin)
        ;
    */
    
    // Register the from-python converters for ODF functions
    bp::def("get_densities_ODF", get_densities_ODF);
    bp::def("ODF_discretization", ODF_discretization);
    
    //Register the from-python converters for rotation
    bp::def("rotate_vec_R", rotate_vec_R);
    bp::def("rotate_vec_angle", rotate_vec_angle);
    bp::def("rotate_mat_R", rotate_mat_R);
    bp::def("rotate_mat_angle", rotate_mat_angle);
    bp::def("fillR_angle", fillR_angle);
    bp::def("fillR_euler", fillR_euler);
    bp::def("fillQS_angle", fillQS_angle);
    bp::def("fillQS_R", fillQS_R);
    bp::def("fillQE_angle", fillQE_angle);
    bp::def("fillQE_R", fillQE_R);

    bp::def("rotateL_angle", rotateL_angle);
    bp::def("rotateL_R", rotateL_R);
    bp::def("rotate_l2g_L", rotate_l2g_L);
    bp::def("rotate_g2l_L", rotate_g2l_L);
    bp::def("rotateM_angle", rotateM_angle);
    bp::def("rotateM_angle", rotateM_angle);
    bp::def("rotate_l2g_M", rotate_l2g_M);
    bp::def("rotate_g2l_M", rotate_g2l_M);
    bp::def("rotateA_angle", rotateA_angle);
    bp::def("rotateA_R", rotateA_R);
    bp::def("rotate_l2g_A", rotate_l2g_A);
    bp::def("rotate_g2l_A", rotate_g2l_A);
    bp::def("rotateB_angle", rotateB_angle);
    bp::def("rotateB_R", rotateB_R);
    bp::def("rotate_l2g_B", rotate_l2g_B);
    bp::def("rotate_g2l_B", rotate_g2l_B);
    bp::def("rotate_strain_angle", rotate_strain_angle);
    bp::def("rotate_strain_R", rotate_strain_R);
    bp::def("rotate_l2g_strain", rotate_l2g_strain);
    bp::def("rotate_g2l_strain", rotate_g2l_strain);
    bp::def("rotate_stress_angle", rotate_stress_angle);
    bp::def("rotate_stress_R", rotate_stress_R);
    bp::def("rotate_l2g_stress", rotate_l2g_stress);
    bp::def("rotate_g2l_stress", rotate_g2l_stress);
    
    //Register the from-python converters for lagrange
    bp::def("lagrange_exp", lagrange_exp);
    bp::def("dlagrange_exp", dlagrange_exp);
    bp::def("lagrange_pow_0", lagrange_pow_0);
    bp::def("dlagrange_pow_0", dlagrange_pow_0);
    bp::def("lagrange_pow_1", lagrange_pow_1);
    bp::def("dlagrange_pow_1", dlagrange_pow_1);
    bp::def("d2lagrange_pow_1", d2lagrange_pow_1);
    
    ///////##### Module for identification ############///////////////////
    // Generation of the constant class
    bp::class_<simcoon::constants>("constants")
    .def(bp::init<const int &, const int &>())
    .def("__init__", build_constants_full)
    .def_readwrite("number", &simcoon::constants::number)
    .def_readwrite("value", &simcoon::constants::value)
    .def_readwrite("key", &simcoon::constants::key)
    .def_readwrite("ninput_files", &simcoon::constants::ninput_files)
    .add_property("input_values", constants_get_input_values, constants_set_input_values)
    .add_property("input_files", constants_get_input_files, constants_set_input_files)
    ;
    
    // Generation of the parameters class
    bp::class_<simcoon::parameters>("parameters")
    .def(bp::init<const int &, const double &, const double &>())
    .def("__init__", build_parameters_full)
    .def_readwrite("number", &simcoon::parameters::number)
    .def_readwrite("value", &simcoon::parameters::value)
    .def_readwrite("min_value", &simcoon::parameters::min_value)
    .def_readwrite("max_value", &simcoon::parameters::max_value)
    .def_readwrite("key", &simcoon::parameters::key)
    .def_readwrite("ninput_files", &simcoon::parameters::ninput_files)
    .add_property("input_files", parameters_get_input_files, parameters_set_input_files)
    ;
    
    // Register the identification solver
    bp::def("identification", identification);
    
    // Register the calc cost
    bp::def("calc_cost", calc_cost);
    
    //Functions related to constants
    bp::def("read_constants", read_constants_py);
    bp::def("copy_constants", copy_constants_py);
    bp::def("apply_constants", apply_constants_py);
    
    //Functions related to parameters
    bp::def("read_parameters", read_parameters_py);
    bp::def("copy_parameters", copy_parameters_py);
    bp::def("apply_parameters", apply_parameters_py);
    
    // Register the function specific for the solver
    bp::def("cost_solver", cost_solver);

    ///////##### Module for Unit_cell ############///////////////////////////////
    // Generation of the Node class
    bp::class_<simcoon::Node>("Node")
    .def("__init__", build_node)
    .def_readwrite("number", &simcoon::Node::number)
    .add_property("coords", Node_get_input_coords, Node_set_input_coords)
    ;
    
    bp::def("nonperioMPC", build_MPC_from_cubic_mesh);
    bp::def("testPerioMesh", test_mesh);
    
    // Generation of the cubic_mesh class
    bp::class_<simcoon::cubic_mesh>("cubic_mesh")
    .def("__init__", build_cubic_mesh)
    .def_readwrite("is_perio", &simcoon::cubic_mesh::is_perio)
    .def_readwrite("Node_list", &simcoon::cubic_mesh::Node_list)
    .def_readwrite("Node_list_name", &simcoon::cubic_mesh::Node_list_name)
    .def_readwrite("volume", &simcoon::cubic_mesh::volume)
    .def_readwrite("Dx", &simcoon::cubic_mesh::Dx)
    .def_readwrite("Dy", &simcoon::cubic_mesh::Dy)
    .def_readwrite("Dz", &simcoon::cubic_mesh::Dz)
    .def_readwrite("size_box", &simcoon::cubic_mesh::size_box)
    .def("get_domain", get_domain)
    .def("construct_lists", construct_lists)
    ;
    
    
}
