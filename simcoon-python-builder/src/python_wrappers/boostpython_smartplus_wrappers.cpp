
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <CGAL/Simple_cartesian.h>

#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/numpy_cgal.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/constitutive.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/contimech.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/criteria.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/recovery_props.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Leff.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/objective_rates.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Umat_fedoo.hpp>

#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/lagrange.hpp>
#include <simcoon/python_wrappers/Libraries/Material/ODF.hpp>
#include <simcoon/python_wrappers/Libraries/Homogenization/eshelby.hpp>

#include <simcoon/python_wrappers/Libraries/Solver/read.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver.hpp>

#include <simcoon/python_wrappers/Libraries/Identification/identification.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/constants.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/parameters.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/optimize.hpp>

#include <simcoon/python_wrappers/Libraries/Unit_cell/unit_cell.hpp>

//#include <simcoon/python_wrappers/Libraries/Abaqus/write.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace simpy;

BOOST_PYTHON_MODULE(simmit) {

    Py_Initialize();
    bn::initialize();
    
    // Register the from-python converters for constitutive
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
    
    // Register the from-python converters for contimech
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
    
    // Register the from-python converters for criteria
    bp::def("Prager_stress", Prager_stress);
    bp::def("dPrager_stress", dPrager_stress);
    bp::def("Tresca_stress", Tresca_stress);
    bp::def("dTresca_stress", dTresca_stress);
    bp::def("Eq_stress", Eq_stress);
    bp::def("dEq_stress", dEq_stress);
    
    // Register the from-python converters for recovery_props
    bp::def("check_symetries", check_symetries);
    bp::def("L_iso_props", L_iso_props);
    bp::def("M_iso_props", M_iso_props);
    bp::def("L_isotrans_props", L_isotrans_props);
    bp::def("M_isotrans_props", M_isotrans_props);
    bp::def("L_cubic_props", L_cubic_props);
    bp::def("M_cubic_props", M_cubic_props);
    bp::def("L_ortho_props", L_ortho_props);
    bp::def("M_ortho_props", M_ortho_props);
    bp::def("M_aniso_props", M_aniso_props);
    
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
    bp::def("finite_W", finite_W);
    bp::def("finite_Omega", finite_Omega);
    bp::def("finite_D", finite_D);
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
    bp::def("solver", solver);
	
	//Wrapper fedoo
	bp::class_<Umat_fedoo>("Umat_fedoo", bp::init < std::string, bn::ndarray, int, int, int> ())
		.def("compute_Detot", &Umat_fedoo::compute_Detot)
		.def("Run", &Umat_fedoo::Run)
		.def("Initialize", &Umat_fedoo::Initialize)
		.def("to_start", &Umat_fedoo::to_start)
		.def("set_start", &Umat_fedoo::set_start)
		.def_readwrite("corate", &Umat_fedoo::corate)
		.def_readonly("Time", &Umat_fedoo::Time)
		.def_readonly("DTime", &Umat_fedoo::DTime)
		.def_readonly("nb_points", &Umat_fedoo::nb_points)
		.def_readonly("nlgeom", &Umat_fedoo::nlgeom)
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
