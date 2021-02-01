
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_L_elastic.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Leff.hpp>


namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//Check the material symetries and the type of elastic response for a given stiffness tensor
bn::ndarray L_eff(const std::string &umat_name_py, const bn::ndarray &props_py, const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve) {

    vec props = array2vec(props_py);
    
//    std::string umat_name = bp::extract<std::string>(umat_name_py);
//    std::string path_data = bp::extract<std::string>(path_data_py);
    double T_init = 273.15;
    
    simcoon::phase_characteristics rve;
    rve.sptr_matprops->update(0, umat_name_py, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve.construct(0,1);
    simcoon::natural_basis nb;
    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(3,3), zeros(3,3), eye(3,3), eye(3,3),T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
    
    auto sv_M = std::dynamic_pointer_cast<simcoon::state_variables_M>(rve.sptr_sv_global);
    
    //Second we call a recursive method that find all the elastic moduli iof the phases
    simcoon::get_L_elastic(rve);
    
    return mat2array(sv_M->Lt);
}

} //namepsace simpy
