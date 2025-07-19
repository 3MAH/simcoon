#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>


#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_L_elastic.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Leff.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//Return the elastic stiffness tensor of a composite material
py::array_t<double> L_eff(const std::string &umat_name, const py::array_t<double> &props, const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve) {

    vec props_cpp = carma::arr_to_col(props);
    
    double T_init = 273.15;
    
    simcoon::phase_characteristics rve;
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props_cpp.n_elem, props_cpp);
    rve.construct(0,1);
    simcoon::natural_basis nb;
    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3),T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
    
    auto sv_M = std::dynamic_pointer_cast<simcoon::state_variables_M>(rve.sptr_sv_global);
    
    //Second we call a recursive method that find all the elastic moduli iof the phases
    simcoon::get_L_elastic(rve);
    
    return carma::mat_to_arr(sv_M->Lt);
}

} //namepsace simpy
