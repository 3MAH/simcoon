
#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/recovery_props.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/recovery_props.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//Check the material symetries and the type of elastic response for a given stiffness tensor
py::dict check_symetries(const py::array_t<double> &input) {
    
    mat L = carma::arr_to_mat(input);
    int axis = 0;
    string umat_type;
    vec props;
    int maj_sym = 0;
    simcoon::check_symetries(L, umat_type, axis, props, maj_sym);
    py::dict d;
    d["umat_type"] = umat_type;
    d["axis"] = axis;
    d["maj_sym"] = maj_sym;
    d["props"] = carma::col_to_arr(props);
    return d;
}

//return a list of elastic properties for the isotropic case (E,nu) from a stiffness tensor
py::array_t<double> L_iso_props(const py::array_t<double> &input) {

    mat Lt = carma::arr_to_mat(input);
    vec props = simcoon::L_iso_props(Lt);
    return carma::col_to_arr(props);
}
    
//return a list of elastic properties for the isotropic case (E,nu) from a compliance tensor
py::array_t<double> M_iso_props(const py::array_t<double> &input) {
    
    mat Mt = carma::arr_to_mat(input);
    vec props = simcoon::M_iso_props(Mt);
    return carma::col_to_arr(props);
}

//return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a stiffness tensor
py::array_t<double> L_isotrans_props(const py::array_t<double> &input, const int &axis) {
    mat Lt = carma::arr_to_mat(input);
    vec props = simcoon::L_isotrans_props(Lt, axis);
    return carma::col_to_arr(props);
}
    
//return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a compliance tensor
py::array_t<double> M_isotrans_props(const py::array_t<double> &input, const int &axis) {
    mat Mt = carma::arr_to_mat(input);
    vec props = simcoon::M_isotrans_props(Mt, axis);
    return carma::col_to_arr(props);
}

//return a list of elastic properties for the cubic case (E,nu,G) from a stiffness tensor
py::array_t<double> L_cubic_props(const py::array_t<double> &input) {
    
    mat Lt = carma::arr_to_mat(input);
    vec props = simcoon::L_cubic_props(Lt);
    return carma::col_to_arr(props);
}
    
//return a list of elastic properties for the cubic case (E,nu,G) from a compliance tensor
py::array_t<double> M_cubic_props(const py::array_t<double> &input) {
    
    mat Mt = carma::arr_to_mat(input);
    vec props = simcoon::M_cubic_props(Mt);
    return carma::col_to_arr(props);
}
 
//return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a stiffness tensor
py::array_t<double> L_ortho_props(const py::array_t<double> &input) {
    
    mat Lt = carma::arr_to_mat(input);
    vec props = simcoon::L_ortho_props(Lt);
    return carma::col_to_arr(props);
}

//return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a compliance tensor
py::array_t<double> M_ortho_props(const py::array_t<double> &input) {
    
    mat Mt = carma::arr_to_mat(input);
    vec props = simcoon::M_ortho_props(Mt);
    return carma::col_to_arr(props);
}
    
//return a list of elastic properties for the anisotropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,deviations) from a compliance tensor
py::array_t<double> M_aniso_props(const py::array_t<double> &input) {
    
    mat Mt = carma::arr_to_mat(input);
    vec props = simcoon::M_aniso_props(Mt);
    return carma::col_to_arr(props);
}


} //namepsace simpy