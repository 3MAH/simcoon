#include <assert.h>
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/numpy_cgal.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>
#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

#include <simcoon/python_wrappers/Libraries/Phase/state_variables.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

using namespace std;
using namespace arma;
using namespace arma2numpy;
using namespace cgal2numpy;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

namespace simpy {

state_variables_T_py::state_variables_T_py() : simcoon::state_variables_T() {}

//-------------------------------------------------------------
/*state_variables_py::state_variables_py(simcoon::state_variables &self, const bn::ndarray& etot_py, const bn::ndarray& F0_py, const bn::ndarray& F1_py, const bn::ndarray& sigma_py, const bn::ndarray& statev_py, const double &T_py, const double &DT_py) {
    
    self.etot = array2vec(etot_py, false);
    self.Detot = zeros(6);

    self.Etot = zeros(6);
    self.DEtot = zeros(6);
    
    self.sigma_start = array2vec(sigma_py, false);
    self.sigma = self.sigma_start;

    self.tau = zeros(6);
    self.tau_start = zeros(6);
    
    self.PKII = zeros(6);
    self.PKII_start = zeros(6);
    
    self.F0 = arrayT2mat_inplace(F0_py);
    self.F1 = arrayT2mat_inplace(F1_py);

    self.R = zeros(3,3);
    self.DR = zeros(3,3);
    
    self.T = T_py;
    self.DT = DT_py;
        
    self.statev = array2vec(statev_py, false);
    self.statev_start = self.statev;
        
    self.nstatev= self.statev.n_elem;    
}*/

//-------------------------------------------------------------
state_variables_T_py::state_variables_T_py(const bn::ndarray& etot_py, const bn::ndarray& F0_py, const bn::ndarray& F1_py, const bn::ndarray& sigma_py, const bn::ndarray& statev_py, const double &T_py, const double &DT_py) : simcoon::state_variables() {
    
    etot = array2vec(etot_py, false);
    Detot = zeros(6);

    sigma_start = array2vec(sigma_py, false);
    sigma = sigma_start;

    F0 = arrayT2mat_inplace(F0_py);
    F1 = arrayT2mat_inplace(F1_py);

    T = T_py;
    DT = DT_py;
        
    statev = array2vec(statev_py, false);
    statev_start = statev;
        
    nstatev= statev.n_elem;
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_F0() {
//-------------------------------------------------------------
    return matT2array_inplace(F0);
}

//------------------------------------------------------
void state_variables_T_py::Set_F0(const bn::ndarray &F0_py) {
    F0 = arrayT2mat_inplace(F0_py);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_F1() {
//-------------------------------------------------------------
    return matT2array_inplace(F1);
}

//------------------------------------------------------
void state_variables_T_py::Set_F1(const bn::ndarray &F1_py) {
    F1 = arrayT2mat_inplace(F1_py);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_etot() {
//-------------------------------------------------------------
    return vec2array(etot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_Detot() {
//-------------------------------------------------------------
    return vec2array(Detot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_Etot() {
//-------------------------------------------------------------
    return vec2array(Etot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_DEtot() {
//-------------------------------------------------------------
    return vec2array(DEtot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_statev() {
//-------------------------------------------------------------
    return vec2array(statev, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_R() {
//-------------------------------------------------------------
    return mat2array_inplace(R);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_DR() {
//-------------------------------------------------------------
    return mat2array_inplace(DR);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_Wm() {
//-------------------------------------------------------------
    return vec2array(Wm, false)
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_Wt() {
//-------------------------------------------------------------
    return vec2array(Wt, false)
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_dSdE() {
//-------------------------------------------------------------
    return matT2array_inplace(dSdE);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_dSdEt() {
//-------------------------------------------------------------
    return matT2array_inplace(dSdEt);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_Get_dSdT() {
//-------------------------------------------------------------
    return matT2array_inplace(Get_dSdT);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_drdE() {
//-------------------------------------------------------------
    return matT2array_inplace(drdE);
}

//-------------------------------------------------------------
bn::ndarray state_variables_T_py::Get_drdT() {
//-------------------------------------------------------------
    return matT2array_inplace(drdT);
}

//----------------------------------------------------------------------
state_variables_T_py& state_variables_T_py::rotate_l2g(const state_variables_T_py& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{
    simcoon::state_variables_T::rotate_l2g(sv, psi, theta, phi);
}

//----------------------------------------------------------------------
state_variables_T_py& state_variables_T_py::rotate_g2l(const state_variables_T_py& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{
    simcoon::state_variables_T::rotate_g2l(sv, psi, theta, phi);
}

}
