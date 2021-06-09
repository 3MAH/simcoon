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
#include <simcoon/Simulation/Maths/rotation.hpp>
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

state_variables_py::state_variables_py() : simcoon::state_variables() {}

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
    
    //self.F0 = arrayT2mat_inplace(F0_py);
    //self.F1 = arrayT2mat_inplace(F1_py);
	self.F0 = array2mat(F0_py, false); //without copy 
	self.F1 = array2mat(F1_py, false); //without copy 

    self.R = zeros(3,3);
    self.DR = zeros(3,3);
    
    self.T = T_py;
    self.DT = DT_py;
        
    self.statev = array2vec(statev_py, false);
    self.statev_start = self.statev;
        
    self.nstatev= self.statev.n_elem;    
}*/

//-------------------------------------------------------------
state_variables_py::state_variables_py(const bn::ndarray& etot_py, const bn::ndarray& F0_py, const bn::ndarray& F1_py, const bn::ndarray& sigma_py, const bn::ndarray& statev_py, const double &T_py, const double &DT_py) : simcoon::state_variables() {
    
    etot = array2vec(etot_py, false);
    Detot = zeros(6);

    sigma_start = array2vec(sigma_py, false);
    sigma = sigma_start;

    //F0 = arrayT2mat_inplace(F0_py);
    //F1 = arrayT2mat_inplace(F1_py);
	F0 = array2mat(F0_py, false); //without copy
	F1 = array2mat(F1_py, false); //without copy
    
	T = T_py;
    DT = DT_py;
        
    statev = array2vec(statev_py, false);
    statev_start = statev;
        
    nstatev= statev.n_elem;
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_F0() {
//-------------------------------------------------------------
    //return matT2array_inplace(F0);
	return mat2array(F0, false); //without copy
}

//------------------------------------------------------
void state_variables_py::Set_F0(const bn::ndarray &F0_py) {
    //F0 = arrayT2mat_inplace(F0_py);
	F0 = array2mat(F0_py, false); //without copy
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_F1() {
//-------------------------------------------------------------
    //return matT2array_inplace(F1);
	return mat2array(F1, false); //without copy
}

//------------------------------------------------------
void state_variables_py::Set_F1(const bn::ndarray &F1_py) {
    //F1 = arrayT2mat_inplace(F1_py);
	F1 = array2mat(F1_py, false); //without copy
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_etot() {
//-------------------------------------------------------------
    return vec2array(etot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_Detot() {
//-------------------------------------------------------------
    return vec2array(Detot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_Etot() {
//-------------------------------------------------------------
    return vec2array(Etot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_DEtot() {
//-------------------------------------------------------------
    return vec2array(DEtot, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_statev() {
//-------------------------------------------------------------
    return vec2array(statev, false);
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_R() {
//-------------------------------------------------------------
    //return mat2array_inplace(R);
	return mat2array(R, false, "F");
}

//-------------------------------------------------------------
bn::ndarray state_variables_py::Get_DR() {
//-------------------------------------------------------------
    //return mat2array_inplace(DR);
	return mat2array(DR, false, "F");
}

//----------------------------------------------------------------------
void state_variables_py::rotate_l2g(const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{
      if(fabs(phi) > sim_iota) {
        Etot = simcoon::rotate_strain(Etot, -phi, axis_phi);
        DEtot = simcoon::rotate_strain(DEtot, -phi, axis_phi);
        etot = simcoon::rotate_strain(etot, -phi, axis_phi);
        Detot = simcoon::rotate_strain(Detot, -phi, axis_phi);
        PKII = simcoon::rotate_stress(PKII, -phi, axis_phi);
        PKII_start = simcoon::rotate_stress(PKII_start, -phi, axis_phi);
        tau = simcoon::rotate_stress(tau, -phi, axis_phi);
        tau_start = simcoon::rotate_stress(tau_start, -phi, axis_phi);
        sigma = simcoon::rotate_stress(sigma, -phi, axis_phi);
        sigma_start = simcoon::rotate_stress(sigma_start, -phi, axis_phi);
        F0 = simcoon::rotate_mat(F0, -phi, axis_phi);
        F1 = simcoon::rotate_mat(F1, -phi, axis_phi);
        R = simcoon::rotate_mat(R, -phi, axis_phi);
        DR = simcoon::rotate_mat(DR, -phi, axis_phi);
    }
      if(fabs(theta) > sim_iota) {
        Etot = simcoon::rotate_strain(Etot, -theta, axis_theta);
        DEtot = simcoon::rotate_strain(DEtot, -theta, axis_theta);
        etot = simcoon::rotate_strain(etot, -theta, axis_theta);
        Detot = simcoon::rotate_strain(Detot, -theta, axis_theta);
        PKII = simcoon::rotate_stress(PKII, -theta, axis_theta);
        PKII_start = simcoon::rotate_stress(PKII_start, -theta, axis_theta);
        tau = simcoon::rotate_stress(tau, -theta, axis_theta);
        tau_start = simcoon::rotate_stress(tau_start, -theta, axis_theta);
        sigma = simcoon::rotate_stress(sigma, -theta, axis_theta);
        sigma_start = simcoon::rotate_stress(sigma_start, -theta, axis_theta);
        F0 = simcoon::rotate_mat(F0, -theta, axis_theta);
        F1 = simcoon::rotate_mat(F1, -theta, axis_theta);
        R = simcoon::rotate_mat(R, -theta, axis_theta);
        DR = simcoon::rotate_mat(DR, -theta, axis_theta);
    }
    if(fabs(psi) > sim_iota) {
        Etot = simcoon::rotate_strain(Etot, -psi, axis_psi);
        DEtot = simcoon::rotate_strain(DEtot, -psi, axis_psi);
        etot = simcoon::rotate_strain(etot, -psi, axis_psi);
        Detot = simcoon::rotate_strain(Detot, -psi, axis_psi);
        PKII = simcoon::rotate_stress(PKII, -psi, axis_psi);
        PKII_start = simcoon::rotate_stress(PKII_start, -psi, axis_psi);
        tau = simcoon::rotate_stress(tau, -psi, axis_psi);
        tau_start = simcoon::rotate_stress(tau_start, -psi, axis_psi);
        sigma = simcoon::rotate_stress(sigma, -psi, axis_psi);
        sigma_start = simcoon::rotate_stress(sigma_start, -psi, axis_psi);
        F0 = simcoon::rotate_mat(F0, -psi, axis_psi);
        F1 = simcoon::rotate_mat(F1, -psi, axis_psi);
        R = simcoon::rotate_mat(R, -psi, axis_psi);
        DR = simcoon::rotate_mat(DR, -psi, axis_psi);
    }
}

//----------------------------------------------------------------------
void state_variables_py::rotate_g2l(const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{
      if(fabs(phi) > sim_iota) {
        Etot = simcoon::rotate_strain(Etot, -phi, axis_phi);
        DEtot = simcoon::rotate_strain(DEtot, -phi, axis_phi);
        etot = simcoon::rotate_strain(etot, -phi, axis_phi);
        Detot = simcoon::rotate_strain(Detot, -phi, axis_phi);
        PKII = simcoon::rotate_stress(PKII, -phi, axis_phi);
        PKII_start = simcoon::rotate_stress(PKII_start, -phi, axis_phi);
        tau = simcoon::rotate_stress(tau, -phi, axis_phi);
        tau_start = simcoon::rotate_stress(tau_start, -phi, axis_phi);
        sigma = simcoon::rotate_stress(sigma, -phi, axis_phi);
        sigma_start = simcoon::rotate_stress(sigma_start, -phi, axis_phi);
        F0 = simcoon::rotate_mat(F0, -phi, axis_phi);
        F1 = simcoon::rotate_mat(F1, -phi, axis_phi);
        R = simcoon::rotate_mat(R, -phi, axis_phi);
        DR = simcoon::rotate_mat(DR, -phi, axis_phi);
    }
      if(fabs(theta) > sim_iota) {
        Etot = simcoon::rotate_strain(Etot, -theta, axis_theta);
        DEtot = simcoon::rotate_strain(DEtot, -theta, axis_theta);
        etot = simcoon::rotate_strain(etot, -theta, axis_theta);
        Detot = simcoon::rotate_strain(Detot, -theta, axis_theta);
        PKII = simcoon::rotate_stress(PKII, -theta, axis_theta);
        PKII_start = simcoon::rotate_stress(PKII_start, -theta, axis_theta);
        tau = simcoon::rotate_stress(tau, -theta, axis_theta);
        tau_start = simcoon::rotate_stress(tau_start, -theta, axis_theta);
        sigma = simcoon::rotate_stress(sigma, -theta, axis_theta);
        sigma_start = simcoon::rotate_stress(sigma_start, -theta, axis_theta);
        F0 = simcoon::rotate_mat(F0, -theta, axis_theta);
        F1 = simcoon::rotate_mat(F1, -theta, axis_theta);
        R = simcoon::rotate_mat(R, -theta, axis_theta);
        DR = simcoon::rotate_mat(DR, -theta, axis_theta);
    }
    if(fabs(psi) > sim_iota) {
        Etot = simcoon::rotate_strain(Etot, -psi, axis_psi);
        DEtot = simcoon::rotate_strain(DEtot, -psi, axis_psi);
        etot = simcoon::rotate_strain(etot, -psi, axis_psi);
        Detot = simcoon::rotate_strain(Detot, -psi, axis_psi);
        PKII = simcoon::rotate_stress(PKII, -psi, axis_psi);
        PKII_start = simcoon::rotate_stress(PKII_start, -psi, axis_psi);
        tau = simcoon::rotate_stress(tau, -psi, axis_psi);
        tau_start = simcoon::rotate_stress(tau_start, -psi, axis_psi);
        sigma = simcoon::rotate_stress(sigma, -psi, axis_psi);
        sigma_start = simcoon::rotate_stress(sigma_start, -psi, axis_psi);
        F0 = simcoon::rotate_mat(F0, -psi, axis_psi);
        F1 = simcoon::rotate_mat(F1, -psi, axis_psi);
        R = simcoon::rotate_mat(R, -psi, axis_psi);
        DR = simcoon::rotate_mat(DR, -psi, axis_psi);
    }
}

}
