#include <assert.h>
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/numpy_cgal.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>
#include <iostream>
#include <fstream>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>

#include <simcoon/python_wrappers/Libraries/Solver/step_thermomeca.hpp>

using namespace std;
using namespace arma;
using namespace arma2numpy;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

namespace simpy {

step_thermomeca_py::step_thermomeca_py() : simcoon::step_thermomeca() {}

//-------------------------------------------------------------
step_thermomeca_py::step_thermomeca_py(const int &mnumber, const double &mDn_init, const double &mDn_mini, const double &mDn_inc, const int &mmode, const unsigned int &mcontrol_type, const bn::ndarray &mcBC_meca, const bn::ndarray &mBC_meca, const bn::ndarray &mmecas, const double &mBC_T, const int &mcBC_T, const bn::ndarray &mTs, const bn::ndarray &mBC_w, const bn::ndarray &mBC_R)
//-------------------------------------------------------------
{
    assert(mnumber>=0);
    assert(mDn_inc<=1.);
    assert(mDn_init<=mDn_inc);
    assert(mDn_mini<=mDn_init);
    assert(fmod(1.,mDn_inc)==0.);
    
    number = mnumber;
    Dn_init = mDn_init;
    Dn_mini = mDn_mini;
    Dn_inc = mDn_inc;
    ninc = std::round(1./mDn_inc);
    mode = mmode;
    control_type = mcontrol_type;
    
    times = zeros(ninc);
    BC_Time = 0.;
    
    file = "";
        
    cBC_meca = array2Col_int(mcBC_meca, false);
    BC_meca = array2vec(mBC_meca, false);
    mecas = array2mat(mmecas, false);
    BC_T = mBC_T;
    cBC_T = mcBC_T;
    Ts = array2vec(mTs, false);
    BC_w = array2mat(mBC_w, false);
    BC_R = array2mat(mBC_R, false);
}

step_thermomeca_py::step_thermomeca_py(const simcoon::step_thermomeca &sttm) : simcoon::step_thermomeca(sttm) { }

//----------------------------------------------------------------------
void step_thermomeca_py::generate(const double& mTime, const bn::ndarray& mEtot_py, const bn::ndarray& msigma_py, const double& mT)
//----------------------------------------------------------------------
{
    vec mEtot = array2vec(mEtot_py, true);
    vec msigma = array2vec(msigma_py, true);
    simcoon::step_thermomeca::generate(mTime, mEtot, msigma, mT);
}

//----------------------------------------------------------------------
void step_thermomeca_py::generate_kin(const double& mTime, const bn::ndarray& mF_py, const double& mT)
//----------------------------------------------------------------------
{
    mat mF = array2mat(mF_py, true);
    simcoon::step_thermomeca::generate_kin(mTime, mF, mT);
}

//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_times() {
//-------------------------------------------------------------
    return vec2array(times, false);
}

//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_cBC_meca() {
//-------------------------------------------------------------
    return Col_int2array(cBC_meca);
}

//------------------------------------------------------
void step_thermomeca_py::Set_cBC_meca(const bn::ndarray &cBC_meca_py) {
    cBC_meca = array2Col_int(cBC_meca_py);
}

//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_BC_meca() {
//-------------------------------------------------------------
    return vec2array(BC_meca, false);
}

//------------------------------------------------------
void step_thermomeca_py::Set_BC_meca(const bn::ndarray &BC_meca_py) {
    BC_meca = array2vec(BC_meca_py, false);
}


//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_mecas() {
//-------------------------------------------------------------
    return mat2array(mecas, false);
}

//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_BC_w() {
//-------------------------------------------------------------
    return mat2array(BC_w, false);
}

//------------------------------------------------------
void step_thermomeca_py::Set_BC_w(const bn::ndarray &BC_w_py) {
    BC_w = array2vec(BC_w_py, false);
}

//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_BC_R() {
//-------------------------------------------------------------
    return mat2array(BC_R, false);
}

//------------------------------------------------------
void step_thermomeca_py::Set_BC_R(const bn::ndarray &BC_R_py) {
    BC_R = array2vec(BC_R_py, false);
}

//-------------------------------------------------------------
bn::ndarray step_thermomeca_py::Get_Ts() {
//-------------------------------------------------------------
    return vec2array(Ts, false);
}

}

