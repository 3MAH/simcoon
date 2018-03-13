
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_Mechanics/Functions/recovery_props.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_Mechanics/recovery_props.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//Check the material symetries and the type of elastic response for a given stiffness tensor
bp::dict check_symetries(const bn::ndarray &nL) {
    
    mat L = array2mat(nL);
    int axis = 0;
    string umat_type;
    vec props;
    int maj_sym = 0;
    simcoon::check_symetries(L, umat_type, axis, props, maj_sym);
    bp::dict d;
    d["umat_type"]  = umat_type;
    d["axis"]  = axis;
    d["maj_sym"]  = maj_sym;
    d["props"]  = vec2array(props);
    return d;
}

//return a list of elastic properties for the isotropic case (E,nu) from a stiffness tensor
bn::ndarray L_iso_props(const bn::ndarray &nLt) {

    mat Lt = array2mat(nLt);
    vec props = simcoon::L_iso_props(Lt);
    return vec2array(props);
}
    
//return a list of elastic properties for the isotropic case (E,nu) from a compliance tensor
bn::ndarray M_iso_props(const bn::ndarray &nMt) {
    
    mat Mt = array2mat(nMt);
    vec props = simcoon::M_iso_props(Mt);
    return vec2array(props);
}

//return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a stiffness tensor
bn::ndarray L_isotrans_props(const bn::ndarray &nLt, const int &axis) {
    mat Lt = array2mat(nLt);
    vec props = simcoon::L_isotrans_props(Lt, axis);
    return vec2array(props);
}
    
//return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a compliance tensor
bn::ndarray M_isotrans_props(const bn::ndarray &nMt, const int &axis) {
    mat Mt = array2mat(nMt);
    vec props = simcoon::M_isotrans_props(Mt, axis);
    return vec2array(props);
}

//return a list of elastic properties for the cubic case (E,nu,G) from a stiffness tensor
bn::ndarray L_cubic_props(const bn::ndarray &nLt) {
    
    mat Lt = array2mat(nLt);
    vec props = simcoon::L_cubic_props(Lt);
    return vec2array(props);
}
    
//return a list of elastic properties for the cubic case (E,nu,G) from a compliance tensor
bn::ndarray M_cubic_props(const bn::ndarray &nMt) {
    
    mat Mt = array2mat(nMt);
    vec props = simcoon::M_cubic_props(Mt);
    return vec2array(props);
}
 
//return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a stiffness tensor
bn::ndarray L_ortho_props(const bn::ndarray &nLt) {
    
    mat Lt = array2mat(nLt);
    vec props = simcoon::L_ortho_props(Lt);
    return vec2array(props);
}

//return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a compliance tensor
bn::ndarray M_ortho_props(const bn::ndarray &nMt) {
    
    mat Mt = array2mat(nMt);
    vec props = simcoon::M_ortho_props(Mt);
    return vec2array(props);
}
    
//return a list of elastic properties for the anisotropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,deviations) from a compliance tensor
bn::ndarray M_aniso_props(const bn::ndarray &nMt) {
    
    mat Mt = array2mat(nMt);
    vec props = simcoon::M_aniso_props(Mt);
    return vec2array(props);
}


} //namepsace simpy