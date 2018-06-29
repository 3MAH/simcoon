
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/transfer.hpp>

namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function transforms the strain Voigt vector into a 3*3 strain matrix
bn::ndarray v2t_strain(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return mat2array(simcoon::v2t_strain(v));
}

//This function transforms a 3*3 strain matrix into a strain Voigt vector
bn::ndarray t2v_strain (const bn::ndarray &nd) {
    mat m = array2mat(nd);
    return vec2array(simcoon::t2v_strain(m));
}

//This function transforms the stress Voigt vector into a 3*3 stress matrix
bn::ndarray v2t_stress(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return mat2array(simcoon::v2t_stress(v));
}

//This function transforms a 3*3 stress matrix into a stress Voigt vector
bn::ndarray t2v_stress (const bn::ndarray &nd) {
    mat m = array2mat(nd);
    return vec2array(simcoon::t2v_stress(m));
}

} //namepsace simpy