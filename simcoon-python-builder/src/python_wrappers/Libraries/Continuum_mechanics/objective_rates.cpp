
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/objective_rates.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
bp::tuple logarithmic(const bn::ndarray &nF0, const bn::ndarray &nF1, const double &DTime) {
    mat F0 = array2mat(nF0);
    mat F1 = array2mat(nF1);
    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    simcoon::logarithmic(DR, D, Omega, DTime, F0, F1);
    bn::ndarray nD = mat2array(D);
    bn::ndarray nDR = mat2array(DR);
    bn::ndarray nOmega = mat2array(Omega);
    return boost::python::make_tuple(nD, nDR, nOmega);
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
bn::ndarray Delta_log_strain(const bn::ndarray &nD, const bn::ndarray &nOmega, const double &DTime) {
    mat D = array2mat(nD);
    mat Omega = array2mat(nOmega);
    return mat2array(simcoon::Delta_log_strain(D, Omega, DTime));
}
    
} //namepsace simpy
