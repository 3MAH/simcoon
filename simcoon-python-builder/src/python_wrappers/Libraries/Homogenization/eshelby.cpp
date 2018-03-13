
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_Mechanics/Homogenization/eshelby.hpp>
#include <simcoon/python_wrappers/Libraries/Homogenization/eshelby.hpp>

namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//Eshelby tensor for a sphere
bn::ndarray Eshelby_sphere(const double &nu) {
    return mat2array(simcoon::Eshelby_sphere(nu));
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
bn::ndarray Eshelby_cylinder(const double &nu) {
    return mat2array(simcoon::Eshelby_cylinder(nu));
}

//	Eshelby tensor determination. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here
bn::ndarray Eshelby_prolate(const double &nu, const double &ar) {
    return mat2array(simcoon::Eshelby_prolate(nu,ar));
}

//	Eshelby tensor determination. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here
bn::ndarray Eshelby_oblate(const double &nu, const double &ar) {
    return mat2array(simcoon::Eshelby_oblate(nu,ar));
}

//Numerical Eshelby tensor determination
bn::ndarray Eshelby(const bn::ndarray &nd, const double &a1, const double &a2, const double &a3, const int &mp, const int &np) {
    mat Lt = array2mat(nd);
    return mat2array(simcoon::Eshelby(Lt, a1, a2, a3, mp, np));
}
    
//Numerical Hill Interaction tensor determination
bn::ndarray T_II(const bn::ndarray &nd, const double &a1, const double &a2, const double &a3, const int &mp, const int &np) {
    mat Lt = array2mat(nd);
    return mat2array(simcoon::T_II(Lt, a1, a2, a3, mp, np));
}

} //namepsace simpy
