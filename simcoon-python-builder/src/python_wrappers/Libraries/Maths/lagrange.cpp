#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>

#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/lagrange.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function is used to determine an exponential Lagrange Multiplier (like contact in Abaqus)
double lagrange_exp(const double &h, const double &c, const double &p0) {
    return simcoon::lagrange_exp(h,c,p0);
}

//This function is used to determine the first derivative of an exponential Lagrange Multiplier
double dlagrange_exp(const double &h, const double &c, const double &p0) {
    return simcoon::dlagrange_exp(h,c,p0);
}
    
//This function is used to determine a power-law Lagrange Multiplier for problem such x >= 0
double lagrange_pow_0(const double &x, const double &c, const double &p0, const double &n, const double &alpha) {
    return simcoon::lagrange_pow_0(x,c,p0,n,alpha);
}

//This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x >= 0
double dlagrange_pow_0(const double &x, const double &c, const double &p0, const double &n, const double &alpha) {
    return simcoon::dlagrange_pow_0(x,c,p0,n,alpha);
}

//This function is used to determine a power-law Lagrange Multiplier for problem such x <= 1
double lagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha) {
    return simcoon::lagrange_pow_1(x,c,p0,n,alpha);
}

//This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x <= 1
double dlagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha) {
    return simcoon::dlagrange_pow_1(x,c,p0,n,alpha);
}

//This function is used to determine the SECOND derivative of a power-law Lagrange Multiplier for problem such x <= 1
double d2lagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha) {
    return simcoon::d2lagrange_pow_1(x,c,p0,n,alpha);
}
    
    
} //namepsace simpy