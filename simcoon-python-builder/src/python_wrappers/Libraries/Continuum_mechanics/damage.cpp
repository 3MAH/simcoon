
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>

#include <simcoon/Continuum_mechanics/Functions/damage.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/damage.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function returns damage evolution (/dt) considering a Weibull damage law
double damage_weibull(const py::array_t<double> &stress, const double &damage, const double &alpha, const double &beta, const double &DTime, const std::string &criterion) {
    vec stress_cpp = simpy::arr_to_col(stress);
    return simcoon::damage_weibull(stress_cpp,damage,alpha,beta,DTime,criterion);
}

//This function returns damage evolution (/dt) considering Kachanov's creep damage law
double damage_kachanov(const py::array_t<double> &stress, const py::array_t<double> &strain, const double &damage, const double &A0, const double &r, const std::string &criterion) {
    vec stress_cpp = simpy::arr_to_col(stress);
    vec strain_cpp = simpy::arr_to_col(strain);
    return simcoon::damage_kachanov(stress_cpp,strain_cpp,damage,A0,r,criterion);
}

//This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law
double damage_miner(const double &S_max, const double &S_mean, const double &S_ult, const double &b, const double &B0, const double &beta, const double &Sl_0) {
    return simcoon::damage_miner(S_max,S_mean,S_ult,b,B0,beta,Sl_0);
}

//This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law
double damage_manson(const double &S_amp, const double &C2, const double &gamma2) {
    return simcoon::damage_manson(S_amp,C2,gamma2);
}

} //namepsace simpy