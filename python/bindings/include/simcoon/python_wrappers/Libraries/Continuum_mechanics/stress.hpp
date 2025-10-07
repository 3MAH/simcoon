#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//Provides the first Piola Kirchoff stress tensor from the Cauchy stress tensor
pybind11::array_t<double> stress_convert(const pybind11::array_t<double> &sigma, const pybind11::array_t<double> &F, const std::string &converter_key, const double &J = 0., const bool & copy=true);

/*
//Provides the second Piola Kirchoff stress tensor from the Cauchy stress tensor
pybind11::array_t<double> Cauchy2PKII(const pybind11::array_t<double> &sigma, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the Kirchoff stress tensor from the Cauchy stress tensor
pybind11::array_t<double> Cauchy2Kirchoff(const pybind11::array_t<double> &sigma, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the Cauchy stress tensor from the Kirchoff stress tensor
pybind11::array_t<double> Kirchoff2Cauchy(const pybind11::array_t<double> &tau, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the first Piola Kirchoff stress tensor from the Kirchoff stress tensor
pybind11::array_t<double> Kirchoff2PKI(const pybind11::array_t<double> &tau, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the second Piola Kirchoff stress tensor from the Kirchoff stress tensor
pybind11::array_t<double> Kirchoff2PKII(const pybind11::array_t<double> &tau, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the Kirchoff stress tensor from the first Piola Kirchoff stress tensor
pybind11::array_t<double> PKI2Kirchoff(const pybind11::array_t<double> &Sigma, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the Kirchoff stress tensor from the second Piola Kirchoff stress tensor
pybind11::array_t<double> PKII2Kirchoff(const pybind11::array_t<double> &S, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the Cauchy stress tensor from the first Piola Kirchoff stress tensor
pybind11::array_t<double> PKI2Cauchy(const pybind11::array_t<double> &Sigma, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);

//Provides the Cauchy stress tensor from the second Piola Kirchoff stress tensor
pybind11::array_t<double> PKII2Cauchy(const pybind11::array_t<double> &S, const pybind11::array_t<double> &F, const double &J = 0., const bool & copy=true);
*/

} //namespace simpy