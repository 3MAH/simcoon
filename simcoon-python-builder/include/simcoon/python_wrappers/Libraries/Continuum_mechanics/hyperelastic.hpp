#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function returns the isochoric invariants from V or b
pybind11::array_t<double> isochoric_invariants(const pybind11::array_t<double> &input, const double &J=0., const bool &copy=true);

//This function returns the isochoric principale stretches from V or b
pybind11::array_t<double> isochoric_pstretch(const pybind11::array_t<double> &input, const std::string &input_tensor="V", const double &J=0., const bool &copy=true);

//This function returns the principal stretches from V or b, and also the principal vectors and orthogonal projectors
pybind11::tuple pstretch(const pybind11::array_t<double> &input, const std::string &input_tensor="V", const bool &N_projectors=true, const double &J=0., const bool &copy=true);

//pybind11::array_t<double> tau_iso_hyper_pstretch(const pybind11::array_t<double> &dWdlambda_bar, const pybind11::array_t<double> &lambda_bar, const pybind11::array_t<double> &N_projectors);

//pybind11::array_t<double> sigma_iso_hyper_pstretch(const pybind11::array_t<double> &dWdlambda_bar, const pybind11::array_t<double> &lambda_bar, const pybind11::array_t<double> &N_projectors);

pybind11::array_t<double> tau_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const pybind11::array_t<double> &input, const double &mJ, const bool &copy=true);

pybind11::array_t<double> sigma_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const pybind11::array_t<double> &input, const double &mJ, const bool &copy=true);

pybind11::array_t<double> tau_vol_hyper(const double &dUdJ, const pybind11::array_t<double> &b, const double &mJ, const bool &copy=true);

pybind11::array_t<double> sigma_vol_hyper(const double &dUdJ, const pybind11::array_t<double> &b, const double &mJ, const bool &copy=true);

//pybind11::array_t<double> L_iso_hyper_pstretch(const pybind11::array_t<double> &dWdlambda_bar, const pybind11::array_t<double> &dW2dlambda_bar2, const pybind11::array_t<double> &lambda_bar, const pybind11::array_t<double> &n_pvectors, const double &J);

//pybind11::array_t<double> L_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const pybind11::array_t<double> &b, const double &mJ);

//pybind11::array_t<double> L_vol_hyper(const double &dUdJ, const pybind11::array_t<double> &b, const double &mJ);


} //namespace simpy