#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Densities of an orientation distribution function at the angles x. `peaks` is a sequence
//of dicts {number, method, mean, s_dev, width, ampl, params}; with radian = false, x and
//the peak angles are degrees.
pybind11::array_t<double> get_densities_ODF(const pybind11::array_t<double> &x, const pybind11::object &peaks, const bool &radian);

//Discretise the sub-phase `num_phase_disc` of a mean-field RVE into `nphases_disc` phases
//spread over [angle_min, angle_max] degrees according to an ODF. `phases` and `peaks` come
//in memory (the dicts of simcoon.solver.micromechanics); the discretised RVE goes back the
//same way, as a list of phase dicts. `angle` selects the Euler angle (0 psi, 1 theta,
//2 phi); `angles_mat` also rotates the material frame of the phases.
pybind11::list ODF_discretization(const pybind11::object &phases, const pybind11::object &peaks,
                                  const std::string &umat_name, const pybind11::array_t<double> &props,
                                  const int &num_phase_disc, const int &nphases_disc,
                                  const double &angle_min, const double &angle_max,
                                  const bool &angles_mat, const int &angle);

} //namespace simpy
