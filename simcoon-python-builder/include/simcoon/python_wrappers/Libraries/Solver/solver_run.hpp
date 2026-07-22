#pragma once

#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace simpy {

/**
 * @brief In-memory solver entry point for the Python bindings (no file I/O).
 *
 * Runs the same incremental solver loop as the file-based binding, but takes the
 * loading programme as Python data and returns the full converged history in memory.
 *
 * @param[in] blocks_py list of block dicts (control_type, ncycle, steps[] with
 *            cBC_meca/BC_meca/mode/times/tabular... — see simcoon.solver.Block.to_dict)
 * @param[in] T_init initial temperature
 * @param[in] umat_name 5-char UMAT name (dispatch key)
 * @param[in] props_py material properties (units MPa)
 * @param[in] nstatev number of internal state variables
 * @param[in] psi_rve,theta_rve,phi_rve material orientation Euler angles
 * @param[in] solver_type solver scheme (0 = Newton)
 * @param[in] corate_type objective rate (0 Jaumann, 1 Green-Naghdi, 2 log, 3 log_R, 4 Truesdell, 5 log_F)
 * @param[in] params_py numeric controls (div/mul_tnew_dt, miniter, maxiter, inforce,
 *            precision, lambda_solver, tangent_mode)
 * @param[in] record_tangent whether to record the per-increment tangent history
 * @return dict of numpy arrays keyed "Strain", "Stress", "F", "R", "Statev", "Wm", "Temp",
 *         "Time", ... plus "status" (0 = completed) — consumed by simcoon.solver.SolverResults
 */
py::dict solver_run(const py::list &blocks_py, const double &T_init,
                    const std::string &umat_name, const py::array_t<double> &props_py,
                    const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve,
                    const int &solver_type, const int &corate_type,
                    const py::dict &params_py, const bool &record_tangent);

} //namespace simpy
