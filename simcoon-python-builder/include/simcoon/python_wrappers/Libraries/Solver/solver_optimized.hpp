#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace simpy {

/**
 * @brief Run the optimized C++ solver with Block/Step objects from Python.
 *
 * This function accepts Python Block and Step objects directly, extracts
 * their configuration, and runs the optimized C++ solver.
 *
 * @param blocks_py List of Block objects from Python
 * @param max_iter Maximum Newton-Raphson iterations (default: 10)
 * @param tol Convergence tolerance (default: 1e-9)
 * @param lambda_solver Penalty stiffness for strain control (default: 10000.0)
 * @return List of HistoryPoint objects matching Python solver output format
 */
py::list solver_optimized(
    const py::list& blocks_py,
    int max_iter = 10,
    double tol = 1e-9,
    double lambda_solver = 10000.0
);

} // namespace simpy
