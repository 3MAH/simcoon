/* This file is part of simcoon.

 simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 simcoon is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with simcoon.  If not, see <http://www.gnu.org/licenses/>.

 */

///@file solver_optimized.cpp
///@brief Python bindings for optimized C++ solver
///@version 1.0

#include <carma>
#include <armadillo>
#include <simcoon/python_wrappers/Libraries/Solver/solver_optimized.hpp>
#include <simcoon/Simulation/Solver_optimized/solver_engine.hpp>

namespace py = pybind11;

namespace simpy {

// Control type mapping (matching Python CONTROL_TYPES)
int get_control_type_int(const std::string& control_type) {
    if (control_type == "small_strain") return 1;
    if (control_type == "green_lagrange") return 2;
    if (control_type == "logarithmic") return 3;
    if (control_type == "biot") return 4;
    if (control_type == "F") return 5;
    if (control_type == "gradU") return 6;
    return 1;  // default: small_strain
}

// Corate type mapping (matching Python CORATE_TYPES)
int get_corate_type_int(const std::string& corate_type) {
    if (corate_type == "jaumann") return 0;
    if (corate_type == "green_naghdi") return 1;
    if (corate_type == "logarithmic") return 2;
    if (corate_type == "logarithmic_R") return 3;
    if (corate_type == "truesdell") return 4;
    if (corate_type == "logarithmic_F") return 5;
    return 0;  // default: jaumann
}

// Extract StepConfig from Python Step object
simcoon::SolverEngine::StepConfig extract_step(const py::object& step) {
    simcoon::SolverEngine::StepConfig sc;

    sc.Dn_init = step.attr("Dn_init").cast<int>();
    sc.Dn_mini = step.attr("Dn_mini").cast<int>();
    sc.Dn_inc = step.attr("Dn_inc").cast<int>();
    sc.time = step.attr("time").cast<double>();

    // Get cBC_meca from get_cBC_meca() method
    py::array_t<int> cBC = step.attr("get_cBC_meca")().cast<py::array_t<int>>();
    sc.cBC_meca = carma::arr_to_col<int>(cBC);

    // Get strain/stress targets
    py::array_t<double> DEtot = step.attr("DEtot_end").cast<py::array_t<double>>();
    sc.DEtot_end = carma::arr_to_col(DEtot);

    py::array_t<double> Dsigma = step.attr("Dsigma_end").cast<py::array_t<double>>();
    sc.Dsigma_end = carma::arr_to_col(Dsigma);

    // Optional thermal fields (StepThermomeca has these)
    if (py::hasattr(step, "DT_end")) {
        sc.DT_end = step.attr("DT_end").cast<double>();
    } else {
        sc.DT_end = 0.0;
    }

    if (py::hasattr(step, "get_cBC_T")) {
        sc.cBC_T = step.attr("get_cBC_T")().cast<int>();
    } else {
        sc.cBC_T = 0;
    }

    return sc;
}

// Extract BlockConfig from Python Block object
simcoon::SolverEngine::BlockConfig extract_block(const py::object& block) {
    simcoon::SolverEngine::BlockConfig bc;

    bc.umat_name = block.attr("umat_name").cast<std::string>();

    py::array_t<double> props = block.attr("props").cast<py::array_t<double>>();
    bc.props = carma::arr_to_col(props);

    bc.nstatev = block.attr("nstatev").cast<int>();

    // Get control and corate types
    std::string control_type = block.attr("control_type").cast<std::string>();
    std::string corate_type = block.attr("corate_type").cast<std::string>();
    bc.control_type = get_control_type_int(control_type);
    bc.corate_type = get_corate_type_int(corate_type);

    bc.ncycle = block.attr("ncycle").cast<int>();

    // Extract steps
    py::list steps = block.attr("steps");
    for (auto step : steps) {
        bc.steps.push_back(extract_step(step.cast<py::object>()));
    }

    return bc;
}

py::list solver_optimized(
    const py::list& blocks_py,
    int max_iter,
    double tol,
    double lambda_solver
) {
    // Extract blocks from Python
    std::vector<simcoon::SolverEngine::BlockConfig> blocks;
    for (auto block : blocks_py) {
        blocks.push_back(extract_block(block.cast<py::object>()));
    }

    // Create solver and run
    simcoon::SolverEngine::SolverParams params;
    params.max_iter = max_iter;
    params.tol = tol;
    params.lambda_solver = lambda_solver;

    simcoon::SolverEngine solver(blocks, params);
    auto history = solver.solve();

    // Convert to Python HistoryPoint objects
    // Import the Python HistoryPoint class from simcoon.solver
    py::module_ solver_module = py::module_::import("simcoon.solver");
    py::object HistoryPoint = solver_module.attr("HistoryPoint");

    py::list result;
    for (const auto& hp : history) {
        result.append(HistoryPoint(
            carma::col_to_arr(hp.Etot),
            carma::col_to_arr(hp.sigma),
            carma::col_to_arr(hp.Wm),
            carma::col_to_arr(hp.statev),
            carma::mat_to_arr(hp.R),
            hp.T
        ));
    }

    return result;
}

} // namespace simpy
