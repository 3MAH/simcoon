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

///@file solver_engine.hpp
///@brief Optimized C++ solver engine with pre-allocated buffers
///@version 1.0

#pragma once

#include <string>
#include <vector>
#include <armadillo>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>

namespace simcoon {

// Control type constants (matching Python CONTROL_TYPES)
constexpr int CTRL_SMALL_STRAIN = 1;
constexpr int CTRL_GREEN_LAGRANGE = 2;
constexpr int CTRL_LOGARITHMIC = 3;
constexpr int CTRL_BIOT = 4;
constexpr int CTRL_F = 5;
constexpr int CTRL_GRADU = 6;

// Corate type constants (matching Python CORATE_TYPES)
constexpr int CORATE_JAUMANN = 0;
constexpr int CORATE_GREEN_NAGHDI = 1;
constexpr int CORATE_LOGARITHMIC = 2;
constexpr int CORATE_LOGARITHMIC_R = 3;
constexpr int CORATE_TRUESDELL = 4;
constexpr int CORATE_LOGARITHMIC_F = 5;

// Forward declaration
class UmatDispatch;

/**
 * @brief History point storing state at each increment.
 *
 * Mirrors the Python HistoryPoint dataclass for compatible output format.
 */
struct HistoryPointCpp {
    arma::vec Etot;      ///< Total strain (6)
    arma::vec sigma;     ///< Cauchy stress (6)
    arma::vec Wm;        ///< Work measures (4): [Wm, Wm_r, Wm_ir, Wm_d]
    arma::vec statev;    ///< State variables (nstatev)
    arma::mat R;         ///< Rotation matrix (3,3)
    double T;            ///< Temperature

    HistoryPointCpp() : Etot(6, arma::fill::zeros),
                        sigma(6, arma::fill::zeros),
                        Wm(4, arma::fill::zeros),
                        statev(),
                        R(3, 3, arma::fill::eye),
                        T(293.15) {}
};

/**
 * @brief Optimized solver engine with pre-allocated Newton-Raphson buffers.
 *
 * This class implements the same algorithm as the Python Solver but with:
 * - Static UMAT dispatch (via UmatDispatch singleton)
 * - Pre-allocated Newton-Raphson buffers (K, residual, Delta)
 * - No Python interpreter overhead
 *
 * Accepts the same Block/Step structure as Python (extracted by bindings).
 */
class SolverEngine {
public:
    /**
     * @brief Step configuration extracted from Python Step object.
     */
    struct StepConfig {
        int Dn_init = 100;        ///< Initial increment count
        int Dn_mini = 10;         ///< Minimum increment count
        int Dn_inc = 200;         ///< Maximum increment count
        double time = 1.0;        ///< Step time
        arma::Col<int> cBC_meca;  ///< Boundary conditions: 0=strain, 1=stress
        arma::vec DEtot_end;      ///< Target strain increment
        arma::vec Dsigma_end;     ///< Target stress increment (for stress control)
        double DT_end = 0.0;      ///< Temperature increment
        int cBC_T = 0;            ///< Temperature BC: 0=prescribed, 1=adiabatic
    };

    /**
     * @brief Block configuration extracted from Python Block object.
     */
    struct BlockConfig {
        std::vector<StepConfig> steps;
        std::string umat_name;
        arma::vec props;
        int nstatev = 1;
        int control_type = 0;     ///< 0=small_strain, 1=finite_strain, 2=logarithmic
        int corate_type = 0;      ///< 0=none, 1=jaumann, 2=green_naghdi
        int ncycle = 1;           ///< Number of cycles
    };

    /**
     * @brief Solver parameters.
     */
    struct SolverParams {
        int max_iter;              ///< Max Newton-Raphson iterations
        double tol;                ///< Convergence tolerance
        double lambda_solver;      ///< Penalty stiffness for strain control

        SolverParams() : max_iter(10), tol(1e-9), lambda_solver(10000.0) {}
        SolverParams(int mi, double t, double ls) : max_iter(mi), tol(t), lambda_solver(ls) {}
    };

    /**
     * @brief Construct solver engine with blocks and parameters.
     *
     * @param blocks Vector of block configurations
     * @param params Solver parameters
     */
    SolverEngine(const std::vector<BlockConfig>& blocks, const SolverParams& params = SolverParams());

    /**
     * @brief Run the simulation.
     *
     * @return Vector of history points (state at each increment)
     */
    std::vector<HistoryPointCpp> solve();

private:
    std::vector<BlockConfig> blocks_;
    SolverParams params_;

    // Pre-allocated Newton-Raphson buffers (reused across ALL increments)
    arma::mat K_;           ///< (6,6) Jacobian matrix
    arma::vec residual_;    ///< (6) Residual vector
    arma::vec Delta_;       ///< (6) Increment vector

    // Working state variables
    state_variables_M sv_;

    // Backup state for rollback on failed increments (legacy C++ pattern)
    // Stores the complete state before attempting an increment
    state_variables_M sv_start_;

    /**
     * @brief Initialize UMAT for a block.
     *
     * Calls UMAT with start=true to initialize state variables.
     */
    void initialize_umat(const BlockConfig& block, double Time);

    /**
     * @brief Solve a single step.
     *
     * @param block Block configuration
     * @param step Step configuration
     * @param Time Current time (updated)
     * @param history Output history vector
     * @return New time after step
     */
    double solve_step(const BlockConfig& block, const StepConfig& step,
                      double Time, std::vector<HistoryPointCpp>& history);

    /**
     * @brief Solve a single increment using Newton-Raphson.
     *
     * @param block Block configuration
     * @param Time Current time
     * @param DTime Time increment
     * @param Dtinc Fraction of increment
     * @param DEtot_target Target strain increment
     * @param Dsigma_target Target stress increment
     * @param cBC_meca Boundary conditions
     * @param nK Increment counter for convergence tracking
     * @return true if converged, false otherwise
     */
    bool solve_increment(const BlockConfig& block,
                         double Time, double DTime, double Dtinc,
                         const arma::vec& DEtot_target, const arma::vec& Dsigma_target,
                         const arma::Col<int>& cBC_meca, int nK);

    /**
     * @brief Compute the residual for Newton-Raphson.
     *
     * @param Dtinc Fraction of step increment
     * @param DEtot_target Target strain increment
     * @param Dsigma_target Target stress increment
     * @param cBC_meca Boundary conditions
     * @param control_type Control type (for selecting stress/strain measures)
     */
    void compute_residual(double Dtinc,
                          const arma::vec& DEtot_target, const arma::vec& Dsigma_target,
                          const arma::Col<int>& cBC_meca, int control_type = CTRL_SMALL_STRAIN);

    /**
     * @brief Build the Jacobian matrix for Newton-Raphson.
     */
    void build_jacobian(const arma::Col<int>& cBC_meca);

    /**
     * @brief Call the UMAT via dispatch singleton.
     */
    void call_umat(const BlockConfig& block, double Time, double DTime, bool start);

    /**
     * @brief Record current state to history.
     */
    void record_history(std::vector<HistoryPointCpp>& history);

    /**
     * @brief Save current state to backup buffers for potential rollback.
     */
    void save_state_for_rollback();

    /**
     * @brief Restore state from backup buffers after failed increment.
     */
    void restore_state_from_backup();

    /**
     * @brief Update kinematic quantities for finite strain.
     *
     * Computes F0, F1, U0, U1, DR from strain and rotation for objective rates.
     *
     * @param control_type Control type (1=small_strain, 2=green_lagrange, 3=logarithmic, etc.)
     * @param corate_type Corotational rate type
     * @param DTime Time increment
     */
    void update_kinematics(int control_type, int corate_type, double DTime);
};

} // namespace simcoon
