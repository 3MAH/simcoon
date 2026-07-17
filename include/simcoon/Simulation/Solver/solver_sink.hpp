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

///@file solver_sink.hpp
///@brief In-memory drivable solver core: numeric controls, results sinks and solver_run()
///@version 1.0

#pragma once

#include <string>
#include <vector>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simcoon{

/** @addtogroup solver
 *  @{
 */

/**
 * @brief Numeric controls of the global Newton-Raphson / adaptive time-stepping loop.
 *
 * Defaults mirror the historical defaults of solver() (see solver.hpp).
 * Named solver_params because read.hpp already declares a solver_control() function.
 */
struct solver_params {
    double div_tnew_dt = 0.5;   ///< time-step division factor on non-convergence
    double mul_tnew_dt = 2.;    ///< time-step multiplication factor on fast convergence
    int miniter = 10;           ///< iteration count below which the step may grow
    int maxiter = 100;          ///< maximum Newton-Raphson iterations per increment
    int inforce = 1;            ///< 0: stop at Dn_mini, 1: enforce & carry residual, 2: enforce silently
    double precision = 1.E-6;   ///< convergence tolerance on the residual norm
    double lambda = 10000.;     ///< penalty stiffness for strain-controlled components
    int tangent_mode = tangent_default; ///< tangent_* constants (parameter.hpp)
};

/**
 * @brief Observer that receives the solver state at output points.
 *
 * solver_run() drives the cadence (solver_output frequency logic) and calls
 * record() once per output point; implementations decide what to do with the
 * converged state carried by the phase_characteristics.
 */
class solver_results_sink {
public:
    virtual ~solver_results_sink() = default;
    /// Called once at the very start of the simulation (start == true), after the
    /// first UMAT call, before the first increment.
    virtual void init(phase_characteristics &rve) { (void)rve; }
    /// Called at each output point. kblock/kcycle/kstep/kinc follow the same
    /// convention as phase_characteristics::output (0-based, written +1).
    virtual void record(phase_characteristics &rve, const solver_output &so,
                        const int &kblock, const int &kcycle, const int &kstep,
                        const int &kinc, const double &Time) = 0;
};

/**
 * @brief File sink reproducing the historical solver() output behaviour:
 * define_output() + phase_characteristics::output() on "global" and "local" streams.
 */
class solver_file_sink : public solver_results_sink {
public:
    solver_file_sink(const std::string &path_results, const std::string &outputfile_global,
                     const std::string &outputfile_local);
    void init(phase_characteristics &rve) override;
    void record(phase_characteristics &rve, const solver_output &so,
                const int &kblock, const int &kcycle, const int &kstep,
                const int &kinc, const double &Time) override;
private:
    std::string path_results;
    std::string outputfile_global;
    std::string outputfile_local;
};

/**
 * @brief Memory sink capturing the canonical converged state per output point.
 *
 * Stores raw state (all stress/strain measures carried by state_variables) so
 * consumers derive whatever output measure they need; Cauchy stress is formed
 * as tau/det(F1), consistent with phase_characteristics::output (o_stress_type 4).
 * Thermomechanical quantities (Wt, Q, r, dSdE/dSdT/drdE/drdT) are captured when
 * the global state variables are of type state_variables_T (sv_type == 2).
 */
class solver_memory_sink : public solver_results_sink {
public:
    bool record_tangent = true; ///< capture Lt (mechanical) or dSdE/dSdT/drdE/drdT (thermomechanical)
    int sv_type = 0;            ///< 1: mechanical, 2: thermomechanical (set by init())

    std::vector<int> blocks_i, cycles_i, steps_i, incs_i; ///< 0-based indices per record
    std::vector<double> time;
    std::vector<arma::vec> Etot, etot, PKII, tau, sigma;  ///< 6-component Voigt vectors
    std::vector<arma::mat> F1, R, DR;                     ///< 3x3 tensors
    std::vector<double> T, Q, r;                          ///< Q, r: thermomechanical only
    std::vector<arma::vec> Wm;                            ///< [Wm, Wm_r, Wm_ir, Wm_d]
    std::vector<arma::vec> Wt;                            ///< [Wt, Wt_r, Wt_ir], thermomechanical only
    std::vector<arma::vec> statev;
    std::vector<arma::mat> Lt;                            ///< 6x6, mechanical only
    std::vector<arma::mat> dSdE, dSdT, drdE, drdT;        ///< thermomechanical only

    void init(phase_characteristics &rve) override;
    void record(phase_characteristics &rve, const solver_output &so,
                const int &kblock, const int &kcycle, const int &kstep,
                const int &kinc, const double &Time) override;
    /// number of records captured so far
    unsigned int n_records() const { return static_cast<unsigned int>(time.size()); }
};

/**
 * @brief In-memory core of the solver: runs the block/cycle/step/increment loops
 * on already-built loading blocks and streams results to a sink.
 *
 * This is the single numerical engine; solver() (solver.hpp) is its file-driven
 * wrapper (read_path/read_output in, solver_file_sink out).
 *
 * @param blocks Loading blocks with fully-defined steps (mutated during the run:
 *        step generation and the inforce residual carry-over write into the steps)
 * @param T_init Initial temperature
 * @param so Output configuration (only the per-block frequency logic and the
 *        statev/tangent selection are used by sinks)
 * @param umat_name Constitutive model name (5 characters)
 * @param props Material properties
 * @param nstatev Number of internal state variables
 * @param psi_rve,theta_rve,phi_rve Euler angles of the RVE orientation (rad)
 * @param solver_type 0: classic Newton-Raphson, 1: RNL (control_type 1 only)
 * @param corate_type Objective rate choice (see corate_kinematics, objective_rates.hpp)
 * @param ctrl Numeric solver controls
 * @param sink Results observer
 * @return 0 on completion, 1 on early abort (invalid solver type / unrecognized
 *         thermal BC / non-convergence with inforce == 0)
 */
int solver_run(std::vector<block> &blocks, const double &T_init, const solver_output &so,
               const std::string &umat_name, const arma::vec &props, const unsigned int &nstatev,
               const double &psi_rve, const double &theta_rve, const double &phi_rve,
               const int &solver_type, const int &corate_type,
               const solver_params &ctrl, solver_results_sink &sink);

/** @} */ // end of solver group

} //namespace simcoon
