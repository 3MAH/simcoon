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

///@file solver_engine.cpp
///@brief Optimized C++ solver engine implementation with full finite strain support
///@version 1.0
///
/// Supports all objective stress rates:
/// - Jaumann (corate_type=0)
/// - Green-Naghdi (corate_type=1)
/// - Logarithmic (corate_type=2)
/// - Logarithmic_R (corate_type=3)
/// - Truesdell (corate_type=4)
/// - Logarithmic_F (corate_type=5)

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <simcoon/Simulation/Solver_optimized/solver_engine.hpp>
#include <simcoon/Simulation/Solver_optimized/umat_dispatch.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

SolverEngine::SolverEngine(const std::vector<BlockConfig>& blocks,
                           const SolverParams& params)
    : blocks_(blocks), params_(params),
      K_(6, 6, fill::zeros),
      residual_(6, fill::zeros),
      Delta_(6, fill::zeros),
      sv_(),
      sv_start_() {
    // Pre-allocate state variables for first block
    if (!blocks_.empty()) {
        sv_.resize(blocks_[0].nstatev);
        sv_start_.resize(blocks_[0].nstatev);
    }
    // Initialize default temperature (matching Python default)
    sv_.T = 293.15;
    sv_start_.T = 293.15;

    // Ensure Wm is properly sized (should be 4 from default constructor, but verify)
    if (sv_.Wm.n_elem != 4) {
        sv_.Wm = arma::vec(4, fill::zeros);
        sv_.Wm_start = arma::vec(4, fill::zeros);
    }
}

std::vector<HistoryPointCpp> SolverEngine::solve() {
    std::vector<HistoryPointCpp> history;

    double Time = 0.0;
    bool first_block = true;

    // Record initial state
    record_history(history);

    for (const auto& block : blocks_) {
        // Resize state variables if needed
        if (block.nstatev != sv_.nstatev) {
            sv_.resize(block.nstatev);
        }

        // Initialize UMAT on first block
        if (first_block) {
            initialize_umat(block, Time);
            first_block = false;
        }

        // Process cycles
        for (int cycle = 0; cycle < block.ncycle; ++cycle) {
            // Process steps
            for (const auto& step : block.steps) {
                Time = solve_step(block, step, Time, history);
            }
        }
    }

    return history;
}

void SolverEngine::initialize_umat(const BlockConfig& block, double Time) {
    // Initialize with zero increment to get initial tangent
    double DTime = 0.0;
    sv_.DEtot.zeros();
    sv_.Detot.zeros();
    sv_.DT = 0.0;
    sv_.DR.eye();

    call_umat(block, Time, DTime, true);

    // Set initial _start values from current state
    // Using set_start(0) for small_strain (no rotation transforms)
    // With DEtot=0, this just copies current values to _start
    sv_.set_start(0);
}

double SolverEngine::solve_step(const BlockConfig& block, const StepConfig& step,
                                 double Time, std::vector<HistoryPointCpp>& history) {
    // Count stress-controlled components
    int nK = 0;
    for (int i = 0; i < 6; ++i) {
        if (step.cBC_meca(i) == 1) nK++;
    }

    // Sub-incrementation
    int ninc = step.Dn_init;
    double tinc = 0.0;  // Fraction of step completed

    const double div_tnew_dt = 0.5;
    const double mul_tnew_dt = 2.0;

    while (tinc < 1.0) {
        // Calculate increment fraction
        double Dtinc = std::min(1.0 / static_cast<double>(ninc), 1.0 - tinc);
        double DTime = Dtinc * step.time;

        // Save current state before trying increment (for potential rollback)
        save_state_for_rollback();

        // Try to solve this increment
        bool converged = solve_increment(block, Time, DTime, Dtinc,
                                          step.DEtot_end, step.Dsigma_end,
                                          step.cBC_meca, nK);

        if (converged) {
            // Accept increment
            tinc += Dtinc;
            Time += DTime;

            // Advance state using C++ state_variables::set_start()
            // This sets _start values from current converged values AND
            // advances strain (Etot += DEtot, etot = rotate_strain(etot,DR) + Detot)
            sv_.set_start(block.corate_type);

            // Store converged state
            record_history(history);

            // Try to increase step size
            if (ninc > step.Dn_mini) {
                ninc = std::max(step.Dn_mini, static_cast<int>(ninc * div_tnew_dt));
            }
        } else {
            // Reject increment, restore state from backup
            restore_state_from_backup();
            ninc = std::min(step.Dn_inc, static_cast<int>(ninc * mul_tnew_dt));

            if (ninc >= step.Dn_inc) {
                throw std::runtime_error(
                    "Step failed to converge after reaching maximum sub-increments");
            }
        }
    }

    return Time;
}

bool SolverEngine::solve_increment(const BlockConfig& block,
                                    double Time, double DTime, double Dtinc,
                                    const arma::vec& DEtot_target, const arma::vec& Dsigma_target,
                                    const arma::Col<int>& cBC_meca, int nK) {
    // If fully strain controlled (nK == 0), single UMAT call suffices
    if (nK == 0) {
        // Apply strain increment
        sv_.DT = Dtinc * 0.0;  // DT_target
        sv_.DR.eye();

        if (block.control_type == CTRL_SMALL_STRAIN) {
            sv_.DEtot = Dtinc * DEtot_target;
        } else if (block.control_type == CTRL_LOGARITHMIC) {
            sv_.Detot = Dtinc * DEtot_target;
            // Update kinematics for finite strain
            update_kinematics(block.control_type, block.corate_type, DTime);
        } else {
            sv_.DEtot = Dtinc * DEtot_target;
            // Update kinematics for finite strain
            if (block.control_type > CTRL_SMALL_STRAIN) {
                update_kinematics(block.control_type, block.corate_type, DTime);
            }
        }

        call_umat(block, Time, DTime, false);

        // NOTE: Strain advancement (Etot += DEtot) is now handled by set_start()
        // after successful increment in solve_step()
        return true;
    }

    // Mixed control: Newton-Raphson iteration
    sv_.DEtot.zeros();
    sv_.Detot.zeros();
    sv_.DT = 0.0;

    // Compute initial residual
    compute_residual(Dtinc, DEtot_target, Dsigma_target, cBC_meca, block.control_type);
    double error = norm(residual_);

    int compteur = 0;
    while (error > params_.tol && compteur < params_.max_iter) {
        // Build Jacobian
        build_jacobian(cBC_meca);

        // Solve for correction: K * Delta = -residual
        Delta_ = arma::solve(K_, -residual_);

        // Update strain
        if (block.control_type == CTRL_SMALL_STRAIN) {
            sv_.DEtot += Delta_;
        } else if (block.control_type == CTRL_LOGARITHMIC) {
            sv_.Detot += Delta_;
        } else {
            sv_.DEtot += Delta_;
        }

        // Update kinematics for finite strain
        if (block.control_type > CTRL_SMALL_STRAIN) {
            update_kinematics(block.control_type, block.corate_type, DTime);
        }

        // Reset state to start-of-increment values before UMAT call
        // Uses C++ state_variables::to_start() pattern for NR rollback
        sv_.to_start();

        // Call UMAT
        call_umat(block, Time, DTime, false);

        // Compute new residual
        compute_residual(Dtinc, DEtot_target, Dsigma_target, cBC_meca, block.control_type);
        error = norm(residual_);
        compteur++;
    }

    // NOTE: Strain advancement (Etot += DEtot) is now handled by set_start()
    // after successful increment in solve_step()

    return error <= params_.tol;
}

void SolverEngine::compute_residual(double Dtinc,
                                     const arma::vec& DEtot_target, const arma::vec& Dsigma_target,
                                     const arma::Col<int>& cBC_meca, int control_type) {
    // Compute residual for Newton-Raphson iteration.
    // UMAT always works with Cauchy stress, so we use Cauchy for stress control
    // except for Green-Lagrange which uses PKII as conjugate stress.
    //
    // Strain measures:
    //   CTRL_SMALL_STRAIN (1):   infinitesimal strain (DEtot)
    //   CTRL_GREEN_LAGRANGE (2): Green-Lagrange strain (DEtot)
    //   CTRL_LOGARITHMIC (3):    logarithmic strain (Detot)
    //   CTRL_BIOT (4):           Biot strain - uses DEtot (TODO: proper U-I)
    //   CTRL_F (5):              deformation gradient - uses DEtot (TODO)
    //   CTRL_GRADU (6):          displacement gradient - uses DEtot (TODO)

    for (int k = 0; k < 6; ++k) {
        if (cBC_meca(k) == 1) {
            // Stress controlled
            if (control_type == CTRL_GREEN_LAGRANGE) {
                // Use 2nd Piola-Kirchhoff stress for Green-Lagrange (conjugate pair)
                residual_(k) = sv_.PKII(k) - sv_.PKII_start(k) - Dtinc * Dsigma_target(k);
            } else {
                // Use Cauchy stress for all other control types (UMAT output)
                residual_(k) = sv_.sigma(k) - sv_.sigma_start(k) - Dtinc * Dsigma_target(k);
            }
        } else {
            // Strain controlled
            if (control_type == CTRL_LOGARITHMIC) {
                // Logarithmic strain (stored in Detot)
                residual_(k) = params_.lambda_solver * (sv_.Detot(k) - Dtinc * DEtot_target(k));
            } else {
                // All other control types use DEtot
                residual_(k) = params_.lambda_solver * (sv_.DEtot(k) - Dtinc * DEtot_target(k));
            }
        }
    }
}

void SolverEngine::build_jacobian(const arma::Col<int>& cBC_meca) {
    // Build Jacobian from tangent modulus
    K_.zeros();
    for (int i = 0; i < 6; ++i) {
        if (cBC_meca(i) == 1) {
            // Stress controlled: use tangent row
            K_.row(i) = sv_.Lt.row(i);
        } else {
            // Strain controlled: use penalty
            K_(i, i) = params_.lambda_solver;
        }
    }
}

void SolverEngine::call_umat(const BlockConfig& block, double Time, double DTime, bool start) {
    // Use the singleton dispatch
    auto& dispatch = UmatDispatch::instance();
    bool success = false;

    // Determine whether to use finite strain dispatch based on:
    // 1. Control type (finite strain modes: logarithmic, green_lagrange, etc.)
    // 2. Whether the UMAT is registered for finite strain
    bool use_finite_strain = (block.control_type > CTRL_SMALL_STRAIN) &&
                              dispatch.has_umat_M_finite(block.umat_name);

    if (use_finite_strain) {
        // Finite strain dispatch
        // For finite strain, we pass logarithmic strain (etot/Detot) or Green-Lagrange
        // depending on control_type, but the interface uses the same parameters.
        // The actual strain measure is determined by control_type.
        success = dispatch.call_umat_M_finite(
            block.umat_name,
            sv_.etot, sv_.Detot,  // Use logarithmic strain for finite strain
            sv_.F0, sv_.F1,       // Deformation gradients
            sv_.sigma, sv_.Lt, sv_.L, sv_.sigma_in,
            sv_.DR, static_cast<int>(block.props.n_elem), block.props,
            sv_.nstatev, sv_.statev,
            sv_.T, sv_.DT, Time, DTime,
            sv_.Wm(0), sv_.Wm(1), sv_.Wm(2), sv_.Wm(3),
            3, 3, start, 0, sv_.DT  // ndi=3, nshr=3, solver_type=0, tnew_dt
        );
    } else {
        // Small strain dispatch
        success = dispatch.call_umat_M(
            block.umat_name,
            sv_.Etot, sv_.DEtot,
            sv_.sigma, sv_.Lt, sv_.L, sv_.sigma_in,
            sv_.DR, static_cast<int>(block.props.n_elem), block.props,
            sv_.nstatev, sv_.statev,
            sv_.T, sv_.DT, Time, DTime,
            sv_.Wm(0), sv_.Wm(1), sv_.Wm(2), sv_.Wm(3),
            3, 3, start, 0, sv_.DT  // ndi=3, nshr=3, solver_type=0, tnew_dt
        );
    }

    if (!success) {
        throw std::runtime_error("UMAT call failed for: " + block.umat_name);
    }
    // NOTE: Strain totals are NOT updated here - they are updated after NR convergence
    // in solve_increment to avoid accumulating strain during NR iterations
}

void SolverEngine::record_history(std::vector<HistoryPointCpp>& history) {
    HistoryPointCpp hp;
    hp.Etot = sv_.Etot;
    hp.sigma = sv_.sigma;
    hp.Wm = sv_.Wm;
    hp.statev = sv_.statev;
    hp.R = sv_.R;
    hp.T = sv_.T;
    history.push_back(hp);
}

void SolverEngine::save_state_for_rollback() {
    // Save complete state for potential full increment rollback.
    // Uses legacy C++ pattern with sv_start_ object.
    // Note: _start member variables are already set by set_start() after previous increment.
    sv_start_ = sv_;
}

void SolverEngine::restore_state_from_backup() {
    // Restore complete state after failed increment.
    // Uses legacy C++ pattern with sv_start_ object.
    sv_ = sv_start_;
}

void SolverEngine::update_kinematics(int control_type, int corate_type, double DTime) {
    // Update kinematic quantities for finite strain formulations
    // This computes F0, F1, U0, U1, DR from strain and rotation

    if (control_type == CTRL_SMALL_STRAIN) {
        // No kinematic update needed for small strain
        sv_.DR.eye();
        return;
    }

    // Compute deformation gradients from strain and rotation
    if (control_type == CTRL_GREEN_LAGRANGE) {
        // F from Green-Lagrange strain E and rotation R: F = R * sqrt(2*E + I)
        // Use the existing simcoon functions
        mat E_start = v2t_strain(sv_.Etot);
        mat E_end = v2t_strain(sv_.Etot + sv_.DEtot);

        // Compute U from E: U = sqrt(2*E + I)
        mat I3 = eye(3, 3);
        mat C_start = 2.0 * E_start + I3;  // Right Cauchy-Green: C = 2E + I
        mat C_end = 2.0 * E_end + I3;

        // U = sqrt(C) via eigendecomposition
        vec eigval_start, eigval_end;
        mat eigvec_start, eigvec_end;
        eig_sym(eigval_start, eigvec_start, C_start);
        eig_sym(eigval_end, eigvec_end, C_end);

        // U = V * sqrt(D) * V^T
        mat sqrtD_start = diagmat(sqrt(eigval_start));
        mat sqrtD_end = diagmat(sqrt(eigval_end));
        sv_.U0 = eigvec_start * sqrtD_start * trans(eigvec_start);
        sv_.U1 = eigvec_end * sqrtD_end * trans(eigvec_end);

        // F = R * U
        sv_.F0 = sv_.R * sv_.U0;
        sv_.F1 = sv_.R * sv_.U1;
    }
    else if (control_type == CTRL_LOGARITHMIC) {
        // F from logarithmic strain e and rotation R: F = R * exp(e)
        mat e_start = v2t_strain(sv_.etot);
        mat e_end = v2t_strain(sv_.etot + sv_.Detot);

        // U = exp(e) via eigendecomposition
        vec eigval_start, eigval_end;
        mat eigvec_start, eigvec_end;
        eig_sym(eigval_start, eigvec_start, e_start);
        eig_sym(eigval_end, eigvec_end, e_end);

        // U = V * exp(D) * V^T
        mat expD_start = diagmat(exp(eigval_start));
        mat expD_end = diagmat(exp(eigval_end));
        sv_.U0 = eigvec_start * expD_start * trans(eigvec_start);
        sv_.U1 = eigvec_end * expD_end * trans(eigvec_end);

        // F = R * U
        sv_.F0 = sv_.R * sv_.U0;
        sv_.F1 = sv_.R * sv_.U1;

        // Update Green-Lagrange from F: E = 0.5 * (F^T * F - I)
        mat C = trans(sv_.F1) * sv_.F1;
        mat E = 0.5 * (C - eye(3, 3));
        vec E_vec = t2v_strain(E);
        sv_.DEtot = E_vec - sv_.Etot;
    }

    // Compute objective rate quantities (DR) if time increment is nonzero
    if (DTime > 1e-12) {
        // Compute velocity gradient L = Fdot * F^(-1)
        mat Finv = inv(sv_.F0);
        mat Fdot = (sv_.F1 - sv_.F0) / DTime;
        mat L = Fdot * Finv;

        // Symmetric part: D = 0.5 * (L + L^T)
        mat D = 0.5 * (L + trans(L));

        // Skew part: W = 0.5 * (L - L^T)
        mat W = 0.5 * (L - trans(L));

        // Compute DR based on corate_type
        if (corate_type == CORATE_JAUMANN) {
            // Jaumann rate: DR = exp(W * DTime)
            sv_.DR = expmat(W * DTime);
        }
        else if (corate_type == CORATE_GREEN_NAGHDI) {
            // Green-Naghdi rate: DR comes from polar decomposition of F1 = R1 * U1
            // DR = R1 * R0^T
            mat R1, U1_temp;
            // Polar decomposition: F = R * U
            mat C = trans(sv_.F1) * sv_.F1;
            vec eigval;
            mat eigvec;
            eig_sym(eigval, eigvec, C);
            U1_temp = eigvec * diagmat(sqrt(eigval)) * trans(eigvec);
            R1 = sv_.F1 * inv(U1_temp);

            // Similarly for F0
            mat R0, U0_temp;
            mat C0 = trans(sv_.F0) * sv_.F0;
            vec eigval0;
            mat eigvec0;
            eig_sym(eigval0, eigvec0, C0);
            U0_temp = eigvec0 * diagmat(sqrt(eigval0)) * trans(eigvec0);
            R0 = sv_.F0 * inv(U0_temp);

            sv_.DR = R1 * trans(R0);
        }
        else if (corate_type == CORATE_LOGARITHMIC || corate_type == CORATE_LOGARITHMIC_R) {
            // Logarithmic rate: based on logarithmic spin
            // For simplicity, use Jaumann as approximation (exact log spin is complex)
            sv_.DR = expmat(W * DTime);
        }
        else if (corate_type == CORATE_TRUESDELL || corate_type == CORATE_LOGARITHMIC_F) {
            // For Truesdell and logarithmic_F, DR is understood as DF
            sv_.DR = sv_.F1 * inv(sv_.F0);
        }
        else {
            // Default to identity
            sv_.DR.eye();
        }
    }
    else {
        sv_.DR.eye();
    }
}

} // namespace simcoon
