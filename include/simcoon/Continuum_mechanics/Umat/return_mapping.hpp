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

/**
 * @file return_mapping.hpp
 * @brief Closest-point projection (CPP) return mapping for the lead-mechanism framework
 * (doc §"Closest-point projection within the lead-mechanism framework").
 *
 * Solves the fully implicit local system
 * \f[
 *   \b{R}_\sigma = \boldsymbol{\sigma} - \boldsymbol{\sigma}^{tr}
 *     + \sum_j \Delta s^j\,\mathbf{L}:\boldsymbol{\Lambda}_\varepsilon^j(\boldsymbol{\sigma},\mathbf{V}) = \mathbf{0},
 *   \qquad
 *   \Delta s^j \ge 0,\; \Phi^j \le 0,\; \Delta s^j\,\Phi^j = 0,
 * \f]
 * with the internal state \f$ \mathbf{V} \f$ resolved by backward Euler at every iterate
 * (inner-consistent scheme), via a condensed semi-smooth Newton iteration on the reduced
 * \f$ N \times N \f$ multiplier system handed to Fischer_Burmeister_m(). The condensation
 * operator \f$ \mathbf{M} = \mathbf{I} + \mathbf{L}:\sum_j \Delta s^j\,
 * \partial\boldsymbol{\Lambda}_\varepsilon^j/\partial\boldsymbol{\sigma} \f$ is the same
 * operator as in assemble_algorithmic_tangent(), so at convergence the exact consistent
 * tangent is available from the returned ingredients at no extra cost
 * (see cpp_consistent_tangent()).
 *
 * Unlike the legacy convex-cutting-plane (CCP) loops, the inelastic strain uses the flow at
 * the CONVERGED state: \f$ \boldsymbol{\varepsilon}^{in} = \boldsymbol{\varepsilon}^{in}_n
 * + \sum_j \Delta s^j \boldsymbol{\Lambda}_\varepsilon^j(\boldsymbol{\sigma}_{n+1},
 * \mathbf{V}_{n+1}) \f$. For radial flows (von Mises) the two integrators coincide; for
 * anisotropic criteria the converged stresses differ by \f$ O(\|\Delta\varepsilon\|^2) \f$.
 *
 * @note Modular-framework alignment: the callback design maps one-to-one onto the modular
 * StrainMechanism interface (PR feature/modular) — compute_constraints -> Phi,
 * compute_flow_directions -> Lambda / dPhi_dsigma, compute_jacobian_contribution -> K,
 * update -> update_state — so a future ClosestPointIntegrator is a thin adapter over
 * closest_point_return_mapping().
 */

#pragma once
#include <armadillo>
#include <functional>
#include <vector>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>

namespace simcoon {

/**
 * @brief One dissipative mechanism, described by callbacks evaluated at the CURRENT iterate.
 *
 * The internal state \f$ \mathbf{V} \f$ is OWNED BY THE CALLER (captured by reference in the
 * lambdas) and is refreshed through ReturnStateHooks::update_state before any of these
 * callbacks is invoked, so each callback only needs the stress argument.
 */
struct ReturnMechanism {
    /// REQUIRED. Criterion value \f$ \Phi^j(\boldsymbol{\sigma}, \mathbf{V}) \f$.
    std::function<double(const arma::vec &sigma)> Phi;
    /// REQUIRED. Criterion gradient \f$ \partial\Phi^j/\partial\boldsymbol{\sigma} \f$ (6).
    std::function<arma::vec(const arma::vec &sigma)> dPhi_dsigma;
    /// OPTIONAL. Strain-flow direction \f$ \boldsymbol{\Lambda}_\varepsilon^j \f$ (6).
    /// Empty => associated flow (Lambda = dPhi_dsigma).
    std::function<arma::vec(const arma::vec &sigma)> Lambda;
    /// REQUIRED. TOTAL derivative \f$ \mathrm{d}\boldsymbol{\Lambda}_\varepsilon^j/
    /// \mathrm{d}\boldsymbol{\sigma} \f$ (6x6, compliance-like Voigt type): analytic Hessian
    /// (deta_stress / ddHill_stress / ddDFA_stress / ddAni_stress) for quadratic criteria, or a
    /// central finite difference of the inner-consistent map
    /// \f$ \boldsymbol{\sigma}' \mapsto \boldsymbol{\Lambda}(\boldsymbol{\sigma}',
    /// \hat{\mathbf{V}}(\boldsymbol{\sigma}', \Delta s)) \f$ for state-coupled mechanisms.
    std::function<arma::mat(const arma::vec &sigma)> dLambda_dsigma;
};

/**
 * @brief Caller-level hooks coupling the mechanisms through the internal state.
 */
struct ReturnStateHooks {
    /// Backward-Euler state refresh: solve
    /// \f$ \mathbf{V} = \mathbf{V}_n + \sum_j \Delta s^j\,\boldsymbol{\Lambda}_V^j \f$ and leave
    /// the caller-captured state consistent (closed forms preferred; short fixed point otherwise).
    /// May be empty for state-free (perfect-plasticity-like) problems. Return false on inner
    /// non-convergence (aborts the outer iteration with converged = false).
    std::function<bool(const arma::vec &sigma, const arma::vec &Dlambda)> update_state;
    /// REQUIRED. Hardening/state block \f$ \mathbf{K}^{lj} = \partial\Phi^l/\partial\mathbf{V}
    /// \cdot \partial\mathbf{V}/\partial\Delta s^j \f$ (NxN) at the current iterate — the same
    /// rows the CCP loops assemble.
    std::function<arma::mat(const arma::vec &sigma, const arma::vec &Dlambda)> K;
    /// OPTIONAL multiplier-side flow/state chain, N vectors of 6:
    /// \f$ \mathbf{c}^j = \sum_q \Delta s^q\,\mathbf{L}:\partial\boldsymbol{\Lambda}^q/
    /// \partial\mathbf{V}\cdot\partial\hat{\mathbf{V}}/\partial\Delta s^j \f$. Empty => omitted
    /// (superlinear instead of quadratic for state-coupled mechanisms; converged solution
    /// unaffected).
    std::function<std::vector<arma::vec>(const arma::vec &sigma, const arma::vec &Dlambda)> flow_state_coupling;
};

/// Iteration controls; zero-valued members fall back to the simcoon defaults.
struct ReturnMappingControl {
    int    maxiter = 0;             ///< 0 => simcoon::maxiter_umat
    double precision = 0.;          ///< 0 => simcoon::precision_umat
    int    max_backtrack = 5;       ///< halvings when the combined error grows > x2
    double max_dsigma_frac = 0.5;   ///< Newton-step cap as a fraction of sigma_ref
};

/**
 * @brief Converged CPP state plus the exact ingredients of the consistent tangent.
 *
 * Bhat_continuum / kappa_j / dPhidsigma_l / dLambda_dsigma_l are evaluated at the CONVERGED
 * state and are exactly the arguments of assemble_algorithmic_tangent() — see
 * cpp_consistent_tangent().
 */
struct ReturnMappingResult {
    arma::vec sigma;                          ///< converged stress (6)
    arma::vec Dlambda;                        ///< converged multipliers \f$ \Delta s^j \f$ (N)
    bool      converged = false;
    int       niter = 0;
    double    error = 0.;                     ///< final combined error (FB + \f$ \|R_\sigma\|/\sigma_{ref} \f$)
    std::vector<double> error_history;        ///< combined error per iteration (rate diagnostics)
    arma::mat Bhat_continuum;                 ///< \f$ \hat{B}^{lj} = n^l:\kappa^j - K^{lj} \f$ (NxN)
    std::vector<arma::vec> kappa_j;           ///< \f$ \mathbf{L}:\boldsymbol{\Lambda}^j \f$ (each 6)
    std::vector<arma::vec> dPhidsigma_l;      ///< \f$ n^l \f$ (each 6)
    std::vector<arma::mat> dLambda_dsigma_l;  ///< \f$ D^j \f$ (each 6x6)
};

/**
 * @brief Coupled closest-point projection return mapping (N mechanisms).
 *
 * @param sigma_tr Elastic trial stress (6), built by the caller with el_pred().
 * @param L        Elastic stiffness (6x6).
 * @param mechanisms Mechanism callbacks (state refreshed via hooks before each evaluation).
 * @param hooks    State hooks (update_state / K / optional flow_state_coupling).
 * @param Y_crit   Per-mechanism normalisation (N), as in the CCP loops (typically sigmaY).
 * @param control  Iteration controls.
 * @return ReturnMappingResult (converged flag; caller signals non-convergence with the
 *         standard tnew_dt < 1 step-cut — no silent CCP fallback).
 *
 * @details ndi == 3 (full 3D) in this version; callers keep the CCP loop for condensed
 * (plane-stress/1D) states. If all \f$ \Phi^j(\sigma_{tr}, V_n) \le 0 \f$ the step is elastic
 * and the trial state is returned immediately (bit-identical elasticity across integrators).
 */
ReturnMappingResult closest_point_return_mapping(
    const arma::vec &sigma_tr,
    const arma::mat &L,
    const std::vector<ReturnMechanism> &mechanisms,
    const ReturnStateHooks &hooks,
    const arma::vec &Y_crit,
    const ReturnMappingControl &control = {});

/**
 * @brief Single-mechanism convenience overload (mirrors the tangent_assembly overloads).
 *
 * @param flow_state_coupling OPTIONAL multiplier-side state chain
 * \f$ \mathbf{c} = \Delta s\,\mathbf{L}:\partial\boldsymbol{\Lambda}/\partial\mathbf{V}
 * \cdot\partial\hat{\mathbf{V}}/\partial\Delta s \f$ (6). Required for the consistent tangent
 * to be exact when the flow depends on \f$ \Delta s \f$ through the state (e.g. backstress);
 * empty otherwise.
 */
ReturnMappingResult closest_point_return_mapping(
    const arma::vec &sigma_tr,
    const arma::mat &L,
    const ReturnMechanism &mechanism,
    const std::function<bool(const arma::vec &sigma, double Dlambda)> &update_state,
    const std::function<double(const arma::vec &sigma, double Dlambda)> &K_scalar,
    double Y_crit,
    const std::function<arma::vec(const arma::vec &sigma, double Dlambda)> &flow_state_coupling = {},
    const ReturnMappingControl &control = {});

/**
 * @brief Exact consistent tangent of the converged CPP map.
 *
 * Forwards to assemble_algorithmic_tangent(r.Bhat_continuum, r.kappa_j, r.dPhidsigma_l,
 * r.Dlambda, L, r.dLambda_dsigma_l). Because the CPP map IS the implicit update the
 * algorithmic tangent linearises, and dLambda_dsigma carries the total (state-chained)
 * derivative, the result is the exact symmetric Jacobian of the discrete update.
 */
ContinuumTangent cpp_consistent_tangent(const ReturnMappingResult &r, const arma::mat &L);

} // namespace simcoon
