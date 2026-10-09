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
 *   \b{R}_\sigma = \boldsymbol{\sigma} - \mathbf{F}\Big(\boldsymbol{\varepsilon}_e^{tr}
 *     - \sum_j \Delta s^j\,\boldsymbol{\Lambda}_\varepsilon^j(\boldsymbol{\sigma},\mathbf{V})\Big) = \mathbf{0},
 *   \qquad
 *   \Delta s^j \ge 0,\; \Phi^j \le 0,\; \Delta s^j\,\Phi^j = 0,
 * \f]
 * with the internal state \f$ \mathbf{V} \f$ resolved by backward Euler at every iterate
 * (inner-consistent scheme), via a condensed semi-smooth Newton iteration on the reduced
 * \f$ N \times N \f$ multiplier system handed to Fischer_Burmeister_m(). \f$ \mathbf{F} \f$ is the
 * elastic block: by default the linear one, \f$ \mathbf{F}(\boldsymbol{\varepsilon}_e) =
 * \mathbf{L}\boldsymbol{\varepsilon}_e \f$, so that \f$ \b{R}_\sigma = \boldsymbol{\sigma} -
 * \boldsymbol{\sigma}^{tr} + \sum_j \Delta s^j \mathbf{L}\boldsymbol{\Lambda}^j \f$; a nonlinear
 * block (hyperelastic potential) is given through ReturnStateHooks::elastic_response, and
 * \f$ \mathbf{L} \f$ is then its tangent at the iterate. The condensation operator
 * \f$ \mathbf{M} = \mathbf{I} + \mathbf{L}:\sum_j \Delta s^j\,
 * \partial\boldsymbol{\Lambda}_\varepsilon^j/\partial\boldsymbol{\sigma} \f$ is the same
 * operator as in assemble_algorithmic_tangent(), so at convergence the exact consistent
 * tangent is available from the returned ingredients at no extra cost
 * (see cpp_consistent_tangent()).
 *
 * Unlike the convex-cutting-plane (CCP) loops, the inelastic strain uses the flow at the
 * CONVERGED state: \f$ \boldsymbol{\varepsilon}^{in} = \boldsymbol{\varepsilon}^{in}_n
 * + \sum_j \Delta s^j \boldsymbol{\Lambda}_\varepsilon^j(\boldsymbol{\sigma}_{n+1},
 * \mathbf{V}_{n+1}) \f$. For radial flows (von Mises with isotropic hardening) the two
 * integrators coincide; where the flow direction rotates within the increment the converged
 * stresses differ by \f$ O(\|\Delta\varepsilon\|^2) \f$ on one increment.
 *
 * The Newton step is taken over the active set — rows with \f$ \Phi^l \ge 0 \f$, or a
 * multiplier, or a non-vanishing diagonal \f$ B^{ll} \f$; a row failing all three (a
 * reorientation surface at zero effective stress) has an identically zero Fischer-Burmeister
 * row and no step, as in the cutting-plane loops of the SMA kernels.
 *
 * Robustness. The semi-smooth Newton is globalised by (i) the \f$ \Delta s^j \ge 0 \f$
 * projection of the multipliers, (ii) a backtracking line search (step halved, at most
 * max_backtrack times) on the merit \f$ e = e_{FB}(\Phi, \Delta s) +
 * \|\mathbf{R}_\sigma\|/\sigma_{ref} \f$ — the Fischer-Burmeister residual of the TRUE
 * constraints (Fischer_Burmeister_residual) plus the stress residual — a trial whose inner
 * state solve fails or is not finite counting as a rejected step, and (iii) the deactivation
 * snap: after 20 iterations a row with \f$ \Phi^l < 0 \f$ and \f$ \Delta s^l |B_{ll}| <
 * 10^{-2} Y^l_{crit} \f$ is set exactly inactive (\f$ \Delta s^l = 0 \f$), which is the
 * complementarity solution the semi-smooth iteration otherwise approaches asymptotically
 * (the FB function has no finite-step root on that boundary). The iterate is never clamped
 * and the tolerance never loosened: a cap on the Newton step or an acceptance at a looser
 * residual would commit a state that is not a solution of the system above. Non-convergence
 * (maxiter, singular condensation, LAPACK failure, NaN) returns `converged = false`; the
 * caller requests the standard step cut (tnew_dt = 0.5), which is the robust answer at the
 * level where the increment can be changed, and never falls back to a cutting-plane update.
 * Alternatives for the local loop, not implemented: a smoothed Fischer-Burmeister function
 * with continuation \f$ \mu \to 0 \f$ (Chen-Harker-Kanzow-Smale), which trades the
 * non-smooth boundary for an outer loop; an active-set Newton (working set from the sign of
 * \f$ \Phi \f$, plain Newton on it, reactivation checks), which is faster per iterate but
 * needs its own combinatorial safeguard. Both would be variants of this function, not of
 * its callers.
 *
 * @note The callbacks map one-to-one onto the modular StrainMechanism interface —
 * compute_constraints -> Phi, closest_point_ingredients -> dPhi_dsigma / dLambda_dsigma / K /
 * flow_state_coupling, refresh_state -> update_state — the overload over StrainMechanism rows
 * (strain_mechanism.hpp) is the adapter, used by ModularUMAT::return_mapping_cpp and by the
 * dedicated plasticity kernels (EPICP, EPCHA, EPICP_T, EPKCP_T), which run the matching
 * PlasticityMechanism on their own props and state; the SMA kernels (unified_T / unified_TR)
 * call this function directly with finite_difference_total_derivatives().
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
    /// OPTIONAL multiplier-side flow/state chain, N vectors of 6 (strain-typed):
    /// \f$ \Delta s^j\,\partial\boldsymbol{\Lambda}^j/\partial\Delta s^j|_\sigma \f$ — the flow's
    /// dependence on its own multiplier through the state (backstress relaxation). The function
    /// forms the flux chain \f$ \mathbf{c}^j = \mathbf{L}\,\cdot \f$ that with the elastic
    /// tangent at the iterate and adds it to \f$ \boldsymbol{\kappa}^j \f$ in the Newton and in
    /// the tangent ingredients. Empty => omitted (superlinear instead of quadratic for
    /// state-coupled mechanisms; converged solution unaffected).
    std::function<std::vector<arma::vec>(const arma::vec &sigma, const arma::vec &Dlambda)> flow_state_coupling;
    /// OPTIONAL nonlinear elastic block: given the elastic strain (6, engineering Voigt),
    /// return the stress (6) and the tangent (6x6). Empty => linear block, \f$ \mathbf{F} =
    /// \mathbf{L}\boldsymbol{\varepsilon}_e \f$ with the L of the call. When set, the caller
    /// also provides the trial elastic strain (ReturnStateHooks::eps_el_tr), and the L argument
    /// of closest_point_return_mapping() is the block's tangent at the trial state.
    std::function<void(const arma::vec &eps_el, arma::vec &sigma, arma::mat &L)> elastic_response;
    arma::vec eps_el_tr;   ///< trial elastic strain (6), required with elastic_response
};

/// Iteration controls; zero-valued members fall back to the simcoon defaults.
struct ReturnMappingControl {
    int    maxiter = 0;             ///< 0 => simcoon::maxiter_umat
    double precision = 0.;          ///< 0 => simcoon::precision_umat
    int    max_backtrack = 5;       ///< halvings of the Newton step when the merit grows by more than x2
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
    double    error = 0.;                     ///< final merit (FB + \f$ \|R_\sigma\|/\sigma_{ref} \f$)
    std::vector<double> error_history;        ///< merit per iteration (rate diagnostics)
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
 * \f$ \Delta s\,\partial\boldsymbol{\Lambda}/\partial\Delta s|_\sigma \f$ (6, strain-typed;
 * the function multiplies by the elastic tangent). Required for the consistent tangent to be
 * exact when the flow depends on \f$ \Delta s \f$ through the state (e.g. backstress); empty
 * otherwise.
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
 * @brief Hooks of the convex-cutting-plane (CCP) loop: the state lives with the caller and
 *        is updated INCREMENTALLY along the flow of each iterate.
 *
 * The constraints of the first iterate are evaluated by the caller (shared with its elastic
 * guard) and handed in; the loop then repeats
 * jacobian -> Fischer_Burmeister_m -> update -> refresh -> consistency -> constraints
 * until the Fischer-Burmeister error plus the consistency residual is below the precision or
 * maxiter is reached. This is the path-dependent integrator the algorithmic tangent is only
 * approximately consistent with (closest_point_return_mapping() for the exact one); it is kept
 * because it is the reference of every result before mode 3 and the only one for criteria
 * without a flow Hessian.
 */
struct CuttingPlaneHooks {
    /// REQUIRED. Constraints and normalisations at the current state and stress (phase 1 of an
    /// iterate; populates whatever caches jacobian() reads).
    std::function<void(arma::vec &Phi, arma::vec &Y_crit)> constraints;
    /// REQUIRED. Local multiplier Jacobian \f$ B^{lj} = -\partial\Phi^l/\partial\boldsymbol{\sigma}
    /// \cdot\boldsymbol{\kappa}^j + K^{lj} \f$ (N x N) from those caches.
    std::function<void(arma::mat &B)> jacobian;
    /// REQUIRED. Incremental state update along the current flow for the multiplier step @p ds.
    std::function<void(const arma::vec &ds)> update;
    /// REQUIRED. Stress from the updated state.
    std::function<void()> refresh;
    /// OPTIONAL. Residual of the states updated outside the FB rows (damage fixed point);
    /// added to the FB error. Empty => 0.
    std::function<double()> consistency;
};

/**
 * @brief Convex-cutting-plane return mapping over caller-owned state.
 *
 * @param[in,out] Phi      constraints at the entering state (N), updated at every iterate
 * @param[in,out] Y_crit   their normalisations (N), updated with them
 * @param[in,out] Ds_total total multipliers (N), accumulated by Fischer_Burmeister_m
 * @param hooks            the caller's state operations
 * @param control          maxiter / precision (zero => simcoon defaults)
 * @param iter0            iterations already spent by the caller on this increment (its
 *                         elastic pass), counted against maxiter
 * @return whether the error fell below the precision (the modular UMAT commits the state either
 *         way, by the reference cutting-plane convention)
 */
bool cutting_plane_return_mapping(
    arma::vec &Phi,
    arma::vec &Y_crit,
    arma::vec &Ds_total,
    const CuttingPlaneHooks &hooks,
    const ReturnMappingControl &control = {},
    int iter0 = 0);

/**
 * @brief Total derivatives by central differences of the composed map, for a kernel whose state
 *        is not an analytic function of \f$ (\boldsymbol{\sigma}, \Delta s) \f$.
 *
 * Fills the derivative callbacks of @p mechanisms (dPhi_dsigma, dLambda_dsigma) and the state
 * hooks (update_state, K, flow_state_coupling) from the kernel's state refresh and its
 * \f$ \Phi^j(\boldsymbol{\sigma}) \f$, \f$ \boldsymbol{\Lambda}^j(\boldsymbol{\sigma}) \f$ (which
 * must already be set on @p mechanisms, evaluated at the refreshed state): every probe
 * re-solves the state at the perturbed \f$ (\boldsymbol{\sigma}, \Delta s) \f$ between
 * @p save and @p restore, so the differences are TOTAL derivatives — the quantities the
 * closest-point Newton and its exact tangent need (the lesson of the state-coupled kernels:
 * partial derivatives at frozen state leave the operator inexact). Steps: relative
 * \f$ h_\sigma = \mathrm{h\_sigma\_rel}\,(\|\boldsymbol{\sigma}\| + 1) \f$ on the stress, absolute
 * \f$ h_{\Delta s} \f$ on the multipliers (one-sided at \f$ \Delta s = 0 \f$).
 *
 * Cost: 12 state refreshes per Newton iterate for the stress derivatives of every row and
 * quantity (one probe pass, cached per \f$ (\boldsymbol{\sigma}, \Delta s) \f$) plus
 * \f$ 2 N \f$ for the multiplier ones — the SMA kernels on this builder stay several times
 * slower than their cutting-plane loop; analytic derivatives (PlasticityMechanism) are the remedy.
 *
 * @param mechanisms  N mechanisms with Phi and Lambda set; the derivatives are filled here
 * @param hooks       filled: update_state (wrapping @p refresh), K, flow_state_coupling
 * @param refresh     the kernel's backward-Euler state refresh at \f$ (\boldsymbol{\sigma}, \Delta s) \f$
 * @param save        copies the kernel's state aside (a lambda capturing a local struct)
 * @param restore     puts it back
 */
void finite_difference_total_derivatives(
    std::vector<ReturnMechanism> &mechanisms,
    ReturnStateHooks &hooks,
    const std::function<bool(const arma::vec &sigma, const arma::vec &Dlambda)> &refresh,
    const std::function<void()> &save,
    const std::function<void()> &restore,
    double h_sigma_rel = 1.e-5,
    double h_Dlambda = 1.e-7);

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
