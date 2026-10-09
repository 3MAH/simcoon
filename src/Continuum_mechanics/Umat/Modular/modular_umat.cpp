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
 * @file modular_umat.cpp
 * @brief Implementation of ModularUMAT class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/viscoelastic_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/damage_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/yield_criterion.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/parameter.hpp>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace simcoon {

// ========== Constructor ==========

ModularUMAT::ModularUMAT()
    : elasticity_()
    , mechanisms_()
    , T_init_(0.0)
    , sigma_start_(arma::zeros(6))
    , L_cur_(arma::zeros(6, 6))
    , initialized_(false)
    , maxiter_(100)
    , precision_(1e-9)
{
}

// ========== Configuration ==========

void ModularUMAT::set_elasticity(ElasticityType type, const arma::vec& props, int& offset) {
    elasticity_.configure(type, props, offset);
}

PlasticityMechanism& ModularUMAT::add_plasticity(
    YieldType yield_type,
    IsoHardType iso_type,
    KinHardType kin_type,
    int N_iso,
    int N_kin,
    const arma::vec& props,
    int& offset
) {
    auto mech = std::make_unique<PlasticityMechanism>(yield_type, iso_type, kin_type, N_iso, N_kin);
    mech->configure(props, offset);
    mechanisms_.push_back(std::move(mech));
    return static_cast<PlasticityMechanism&>(*mechanisms_.back());
}

ViscoelasticMechanism& ModularUMAT::add_viscoelasticity(
    int N_prony,
    const arma::vec& props,
    int& offset
) {
    auto mech = std::make_unique<ViscoelasticMechanism>(N_prony);
    mech->configure(props, offset);
    mech->set_reference_stiffness(elasticity_.L0());
    mechanisms_.push_back(std::move(mech));
    return static_cast<ViscoelasticMechanism&>(*mechanisms_.back());
}

DamageMechanism& ModularUMAT::add_damage(
    DamageType damage_type,
    const arma::vec& props,
    int& offset
) {
    auto mech = std::make_unique<DamageMechanism>(damage_type);
    mech->configure(props, offset);
    mechanisms_.push_back(std::move(mech));
    return static_cast<DamageMechanism&>(*mechanisms_.back());
}

void ModularUMAT::configure_from_props(const arma::vec& props, int offset) {
    // Read elasticity type
    int el_type = static_cast<int>(props(offset));
    offset += 1;

    // Configure elasticity
    set_elasticity(static_cast<ElasticityType>(el_type), props, offset);

    // Read number of mechanisms
    int num_mech = static_cast<int>(props(offset));
    offset += 1;

    // Configure each mechanism
    for (int i = 0; i < num_mech; ++i) {
        int mech_type = static_cast<int>(props(offset));
        offset += 1;

        switch (static_cast<MechanismType>(mech_type)) {
            case MechanismType::PLASTICITY: {
                // Read plasticity configuration
                int yield_type = static_cast<int>(props(offset));
                int iso_type = static_cast<int>(props(offset + 1));
                int kin_type = static_cast<int>(props(offset + 2));
                int N_iso = static_cast<int>(props(offset + 3));
                int N_kin = static_cast<int>(props(offset + 4));
                offset += 5;

                add_plasticity(
                    static_cast<YieldType>(yield_type),
                    static_cast<IsoHardType>(iso_type),
                    static_cast<KinHardType>(kin_type),
                    N_iso,
                    N_kin,
                    props,
                    offset
                );
                break;
            }
            case MechanismType::VISCOELASTICITY: {
                int N_prony = static_cast<int>(props(offset));
                offset += 1;
                add_viscoelasticity(N_prony, props, offset);
                break;
            }
            case MechanismType::DAMAGE: {
                int dmg_type = static_cast<int>(props(offset));
                offset += 1;
                add_damage(static_cast<DamageType>(dmg_type), props, offset);
                break;
            }
            default:
                throw std::runtime_error("ModularUMAT: unknown mechanism type " +
                                        std::to_string(mech_type));
        }
    }

    // An anisotropic elastic potential composed with an inelastic mechanism is accepted on
    // purpose; the fibre-convection approximation that entails is documented on
    // structure_tensors_push_forward (hyperelastic.hpp).
}

void ModularUMAT::initialize(int nstatev, arma::vec& statev) {
    // Each mechanism registers its variables into its OWN collection, then
    // receives a base offset into the shared statev layout:
    //   statev = [T_init | mech 0 | mech 1 | ...]  (composition order)
    for (auto& mech : mechanisms_) {
        mech->register_variables();
    }
    unsigned int base = 1;  // statev(0) = T_init (orchestrator-owned)
    for (auto& mech : mechanisms_) {
        mech->compute_offsets(base);
        base += mech->statev_size();
    }

    // Check that we have enough state variables
    int required = static_cast<int>(base);
    if (nstatev < required) {
        throw std::runtime_error("ModularUMAT: nstatev (" + std::to_string(nstatev) +
                                ") < required (" + std::to_string(required) + ")");
    }

    // required_nstatev() is a hand-counted pre-initialize estimate (callers
    // size their statev array with it). Guard it against silent drift from the
    // authoritative count that register_variables() just produced — any
    // mechanism that gains/loses a variable without updating required_nstatev()
    // fails here on the first initialize, not with a corrupted statev later.
    if (required_nstatev() != required) {
        throw std::runtime_error(
            "ModularUMAT: required_nstatev() (" + std::to_string(required_nstatev()) +
            ") disagrees with the registered statev size (" + std::to_string(required) +
            ") — the hand-counted estimate has drifted from register_variables()");
    }

    // Unpack initial values from statev
    T_init_ = statev(0);
    for (auto& mech : mechanisms_) {
        mech->unpack(statev);
    }

    // Cache the composition invariants (constant for the rest of the
    // object's lifetime): per-mechanism constraint-row offsets and the
    // drift-guard arming decision (see the member doc).
    mech_offset_.assign(mechanisms_.size(), 0);
    int acc = 0;
    int n_guarded_rows = 0;
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        mech_offset_[m] = acc;
        acc += mechanisms_[m]->num_constraints();
        if (mechanisms_[m]->guarded_constraints()) {
            // Count ROWS, not mechanisms: a single future mechanism carrying
            // several coupled surfaces (SMA forward/reverse transformation)
            // is already the multi-surface configuration the guard watches.
            n_guarded_rows += mechanisms_[m]->num_constraints();
        }
    }
    drift_guard_armed_ = (n_guarded_rows >= 2);

    initialized_ = true;
}

int ModularUMAT::required_nstatev() const {
    // T_init + all mechanism variables
    int count = 1;  // T_init
    for (const auto& mech : mechanisms_) {
        switch (mech->type()) {
            case MechanismType::PLASTICITY: {
                count += 7;  // p(1) + EP(6)
                auto* pm = dynamic_cast<const PlasticityMechanism*>(mech.get());
                if (pm) {
                    count += 6 * pm->kinematic_hardening().num_backstresses();
                }
                break;
            }
            case MechanismType::VISCOELASTICITY: {
                auto* vm = dynamic_cast<const ViscoelasticMechanism*>(mech.get());
                if (vm) {
                    // v_i (1 scalar) + EV_i (6 Voigt) per Prony branch
                    count += 7 * vm->num_prony_terms();
                }
                break;
            }
            case MechanismType::DAMAGE: {
                count += 2;  // D(1) + Y_max(1)
                break;
            }
        }
    }
    return count;
}

// ========== Main UMAT Entry Point ==========

void ModularUMAT::run(
    const std::string& umat_name,
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    arma::mat& Lt,
    arma::mat& L,
    const arma::mat& DR,
    int nprops,
    const arma::vec& props,
    int nstatev,
    arma::vec& statev,
    double T,
    double DT,
    double Time,
    double DTime,
    double& Wm,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d,
    int ndi,
    int nshr,
    bool start,
    double& tnew_dt,
    int tangent_mode
) {
    // Initialize if first call
    if (!initialized_ || start) {
        // On first call, set up internal variables from statev
        if (!initialized_) {
            // Configure from props if not already done
            if (!elasticity_.is_configured()) {
                int offset = 0;
                configure_from_props(props, offset);
            }
            initialize(nstatev, statev);
        }

        // Store initial temperature
        T_init_ = statev(0);
        if (start) {
            T_init_ = T;
            // Legacy contract: the cumulative work accumulators are RESET on
            // the start increment (every legacy kernel zeroes them in its
            // if(start) block); stale caller-supplied values must not leak in.
            Wm = 0.0;
            Wm_r = 0.0;
            Wm_ir = 0.0;
            Wm_d = 0.0;
        }
    } else {
        // Unpack state variables
        T_init_ = statev(0);
        for (auto& mech : mechanisms_) {
            mech->unpack(statev);
        }
    }

    // Apply rotation for objectivity (DR = I under small strain: nothing to do)
    const bool rotate = !arma::approx_equal(DR, arma::eye(3, 3), "absdiff", 0.0);
    if (rotate) {
        for (auto& mech : mechanisms_) {
            mech->rotate(DR);
        }
    }

    // Save start values
    sigma_start_ = sigma;
    for (auto& mech : mechanisms_) {
        mech->set_start();
    }

    // Set elastic stiffness
    L = elasticity_.L0();

    // Total number of constraints
    int n_total = 0;
    for (const auto& mech : mechanisms_) {
        n_total += mech->num_constraints();
    }

    // Perform return mapping (unconverged-at-maxiter commit semantics: see
    // the return_mapping doc).
    arma::vec Ds_total = arma::zeros(n_total);
    const bool converged = return_mapping(Etot, DEtot, sigma, T_init_, T, DT, DTime, ndi, Ds_total, tangent_mode);

    // Reject unusable committed states with a step cut: non-finite stress
    // (true divergence), a runaway COMMITTED multiplier (multiplier_cap —
    // checked on the state, not the raw FB row, which can be smaller after
    // clipped excursions), or a state that drifted from the flow rule
    // (state_drift) — a spurious CONVERGED root of an oscillating FB Newton,
    // invisible to every residual check (arming policy: see
    // drift_guard_armed_). Thresholds are per-mechanism, infinity = opt out.
    // (Negative multipliers are handled at the source: PlasticityMechanism
    // projects the state update onto p >= p_start, and PowerLawHardening is
    // C1-regularized at onset, so no negative-Ds check is needed here.)
    // The closest-point branch reports its own non-convergence (no commit-at-maxiter
    // semantics there) and satisfies the flow rule by construction (drift guard moot).
    bool reject = !converged || !sigma.is_finite();
    for (size_t m = 0; !reject && m < mechanisms_.size(); ++m) {
        const double drift_tol = mechanisms_[m]->drift_tolerance();
        reject = (mechanisms_[m]->committed_multiplier() >
                      mechanisms_[m]->multiplier_cap())
              || (!cpp_result_.converged && drift_guard_armed_ && std::isfinite(drift_tol) &&
                  mechanisms_[m]->state_drift(sigma_eff_) > drift_tol);
    }
    if (reject) {
        // Ask the global solver to halve the increment and retry. statev is
        // left at its incoming values (pack_all is skipped) so the retry
        // restarts from the correct state; sigma is reset to its incoming
        // value and Lt to elastic so the rejected-step output stays a
        // self-consistent triple — the solver's inforce-at-Dn_mini branch
        // commits this output WITHOUT re-running the UMAT, and a rejected
        // stress paired with the un-updated statev would silently corrupt
        // the rest of the history.
        sigma = sigma_start_;
        tnew_dt = 0.5;
        // Lt is the START state's elastic tangent too: restore the mechanisms
        // from the untouched statev (with the same rotation as above) and
        // evaluate the elastic block there, so neither the rejected iterate's
        // damage nor its elastic strain leaks into the committed tangent.
        for (auto& mech : mechanisms_) {
            mech->unpack(statev);
            if (rotate) mech->rotate(DR);
        }
        arma::vec sigma_at_start;
        refresh_stress(Etot, T - T_init_, ndi, sigma_at_start);
        Lt = stiffness_reduction() * L_cur_;
        return;
    }

    // Compute consistent tangent (it seeds Lt itself)
    compute_tangent(sigma, Ds_total, Lt, tangent_mode);

    // Work quantities — CUMULATIVE in/out, the legacy UMAT contract
    // (Wm += increment each call; the solver reports the path integral).
    // Total mechanical work increment: trapezoidal sigma:dE on the total
    // strain increment, exactly as the legacy kernels.
    const double Wm_inc = 0.5 * arma::dot(sigma_start_ + sigma, DEtot);

    // Stored/dissipated increments from the mechanisms; the recoverable part
    // closes the energy balance (guarantees Wm = Wm_r + Wm_ir + Wm_d and
    // reduces to the legacy elastic trapezoid when no mechanism is present).
    double Wm_ir_inc = 0.0, Wm_d_inc = 0.0;
    for (const auto& mech : mechanisms_) {
        double Wm_r_m = 0.0, Wm_ir_m = 0.0, Wm_d_m = 0.0;
        mech->compute_work(sigma_start_, sigma, Etot, Etot + DEtot, Wm_r_m, Wm_ir_m, Wm_d_m);
        Wm_ir_inc += Wm_ir_m;
        Wm_d_inc += Wm_d_m;
    }

    Wm += Wm_inc;
    Wm_ir += Wm_ir_inc;
    Wm_d += Wm_d_inc;
    Wm_r += Wm_inc - Wm_ir_inc - Wm_d_inc;

    // Pack state variables
    statev(0) = T_init_;
    for (const auto& mech : mechanisms_) {
        mech->pack(statev);
    }
}

void ModularUMAT::refresh_stress(const arma::vec& Etot_end, double DT_init,
                                 int ndi, arma::vec& sigma) {
    arma::vec E_inel = arma::zeros(6);
    for (const auto& mech : mechanisms_) {
        E_inel += mech->inelastic_strain();
    }
    const arma::vec Eel = Etot_end - elasticity_.thermal_strain(DT_init) - E_inel;
    elasticity_.evaluate(Eel, ndi, sigma_eff_, L_cur_);
    // strain equivalence (Lemaitre): the mechanisms work on sigma_eff_, damage scales it
    sigma = stiffness_reduction() * sigma_eff_;
}

bool ModularUMAT::return_mapping(
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    double T_init,
    double T,
    double DT,
    double DTime,
    int ndi,
    arma::vec& Ds_total,
    int tangent_mode
) {
    cpp_result_.converged = false;   // no closest-point solve yet for this increment

    // Total number of constraints
    int n_total = 0;
    for (const auto& mech : mechanisms_) {
        n_total += mech->num_constraints();
    }

    // Elastic prediction. With no constraints there is no inelastic strain and
    // no damage, so this same call IS the whole elastic response — only the
    // constraint machinery below is skipped.
    const arma::vec Etot_end = Etot + DEtot;
    for (auto& mech : mechanisms_) {
        mech->predict(Etot_end, DTime);   // closed-form parts first (viscoelastic branches)
    }
    refresh_stress(Etot_end, T + DT - T_init, ndi, sigma);
    if (n_total == 0) {
        return true;
    }

    // Constraints at the trial state (phase 1 of the first cutting-plane iterate; also what the
    // closest-point branch starts from).
    arma::vec Phi = arma::zeros(n_total);
    arma::vec Y_crit = arma::zeros(n_total);
    evaluate_constraints(Etot_end, DTime, Phi, Y_crit);
    Ds_total.zeros(n_total);

    // Admissible trial state: the whole answer for every integrator. The FB solve would return
    // ds = 0 and error = 0 exactly (phi(Phi, 0) = |Phi| + Phi = 0 for Phi <= 0) and the zero
    // update is the identity for a multiplier row; only the rows without multipliers (damage)
    // commit their fixed point here, after which the stress is refreshed. With no closest-point
    // solve, cpp_consistent_tangent() returns the elastic operator.
    bool elastic_pass = false;   // the admissible-trial pass counted as the loop's first iterate
    if (Phi.max() <= 0.0) {
        bool inert_rows = false;
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (mechanisms_[m]->carries_multipliers()) continue;
            mechanisms_[m]->update(Ds_total, mech_offset_[m]);
            inert_rows = true;
        }
        if (inert_rows) {
            refresh_stress(Etot_end, T + DT - T_init, ndi, sigma);
        }
        double residual = 0.0;
        for (const auto& mech : mechanisms_) {
            residual += mech->consistency_residual(sigma_eff_);
        }
        if (residual <= precision_) {
            if (tangent_mode == tangent_closest_point) {
                cpp_result_ = ReturnMappingResult{};   // the empty solve: the elastic operator
                cpp_result_.converged = true;
            }
            return true;
        }
        evaluate_constraints(Etot_end, DTime, Phi, Y_crit);   // a damage row moved the state: iterate as before
        elastic_pass = true;
    }

    // Dispatch: closest-point when every multiplier row has that form, else cutting plane.
    bool all_cpp = (tangent_mode == tangent_closest_point && ndi == 3);
    for (size_t m = 0; all_cpp && m < mechanisms_.size(); ++m) {
        all_cpp = !mechanisms_[m]->carries_multipliers() || mechanisms_[m]->supports_closest_point();
    }
    if (all_cpp) {
        return return_mapping_cpp(Etot_end, T + DT - T_init, DTime, ndi, Y_crit, sigma, Ds_total);
    }
    return_mapping_ccp(Etot_end, T + DT - T_init, DT, DTime, ndi, Phi, Y_crit, sigma, Ds_total,
                       elastic_pass ? 1 : 0);
    return true;   // the reference CCP convention commits an unconverged-at-maxiter state
}

void ModularUMAT::evaluate_constraints(const arma::vec& Etot_end, double DTime,
                                       arma::vec& Phi, arma::vec& Y_crit) {
    arma::vec Phi_m, Y_crit_m;
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        const int n = mechanisms_[m]->num_constraints();
        mechanisms_[m]->compute_constraints(sigma_eff_, Etot_end, L_cur_, DTime, Phi_m, Y_crit_m);
        Phi.subvec(mech_offset_[m], mech_offset_[m] + n - 1) = Phi_m;
        Y_crit.subvec(mech_offset_[m], mech_offset_[m] + n - 1) = Y_crit_m;
    }
}

void ModularUMAT::return_mapping_ccp(
    const arma::vec& Etot_end,
    double DT_init,
    double DT,
    double DTime,
    int ndi,
    arma::vec& Phi,
    arma::vec& Y_crit,
    arma::vec& sigma,
    arma::vec& Ds_total,
    int iter0
) {
    CuttingPlaneHooks hooks;
    hooks.constraints = [&](arma::vec& Phi_all, arma::vec& Y_all) {
        evaluate_constraints(Etot_end, DTime, Phi_all, Y_all);
    };
    hooks.jacobian = [&](arma::mat& B) { assemble_jacobian(sigma, DT, B); };
    // Incremental CCP update (see the header for why this is not a total-multiplier refresh)
    hooks.update = [&](const arma::vec& ds) {
        for (size_t m = 0; m < mechanisms_.size(); ++m) mechanisms_[m]->update(ds, mech_offset_[m]);
    };
    // D may have evolved in update: the reduction factor is re-evaluated with the stress
    hooks.refresh = [&]() { refresh_stress(Etot_end, DT_init, ndi, sigma); };
    hooks.consistency = [&]() {
        double res = 0.0;
        for (const auto& mech : mechanisms_) res += mech->consistency_residual(sigma_eff_);
        return res;
    };
    ReturnMappingControl control;
    control.maxiter = maxiter_;
    control.precision = precision_;
    cutting_plane_return_mapping(Phi, Y_crit, Ds_total, hooks, control, iter0);
}

ReturnMappingResult closest_point_return_mapping(
    const arma::vec& sigma_tr,
    const arma::mat& L,
    const std::vector<StrainMechanism*>& mechanisms,
    const std::vector<int>& offsets,
    arma::vec& Ds_total,
    const arma::vec& Y_crit_all,
    ReturnStateHooks hooks,
    const ReturnMappingControl& control
) {
    // Helper rows = the multiplier-carrying mechanisms' rows in order; glob maps them into
    // Ds_total.
    std::vector<size_t> active, row_mech;
    std::vector<int> row_c;
    for (size_t m = 0; m < mechanisms.size(); ++m) {
        if (!mechanisms[m]->carries_multipliers()) continue;
        active.push_back(m);
        for (int c = 0; c < mechanisms[m]->num_constraints(); ++c) {
            row_mech.push_back(m);
            row_c.push_back(c);
        }
    }
    const int N = static_cast<int>(row_mech.size());
    arma::uvec glob(N);
    for (int k = 0; k < N; ++k) glob(k) = offsets[row_mech[k]] + row_c[k];

    // Per-row ingredients, refreshed by the state hook (Phi included); the row callbacks
    // read them.
    std::vector<const ClosestPointIngredients*> ing(N, nullptr);
    arma::vec Ds_local = arma::zeros(Ds_total.n_elem);

    hooks.update_state = [&](const arma::vec& sig, const arma::vec& Dl) -> bool {
        Ds_local.elem(glob) = Dl;
        for (size_t m : active) {
            auto& mech = *mechanisms[m];
            if (!mech.refresh_state(sig, Ds_local, offsets[m])) return false;
            const auto& ing_m = *mech.closest_point_ingredients();
            for (int k = 0; k < N; ++k) {
                if (row_mech[k] == m) ing[k] = &ing_m[row_c[k]];
            }
        }
        return true;
    };
    hooks.K = [&](const arma::vec&, const arma::vec&) {
        arma::mat K(N, N);
        for (int l = 0; l < N; ++l) {
            for (int j = 0; j < N; ++j) {
                K(l, j) = (l == j) ? ing[l]->K
                                   : mechanisms[row_mech[l]]->K_cross(row_c[l], *mechanisms[row_mech[j]], row_c[j]);
            }
        }
        return K;
    };
    hooks.flow_state_coupling = [&](const arma::vec&, const arma::vec& Dl) {
        std::vector<arma::vec> c(N);
        for (int k = 0; k < N; ++k) c[k] = Dl(k) * ing[k]->dLambda_dDs;
        return c;
    };

    std::vector<ReturnMechanism> mechs(N);
    for (int k = 0; k < N; ++k) {
        mechs[k].Phi = [&, k](const arma::vec&) { return ing[k]->Phi; };
        mechs[k].dPhi_dsigma = [&, k](const arma::vec&) { return ing[k]->dPhi_dsigma; };
        mechs[k].Lambda = [&, k](const arma::vec&) { return ing[k]->Lambda; };
        mechs[k].dLambda_dsigma = [&, k](const arma::vec&) { return ing[k]->dLambda_dsigma; };
    }

    ReturnMappingResult r = closest_point_return_mapping(sigma_tr, L, mechs, hooks,
                                                         arma::vec(Y_crit_all.elem(glob)), control);
    if (r.converged) Ds_total.elem(glob) = r.Dlambda;
    return r;
}

bool ModularUMAT::return_mapping_cpp(
    const arma::vec& Etot_end,
    double DT_init,
    double DTime,
    int ndi,
    const arma::vec& Y_crit_all,
    arma::vec& sigma,
    arma::vec& Ds_total
) {
    ReturnStateHooks hooks;
    if (!elasticity_.has_constant_stiffness()) {
        // Hyperelastic block: the helper evaluates it at every iterate; the trial elastic strain
        // is the one of the elastic prediction (mechanisms at their start state).
        arma::vec E_inel = arma::zeros(6);
        for (const auto& mech : mechanisms_) E_inel += mech->inelastic_strain();
        hooks.eps_el_tr = Etot_end - elasticity_.thermal_strain(DT_init) - E_inel;
        hooks.elastic_response = [&, ndi](const arma::vec& eps_el, arma::vec& sig, arma::mat& Lt) {
            elasticity_.evaluate(eps_el, ndi, sig, Lt);
        };
    }
    std::vector<StrainMechanism*> rows(mechanisms_.size());
    for (size_t m = 0; m < mechanisms_.size(); ++m) rows[m] = mechanisms_[m].get();

    ReturnMappingControl control;
    control.maxiter = maxiter_;
    control.precision = precision_;
    // sigma_eff_ is the trial effective stress, L_cur_ the elastic tangent there.
    ReturnMappingResult r = closest_point_return_mapping(sigma_eff_, L_cur_, rows, mech_offset_,
                                                         Ds_total, Y_crit_all, std::move(hooks), control);
    if (!r.converged) {
        return false;
    }

    // Rows without multipliers: one evaluation at the converged effective stress (damage
    // fixed point under strain equivalence, viscoelastic rows already taken in predict).
    arma::vec Phi_scratch, Y_scratch;
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        if (mechanisms_[m]->carries_multipliers()) continue;
        mechanisms_[m]->compute_constraints(r.sigma, Etot_end, L_cur_, DTime, Phi_scratch, Y_scratch);
        mechanisms_[m]->update(Ds_total, mech_offset_[m]);
    }

    // Commit: the stress from the refreshed state (self-consistent with it), the solve for
    // the tangent.
    refresh_stress(Etot_end, DT_init, ndi, sigma);
    cpp_result_ = std::move(r);
    return true;
}

void ModularUMAT::assemble_jacobian(
    const arma::vec& sigma,
    double DT,
    arma::mat& B
) {
    // Phase 2: cross-mechanism off-diagonal Jacobian entries per theory eq
    //          B_{lj} = -dPhi^l/dsigma · kappa^j + K^{lj}.
    //          Mechanisms whose Phi is strain-form (viscoelastic) return
    //          empty dPhi_dsigma → no rows contributed here, matching the
    //          decoupled-from-stress structure.
    //          The dot is taken on engineering Voigt components — the
    //          work-conjugate pairing for strain-typed dPhi with
    //          stress-typed kappa (and the established numerics for the
    //          damage strain-typed kappa; see DamageMechanism::kappa).
    B.zeros();
    // Fetch each mechanism's kappa list once — the const-refs stay valid
    // for the whole phase (no compute_constraints call in between).
    std::vector<const std::vector<tensor2>*> kappa_all(mechanisms_.size());
    for (size_t jm = 0; jm < mechanisms_.size(); ++jm) {
        kappa_all[jm] = &mechanisms_[jm]->kappa(sigma_eff_, DT, L_cur_);
    }
    for (size_t lm = 0; lm < mechanisms_.size(); ++lm) {
        const auto& dPhi_l_all = mechanisms_[lm]->dPhi_dsigma(sigma_eff_);
        if (dPhi_l_all.empty()) continue;
        for (size_t l_c = 0; l_c < dPhi_l_all.size(); ++l_c) {
            const int row = mech_offset_[lm] + static_cast<int>(l_c);
            const arma::vec::fixed<6> dPhi_l = dPhi_l_all[l_c].voigt();
            for (size_t jm = 0; jm < mechanisms_.size(); ++jm) {
                const auto& kappa_j_all = *kappa_all[jm];
                for (size_t j_c = 0; j_c < kappa_j_all.size(); ++j_c) {
                    if (lm == jm && l_c == j_c) continue;  // diagonal last
                    const int col = mech_offset_[jm] + static_cast<int>(j_c);
                    B(row, col) = -arma::dot(dPhi_l, kappa_j_all[j_c].voigt())
                                + mechanisms_[lm]->K_cross(
                                      static_cast<int>(l_c),
                                      *mechanisms_[jm],
                                      static_cast<int>(j_c));
                }
            }
        }
    }

    // Phase 3: each mechanism fills its own diagonal (self-stress + K^{ll}).
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        mechanisms_[m]->compute_jacobian_contribution(
            sigma_eff_, L_cur_, B, mech_offset_[m]);
    }
}

double ModularUMAT::stiffness_reduction() const {
    double f = 1.0;
    for (const auto& mech : mechanisms_) {
        f *= mech->stiffness_reduction();
    }
    return f;
}

void ModularUMAT::compute_tangent(
    const arma::vec& sigma,
    const arma::vec& Ds_total,
    arma::mat& Lt,
    int tangent_mode
) {
    Lt = L_cur_;

    if (tangent_mode == tangent_none) {
        // Explicit integration: elastic operator, no assembly.
        return;
    }
    if (tangent_mode == tangent_closest_point && !cpp_result_.converged) {
        tangent_mode = tangent_algorithmic;   // closest-point branch did not run: documented degradation
    }
    if (tangent_mode != tangent_continuum && tangent_mode != tangent_algorithmic
        && tangent_mode != tangent_closest_point) {
        throw std::invalid_argument(
            "ModularUMAT::compute_tangent: unknown tangent_mode "
            + std::to_string(tangent_mode));
    }

    // The assembly below is the tangent of the effective (undamaged) composition: the
    // mechanisms see sigma_eff_; damage enters last through its stress_map().
    const arma::vec& sigma_eff = sigma_eff_;

    // Mechanisms opting into the algorithmic assembly: stress-dependent Phi
    // AND an analytic flow Hessian. Others (viscoelastic, damage, Hessian-less
    // criteria) keep their continuum tangent_contribution in every mode.
    std::vector<size_t> algo;
    if (tangent_mode == tangent_closest_point) {
        // Exact operator of the closest-point solve; every multiplier-carrying row is inside it.
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (mechanisms_[m]->carries_multipliers()) algo.push_back(m);
        }
        Lt = cpp_consistent_tangent(cpp_result_, L_cur_).Lt;
    } else if (tangent_mode == tangent_algorithmic && arma::any(Ds_total > simcoon::iota)) {
        // With no active multiplier the assembly masks every mechanism (same Ds > iota
        // threshold, tangent_assembly.cpp) and returns L: skip the Hessians, which the
        // opt-in predicate below would evaluate before the assembly could; the continuum
        // contributions are no-ops at Ds <= iota.
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (mechanisms_[m]->dLambda_dsigma(sigma_eff) != nullptr &&
                !mechanisms_[m]->dPhi_dsigma(sigma_eff).empty()) {
                algo.push_back(m);
            }
        }
    }

    if (tangent_mode == tangent_algorithmic && !algo.empty()) {
        // Rebuild the local Jacobian at the converged state (the mechanism
        // caches were refreshed by the last compute_constraints call) and
        // hand the opted-in sub-block to the Simo-Hughes assembly.
        // Sign: modular B = -dPhi·kappa + K, tangent_assembly Bhat = dPhi·kappa - K,
        // hence Bhat = -B.
        int n_total = 0;
        for (const auto& mech : mechanisms_) {
            n_total += mech->num_constraints();
        }
        arma::mat B(n_total, n_total);
        assemble_jacobian(sigma, 0.0, B);

        std::vector<int> rows;
        std::vector<arma::vec> kappa_j;
        std::vector<arma::vec> dPhi_l;
        std::vector<arma::mat> dLambda_l;
        for (size_t m : algo) {
            const auto& dPhi_all = mechanisms_[m]->dPhi_dsigma(sigma_eff);
            const auto& kappa_all = mechanisms_[m]->kappa(sigma_eff, 0.0, L_cur_);
            const auto& hess_all = *mechanisms_[m]->dLambda_dsigma(sigma_eff);
            for (size_t c = 0; c < dPhi_all.size(); ++c) {
                rows.push_back(mech_offset_[m] + static_cast<int>(c));
                dPhi_l.emplace_back(dPhi_all[c].voigt());
                kappa_j.emplace_back(kappa_all[c].voigt());
                dLambda_l.emplace_back(hess_all[c].mat());
            }
        }
        const size_t nb = rows.size();
        arma::mat Bhat(nb, nb);
        arma::vec Ds_sub(nb);
        for (size_t i = 0; i < nb; ++i) {
            Ds_sub(i) = Ds_total(rows[i]);
            for (size_t j = 0; j < nb; ++j) {
                Bhat(i, j) = -B(rows[i], rows[j]);
            }
        }
        const ContinuumTangent ct =
            assemble_algorithmic_tangent(Bhat, kappa_j, dPhi_l, Ds_sub, L_cur_, dLambda_l);
        Lt = ct.Lt;
    }

    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        if (std::find(algo.begin(), algo.end(), m) != algo.end()) {
            continue;
        }
        mechanisms_[m]->tangent_contribution(
            sigma_eff, L_cur_, Ds_total, mech_offset_[m], Lt);
    }

    // Left and right maps: d sigma = Q d sigma_eff (damage), eps_in(eps) (viscoelastic branches)
    for (const auto& mech : mechanisms_) {
        const arma::mat Q = mech->stress_map(sigma_eff_);
        if (!Q.is_empty()) {
            Lt = Q * Lt;
        }
    }
    apply_total_strain_maps(Lt);
}

void ModularUMAT::apply_total_strain_maps(arma::mat& Lt) const {
    // Inelastic strains driven by the total strain alone (viscoelastic branches): chain rule,
    // applied last, sigma = F(eps - eps_in(eps)).
    for (const auto& mech : mechanisms_) {
        const arma::mat map = mech->total_strain_map();
        if (!map.is_empty()) {
            Lt = Lt * map;
        }
    }
}

// ========== Standalone UMAT Function ==========

void umat_modular(
    const std::string& umat_name,
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    arma::mat& Lt,
    arma::mat& L,
    const arma::mat& DR,
    const int& nprops,
    const arma::vec& props,
    const int& nstatev,
    arma::vec& statev,
    const double& T,
    const double& DT,
    const double& Time,
    const double& DTime,
    double& Wm,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d,
    const int& ndi,
    const int& nshr,
    const bool& start,
    double& tnew_dt,
    const int& tangent_mode
) {
    // Create a fresh instance each call. Configuration is re-parsed from props
    // and state is restored from statev, ensuring correctness for multi-point
    // simulations. The cost of re-parsing props is negligible compared to the
    // Newton iteration in return_mapping.
    ModularUMAT mumat;

    mumat.run(
        umat_name, Etot, DEtot, sigma, Lt, L, DR,
        nprops, props, nstatev, statev,
        T, DT, Time, DTime,
        Wm, Wm_r, Wm_ir, Wm_d,
        ndi, nshr, start, tnew_dt, tangent_mode
    );
}

} // namespace simcoon
