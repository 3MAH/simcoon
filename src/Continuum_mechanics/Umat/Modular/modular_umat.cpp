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
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/Simulation/Maths/small_lu.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/parameter.hpp>
#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <utility>
#include <vector>

namespace simcoon {


// ========== Constructor ==========

ModularUMAT::ModularUMAT()
    : elasticity_()
    , mechanisms_()
    , T_init_(0.0)
    , sigma_start_(arma::fill::zeros)
    , L_cur_(arma::fill::zeros)
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
    mech->set_thermal_expansion(elasticity_.alpha());
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

void ModularUMAT::set_ncp(NCPMethod ncp) {
    if (ncp == NCPMethod::min_fb_fallback && nonlinear_kinematic_) {
        throw std::invalid_argument("MODUL: ncp='min' serves isotropic and linear (Prager) kinematic "
                                    "hardening only, not Armstrong-Frederick back stresses");
    }
    ncp_ = ncp;
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
    nonlinear_kinematic_ = false;   // an Armstrong-Frederick back stress (see set_ncp)

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
                nonlinear_kinematic_ = nonlinear_kinematic_
                    || (kin_type != static_cast<int>(KinHardType::NONE) && kin_type != static_cast<int>(KinHardType::PRAGER));

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
                throw std::invalid_argument("ModularUMAT: unknown mechanism type " +
                                        std::to_string(mech_type));
        }
    }

    // anisotropic potential + inelastic mechanism accepted: fibre convection approximate (structure_tensors_push_forward)
    set_ncp(ncp_);   // the complementarity function asked before configuration, checked against the composition
}

void ModularUMAT::initialize(int nstatev, arma::vec& statev) {
    // statev = [T_init | mech 0 | mech 1 | ...] (composition order)
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
        throw std::invalid_argument("ModularUMAT: nstatev (" + std::to_string(nstatev) +
                                ") < required (" + std::to_string(required) + ")");
    }

    // required_nstatev() is hand-counted: catch its drift from register_variables() here, not as corrupted statev
    if (required_nstatev() != required) {
        throw std::logic_error(
            "ModularUMAT: required_nstatev() (" + std::to_string(required_nstatev()) +
            ") disagrees with the registered statev size (" + std::to_string(required) +
            ") — the hand-counted estimate has drifted from register_variables()");
    }

    // Unpack initial values from statev
    T_init_ = statev(0);
    for (auto& mech : mechanisms_) {
        mech->unpack(statev);
    }

    // composition invariants: per-mechanism row offsets and the drift-guard arming
    mech_offset_.assign(mechanisms_.size(), 0);
    int acc = 0;
    int n_guarded_rows = 0;
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        mech_offset_[m] = acc;
        acc += mechanisms_[m]->num_constraints();
        if (mechanisms_[m]->guarded_constraints()) {
            // count rows, not mechanisms: one mechanism with coupled surfaces is already multi-surface
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
        else {
            // Reused instance (umat_modular keeps one per configuration): the state
            // comes from statev, exactly as initialize() unpacks it on a fresh one.
            for (auto& mech : mechanisms_) {
                mech->unpack(statev);
            }
        }

        // Store initial temperature
        T_init_ = statev(0);
        if (start) {
            T_init_ = T;
            // legacy contract: the cumulative work accumulators are reset on the start increment
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
    const bool rotate = !is_identity_rotation(DR);
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
    const bool converged = return_mapping(Etot, DEtot, sigma, T_init_, T, DT, DTime, ndi, nshr, Ds_total, tangent_mode);

    // rejection policy: see run() in modular_umat.hpp
    bool reject = !converged || !sigma.is_finite();
    for (size_t m = 0; !reject && m < mechanisms_.size(); ++m) {
        const double drift_tol = mechanisms_[m]->drift_tolerance();
        reject = (mechanisms_[m]->committed_multiplier() >
                      mechanisms_[m]->multiplier_cap())
              || (!cpp_result_.converged && drift_guard_armed_ && std::isfinite(drift_tol) &&
                  mechanisms_[m]->state_drift(sigma_eff_) > drift_tol);
    }
    if (reject) {
        // statev untouched (no pack_all), sigma and Lt back to the start state: a self-consistent triple
        sigma = sigma_start_;
        tnew_dt = 0.5;
        // elastic tangent of the start state: neither the rejected damage nor its elastic strain leaks in
        for (auto& mech : mechanisms_) {
            mech->unpack(statev);
            if (rotate) mech->rotate(DR);
        }
        arma::vec sigma_at_start;
        refresh_stress(Etot, T - T_init_, ndi, nshr, sigma_at_start);
        Lt = stiffness_reduction() * L_cur_;
        return;
    }

    // Compute consistent tangent (it seeds Lt itself)
    compute_tangent(sigma, Ds_total, Lt, tangent_mode);

    // cumulative work (legacy contract): trapezoidal sigma:dE on the total strain increment
    const double Wm_inc = 0.5 * arma::dot(sigma_start_ + sigma, DEtot);

    // Wm_r closes the balance Wm = Wm_r + Wm_ir + Wm_d
    double Wm_ir_inc = 0.0, Wm_d_inc = 0.0;
    if (!mechanisms_.empty()) {
        // the strains the mechanisms see: thermal part removed (the branch forces of PRONK)
        const arma::vec::fixed<6> Eth_start = Etot - elasticity_.alpha() * (T - T_init_);
        const arma::vec::fixed<6> Eth_end = Etot + DEtot - elasticity_.alpha() * (T + DT - T_init_);
        for (const auto& mech : mechanisms_) {
            double Wm_r_m = 0.0, Wm_ir_m = 0.0, Wm_d_m = 0.0;
            mech->compute_work(sigma_start_, sigma, Eth_start, Eth_end, Wm_r_m, Wm_ir_m, Wm_d_m);
            Wm_ir_inc += Wm_ir_m;
            Wm_d_inc += Wm_d_m;
        }
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
                                 int ndi, int nshr, arma::vec& sigma) {
    // no mechanism: the inelastic sum is +0 and x - (+0) = x bit for bit, so it is skipped
    if (!elasticity_.is_configured()) {
        throw std::runtime_error("ElasticityModule: not configured");   // thermal_strain's guard
    }
    arma::vec::fixed<6> Eel = Etot_end - elasticity_.alpha() * DT_init;
    if (elasticity_.omega() != 0.) {   // entropic hyperelasticity: f(theta) = 1 + omega (theta - theta_0)/theta_0
        if (T_init_ <= 0.) {
            throw std::invalid_argument("ModularUMAT: an entropic hyperelastic block (omega != 0) needs an absolute "
                                        "reference temperature (T_init > 0)");
        }
        elasticity_.set_entropic_factor(1. + elasticity_.omega() * DT_init / T_init_);
    }
    if (!mechanisms_.empty()) {
        arma::vec::fixed<6> E_inel(arma::fill::zeros);
        for (const auto& mech : mechanisms_) {
            E_inel += mech->inelastic_strain();
        }
        Eel -= E_inel;
    }
    elasticity_.evaluate(Eel, ndi, nshr, sigma_eff_, L_cur_);
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
    int nshr,
    arma::vec& Ds_total,
    int tangent_mode
) {
    cpp_result_.converged = false;   // no closest-point solve yet for this increment

    // Total number of constraints
    int n_total = 0;
    for (const auto& mech : mechanisms_) {
        n_total += mech->num_constraints();
    }

    // with no constraints this elastic call is the whole response
    const arma::vec::fixed<6> Etot_end = Etot + DEtot;
    const arma::vec::fixed<6> Eth_end = Etot_end - elasticity_.alpha() * (T + DT - T_init);
    for (auto& mech : mechanisms_) {
        mech->predict(Eth_end, DTime);   // closed-form parts first (viscoelastic branches), thermal part removed
    }
    refresh_stress(Etot_end, T + DT - T_init, ndi, nshr, sigma);
    if (n_total == 0) {
        return true;
    }

    // Constraints at the trial state (phase 1 of the first cutting-plane iterate; also what the
    // closest-point branch starts from).
    arma::vec Phi = arma::zeros(n_total);
    arma::vec Y_crit = arma::zeros(n_total);
    evaluate_constraints(Etot_end, DTime, Phi, Y_crit);
    Ds_total.zeros(n_total);

    // admissible trial: the answer for every integrator (FB would return ds = 0); only the rows without
    // multipliers (damage) commit their fixed point here
    bool elastic_pass = false;   // the admissible-trial pass counted as the loop's first iterate
    if (Phi.max() <= 0.0) {
        bool inert_rows = false;
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (mechanisms_[m]->carries_multipliers()) continue;
            mechanisms_[m]->update(Ds_total, mech_offset_[m]);
            inert_rows = true;
        }
        if (inert_rows) {
            refresh_stress(Etot_end, T + DT - T_init, ndi, nshr, sigma);
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
    bool all_cpp = (tangent_mode == tangent_closest_point);
    for (size_t m = 0; all_cpp && m < mechanisms_.size(); ++m) {
        all_cpp = !mechanisms_[m]->carries_multipliers() || mechanisms_[m]->supports_closest_point();
    }
    if (all_cpp) {
        return return_mapping_cpp(Etot_end, T + DT - T_init, DTime, ndi, nshr, Y_crit, sigma, Ds_total);
    }
    return_mapping_ccp(Etot_end, T + DT - T_init, DT, DTime, ndi, nshr, Phi, Y_crit, sigma, Ds_total,
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
    int nshr,
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
    // D may have evolved in update: the reduction factor is re-evaluated with the stress.
    // one captured reference keeps the closure in std::function's small buffer (no heap per call)
    struct refresh_args {
        const arma::vec& Etot_end;
        double DT_init;
        int ndi, nshr;
        arma::vec& sigma;
    } const ra{Etot_end, DT_init, ndi, nshr, sigma};
    hooks.refresh = [this, &ra]() { refresh_stress(ra.Etot_end, ra.DT_init, ra.ndi, ra.nshr, ra.sigma); };
    hooks.consistency = [&]() {
        double res = 0.0;
        for (const auto& mech : mechanisms_) res += mech->consistency_residual(sigma_eff_);
        return res;
    };
    ReturnMappingControl control;
    control.maxiter = maxiter_;
    control.precision = precision_;
    control.ncp = ncp_;
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
    std::vector<int> row_c, row0(mechanisms.size(), -1);   // row0[m] = first helper row of mechanism m
    for (size_t m = 0; m < mechanisms.size(); ++m) {
        if (!mechanisms[m]->carries_multipliers()) continue;
        active.push_back(m);
        row0[m] = static_cast<int>(row_mech.size());
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
            for (int c = 0; c < mech.num_constraints(); ++c) ing[row0[m] + c] = &ing_m[c];
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
    int nshr,
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
        hooks.elastic_response = [&, ndi, nshr](const arma::vec& eps_el, arma::vec& sig, arma::mat& Lt) {
            elasticity_.evaluate(eps_el, ndi, nshr, sig, Lt);
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
    refresh_stress(Etot_end, DT_init, ndi, nshr, sigma);
    cpp_result_ = std::move(r);
    return true;
}

void ModularUMAT::assemble_jacobian(
    const arma::vec& sigma,
    double DT,
    arma::mat& B
) {
    // Phase 2: B_lj = -dPhi^l/dsigma . kappa^j + K^lj (engineering Voigt dot); strain-form Phi
    // (viscoelastic) returns an empty dPhi_dsigma and contributes no row
    B.zeros();
    // single row: no off-diagonal entry, phase 2 skipped (compute_tangent fetches dPhi and kappa itself)
    if (B.n_rows > 1) {
        // kappa lists fetched once: the refs stay valid, nothing below calls compute_constraints
        std::vector<const arma::mat*> kappa_all(mechanisms_.size());
        for (size_t jm = 0; jm < mechanisms_.size(); ++jm) {
            kappa_all[jm] = &mechanisms_[jm]->kappa(sigma_eff_, DT, L_cur_);
        }
        for (size_t lm = 0; lm < mechanisms_.size(); ++lm) {
            const arma::mat& dPhi_l_all = mechanisms_[lm]->dPhi_dsigma(sigma_eff_);
            if (dPhi_l_all.is_empty()) continue;
            for (arma::uword l_c = 0; l_c < dPhi_l_all.n_cols; ++l_c) {
                const int row = mech_offset_[lm] + static_cast<int>(l_c);
                const arma::vec::fixed<6> dPhi_l = dPhi_l_all.col(l_c);
                for (size_t jm = 0; jm < mechanisms_.size(); ++jm) {
                    const arma::mat& kappa_j_all = *kappa_all[jm];
                    for (arma::uword j_c = 0; j_c < kappa_j_all.n_cols; ++j_c) {
                        if (lm == jm && l_c == j_c) continue;  // diagonal last
                        const int col = mech_offset_[jm] + static_cast<int>(j_c);
                        const arma::vec::fixed<6> kappa_j = kappa_j_all.col(j_c);
                        B(row, col) = -arma::dot(dPhi_l, kappa_j)
                                    + mechanisms_[lm]->K_cross(
                                          static_cast<int>(l_c),
                                          *mechanisms_[jm],
                                          static_cast<int>(j_c));
                    }
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
    stress_mapped_ = false;
    ct_cached_ = false;

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

    // algorithmic assembly only for stress-dependent Phi with an analytic flow Hessian (the Hessians are kept)
    std::vector<std::pair<size_t, const arma::mat*>>& algo = algo_;
    algo.clear();
    if (tangent_mode == tangent_closest_point) {
        // Exact operator of the closest-point solve; every multiplier-carrying row is inside it.
        ct_rows_.clear();
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (!mechanisms_[m]->carries_multipliers()) continue;
            algo.emplace_back(m, nullptr);
            for (int c = 0; c < mechanisms_[m]->num_constraints(); ++c) ct_rows_.emplace_back(m, c);
        }
        ct_ = cpp_consistent_tangent(cpp_result_, L_cur_);
        ct_kappa_ = cpp_result_.kappa_j;
        ct_dPhi_ = cpp_result_.dPhidsigma_l;
        ct_cached_ = true;
        Lt = ct_.Lt;
    } else if (tangent_mode == tangent_algorithmic && arma::any(Ds_total > simcoon::iota)) {
        // no active multiplier: the assembly masks everything and returns L, skip the Hessians
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            const arma::mat* hess = mechanisms_[m]->dLambda_dsigma(sigma_eff);
            if (hess != nullptr && !mechanisms_[m]->dPhi_dsigma(sigma_eff).is_empty()) {
                algo.emplace_back(m, hess);
            }
        }
    }

    if (tangent_mode == tangent_algorithmic && !algo.empty()) {
        std::vector<int>& rows = rows_;
        rows.clear();
        ct_rows_.clear();
        size_t nb = 0;
        for (const auto& [m, hess] : algo) nb += mechanisms_[m]->dPhi_dsigma(sigma_eff).n_cols;
        ct_kappa_.set_size(6, nb);
        ct_dPhi_.set_size(6, nb);
        dLambda_l_.set_size(6, 6 * nb);
        size_t i = 0;
        for (const auto& [m, hess] : algo) {
            const arma::mat& dPhi = mechanisms_[m]->dPhi_dsigma(sigma_eff);
            const arma::mat& kappa = mechanisms_[m]->kappa(sigma_eff, 0.0, L_cur_);
            for (int c = 0; c < static_cast<int>(dPhi.n_cols); ++c, ++i) {
                rows.push_back(mech_offset_[m] + c);
                ct_rows_.emplace_back(m, c);
                ct_dPhi_.col(i) = dPhi.col(c);
                ct_kappa_.col(i) = kappa.col(c);
                dLambda_l_.cols(6 * i, 6 * i + 5) = hess->cols(6 * c, 6 * c + 5);
            }
        }
        arma::mat Bhat;
        multiplier_jacobian(sigma, rows, Bhat);
        arma::vec Ds_sub(nb);
        for (size_t k = 0; k < nb; ++k) Ds_sub(k) = Ds_total(rows[k]);
        ct_ = assemble_algorithmic_tangent(Bhat, ct_kappa_, ct_dPhi_, Ds_sub, L_cur_, dLambda_l_);
        ct_cached_ = true;
        Lt = ct_.Lt;
    }

    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        if (std::find_if(algo.begin(), algo.end(), [m](const auto& a) { return a.first == m; }) != algo.end()) {
            continue;
        }
        mechanisms_[m]->tangent_contribution(
            sigma_eff, L_cur_, Ds_total, mech_offset_[m], Lt);
    }

    // Left and right maps: d sigma = Q d sigma_eff (damage), eps_in(eps) (viscoelastic branches);
    // the effective tangent is kept for the heat of the thermomechanical face
    for (const auto& mech : mechanisms_) {
        const arma::mat Q = mech->stress_map(sigma_eff_);
        if (!Q.is_empty()) {
            if (!stress_mapped_) {
                Lt_eff_ = Lt;
                stress_mapped_ = true;
            }
            Lt = Q * Lt;
        }
    }
    apply_total_strain_maps(Lt);
    if (stress_mapped_) apply_total_strain_maps(Lt_eff_);
}

void ModularUMAT::multiplier_jacobian(const arma::vec& sigma, const std::vector<int>& rows, arma::mat& Bhat) {
    // sign: modular B = -dPhi.kappa + K, tangent_assembly Bhat = dPhi.kappa - K, hence Bhat = -B
    int n_total = 0;
    for (const auto& mech : mechanisms_) n_total += mech->num_constraints();
    arma::mat B(n_total, n_total);
    assemble_jacobian(sigma, 0.0, B);
    const size_t nb = rows.size();
    Bhat.set_size(nb, nb);
    for (size_t i = 0; i < nb; ++i) {
        for (size_t j = 0; j < nb; ++j) Bhat(i, j) = -B(rows[i], rows[j]);
    }
}

void ModularUMAT::apply_total_strain_maps(arma::mat& Lt) const {
    // Inelastic strains driven by the total strain alone (viscoelastic branches): chain rule,
    // applied last, sigma = F(eps - eps_in(eps)).
    for (const auto& mech : mechanisms_) {
        const arma::mat& map = mech->total_strain_map();
        if (!map.is_empty()) {
            Lt = mul66(Lt, map);
        }
    }
}

// ========== Standalone UMAT Function ==========

namespace {

// configured instances of this thread, most recent first, keyed on (props bitwise, nstatev, kinematic form);
// the state is unpacked from statev at every call, so a reused instance returns a fresh one's result
class modular_umat_cache {
public:
    static constexpr std::size_t capacity = 8;

    /// the instance of (props, nstatev, kinematic_form), created on a miss; least recently used evicted
    ModularUMAT& acquire(const arma::vec& props, int nstatev, KinematicForm kinematic_form) {
        for (auto it = entries_.begin(); it != entries_.end(); ++it) {
            if (it->matches(props, nstatev, kinematic_form)) {
                std::rotate(entries_.begin(), it, it + 1);
                return *entries_.front().umat;
            }
        }
        if (entries_.size() == capacity) {
            entries_.pop_back();
        }
        entry e;
        e.props.assign(props.begin(), props.end());
        e.nstatev = nstatev;
        e.kinematic_form = kinematic_form;
        e.umat = std::make_unique<ModularUMAT>();
        e.umat->set_kinematic_form(kinematic_form);
        entries_.insert(entries_.begin(), std::move(e));
        return *entries_.front().umat;
    }

    /// Forget the instance returned by the last acquire().
    void drop_most_recent() {
        if (!entries_.empty()) {
            entries_.erase(entries_.begin());
        }
    }

private:
    struct entry {
        std::vector<double> props;
        int nstatev = 0;
        KinematicForm kinematic_form = KinematicForm::rate;
        std::unique_ptr<ModularUMAT> umat;

        bool matches(const arma::vec& p, int n, KinematicForm kf) const {
            return nstatev == n && kinematic_form == kf && props.size() == p.n_elem
                && (props.empty() || std::memcmp(props.data(), p.memptr(), props.size() * sizeof(double)) == 0);
        }
    };
    std::vector<entry> entries_;
};

}  // namespace

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
    const int& tangent_mode,
    const KinematicForm& kinematic_form,
    const int& ncp
) {
    thread_local modular_umat_cache cache;
    ModularUMAT& mumat = cache.acquire(props, nstatev, kinematic_form);
    try {
        mumat.set_ncp(ncp_method_of(ncp, "MODUL"));
        mumat.run(
            umat_name, Etot, DEtot, sigma, Lt, L, DR,
            nprops, props, nstatev, statev,
            T, DT, Time, DTime,
            Wm, Wm_r, Wm_ir, Wm_d,
            ndi, nshr, start, tnew_dt, tangent_mode
        );
    } catch (...) {
        // a failed configuration or initialization must not be reused
        cache.drop_most_recent();
        throw;
    }
}

// ========== Thermomechanical face ==========

namespace {
// the memo key of start_iso_energy: Etot, DR (identity when null), T, statev
void iso_memo_key(std::vector<double>& key, const arma::vec& Etot, const double* DR, double T, const arma::vec& statev) {
    key.resize(16 + statev.n_elem);
    std::memcpy(key.data(), Etot.memptr(), 6 * sizeof(double));
    if (DR) {
        std::memcpy(key.data() + 6, DR, 9 * sizeof(double));
    } else {
        for (int k = 0; k < 9; ++k) key[6 + k] = (k % 4 == 0) ? 1. : 0.;
    }
    key[15] = T;
    std::memcpy(key.data() + 16, statev.memptr(), statev.n_elem * sizeof(double));
}
}  // namespace

void ModularUMAT::run_T(
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    double& r,
    arma::mat& dSdE,
    arma::vec& dSdT,
    arma::vec& drdE,
    double& drdT,
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
    double& Wt,
    double& Wt_r,
    double& Wt_ir,
    int ndi,
    int nshr,
    bool start,
    double& tnew_dt,
    int tangent_mode
) {
    if (!elasticity_.is_configured()) {
        int offset = 2;   // [rho, c_p, the MODUL stream]
        configure_from_props(props, offset);
    }
    const bool entropic = (elasticity_.omega() != 0.);
    const double W_iso_start = entropic ? start_iso_energy(Etot, DR, T, nstatev, statev, start) : 0.;
    const double c_0 = props(0) * props(1);
    const double T1 = T + DT;
    if (start) {
        Wt = 0.0;
        Wt_r = 0.0;
        Wt_ir = 0.0;
    }
    const double Wm_d_before = start ? 0.0 : Wm_d;

    arma::mat L;
    run("MODUL", Etot, DEtot, sigma, dSdE, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime,
        Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt, tangent_mode);
    const arma::vec& alpha = elasticity_.alpha();
    if (tnew_dt < 1.0) {   // rejected by run(): the start state, no heat
        dSdT = -dSdE * alpha;
        r = 0.0;
        drdE.zeros(6);
        drdT = 0.0;
        return;
    }

    // multiplier rows (mechanism, constraint) of the committed state; the viscoelastic rows are inert
    const arma::vec::fixed<6> Etot_end = Etot + DEtot;
    const arma::vec::fixed<6> Eth_start = Etot - elasticity_.thermal_strain(T - T_init_);
    const arma::vec::fixed<6> Eth_end = Etot_end - elasticity_.thermal_strain(T1 - T_init_);
    std::vector<std::pair<size_t, int>>& rows_T = rows_T_;
    rows_T.clear();
    bool active = false;
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        if (!mechanisms_[m]->carries_multipliers()) continue;
        const int nc = static_cast<int>(mechanisms_[m]->dPhi_dsigma(sigma_eff_).n_cols);
        for (int c = 0; c < nc; ++c) rows_T.emplace_back(m, c);
        active = active || (nc > 0 && mechanisms_[m]->committed_multiplier() > simcoon::iota);
    }
    const size_t n_rows = rows_T.size();
    // the assembly of compute_tangent serves when it covers these rows (modes 2 and 3); the mechanisms
    // without multipliers refresh their caches at the committed state either way
    const bool reuse = (tangent_mode == tangent_algorithmic || tangent_mode == tangent_closest_point)
                       && (!active || (ct_cached_ && ct_rows_ == rows_T));
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        if (reuse && mechanisms_[m]->carries_multipliers()) continue;
        arma::vec Phi_m, Y_m;
        mechanisms_[m]->compute_constraints(sigma_eff_, Eth_end, L_cur_, DTime, Phi_m, Y_m);
    }
    ContinuumTangent ct_local;
    arma::mat kappa_local, dPhi_local;
    if (n_rows == 0 || (reuse && !active)) {   // no row, or every row masked: zero sensitivities, no flux needed
        ct_local.L_tilde = L_cur_;
        ct_local.P_epsilon.zeros(6, n_rows);
        ct_local.invBhat.zeros(n_rows, n_rows);
        kappa_local.zeros(6, n_rows);
        dPhi_local.zeros(6, n_rows);
    } else if (!reuse) {   // sensitivities assembled at the committed state (continuum for tangent_none)
        std::vector<int>& rows = rows_;
        rows.clear();
        std::vector<std::pair<size_t, const arma::mat*>>& hess = algo_;   // (mechanism, flow Hessian) per row
        hess.clear();
        kappa_local.set_size(6, n_rows);
        dPhi_local.set_size(6, n_rows);
        arma::vec Ds_rows(n_rows);
        int mode = (tangent_mode == tangent_none) ? tangent_continuum : (tangent_mode == tangent_closest_point) ? tangent_algorithmic : tangent_mode;
        for (size_t l = 0; l < n_rows; ++l) {
            const auto [m, c] = rows_T[l];
            rows.push_back(mech_offset_[m] + c);
            dPhi_local.col(l) = mechanisms_[m]->dPhi_dsigma(sigma_eff_).col(c);
            kappa_local.col(l) = mechanisms_[m]->kappa(sigma_eff_, 0.0, L_cur_).col(c);
            Ds_rows(l) = mechanisms_[m]->committed_multiplier();
            hess.emplace_back(m, mechanisms_[m]->dLambda_dsigma(sigma_eff_));
            if (hess.back().second == nullptr) mode = tangent_continuum;   // no flow Hessian: the continuum sensitivities
        }
        arma::mat Bhat;
        multiplier_jacobian(sigma, rows, Bhat);
        ct_local = compute_tangent_operator(mode, Bhat, kappa_local, dPhi_local, Ds_rows, L_cur_,
            [&]() -> arma::mat {
                arma::mat H(6, 6 * n_rows);
                for (size_t l = 0; l < n_rows; ++l) {
                    const int c = rows_T[l].second;
                    H.cols(6 * l, 6 * l + 5) = hess[l].second->cols(6 * c, 6 * c + 5);
                }
                return H;
            });
    }
    const bool shared = reuse && active;
    const ContinuumTangent& ct = shared ? ct_ : ct_local;
    const arma::mat& kappa_j = shared ? ct_kappa_ : kappa_local;
    const arma::mat& dPhi_l = shared ? ct_dPhi_ : dPhi_local;

    // trial stress rate q of a temperature move and P_theta (run_T in the header)
    arma::vec alpha_eff = alpha;
    if (tangent_mode != tangent_none) {
        for (auto it = mechanisms_.rbegin(); it != mechanisms_.rend(); ++it) {
            const arma::mat& map = (*it)->total_strain_map();
            if (!map.is_empty()) mul6(map, alpha_eff, alpha_eff);
        }
    }
    arma::vec T_iso(6, arma::fill::zeros);
    const double w0 = entropic ? elasticity_.omega() / T_init_ : 0.;
    if (entropic) {
        T_iso = sigma_eff_;
        T_iso.subvec(0, 2) -= elasticity_.vol_pressure();
        T_iso /= elasticity_.entropic_factor();
    }
    arma::vec::fixed<6> q;
    mul6(ct.L_tilde, alpha_eff, q);
    q = -q;
    // L^-1 T_iso of the entropic term, and L^-T T_iso for the entropy chain below (entropic only)
    arma::vec::fixed<6> L_inv_T, Lt_inv_T;
    if (entropic) {
        small_lu<6> lu;
        if (!lu.factor(L_cur_) ) throw simcoon::exception_inv("ModularUMAT::run_T: singular elastic tangent");
        L_inv_T = T_iso;
        lu.solve_in_place(L_inv_T.memptr());
        arma::vec::fixed<6> Lt_L_inv_T;
        mul6(ct.L_tilde, L_inv_T, Lt_L_inv_T);
        q += w0 * Lt_L_inv_T;
        const arma::mat::fixed<6, 6> L_t = L_cur_.t();
        if (!lu.factor(L_t)) throw simcoon::exception_inv("ModularUMAT::run_T: singular elastic tangent");
        Lt_inv_T = T_iso;
        lu.solve_in_place(Lt_inv_T.memptr());
    }
    arma::vec g(n_rows);
    for (size_t l = 0; l < n_rows; ++l) {
        const arma::vec::fixed<6> dPhi_m = dPhi_l.col(l);
        g(l) = arma::dot(dPhi_m, q);
    }
    const arma::vec P_theta = ct.invBhat * g;
    arma::vec dSdT_eff = q;   // of the effective stress
    const arma::mat& kappa_t = ct.kappa_tilde.is_empty() ? kappa_j : ct.kappa_tilde;   // M^-1 kappa (kappa: no correction)
    for (size_t l = 0; l < n_rows; ++l) {
        const arma::vec::fixed<6> kappa_l = kappa_t.col(l);
        dSdT_eff -= kappa_l * P_theta(l);
    }
    dSdT = dSdT_eff;
    if (stress_mapped_) {   // damage: d sigma = Q d sigma_eff, as the tangent
        for (const auto& mech : mechanisms_) {
            const arma::mat Q = mech->stress_map(sigma_eff_);
            if (!Q.is_empty()) dSdT = Q * dSdT;
        }
    }
    const arma::mat& dSdE_eff = stress_mapped_ ? Lt_eff_ : dSdE;

    // entropy, dissipation and the heat of the increment, per reference volume
    const double eta = c_0 * std::log(T1 / T_init_) + arma::dot(alpha, sigma) - w0 * elasticity_.W_iso();
    const double eta_start = c_0 * std::log(T / T_init_) + arma::dot(alpha, sigma_start_) - w0 * W_iso_start;
    const double Deta = eta - eta_start;
    const double Dgamma = Wm_d - Wm_d_before;
    arma::vec dDgamma_dE(6, arma::fill::zeros);
    double dDgamma_dT = 0.0;
    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        size_t first = n_rows, n_m = 0;   // the rows of mechanism m are consecutive
        for (size_t l = 0; l < n_rows; ++l) {
            if (rows_T[l].first == m) {
                first = std::min(first, l);
                ++n_m;
            }
        }
        const arma::mat P_eps_m = (n_m > 0) ? arma::mat(ct.P_epsilon.cols(first, first + n_m - 1)) : arma::mat(6, 0);
        const arma::vec P_th_m = (n_m > 0) ? arma::vec(P_theta.subvec(first, first + n_m - 1)) : arma::vec();
        const StrainMechanism::DissipationInputs in{sigma_start_, sigma, sigma_eff_, Eth_start, Eth_end,
                                                    dSdE, dSdT, dSdE_eff, dSdT_eff, P_eps_m, P_th_m};
        arma::vec dE;
        double dT_m;
        mechanisms_[m]->dissipation_sensitivity(in, dE, dT_m);
        dDgamma_dE += dE;
        dDgamma_dT += dT_m;
    }
    // dW_iso = T_iso : dH_e
    arma::vec dW_dE = T_iso;
    double dW_dT = -arma::dot(T_iso, alpha);
    if (entropic && !mechanisms_.empty()) {
        dW_dE = mul6_At(dSdE_eff, Lt_inv_T);
        dW_dT = arma::dot(Lt_inv_T, dSdT_eff - w0 * T_iso);
    }
    heat_source_increment(Dgamma, dDgamma_dE, dDgamma_dT, Deta, arma::vec(mul6_At(dSdE, alpha) - w0 * dW_dE),
                          c_0 / T1 + arma::dot(dSdT, alpha) - w0 * dW_dT, T1, r, drdE, drdT);
    Wt += T1 * Deta;
    Wt_r += T1 * Deta;
    if (entropic) {   // the end state is the next start of a sequential path: its W_iso is at hand
        iso_memo_key(iso_key_[1], Etot_end, nullptr, T1, statev);
        iso_W_[1] = elasticity_.W_iso();
    }
}


double ModularUMAT::start_iso_energy(const arma::vec& Etot, const arma::mat& DR, double T, int nstatev,
                                     const arma::vec& statev, bool start) {
    if (!initialized_) {
        arma::vec sv = statev;
        initialize(nstatev, sv);
    }
    std::vector<double>& key = iso_key_scratch_;
    iso_memo_key(key, Etot, DR.memptr(), T, statev);
    for (int k = 0; k < 2; ++k) {
        if (iso_key_[k] == key) return iso_W_[k];
    }
    T_init_ = start ? T : statev(0);
    const bool rotate = !is_identity_rotation(DR);   // Etot arrives transported with DR, the state follows as in run()
    for (auto& mech : mechanisms_) {
        mech->unpack(statev);
        if (rotate) mech->rotate(DR);
    }
    arma::vec::fixed<6> Eel = Etot - elasticity_.alpha() * (T - T_init_);
    if (!mechanisms_.empty()) {   // as refresh_stress
        arma::vec::fixed<6> E_inel(arma::fill::zeros);
        for (const auto& mech : mechanisms_) {
            E_inel += mech->inelastic_strain();
        }
        Eel -= E_inel;
    }
    iso_W_[0] = elasticity_.iso_energy(Eel);
    iso_key_[0].swap(key);
    return iso_W_[0];
}

void umat_modular_T(
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    double& r,
    arma::mat& dSdE,
    arma::vec& dSdT,
    arma::vec& drdE,
    double& drdT,
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
    double& Wt,
    double& Wt_r,
    double& Wt_ir,
    const int& ndi,
    const int& nshr,
    const bool& start,
    double& tnew_dt,
    const int& tangent_mode,
    const KinematicForm& kinematic_form,
    const int& ncp
) {
    thread_local modular_umat_cache cache;   // its own: these props carry [rho, c_p] in front
    ModularUMAT& mumat = cache.acquire(props, nstatev, kinematic_form);
    try {
        mumat.set_ncp(ncp_method_of(ncp, "MODUL"));
        mumat.run_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR, nprops, props, nstatev, statev,
                    T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir, ndi, nshr, start, tnew_dt, tangent_mode);
    } catch (...) {
        cache.drop_most_recent();
        throw;
    }
}

} // namespace simcoon
