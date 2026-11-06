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
#include <cstring>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <utility>
#include <vector>

namespace simcoon {

namespace {

/// approx_equal(DR, eye(3, 3), "absdiff", 0.0) without building the identity (NaN or wrong size: false)
bool is_identity_rotation(const arma::mat& DR) {
    if (DR.n_rows != 3 || DR.n_cols != 3) {
        return false;
    }
    for (arma::uword j = 0; j < 3; ++j) {
        for (arma::uword i = 0; i < 3; ++i) {
            if (DR(i, j) != ((i == j) ? 1.0 : 0.0)) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

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
        const arma::vec::fixed<6> Etot_end = Etot + DEtot;   // once, not per mechanism
        for (const auto& mech : mechanisms_) {
            double Wm_r_m = 0.0, Wm_ir_m = 0.0, Wm_d_m = 0.0;
            mech->compute_work(sigma_start_, sigma, Etot, Etot_end, Wm_r_m, Wm_ir_m, Wm_d_m);
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
    for (auto& mech : mechanisms_) {
        mech->predict(Etot_end, DTime);   // closed-form parts first (viscoelastic branches)
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

    // algorithmic assembly only for stress-dependent Phi with an analytic flow Hessian (the Hessians are kept)
    std::vector<std::pair<size_t, const std::vector<tensor4>*>> algo;
    if (tangent_mode == tangent_closest_point) {
        // Exact operator of the closest-point solve; every multiplier-carrying row is inside it.
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (mechanisms_[m]->carries_multipliers()) algo.emplace_back(m, nullptr);
        }
        Lt = cpp_consistent_tangent(cpp_result_, L_cur_).Lt;
    } else if (tangent_mode == tangent_algorithmic && arma::any(Ds_total > simcoon::iota)) {
        // no active multiplier: the assembly masks everything and returns L, skip the Hessians
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            const std::vector<tensor4>* hess = mechanisms_[m]->dLambda_dsigma(sigma_eff);
            if (hess != nullptr && !mechanisms_[m]->dPhi_dsigma(sigma_eff).empty()) {
                algo.emplace_back(m, hess);
            }
        }
    }

    if (tangent_mode == tangent_algorithmic && !algo.empty()) {
        // sign: modular B = -dPhi.kappa + K, tangent_assembly Bhat = dPhi.kappa - K, hence Bhat = -B
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
        rows.reserve(n_total);
        kappa_j.reserve(n_total);
        dPhi_l.reserve(n_total);
        dLambda_l.reserve(n_total);
        for (const auto& [m, hess] : algo) {
            const auto& dPhi_all = mechanisms_[m]->dPhi_dsigma(sigma_eff);
            const auto& kappa_all = mechanisms_[m]->kappa(sigma_eff, 0.0, L_cur_);
            const auto& hess_all = *hess;
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
        if (std::find_if(algo.begin(), algo.end(), [m](const auto& a) { return a.first == m; }) != algo.end()) {
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

} // namespace simcoon
