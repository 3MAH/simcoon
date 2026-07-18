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
    mech->set_reference_stiffness(elasticity_.L());
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

    // Cache per-mechanism constraint-row offsets (constant for the rest of
    // the object's lifetime).
    mech_offset_.assign(mechanisms_.size(), 0);
    {
        int acc = 0;
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            mech_offset_[m] = acc;
            acc += mechanisms_[m]->num_constraints();
        }
    }

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
        }
    } else {
        // Unpack state variables
        T_init_ = statev(0);
        for (auto& mech : mechanisms_) {
            mech->unpack(statev);
        }
    }

    // Apply rotation for objectivity
    for (auto& mech : mechanisms_) {
        mech->rotate(DR);
    }

    // Save start values
    sigma_start_ = sigma;
    for (auto& mech : mechanisms_) {
        mech->to_start();
    }

    // Set elastic stiffness
    L = elasticity_.L();

    // Total number of constraints
    int n_total = 0;
    for (const auto& mech : mechanisms_) {
        n_total += mech->num_constraints();
    }

    // Perform return mapping. return_mapping reports FB convergence, but by the
    // reference CCP convention a finite unconverged-at-maxiter state is still
    // committed (the damage row in particular is integrated explicitly and
    // never drives the FB error to precision_). Only a NON-FINITE result
    // (true divergence -> NaN/Inf) is unusable and triggers a step cut.
    arma::vec Ds_total = arma::zeros(n_total);
    return_mapping(Etot, DEtot, sigma, T_init_, T, DT, DTime, ndi, Ds_total);
    if (!sigma.is_finite()) {
        // Ask the global solver to halve the increment and retry. statev is
        // left at its incoming values (pack_all is skipped) so the retry
        // restarts from the correct state; Lt is set elastic to keep the
        // rejected-step output well-defined.
        tnew_dt = 0.5;
        Lt = elasticity_.L();
        return;
    }

    // Compute consistent tangent
    Lt = L;  // Start with elastic stiffness
    compute_tangent(sigma, Ds_total, Lt, tangent_mode);

    // Compute work quantities
    // Elastic strain
    arma::vec E_inel = arma::zeros(6);
    for (const auto& mech : mechanisms_) {
        E_inel += mech->inelastic_strain();
    }
    arma::vec Eel = Etot + DEtot - elasticity_.alpha() * (T - T_init_) - E_inel;

    // Elastic work (recoverable)
    double Wm_r_new = 0.5 * arma::dot(sigma, Eel);
    Wm_r = Wm_r_new - 0.5 * arma::dot(sigma_start_, Eel);

    // Dissipated and stored work from mechanisms
    Wm_ir = 0.0;
    Wm_d = 0.0;
    for (const auto& mech : mechanisms_) {
        double Wm_r_m = 0.0, Wm_ir_m = 0.0, Wm_d_m = 0.0;
        mech->compute_work(sigma_start_, sigma, Wm_r_m, Wm_ir_m, Wm_d_m);
        Wm_ir += Wm_ir_m;
        Wm_d += Wm_d_m;
    }

    // Total mechanical work
    Wm = Wm_r + Wm_ir + Wm_d;

    // Pack state variables
    statev(0) = T_init_;
    for (const auto& mech : mechanisms_) {
        mech->pack(statev);
    }
}

void ModularUMAT::return_mapping(
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    double T_init,
    double T,
    double DT,
    double DTime,
    int ndi,
    arma::vec& Ds_total
) {
    // Total number of constraints
    int n_total = 0;
    for (const auto& mech : mechanisms_) {
        n_total += mech->num_constraints();
    }

    if (n_total == 0) {
        // Pure elastic: just compute stress
        arma::vec Eel = Etot + DEtot - elasticity_.alpha() * (T - T_init);
        sigma = el_pred(elasticity_.L(), Eel, ndi);
        return;
    }

    // Elastic prediction. stiffness_reduction() is the multiplicative CDM
    // (1-D) factor — without it the stress ignores damage entirely while the
    // tangent is softened (inconsistent Newton, unsoftened response).
    arma::vec E_inel = arma::zeros(6);
    for (const auto& mech : mechanisms_) {
        E_inel += mech->inelastic_strain();
    }
    arma::vec Eel = Etot + DEtot - elasticity_.alpha() * (T - T_init) - E_inel;
    sigma = stiffness_reduction() * el_pred(elasticity_.L(), Eel, ndi);

    // Allocate constraint arrays
    arma::vec Phi = arma::zeros(n_total);
    arma::vec Y_crit = arma::zeros(n_total);
    arma::mat B = arma::zeros(n_total, n_total);
    arma::vec ds = arma::zeros(n_total);
    Ds_total.zeros(n_total);

    // Iteration loop
    double error = 1.0;
    int iter = 0;

    // Scratch buffers for per-mechanism constraint returns (reused across FB
    // iterations — avoid reallocating arma::vecs inside the loop).
    arma::vec Phi_m, Y_crit_m;

    while (iter < maxiter_ && error > precision_) {
        // Phase 1: evaluate constraints per mechanism (populates mechanism
        //          caches used by dPhi_dsigma / kappa / K_cross in phase 2).
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            const int n = mechanisms_[m]->num_constraints();
            mechanisms_[m]->compute_constraints(
                sigma, Etot + DEtot, elasticity_.L(), DTime, Phi_m, Y_crit_m);
            Phi.subvec(mech_offset_[m], mech_offset_[m] + n - 1) = Phi_m;
            Y_crit.subvec(mech_offset_[m], mech_offset_[m] + n - 1) = Y_crit_m;
        }

        // Phases 2 and 3: local multiplier Jacobian from the mechanism caches.
        assemble_jacobian(sigma, DT, B);

        // Solve using Fischer-Burmeister
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_total, ds, error);

        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            mechanisms_[m]->update(ds, mech_offset_[m]);
        }

        // Recompute stress (D may have evolved in update, so re-evaluate the
        // reduction factor)
        E_inel.zeros();
        for (const auto& mech : mechanisms_) {
            E_inel += mech->inelastic_strain();
        }
        Eel = Etot + DEtot - elasticity_.alpha() * (T - T_init) - E_inel;
        sigma = stiffness_reduction() * el_pred(elasticity_.L(), Eel, ndi);

        ++iter;
    }
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
        kappa_all[jm] = &mechanisms_[jm]->kappa(sigma, DT, elasticity_.L());
    }
    for (size_t lm = 0; lm < mechanisms_.size(); ++lm) {
        const auto& dPhi_l_all = mechanisms_[lm]->dPhi_dsigma(sigma);
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
            sigma, elasticity_.L(), B, mech_offset_[m]);
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
    const arma::mat& L = elasticity_.L();
    Lt = L;

    // Mechanisms opting into the algorithmic assembly: stress-dependent Phi
    // AND an analytic flow Hessian. Others (viscoelastic, damage, Hessian-less
    // criteria) keep their continuum tangent_contribution in every mode.
    std::vector<size_t> algo;
    if (tangent_mode == 1) {
        for (size_t m = 0; m < mechanisms_.size(); ++m) {
            if (mechanisms_[m]->dLambda_dsigma(sigma) != nullptr &&
                !mechanisms_[m]->dPhi_dsigma(sigma).empty()) {
                algo.push_back(m);
            }
        }
    }

    if (!algo.empty()) {
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
            const auto& dPhi_all = mechanisms_[m]->dPhi_dsigma(sigma);
            const auto& kappa_all = mechanisms_[m]->kappa(sigma, 0.0, L);
            const auto& hess_all = *mechanisms_[m]->dLambda_dsigma(sigma);
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
            assemble_algorithmic_tangent(Bhat, kappa_j, dPhi_l, Ds_sub, L, dLambda_l);
        Lt = ct.Lt;
    }

    for (size_t m = 0; m < mechanisms_.size(); ++m) {
        if (std::find(algo.begin(), algo.end(), m) != algo.end()) {
            continue;
        }
        mechanisms_[m]->tangent_contribution(
            sigma, L, Ds_total, mech_offset_[m], Lt);
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
