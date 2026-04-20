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
#include <stdexcept>
#include <cmath>

namespace simcoon {

// ========== Constructor ==========

ModularUMAT::ModularUMAT()
    : elasticity_()
    , mechanisms_()
    , ivc_()
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
    // Register all internal variables
    // First: T_init (temperature reference)
    ivc_.add_scalar("T_init", 0.0, false);

    // Register variables from each mechanism
    for (auto& mech : mechanisms_) {
        mech->register_variables(ivc_);
    }

    // Compute offsets (T_init is at offset 0)
    ivc_.compute_offsets(0);

    // Check that we have enough state variables
    int required = ivc_.total_statev_size();
    if (nstatev < required) {
        throw std::runtime_error("ModularUMAT: nstatev (" + std::to_string(nstatev) +
                                ") < required (" + std::to_string(required) + ")");
    }

    // Unpack initial values from statev
    ivc_.unpack_all(statev);

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
                    count += 6 * vm->num_prony_terms();  // EV_i(6) per Prony term
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
    double& tnew_dt
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
            ivc_.get("T_init").scalar() = T_init_;
        }
    } else {
        // Unpack state variables
        ivc_.unpack_all(statev);
        T_init_ = ivc_.get("T_init").scalar();
    }

    // Apply rotation for objectivity
    ivc_.rotate_all(DR);

    // Save start values
    sigma_start_ = sigma;
    ivc_.to_start_all();

    // Set elastic stiffness
    L = elasticity_.L();

    // Total number of constraints
    int n_total = 0;
    for (const auto& mech : mechanisms_) {
        n_total += mech->num_constraints();
    }

    // Perform return mapping
    arma::vec Ds_total = arma::zeros(n_total);
    return_mapping(Etot, DEtot, sigma, T_init_, T, DT, DTime, ndi, Ds_total);

    // Compute consistent tangent
    Lt = L;  // Start with elastic stiffness
    compute_tangent(sigma, Ds_total, Lt);

    // Compute work quantities
    // Elastic strain
    arma::vec E_inel = arma::zeros(6);
    for (const auto& mech : mechanisms_) {
        E_inel += mech->inelastic_strain(ivc_);
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
        mech->compute_work(sigma_start_, sigma, ivc_, Wm_r_m, Wm_ir_m, Wm_d_m);
        Wm_ir += Wm_ir_m;
        Wm_d += Wm_d_m;
    }

    // Total mechanical work
    Wm = Wm_r + Wm_ir + Wm_d;

    // Pack state variables
    ivc_.pack_all(statev);
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

    // Elastic prediction
    arma::vec E_inel = arma::zeros(6);
    for (const auto& mech : mechanisms_) {
        E_inel += mech->inelastic_strain(ivc_);
    }
    arma::vec Eel = Etot + DEtot - elasticity_.alpha() * (T - T_init) - E_inel;
    sigma = el_pred(elasticity_.L(), Eel, ndi);

    // Allocate constraint arrays
    arma::vec Phi = arma::zeros(n_total);
    arma::vec Y_crit = arma::zeros(n_total);
    arma::mat B = arma::zeros(n_total, n_total);
    arma::vec ds = arma::zeros(n_total);
    Ds_total.zeros(n_total);

    // Iteration loop
    double error = 1.0;
    int iter = 0;

    while (iter < maxiter_ && error > precision_) {
        int offset = 0;

        // Compute constraints and Jacobian from all mechanisms
        for (const auto& mech : mechanisms_) {
            int n = mech->num_constraints();

            // Constraint functions
            arma::vec Phi_m, Y_crit_m;
            mech->compute_constraints(sigma, elasticity_.L(), DTime, ivc_, Phi_m, Y_crit_m);

            Phi.subvec(offset, offset + n - 1) = Phi_m;
            Y_crit.subvec(offset, offset + n - 1) = Y_crit_m;

            // Jacobian contribution
            mech->compute_jacobian_contribution(sigma, elasticity_.L(), ivc_, B, offset);

            offset += n;
        }

        // Solve using Fischer-Burmeister
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_total, ds, error);

        // Update all mechanisms
        offset = 0;
        for (auto& mech : mechanisms_) {
            mech->update(ds, offset, ivc_);
            offset += mech->num_constraints();
        }

        // Recompute stress
        E_inel.zeros();
        for (const auto& mech : mechanisms_) {
            E_inel += mech->inelastic_strain(ivc_);
        }
        Eel = Etot + DEtot - elasticity_.alpha() * (T - T_init) - E_inel;
        sigma = el_pred(elasticity_.L(), Eel, ndi);

        ++iter;
    }
}

void ModularUMAT::compute_tangent(
    const arma::vec& sigma,
    const arma::vec& Ds_total,
    arma::mat& Lt
) {
    // Start with elastic stiffness
    Lt = elasticity_.L();

    // Subtract plastic/viscous/damage contributions
    int offset = 0;
    for (const auto& mech : mechanisms_) {
        mech->tangent_contribution(sigma, elasticity_.L(), Ds_total, offset, ivc_, Lt);
        offset += mech->num_constraints();
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
    double& tnew_dt
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
        ndi, nshr, start, tnew_dt
    );
}

} // namespace simcoon
