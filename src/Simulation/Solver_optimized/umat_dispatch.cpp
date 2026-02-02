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

///@file umat_dispatch.cpp
///@brief Static UMAT dispatch singleton - single source of truth for all UMAT selection
///@version 2.0

#include <iostream>
#include <simcoon/Simulation/Solver_optimized/umat_dispatch.hpp>

// ============================================================================
// Mechanical small strain UMAT headers
// ============================================================================
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Ani_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/DFA_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Generic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/aniso_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono_cubic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>

// ============================================================================
// Finite strain UMAT headers
// ============================================================================
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>

// ============================================================================
// Thermomechanical UMAT headers
// ============================================================================
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Prony_Nfast.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

UmatDispatch::UmatDispatch() {
    // ========================================================================
    // Mechanical small strain dispatch map
    // Single source of truth - matches select_umat_M in umat_smart.cpp
    // ========================================================================
    dispatch_map_M_ = {
        {"UMEXT", 0},
        {"UMABA", 1},
        {"ELISO", 2},
        {"ELIST", 3},
        {"ELORT", 4},
        {"EPICP", 5},
        {"EPKCP", 6},
        {"EPCHA", 7},
        {"SMAUT", 8},
        {"SMANI", 9},
        {"LLDM0", 10},
        {"ZENER", 11},
        {"ZENNK", 12},
        {"PRONK", 13},
        {"EPHIL", 17},
        {"EPHAC", 18},
        {"EPANI", 19},
        {"EPDFA", 20},
        {"EPCHG", 21},
        {"EPHIN", 22},
        {"SMAMO", 23},
        {"SMAMC", 24},
        {"MIHEN", 100},
        {"MIMTN", 101},
        {"MISCN", 103},
        {"MIPLN", 104}
    };

    // ========================================================================
    // Mechanical finite strain dispatch map
    // Single source of truth - matches select_umat_M_finite in umat_smart.cpp
    // ========================================================================
    dispatch_map_M_finite_ = {
        {"UMEXT", 0},
        {"UMABA", 1},
        {"ELISO", 2},
        {"ELIST", 3},
        {"ELORT", 4},
        {"HYPOO", 5},
        {"EPICP", 6},
        {"EPKCP", 7},
        {"SNTVE", 8},
        {"NEOHI", 9},
        {"NEOHC", 10},
        {"MOORI", 11},
        {"YEOHH", 12},
        {"ISHAH", 13},
        {"GETHH", 14},
        {"SWANH", 15}
    };

    // ========================================================================
    // Thermomechanical dispatch map
    // Single source of truth - matches select_umat_T in umat_smart.cpp
    // ========================================================================
    dispatch_map_T_ = {
        {"UMEXT", 0},
        {"ELISO", 1},
        {"ELIST", 2},
        {"ELORT", 3},
        {"EPICP", 4},
        {"EPKCP", 5},
        {"ZENER", 6},
        {"ZENNK", 7},
        {"PRONK", 8},
        {"SMAUT", 9}
    };
}

UmatDispatch& UmatDispatch::instance() {
    // Thread-safe singleton via C++11 magic statics
    static UmatDispatch instance;
    return instance;
}

// ============================================================================
// Query functions
// ============================================================================

bool UmatDispatch::has_umat_M(const std::string& umat_name) const {
    return dispatch_map_M_.find(umat_name) != dispatch_map_M_.end();
}

bool UmatDispatch::has_umat_M_finite(const std::string& umat_name) const {
    return dispatch_map_M_finite_.find(umat_name) != dispatch_map_M_finite_.end();
}

bool UmatDispatch::has_umat_T(const std::string& umat_name) const {
    return dispatch_map_T_.find(umat_name) != dispatch_map_T_.end();
}

int UmatDispatch::get_umat_M_id(const std::string& umat_name) const {
    auto it = dispatch_map_M_.find(umat_name);
    return (it != dispatch_map_M_.end()) ? it->second : -1;
}

int UmatDispatch::get_umat_M_finite_id(const std::string& umat_name) const {
    auto it = dispatch_map_M_finite_.find(umat_name);
    return (it != dispatch_map_M_finite_.end()) ? it->second : -1;
}

int UmatDispatch::get_umat_T_id(const std::string& umat_name) const {
    auto it = dispatch_map_T_.find(umat_name);
    return (it != dispatch_map_T_.end()) ? it->second : -1;
}

// ============================================================================
// Mechanical small strain dispatch
// ============================================================================

bool UmatDispatch::call_umat_M(
    const std::string& umat_name,
    const arma::vec& Etot, const arma::vec& DEtot,
    arma::vec& sigma, arma::mat& Lt,
    arma::mat& L, arma::vec& sigma_in,
    const arma::mat& DR, int nprops, const arma::vec& props,
    int nstatev, arma::vec& statev,
    double T, double DT, double Time, double DTime,
    double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
    int ndi, int nshr, bool start, int solver_type, double& tnew_dt) {

    int umat_id = get_umat_M_id(umat_name);
    if (umat_id < 0) {
        return false;
    }

    switch (umat_id) {
        case 0: {
            // UMEXT - external UMAT (requires plugin, handled separately)
            cerr << "Error: External UMAT (UMEXT) requires plugin loading.\n";
            return false;
        }
        case 1: {
            // UMABA - Abaqus UMAT (requires plugin, handled separately)
            cerr << "Error: Abaqus UMAT (UMABA) requires plugin loading.\n";
            return false;
        }
        case 2: {
            umat_elasticity_iso(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 3: {
            umat_elasticity_trans_iso(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 4: {
            umat_elasticity_ortho(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 5: {
            umat_plasticity_iso_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 6: {
            umat_plasticity_kin_iso_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 7: {
            umat_plasticity_chaboche_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 8: {
            umat_sma_unified_T(Etot, DEtot, sigma, Lt, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 9: {
            umat_sma_aniso_T(Etot, DEtot, sigma, Lt, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 10: {
            umat_damage_LLD_0(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 11: {
            umat_zener_fast(Etot, DEtot, sigma, Lt, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 12: {
            umat_zener_Nfast(Etot, DEtot, sigma, Lt, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 13: {
            umat_prony_Nfast(Etot, DEtot, sigma, Lt, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 17: {
            umat_plasticity_hill_isoh_CCP(Etot, DEtot, sigma, Lt, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 18: {
            umat_hill_chaboche_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 19: {
            umat_ani_chaboche_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 20: {
            umat_dfa_chaboche_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 21: {
            umat_generic_chaboche_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 23: {
            umat_sma_mono(Etot, DEtot, sigma, Lt, L, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 24: {
            umat_sma_mono_cubic(Etot, DEtot, sigma, Lt, L, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 100: case 101: case 103: case 104: {
            // Multiphase UMATs - handled separately via umat_multi
            cerr << "Error: Multiphase UMATs require phase_characteristics context.\n";
            return false;
        }
        default: {
            cerr << "Error: UMAT '" << umat_name << "' not supported.\n";
            return false;
        }
    }

    return true;
}

// ============================================================================
// Mechanical finite strain dispatch
// ============================================================================

bool UmatDispatch::call_umat_M_finite(
    const std::string& umat_name,
    const arma::vec& etot, const arma::vec& Detot,
    const arma::mat& F0, const arma::mat& F1,
    arma::vec& sigma, arma::mat& Lt,
    arma::mat& L, arma::vec& sigma_in,
    const arma::mat& DR, int nprops, const arma::vec& props,
    int nstatev, arma::vec& statev,
    double T, double DT, double Time, double DTime,
    double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
    int ndi, int nshr, bool start, int solver_type, double& tnew_dt) {

    int umat_id = get_umat_M_finite_id(umat_name);
    if (umat_id < 0) {
        return false;
    }

    switch (umat_id) {
        case 0: {
            cerr << "Error: External UMAT (UMEXT) requires plugin loading.\n";
            return false;
        }
        case 1: {
            cerr << "Error: Abaqus UMAT (UMABA) requires plugin loading.\n";
            return false;
        }
        case 2: {
            umat_elasticity_iso(etot, Detot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 3: {
            umat_elasticity_trans_iso(etot, Detot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 4: {
            umat_elasticity_ortho(etot, Detot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 5: {
            umat_hypoelasticity_ortho(etot, Detot, F0, F1, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 6: {
            umat_plasticity_iso_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 7: {
            umat_plasticity_kin_iso_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 8: {
            umat_saint_venant(etot, Detot, F0, F1, sigma, Lt, L, sigma_in, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 9: {
            umat_neo_hookean_incomp(etot, Detot, F0, F1, sigma, Lt, L, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        case 10: case 11: case 12: case 13: case 14: case 15: {
            umat_generic_hyper_invariants(umat_name, etot, Detot, F0, F1, sigma, Lt, L, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
            break;
        }
        default: {
            cerr << "Error: Finite strain UMAT '" << umat_name << "' not supported.\n";
            return false;
        }
    }

    return true;
}

// ============================================================================
// Thermomechanical dispatch
// ============================================================================

bool UmatDispatch::call_umat_T(
    const std::string& umat_name,
    const arma::vec& Etot, const arma::vec& DEtot,
    arma::vec& sigma, double& r,
    arma::mat& dSdE, arma::mat& dSdT,
    arma::mat& drdE, arma::mat& drdT,
    const arma::mat& DR, int nprops, const arma::vec& props,
    int nstatev, arma::vec& statev,
    double T, double DT, double Time, double DTime,
    double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
    double& Wt0, double& Wt1, double& Wt2,
    int ndi, int nshr, bool start, double& tnew_dt) {

    int umat_id = get_umat_T_id(umat_name);
    if (umat_id < 0) {
        return false;
    }

    switch (umat_id) {
        case 0: {
            cerr << "Error: External UMAT (UMEXT) requires plugin loading.\n";
            return false;
        }
        case 1: {
            umat_elasticity_iso_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 2: {
            umat_elasticity_trans_iso_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 3: {
            umat_elasticity_ortho_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 4: {
            umat_plasticity_iso_CCP_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 5: {
            umat_plasticity_kin_iso_CCP_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 6: {
            umat_zener_fast_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 7: {
            umat_zener_Nfast_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 8: {
            umat_prony_Nfast_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        case 9: {
            umat_sma_unified_T_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                nprops, props, nstatev, statev, T, DT, Time, DTime,
                Wm, Wm_r, Wm_ir, Wm_d, Wt0, Wt1, Wt2, ndi, nshr, start, tnew_dt);
            break;
        }
        default: {
            cerr << "Error: Thermomechanical UMAT '" << umat_name << "' not supported.\n";
            return false;
        }
    }

    return true;
}

} // namespace simcoon
