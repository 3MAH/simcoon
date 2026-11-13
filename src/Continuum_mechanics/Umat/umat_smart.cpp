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

///@file umat_smart.cpp
///@brief Selection of constitutive laws and transfer between Abaqus and simcoon formats,
///       implemented in 1D-2D-3D
///@version 1.0

#include <iostream>
#include <map>
#include <string>
#include <stdexcept>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <armadillo>
#include <memory>
#include <dylib.hpp>
#include <filesystem>

#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/condensed_state.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>

#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_pstretch.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_callback.hpp>

#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Prony_Nfast.hpp>

#include <simcoon/Continuum_mechanics//Umat/Thermomechanical/SMA/unified_T.hpp>

#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>

using namespace std;
using namespace arma;
namespace fs = std::filesystem;
namespace simcoon{


void size_statev(phase_characteristics &rve, unsigned int &size) {

    for (auto &r:rve.sub_phases) {   // by reference: a copy would duplicate the whole sub-tree
        size = size + r.sptr_sv_local->nstatev + 57;
        size_statev(r,size);
    }
}

void statev_2_phases(phase_characteristics &rve, unsigned int &pos, const vec &statev) {

    for(auto &r : rve.sub_phases) {
        shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        unsigned int nstatev = umat_phase_M->statev.n_elem;

        vec vide = zeros(6);
        umat_phase_M->Etot = statev.subvec(pos,size(vide));
        umat_phase_M->DEtot = statev.subvec(pos+6,size(vide));
        umat_phase_M->sigma = statev.subvec(pos+12,size(vide));
        umat_phase_M->T = statev(pos+19);
        umat_phase_M->DT = statev(pos+20);

        for (int i=0; i<6; i++) {
            umat_phase_M->Lt.col(i) = statev.subvec(pos+21+i*6,size(vide));
        }
        
        umat_phase_M->statev = statev.subvec(pos+57,size(umat_phase_M->statev));
        pos+=57+nstatev;
        statev_2_phases(r,pos,statev);
    }

}
    
void phases_2_statev(vec &statev, unsigned int &pos, const phase_characteristics &rve) {
    
    for(const auto &r : rve.sub_phases) {
        shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        unsigned int nstatev = umat_phase_M->statev.n_elem;
        
        vec vide = zeros(6);
        statev.subvec(pos,size(vide)) = umat_phase_M->Etot;
        statev.subvec(pos+6,size(vide)) = umat_phase_M->DEtot;
        statev.subvec(pos+12,size(vide)) = umat_phase_M->sigma;
        statev(pos+19) = umat_phase_M->T;
        statev(pos+20) = umat_phase_M->DT;
        for (int i=0; i<6; i++) {
            statev.subvec(pos+21+i*6,size(vide)) = umat_phase_M->Lt.col(i);
        }
        statev.subvec(pos+57,size(umat_phase_M->statev)) = umat_phase_M->statev;
        
        pos+=57+nstatev;
        phases_2_statev(statev,pos,r);
    }
    
}

bool is_external_plugin(const std::string &umat_name)
{
    const umat_entry *entry = umat_entry_of(umat_name);
    if (entry == nullptr) return false;
    const UmatKernel k = entry->mechanical;
    return k == UmatKernel::plugin_ext || k == UmatKernel::plugin_abaqus || k == UmatKernel::plugin_ansys;
}

namespace {

// the row of a name on the finite-strain route, nullptr when that route does not serve it
const umat_entry *finite_entry_of(const std::string &umat_name)
{
    const umat_entry *entry = umat_entry_of(umat_name);
    return (entry != nullptr && entry->finite != UmatKernel::none) ? entry : nullptr;
}

}  // namespace

umat_convention output_convention_of(const std::string &umat_name)
{
    const umat_entry *entry = finite_entry_of(umat_name);
    if (entry == nullptr) {
        throw std::invalid_argument(
            "output_convention_of: the finite-strain route does not serve '" + umat_name + "' (the registry of "
            "umat_registry.cpp has no finite kernel for it): a missing stress measure is an error of exactly J, "
            "a missing tangent rate is a wrong rate, and both are silent.");
    }
    return entry->convention;
}

bool stress_output_is_kirchhoff(const std::string &umat_name)
{
    const umat_entry *entry = finite_entry_of(umat_name);
    return entry != nullptr && entry->convention.stress == StressMeasure::kirchhoff;
}

bool total_formulation_serves(const std::string &umat_name)
{
    return umat_name != "HYPOO" && stress_output_is_kirchhoff(umat_name);
}

bool material_tangent_serves(const std::string &umat_name)
{
    const umat_entry *entry = finite_entry_of(umat_name);
    return entry != nullptr && entry->convention.material_tangent;
}

void check_ncp_serves(const std::string &umat_name, const int &ncp, const bool &thermomechanical)
{
    if (ncp_method_of(ncp, "ncp") == NCPMethod::fischer_burmeister) return;
    const umat_entry *entry = umat_entry_of(umat_name);
    const bool serves = entry != nullptr
        && (thermomechanical ? kernel_reads_ncp(entry->thermomechanical, true)
                             : (kernel_reads_ncp(entry->mechanical, false) || kernel_reads_ncp(entry->finite, false)));
    if (!serves) {
        throw std::invalid_argument("ncp='min' serves the cutting planes of the modular engine (MODUL and the legacy "
                                    "names it serves, mechanical and thermomechanical) and of EPICP_closed; '"
                                    + umat_name + "'" + (thermomechanical ? " (thermomechanical)" : "") + " does not read it");
    }
}

void check_integrator_serves(const std::string &umat_name, const int &tangent_mode, const bool &thermomechanical)
{
    if (tangent_mode != tangent_closest_point) return;
    const umat_entry *entry = umat_entry_of(umat_name);
    const bool serves = entry != nullptr
        && (thermomechanical ? kernel_serves_closest_point(entry->thermomechanical, true)
                             : (kernel_serves_closest_point(entry->mechanical, false)
                                || kernel_serves_closest_point(entry->finite, false)));
    if (!serves) {
        throw std::invalid_argument("integrator='closest_point' serves the laws with a closest-point update only "
                                    "(mechanical: EPICP_closed, EPCHA_closed, EPJCK, MODUL and the legacy names it serves, "
                                    "the SMA kernels; thermomechanical: EPICP_closed, EPKCP_closed, MODUL and the legacy "
                                    "names it serves); '" + umat_name + "'"
                                    + (thermomechanical ? " (thermomechanical)" : "") + " has none");
    }
}

int kernel_tangent_mode(const int &integrator, const int &tangent_mode)
{
    // Migration message of 2.x, removed in 3.0 (CHANGELOG.md, "Planned for 3.0"): 3 then falls in the range check
    if (tangent_mode == tangent_closest_point) {
        throw std::invalid_argument("tangent_mode 3 (closest point) is no longer a tangent choice: give "
                                    "integrator='closest_point' with tangent_mode='algorithmic'");
    }
    if (tangent_mode < tangent_none || tangent_mode > tangent_algorithmic) {
        throw std::invalid_argument("tangent_mode must be 0 (none), 1 (continuum) or 2 (algorithmic); got "
                                    + std::to_string(tangent_mode));
    }
    switch (integrator) {
        case integrator_cutting_plane:
            return tangent_mode;
        case integrator_closest_point:
            if (tangent_mode == tangent_algorithmic) return tangent_closest_point;
            throw std::invalid_argument("integrator='closest_point' returns its exact (algorithmic) operator only: "
                                        "give tangent_mode='algorithmic'");
        default:
            throw std::invalid_argument("integrator must be 'cutting_plane' or 'closest_point'; got code "
                                        + std::to_string(integrator));
    }
}

bool total_formulation_serves_T(const std::string &umat_name)
{
    // every thermomechanical kernel is a log-strain box: the face is the admission
    const umat_entry *entry = umat_entry_of(umat_name);
    return entry != nullptr && entry->thermomechanical != UmatKernel::none;
}

vec total_formulation_tau(const vec &sigma, const vec &Detot, const mat &F1, const bool &box, double &dW)
{
    dW = 0.;
    if (box)
        return t2v_stress(LogStressU2Kirchoff(v2t_stress(sigma), F1));
    const vec T1 = t2v_stress(Kirchoff2LogStressU(v2t_stress(sigma), F1));
    dW = 0.5*sum((T1 - sigma)%Detot);
    return sigma;
}

namespace {

// One call of the thermomechanical kernel, fed the strain pair of the route: (Etot, DEtot) under
// small strain, the box strain (etot, Detot) under finite strain
void call_umat_T(const UmatKernel &kernel, phase_characteristics &rve, state_variables_T &sv, const vec &strain, const vec &Dstrain, const mat &DR, const KinematicForm &form, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, double &tnew_dt)
{
    switch (kernel) {
        case UmatKernel::elasticity_iso: {
            umat_elasticity_iso_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::elasticity_trans_iso: {
            umat_elasticity_trans_iso_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::elasticity_ortho: {
            umat_elasticity_ortho_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::plasticity_iso: {
            umat_plasticity_iso_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::johnson_cook: {
            umat_plasticity_johnson_cook_CCP_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::modular: {
            umat_modular_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode, form, sv.ncp);
            break;
        }
        case UmatKernel::legacy_modular: {   // legacy names served by the modular engine (legacy_adapters.cpp)
            umat_legacy_modular_T(rve.sptr_matprops->umat_name, strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode, form, sv.ncp);
            break;
        }
        case UmatKernel::plasticity_kin_iso: {
            umat_plasticity_kin_iso_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
           break;
        }
        case UmatKernel::zener: {
            umat_zener_fast_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::zener_N: {
            umat_zener_Nfast_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::prony_N: {
            umat_prony_Nfast_T(strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        case UmatKernel::sma_unified: {
            umat_sma_unified_T_T(rve.sptr_matprops->umat_name, strain, Dstrain, sv.sigma, sv.r, sv.dSdE, sv.dSdT, sv.drdE, sv.drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, sv.nstatev, sv.statev, sv.T, sv.DT, Time, DTime, sv.Wm(0), sv.Wm(1), sv.Wm(2), sv.Wm(3), sv.Wt(0), sv.Wt(1), sv.Wt(2), ndi, nshr, start, tnew_dt, sv.tangent_mode);
            break;
        }
        default: {
            throw std::logic_error("call_umat_T: no kernel case for '" + rve.sptr_matprops->umat_name
                                   + "' (the registry and the switch disagree)");
        }
    }
}

}  // namespace

static void select_umat_T_3D(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    UNUSED(solver_type);

    // caller's DR is global: rotated with the state so the local UMAT gets the material-frame increment
    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_local);
    const umat_entry *entry = umat_entry_of(rve.sptr_matprops->umat_name);
    if (entry == nullptr || entry->thermomechanical == UmatKernel::none) {
        throw std::invalid_argument("Unknown umat name in the thermomechanical dispatch: " + rve.sptr_matprops->umat_name);
    }
    call_umat_T(entry->thermomechanical, rve, *umat_T, umat_T->Etot, umat_T->DEtot, umat_T->DR, KinematicForm::rate, Time, DTime, ndi, nshr, start, tnew_dt);
    // small strain (F = I): tau = PKII = sigma, else the Kirchhoff-route consumers read an empty tau
    umat_T->tau = umat_T->sigma;
    umat_T->PKII = umat_T->sigma;
    rve.local2global();
    
}

// The instance of an external plugin API (one library per API type), loaded on first use and kept
// until exit: the library before the object, so the object is destroyed while its code is loaded
template <class Api>
static Api &plugin(const char *path)
{
    static dylib::library lib(path, dylib::decorations::os_default());
    using destroy_fn = void(Api*);
    static std::unique_ptr<Api, destroy_fn*> api(lib.get_function<Api*()>("create_api")(),
                                                 lib.get_function<destroy_fn>("destroy_api"));
    return *api;
}

namespace {

// transport of a statev tensor by the rotation-free M: Truesdell (4) = tensor2 push-forward (no Piola
// factor), log_F (5) = symmetrised M X M^-1; same rules as the transport of tau and etot
vec transport_convected(const vec &v, const mat &M, const int &corate_type, const Tensor2Type &type)
{
    const tensor2 t = tensor2::from_voigt(v, type);
    if (corate_type == 4)
        return t.push_forward(M, false).to_arma_voigt();
    mat M_inv;
    if (!inv(M_inv, M))
        throw simcoon::exception_inv("transport_convected: the frame-relative stretch is not invertible");
    const mat X = M*t.to_arma_mat()*M_inv;
    return tensor2(mat(0.5*(X + X.t())), type).to_arma_voigt();
}


// the work around a box kernel on the finite route (select_umat_M_finite, select_umat_T_finite)
struct finite_box {
    vec etot_stored;     // set_start transports it at commit, with the increment it belongs to
    vec Detot_stored;
    vec tau_start_lab;   // committed, untransported: for the work correction
    vec tau_start_tr;    // the start stress the kernel integrates from
    umat_convention conv{StressMeasure::kirchhoff};
    bool total = false;
    mat R_hat;
    mat DR_kernel;
};

finite_box enter_finite_box(state_variables &sv, const umat_convention &conv, const std::string &umat_name, const int &corate_type)
{
    finite_box box;
    box.conv = conv;
    box.total = (corate_type == corate_type_total);
    const mat &DR = sv.DR;

    // transport the start state with THIS increment's DR: strain, start stress and statev share one configuration
    box.etot_stored = sv.etot;
    box.tau_start_lab = sv.sigma;
    if (corate_type != 4 && corate_type != 5) {
        sv.etot = rotate_strain(box.etot_stored, DR);
        sv.sigma = rotate_stress(sv.sigma, DR);
    }
    else if (corate_type == 4) {   // Truesdell, DR = DF: strain lower-, Kirchhoff stress upper-convected
        mat DR_inv;
        if (!inv(DR_inv, DR))
            throw simcoon::exception_inv("finite route: DF is not invertible");
        sv.etot = t2v_strain(DR_inv.t()*v2t_strain(box.etot_stored)*DR_inv);
        sv.sigma = t2v_stress(DR*v2t_stress(sv.sigma)*DR.t());
    }
    else {   // log_F, DR = DF: similarity transport, as in set_start
        mat DR_inv;
        if (!inv(DR_inv, DR))
            throw simcoon::exception_inv("finite route: DF is not invertible");
        sv.etot = t2v_strain(DR*v2t_strain(box.etot_stored)*DR_inv);
        sv.sigma = t2v_stress(DR*v2t_stress(sv.sigma)*DR_inv);
    }
    box.tau_start_tr = sv.sigma;
    if (conv.stress == StressMeasure::cauchy)   // a host kernel integrates from its own measure
        sv.sigma = Kirchoff2Cauchy(sv.sigma, sv.F0);

    if (box.total)
        sv.sigma = t2v_stress(Kirchoff2LogStressU(v2t_stress(sv.sigma), sv.F0));

    // box kernel frame R_hat: increment I for corates 0-3, stretch M (applied to the declared statev) for 4/5
    box.Detot_stored = sv.Detot;
    if (corate_type == 4) {
        // Truesdell: the kernel's increment is the closed-form Almansi increment of DR = DF,
        // whatever the solver's control variable is (logarithmic control increments ln V).
        sv.Detot = t2v_strain(Euler_Almansi(DR));
    }
    box.R_hat = eye(3,3);
    box.DR_kernel = DR;
    if (conv.material_frame) {
        const bool convected = (corate_type == 4 || corate_type == 5);
        mat R_hat_n;
        if (!convected) {
            R_hat_n = sv.R;
            box.R_hat = DR*sv.R;
        }
        else {
            mat U;
            RU_decomposition(R_hat_n, U, sv.F0);
            RU_decomposition(box.R_hat, U, sv.F1);
        }
        sv.etot = rotate_strain(sv.etot, box.R_hat.t());
        sv.Detot = rotate_strain(sv.Detot, box.R_hat.t());   // the Truesdell override, if any
        sv.sigma = rotate_stress(sv.sigma, box.R_hat.t());
        box.DR_kernel = eye(3,3);
        if (convected) {
            if (!conv.layout_declared) {
                throw simcoon::exception_solver(
                    "corate_type " + std::to_string(corate_type) + " requires the tensorial "
                    "internal variables of '" + umat_name + "' to be declared "
                    "(umat_entry::convention, umat_registry.cpp); use corate_type 3");
            }
            const mat M = box.R_hat.t()*DR*R_hat_n;
            for (const StatevTensor &t : conv.statev_tensors) {
                if (t.offset + 6 > static_cast<int>(sv.statev.n_elem)) {
                    throw simcoon::exception_solver(
                        "'" + umat_name + "' declares a statev tensor at offset "
                        + std::to_string(t.offset) + ", past nstatev = " + std::to_string(sv.statev.n_elem));
                }
                const vec x = sv.statev.subvec(t.offset, t.offset + 5);
                sv.statev.subvec(t.offset, t.offset + 5) = transport_convected(x, M, corate_type, t.type);
            }
        }
    }
    return box;
}

void leave_finite_box(state_variables &sv, const finite_box &box, const int &corate_type, vec &Wm)
{
    // Back from the material frame (statev stays in it); the caller rotates its own tangents
    if (box.conv.material_frame)
        sv.sigma = rotate_stress(sv.sigma, box.R_hat);

    // tau is the route stress: pass-through, converted for a kernel declaring cauchy
    if (box.total) {
        double dW_T;
        sv.tau = total_formulation_tau(sv.sigma, sv.Detot, sv.F1, box.conv.material_frame, dW_T);
        Wm(0) += dW_T;
        Wm(1) += dW_T;
        sv.sigma = sv.tau;
    }
    else if (box.conv.stress == StressMeasure::kirchhoff)
        sv.tau = sv.sigma;                                                      // the kernel's output IS the Kirchhoff stress
    else {
        sv.tau = t2v_stress(Cauchy2Kirchoff(v2t_stress(sv.sigma), sv.F1));  // foreign Cauchy -> Kirchhoff
        sv.sigma = sv.tau;                                                      // the state travels in Kirchhoff
    }
    sv.PKII = t2v_stress(Kirchoff2PKII(v2t_stress(sv.tau), sv.F1));

    // Log corates: replace the kernel's work by the stress power, in one frame.
    const double dW = Delta_work_conjugacy(box.tau_start_lab, box.tau_start_tr, sv.tau, box.Detot_stored,
                                           sv.F0, sv.F1, corate_type);
    Wm(0) += dW;
    Wm(1) += dW;

    sv.etot = box.etot_stored;
    sv.Detot = box.Detot_stored;
}

// The host plugins (Abaqus, Ansys) read the Etot/DEtot slots of the phase: hand them the box strain
// there and put the Green-Lagrange pair back (it travels back to the global state)
void call_host_plugin(const UmatKernel &kernel, phase_characteristics &rve, state_variables_M &sv, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    const vec::fixed<6> Etot_GL = sv.Etot, DEtot_GL = sv.DEtot;
    sv.Etot = sv.etot;
    sv.DEtot = sv.Detot;
    try {
        if (kernel == UmatKernel::plugin_abaqus)
            plugin<umat_plugin_aba_api>("external/umat_plugin_aba").umat_abaqus(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
        else
            plugin<umat_plugin_ans_api>("external/umat_plugin_ans").umat_ansys(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
    } catch (...) {
        sv.Etot = Etot_GL;
        sv.DEtot = DEtot_GL;
        throw;
    }
    sv.Etot = Etot_GL;
    sv.DEtot = DEtot_GL;
}

}  // namespace

void select_umat_M_finite(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const int &corate_type, double &tnew_dt, const bool &material_tangent)
{
    // The finite dispatch runs 3D: a condensed state is solved around it (run_umat_M, select_finite_condensed)
    if (ndi < 3) {
        throw std::invalid_argument("select_umat_M_finite: condensed stress states (ndi < 3) are solved around the finite dispatch (run_umat_M)");
    }
    const std::string &law = rve.sptr_matprops->umat_name;
    const umat_entry *entry = umat_entry_of(law);
    if (entry == nullptr || entry->finite == UmatKernel::none) {
        if (entry != nullptr) {
            throw std::invalid_argument("'" + law + "' serves the small-strain route only (no finite-strain kernel)");
        }
        throw std::invalid_argument("Unknown umat name in the finite-strain dispatch: " + law);
    }
    const UmatKernel kernel = entry->finite;

    // The caller provides DR in global coordinates; rotate it with the other
    // state variables so that the local UMAT receives a consistent increment.
    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);

    // one declaration per kernel: output_convention_of throws for an undeclared one
    const umat_convention conv = output_convention_of(law);
    // The corate the kernel's tangent is expressed in: the solver's, or the material tangent
    // dS/dE when asked for and declared (a kernel built on F: the corate shapes nothing else)
    const int corate_Lt = (material_tangent && conv.material_tangent) ? corate_type_material_tangent : corate_type;

    const std::string &law_name = rve.sptr_matprops->umat_name;
    if (kernel == UmatKernel::modular || (kernel == UmatKernel::legacy_modular && adapter_is_hyperelastic(law_name))) {
        // tau = d psi / d eps_el is a stored-energy law only when the accumulated strain is ln V:
        // corate 3, or 2 in closed form. Only 3 is accepted (policy: the production rate).
        if (corate_type != 3 && corate_type != corate_type_total) {
            const bool closed = umat_entry_of(law_name + "_closed") != nullptr;
            throw simcoon::exception_solver(
                (kernel == UmatKernel::modular ? std::string("MODUL") : "'" + law_name + "' (modular engine)")
                + " under finite strain requires corate_type = 3 (log_R) or the total "
                "formulation: the modular composition is hyperelastic in the logarithmic "
                "strain; got corate_type = " + std::to_string(corate_type)
                + (closed ? "; " + law_name + "_closed, built on F, serves every corate" : ""));
        }
    }

    // Total formulation: admission and stress contract in total_formulation_serves().
    const bool total = (corate_type == corate_type_total);
    if (total && !total_formulation_serves(rve.sptr_matprops->umat_name)) {
        throw simcoon::exception_solver(
            "the total formulation (strain space) does not serve '"
            + rve.sptr_matprops->umat_name + "'; use the rate formulation");
    }

    if (kernel == UmatKernel::plasticity_chaboche && corate_type == 4) {   // EPCHA_closed
        // X_i = 2/3 C_i a_i is stored: Truesdell convects a (strain) and X (stress) differently,
        // so the pair would drift apart. log_F transports both by the same similarity.
        throw simcoon::exception_solver(
            "EPCHA_closed stores its back stresses X_i = 2/3 C_i a_i, which the Truesdell rate "
            "(corate_type 4) cannot keep consistent; use corate_type 3 (or 5)");
    }

    const finite_box box = enter_finite_box(*umat_M, conv, rve.sptr_matprops->umat_name, corate_type);
    const mat &DR_kernel = box.DR_kernel;
    const vec::fixed<4> Wm_before = umat_M->Wm;

    switch (kernel) {

            case UmatKernel::hypoelasticity_ortho: {
                umat_hypoelasticity_ortho(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_type, umat_M->tangent_mode);
                break;
            }
             case UmatKernel::plasticity_iso: {
                umat_plasticity_iso(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode, umat_M->ncp);
                 break;
             }
            case UmatKernel::johnson_cook: {
                // Johnson-Cook on the log-strain/Kirchhoff box route, same as EPICP
                umat_plasticity_johnson_cook_CCP(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                break;
            }
             case UmatKernel::plasticity_chaboche: {
                // Chaboche on the log-strain/Kirchhoff box route, same as EPICP
                umat_plasticity_chaboche(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                 break;
             }
            case UmatKernel::saint_venant: {
                umat_saint_venant(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_Lt, umat_M->tangent_mode);
                break;
            }
            case UmatKernel::legacy_modular: {
                // Legacy names served by the modular engine on the log-strain box, like EPICP:
                // run in the material frame, so their anisotropy axes follow the body.
                umat_legacy_modular(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode, umat_M->ncp, total ? KinematicForm::total : KinematicForm::rate);
                break;
            }
            case UmatKernel::modular: {
                // MODUL under NLGEOM is a Hencky hyperelastic composition:
                umat_modular(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode, total ? KinematicForm::total : KinematicForm::rate, umat_M->ncp);
                break;
            }
            case UmatKernel::neo_hookean_incomp: {
                umat_neo_hookean_incomp(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_Lt, umat_M->tangent_mode);
                break;
            }                         
            case UmatKernel::hyper_invariants: {
                umat_generic_hyper_invariants(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_Lt, umat_M->tangent_mode);
                break;
            }
            case UmatKernel::hyper_pstretch: {
                umat_generic_hyper_pstretch(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_Lt, umat_M->tangent_mode);
                break;
            }
            case UmatKernel::python_callback: {
                // PYEXT: process-wide callback UMAT (umat_callback.hpp; registered by the Python
                // bindings). Small-strain convention on the log-strain / Kirchhoff box, as EPICP.
                umat_callback_M(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                break;
            }
            case UmatKernel::plugin_ext: {
                // UMEXT: the 24-argument convention on the log-strain / Kirchhoff box, as PYEXT
                plugin<umat_plugin_ext_api>("external/umat_plugin_ext").umat_external_M(law, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case UmatKernel::plugin_abaqus:
            case UmatKernel::plugin_ansys: {
                call_host_plugin(kernel, rve, *umat_M, DR_kernel, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
                break;
            }
            default: {
                throw std::logic_error("select_umat_M_finite: no kernel case for '" + law
                                       + "' (the registry and the switch disagree)");
            }
        }
    
        // Back from the material frame: the box tangent, in the kernel's frame, to the lab
        if (conv.material_frame) {
            umat_M->Lt = rotate_stiffness(umat_M->Lt, box.R_hat);
            umat_M->L = rotate_stiffness(umat_M->L, box.R_hat);
        }
        const vec::fixed<4> Wm_host = umat_M->Wm;
        leave_finite_box(*umat_M, box, corate_type, umat_M->Wm);
        if (conv.stress == StressMeasure::cauchy) {
            // A host plugin answers in its own conventions, lab frame here: Cauchy stress; Abaqus DDSDDE =
            // (1/J) d(J sigma)/d eps in the Jaumann rate, Ansys dsdePl = d sigma/d eps in the Green-Naghdi rate
            // (its sigma (x) I term added, symmetric part, as the Abaqus Jacobian carries it). Tangent -> dS/dE ->
            // the box of the solver's corate; the host's work trapezoid (Cauchy) replaced by the Kirchhoff one,
            // the conjugacy correction kept.
            const vec sigma_cv = Kirchoff2Cauchy(umat_M->tau, umat_M->F1);
            const mat sigma_c = v2t_stress(sigma_cv);
            const mat DSDE = (kernel == UmatKernel::plugin_abaqus)
                ? DsigmaDe_JaumannDD_2_DSDE(umat_M->Lt, umat_M->F1, sigma_c)
                : DsigmaDe_GreenNaghdiDD_2_DSDE(abaqus_jacobian(umat_M->Lt, sigma_cv), umat_M->F1, sigma_c);
            umat_M->Lt = DSDE_2_DtauDe_corate(DSDE, corate_type, umat_M->F1, v2t_stress(umat_M->tau));
            const double trap = 0.5*dot(box.tau_start_tr + umat_M->tau, box.Detot_stored);
            umat_M->Wm(0) += trap - (Wm_host(0) - Wm_before(0));
            umat_M->Wm(1) += trap - (Wm_host(1) - Wm_before(1));
        }
        // Otherwise no tangent conversion: every kernel emits Lt in corate_type already (in dS/dE when
        // the material tangent was asked for and served).
        rve.local2global();
}
    
static void select_umat_T_finite(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const int &corate_type, double &tnew_dt)
{
    UNUSED(solver_type);
    if (ndi < 3) {
        throw std::invalid_argument("select_umat_T_finite: condensed stress states (ndi < 3) are solved around the finite dispatch (run_umat_T)");
    }
    // log_R (3) and total only: under 4/5 the stretched start stress shifts eta_start (spurious heat),
    // 0-2 are not the production rate
    if (corate_type != 3 && corate_type != corate_type_total) {
        throw simcoon::exception_solver(
            "thermomechanical finite strain serves corate_type 3 (log_R) or the total "
            "formulation; got corate_type = " + std::to_string(corate_type));
    }
    const umat_entry *entry = umat_entry_of(rve.sptr_matprops->umat_name);
    if (entry == nullptr || entry->thermomechanical == UmatKernel::none) {
        throw std::invalid_argument("Unknown umat name in the thermomechanical dispatch: " + rve.sptr_matprops->umat_name);
    }
    if (corate_type == corate_type_total && !total_formulation_serves_T(rve.sptr_matprops->umat_name)) {
        throw simcoon::exception_solver(
            "the total formulation (strain space) does not serve the thermomechanical '"
            + rve.sptr_matprops->umat_name + "'; use the rate formulation");
    }

    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_local);

    // thermo kernels: Kirchhoff boxes in the material frame, no statev layout (corates 4/5 refused above)
    const umat_convention conv{StressMeasure::kirchhoff, true, false, {}};
    const finite_box box = enter_finite_box(*umat_T, conv, rve.sptr_matprops->umat_name, corate_type);
    call_umat_T(entry->thermomechanical, rve, *umat_T, umat_T->etot, umat_T->Detot, box.DR_kernel, box.total ? KinematicForm::total : KinematicForm::rate, Time, DTime, ndi, nshr, start, tnew_dt);

    // Back from the material frame: dSdT and drdE rotate with the stress operator (as in rotate_l2g)
    umat_T->dSdE = rotate_stiffness(umat_T->dSdE, box.R_hat);
    umat_T->dSdT = rotate_stress(umat_T->dSdT, box.R_hat);
    umat_T->drdE = rotate_stress(umat_T->drdE, box.R_hat);
    leave_finite_box(*umat_T, box, corate_type, umat_T->Wm);
    // total form: dtau/dtheta at fixed F is the map of dT/dtheta (linear at fixed F); dSdE and
    // drdE stay box derivatives, r per reference volume
    if (box.total)
        umat_T->dSdT = t2v_stress(LogStressU2Kirchoff(v2t_stress(umat_T->dSdT), umat_T->F1));
    rve.local2global();
}
    
static void select_umat_M_3D(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{

    const umat_entry *entry = umat_entry_of(rve.sptr_matprops->umat_name);
    if (entry == nullptr || entry->mechanical == UmatKernel::none) {
        throw std::invalid_argument("Unknown umat name in the small-strain dispatch: " + rve.sptr_matprops->umat_name);
    }

    // caller's DR is global: rotated with the state so the local UMAT gets the material-frame increment
    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    const mat &DR = umat_M->DR;

    switch (entry->mechanical) {

        case UmatKernel::plugin_ext: {
            // plugin ABI predates tangent_mode: never forward it (symbol mismatch); bump umat_plugin_ext_api to add it
            plugin<umat_plugin_ext_api>("external/umat_plugin_ext").umat_external_M(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case UmatKernel::plugin_abaqus: {
            plugin<umat_plugin_aba_api>("external/umat_plugin_aba").umat_abaqus(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case UmatKernel::plugin_ansys: {
            plugin<umat_plugin_ans_api>("external/umat_plugin_ans").umat_ansys(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case UmatKernel::plasticity_iso: {
            umat_plasticity_iso(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode, umat_M->ncp);
            break;
        }
        case UmatKernel::johnson_cook: {
            umat_plasticity_johnson_cook_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::plasticity_chaboche: {
            umat_plasticity_chaboche(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::sma_unified: {
            umat_sma_unified_T(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::sma_unified_R: {
            umat_sma_unified_TR(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::damage_LLD0: {
            umat_damage_LLD_0(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::zener: {
            umat_zener_fast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::zener_N: {
            umat_zener_Nfast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::prony_N: {
            umat_prony_Nfast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::sma_mono: {
            // SMAMO (isotropic) and SMAMC (cubic) both use unified umat_sma_mono
            umat_sma_mono(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case UmatKernel::modular: {
            // Modular UMAT - composable constitutive model
            umat_modular(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode, KinematicForm::rate, umat_M->ncp);
            break;
        }
        case UmatKernel::legacy_modular: {
            // legacy names on the modular engine: props translated per name (legacy_adapters.cpp)
            umat_legacy_modular(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode, umat_M->ncp);
            break;
        }
        case UmatKernel::mean_field_homogeneous: case UmatKernel::mean_field_mori_tanaka:
        case UmatKernel::mean_field_self_consistent: case UmatKernel::mean_field_periodic_layer: {
            umat_multi(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt, mean_field_scheme(entry->mechanical));
            break;
        }
        case UmatKernel::python_callback: {
            // PYEXT: process-wide callback UMAT (umat_callback.hpp; registered by the Python
            // bindings). Small-strain convention on the small-strain box, as EPICP.
            umat_callback_M(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        default: {
            throw std::logic_error("select_umat_M_3D: no kernel case for '" + rve.sptr_matprops->umat_name
                                   + "' (the registry and the switch disagree)");
        }
    }
    // small strain: F = I, J = 1, so the box stress is both Cauchy and Kirchhoff
    umat_M->tau = umat_M->sigma;
    if (umat_M->F1.is_diagmat() && arma::all(umat_M->F1.diag() == 1.)) {
        // F1 = I exactly (the small-strain route): F^-1 tau F^-T is tau itself
        umat_M->PKII = umat_M->tau;
    }
    else {
        umat_M->PKII = t2v_stress(Kirchoff2PKII(v2t_stress(umat_M->tau), umat_M->F1));
    }
    rve.local2global();
}
    
namespace {

// State a kernel call overwrites, sub-phases included, restored before each condensed evaluation.
// Not to_start(): the Abaqus entries do not fill the *_start fields (same field list as to_start).
struct phase_snapshot {
    state_variables *sv;
    vec *Wm = nullptr, *Wt = nullptr;
    vec sigma, tau, PKII, statev, Wm0, Wt0;
    std::vector<phase_snapshot> sub;
};

phase_snapshot take_snapshot(phase_characteristics &pc) {
    phase_snapshot s;
    s.sv = pc.sptr_sv_global.get();
    if (auto sv_M = dynamic_cast<state_variables_M *>(s.sv)) s.Wm = &sv_M->Wm;
    if (auto sv_T = dynamic_cast<state_variables_T *>(s.sv)) {
        s.Wm = &sv_T->Wm;
        s.Wt = &sv_T->Wt;
    }
    s.sigma = s.sv->sigma;
    s.tau = s.sv->tau;
    s.PKII = s.sv->PKII;
    s.statev = s.sv->statev;
    if (s.Wm) s.Wm0 = *s.Wm;
    if (s.Wt) s.Wt0 = *s.Wt;
    s.sub.reserve(pc.sub_phases.size());
    for (auto &sub : pc.sub_phases) s.sub.push_back(take_snapshot(sub));
    return s;
}

void restore_snapshot(const phase_snapshot &s) {
    s.sv->sigma = s.sigma;
    s.sv->tau = s.tau;
    s.sv->PKII = s.PKII;
    s.sv->statev = s.statev;
    if (s.Wm) *s.Wm = s.Wm0;
    if (s.Wt) *s.Wt = s.Wt0;
    for (const auto &sub : s.sub) restore_snapshot(sub);
}

// Condensed state around the 3D dispatch run_3D (solve_condensed); no convergence: start state
// restored, step cut requested
template <class Run>
bool select_condensed(phase_characteristics &rve, const condensed_state &cs, bool &start, double &tnew_dt,
                      const mat &Lt, Run &&run_3D) {
    state_variables &sv = *rve.sptr_sv_global;
    const phase_snapshot snap = take_snapshot(rve);
    const bool start0 = start;
    auto restore = [&]() {
        restore_snapshot(snap);
        start = start0;
    };
    if (!solve_condensed(cs, sv.DEtot, sv.sigma, Lt, [&](double &tn) { run_3D(start, tn); }, restore, tnew_dt))
        return false;
    sv.tau.elem(cs.C).zeros();
    sv.PKII.elem(cs.C).zeros();
    return true;
}

}  // namespace

// A statev index past nstatev surfaces as Armadillo's bounds error: name the law and the likely cause
[[noreturn]] static void rethrow_statev_bounds(const std::out_of_range &e, const phase_characteristics &rve)
{
    throw std::invalid_argument("umat " + rve.sptr_matprops->umat_name + ": " + e.what() + " (nstatev = "
                                + std::to_string(rve.sptr_sv_global->nstatev)
                                + " is most likely too small for this law: see the statev layout in its header)");
}

void select_umat_M(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    try {
        if (ndi == 3) {
            select_umat_M_3D(rve, DR_global, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            return;
        }
        const condensed_state cs = make_condensed_state(ndi, nshr);
        auto sv = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
        if (select_condensed(rve, cs, start, tnew_dt, sv->Lt,
                             [&](bool &st, double &tn) { select_umat_M_3D(rve, DR_global, Time, DTime, 3, 3, st, solver_type, tn); }))
            condense_tangent(cs, sv->Lt);
    } catch (const std::out_of_range &e) {
        rethrow_statev_bounds(e, rve);
    }
}

void select_umat_T(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    try {
        if (ndi == 3) {
            select_umat_T_3D(rve, DR_global, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            return;
        }
        const condensed_state cs = make_condensed_state(ndi, nshr);
        auto sv = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
        if (select_condensed(rve, cs, start, tnew_dt, sv->dSdE,
                             [&](bool &st, double &tn) { select_umat_T_3D(rve, DR_global, Time, DTime, 3, 3, st, solver_type, tn); }))
            condense_tangent_T(cs, sv->dSdE, sv->dSdT, sv->drdE, sv->drdT);
    } catch (const std::out_of_range &e) {
        rethrow_statev_bounds(e, rve);
    }
}

namespace {

// Condensed state on the finite route: Newton on the free components of the box strain increment
// (solve_condensed), the stress the law returns vanishing on them; F1 follows the strain state
// (condensed_F, from the F1 of the start free stretches that abaqus2smart_finite built). The box
// tangent is the Jacobian (dT/dH under the total form: a secant there, the residual is exact).
template <class Run>
bool select_finite_condensed(phase_characteristics &rve, const condensed_state &cs, const int &corate_type,
                             bool &start, double &tnew_dt, const mat &Lt, Run &&run_3D) {
    state_variables &sv = *rve.sptr_sv_global;
    const condensed_F F1_of(sv.F1);
    const phase_snapshot snap = take_snapshot(rve);
    const bool start0 = start;
    auto restore = [&]() {
        restore_snapshot(snap);
        start = start0;
    };
    auto run = [&](double &tn) {
        const vec e1 = sv.etot + sv.Detot;
        sv.F1 = (corate_type == corate_type_total) ? F1_of.total(e1) : F1_of(cs, e1);
        run_3D(start, tn);
    };
    if (!solve_condensed(cs, sv.Detot, sv.sigma, Lt, run, restore, tnew_dt))
        return false;
    sv.tau.elem(cs.C).zeros();
    return true;
}

}  // namespace

void run_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, const int &corate_type, double &tnew_dt)
{
    
    tnew_dt = 1.;
    
    if (Time > simcoon::limit) {
        start = false;
    }
    switch (control_type) {
        case 1: {
            select_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 2: case 3: case 4: case 5: case 6: {
            try {
                if (ndi < 3) {
                    const condensed_state cs = make_condensed_state(ndi, nshr);
                    auto sv = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                    if (select_finite_condensed(rve, cs, corate_type, start, tnew_dt, sv->dSdE,
                            [&](bool &st, double &tn) { select_umat_T_finite(rve, DR, Time, DTime, 3, 3, st, solver_type, corate_type, tn); }))
                        condense_tangent_T(cs, sv->dSdE, sv->dSdT, sv->drdE, sv->drdT);
                }
                else
                    select_umat_T_finite(rve, DR, Time, DTime, ndi, nshr, start, solver_type, corate_type, tnew_dt);
            } catch (const std::out_of_range &e) {
                rethrow_statev_bounds(e, rve);
            }
            break;
        }
        default: {
            throw simcoon::exception_solver("run_umat: control type " + std::to_string(control_type)
                                            + " is not supported by this block type");
        }
    }
}

void run_umat_M(phase_characteristics &rve, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, const int &corate_type, double &tnew_dt, const bool &material_tangent)
{

    tnew_dt = 1.;
    if (Time > simcoon::limit) {
        start = false;
    }
    switch (control_type) {
        case 1: {
            select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 2: case 3: case 4: case 5: case 6: {
            try {
                if (ndi < 3) {
                    const condensed_state cs = make_condensed_state(ndi, nshr);
                    auto sv = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                    if (select_finite_condensed(rve, cs, corate_type, start, tnew_dt, sv->Lt,
                            [&](bool &st, double &tn) { select_umat_M_finite(rve, DR, Time, DTime, 3, 3, st, solver_type, corate_type, tn, material_tangent); }))
                        condense_tangent(cs, sv->Lt);
                }
                else
                    select_umat_M_finite(rve, DR, Time, DTime, ndi, nshr, start, solver_type, corate_type, tnew_dt, material_tangent);
            } catch (const std::out_of_range &e) {
                rethrow_statev_bounds(e, rve);
            }
            break;
        }
        default: {
            throw simcoon::exception_solver("run_umat: control type " + std::to_string(control_type)
                                            + " is not supported by this block type");
        }
    }
}
    
} //namespace simcoon
