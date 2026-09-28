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
///@file umat_smart.cpp
///@brief Selection of constitutive laws and transfer between Abaqus and simcoon formats,
///       implemented in 1D-2D-3D
///@version 1.0

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <armadillo>
#include <memory>
#include <dylib.hpp>
#include <filesystem>

#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>

#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_pstretch.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_callback.hpp>

#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso_ccp.hpp>
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

    for (auto r:rve.sub_phases) {
        size = size + r.sptr_sv_local->nstatev + 57;
        size_statev(r,size);
    }
}

void statev_2_phases(phase_characteristics &rve, unsigned int &pos, const vec &statev) {

    for(auto r : rve.sub_phases) {
        //The number of statev is here determined for each phase, then a sub_vector of the statev vector is taken from this
        //vec Etot -> 6X
        //vec DEtot -> 6X 12
        //vec sigma -> 6X 18
        //vec sigma_start -> 6
        //mat F0 -> 9
        //mat F1 -> 9
        //double T -> 1X 19
        //double DT -> 1X 20
        //vec sigma_in -> 6
        //vec sigma_in_start -> 6
    
        //vec Wm    -> 3X 23
        //vec Wm_start  -> 3X 26
        //mat L -> 36X 62
        //mat Lt -> 36
        
        //int nstatev
        //vec statev    nstatevX 62+nstatev
        //vec statev_start

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
    
    for(auto r : rve.sub_phases) {
        //The number of statev is here determined for each phase, then a sub_vector of the statev vector is taken from this
        //vec Etot -> 6X
        //vec DEtot -> 6X 12
        //vec sigma -> 6X 18
        //vec sigma_start -> 6
        //mat F0 -> 9
        //mat F1 -> 9
        //double T -> 1X 19
        //double DT -> 1X 20
        //vec sigma_in -> 6
        //vec sigma_in_start -> 6
        
        //vec Wm    -> 3X 23
        //vec Wm_start  -> 3X 26
        //mat L -> 36X 62
        //mat Lt -> 36
        
        //int nstatev
        //vec statev    nstatevX 62+nstatev
        //vec statev_start
        
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

const std::map<string, int> &finite_umat_names()
{
    // The finite dispatch's name -> id map. A file-scope accessor rather than a
    // function-local static inside select_umat_M_finite so the convention test can iterate it.
    static const std::map<string, int> list_umat = {{"UMEXT",0},{"UMABA",1},{"ELISO",201},{"ELIST",201},{"ELORT",201},{"HYPOO",5},{"EPICP",6},{"EPCHA",7},{"EPKCP",201},{"SNTVE",8},{"NEOHI",9},{"NEOHC",10},{"MOORI",11},{"YEOHH",12},{"ISHAH",13},{"GETHH",14},{"SWANH",15},{"HOLZA",16},{"EPHIL",201},{"EPTRI",201},{"EPHAC",201},{"EPANI",201},{"EPDFA",201},{"EPCHG",201},{"EPHIN",201},{"MODUL",200},{"OGDEN",22},{"PYEXT",300}};
    return list_umat;
}

namespace {

// The single convention table, plain types only: a function-local static of an armadillo type 
// registers a destructor that runs at DLL unload, and on Windows that order is undefined.
const std::map<string, umat_convention> &umat_conventions()
{
    using SM = StressMeasure;

    // Every NATIVE kernel is Kirchhoff. The groups below differ in how the tangent reaches the
    // solver's corate (fed the corotated strain: in-rate for free; built from F: converted by
    // Dtau_LieDD_2_DtauDe_corate) and in the frame they run in (material or lab).
    static const std::map<string, umat_convention> conventions = {
        // --- log-strain boxes: fed the corotated strain, run in the material frame ---
        // {measure, material_frame, layout_declared, tensorial statev {offset, Tensor2Type}}
        {"ELISO", {SM::kirchhoff, true, true, {}}},
        {"ELIST", {SM::kirchhoff, true, true, {}}},
        {"ELORT", {SM::kirchhoff, true, true, {}}},
        {"EPICP", {SM::kirchhoff, true, true, {{2, Tensor2Type::strain}}}},                  // EP
        {"EPCHA", {SM::kirchhoff, true, true, {{2, Tensor2Type::strain},                     // EP
                                               {8, Tensor2Type::strain},                     // a_1
                                               {14, Tensor2Type::strain},                    // a_2
                                               {20, Tensor2Type::stress},                    // X_1
                                               {26, Tensor2Type::stress}}}},                 // X_2
        // 201 adapters (modular engine): layout not declared -> refused under corates 4/5
        {"EPKCP", {SM::kirchhoff, true, false, {}}},
        {"EPHIL", {SM::kirchhoff, true, false, {}}},
        {"EPTRI", {SM::kirchhoff, true, false, {}}},
        {"EPHAC", {SM::kirchhoff, true, false, {}}},
        {"EPANI", {SM::kirchhoff, true, false, {}}},
        {"EPDFA", {SM::kirchhoff, true, false, {}}},
        {"EPCHG", {SM::kirchhoff, true, false, {}}},
        {"EPHIN", {SM::kirchhoff, true, false, {}}},
        // MODUL builds its tangent from V_el, so its log box IS d(tau)/d(eps_el); it enforces
        // corate 3, and in the material frame its fibres rotate with the body.
        {"MODUL", {SM::kirchhoff, true, false, {}}},
        // PYEXT: fed the log strain, returns tau (see umat_callback.hpp)
        {"PYEXT", {SM::kirchhoff, true, false, {}}},
        // HYPOO: the rate counterpart of ELORT, tau_{n+1} = tau_n + L : DEel
        {"HYPOO", {SM::kirchhoff, true, true, {}}},

        // --- finite kernels: build the tangent from F, then convert it to corate_type ---
        {"SNTVE", {SM::kirchhoff}},
        {"NEOHI", {SM::kirchhoff}},
        {"NEOHC", {SM::kirchhoff}},
        {"MOORI", {SM::kirchhoff}},
        {"YEOHH", {SM::kirchhoff}},
        {"ISHAH", {SM::kirchhoff}},
        {"GETHH", {SM::kirchhoff}},
        {"SWANH", {SM::kirchhoff}},
        {"OGDEN", {SM::kirchhoff}},
        {"HOLZA", {SM::kirchhoff}},

        // --- foreign conventions simcoon does not own ---
        // Plugin adapters: the contract is the host code's (Abaqus DDSDDE is Cauchy-based).
        // UMABA is in fact unreachable today -- id 1 has no case in the finite switch, so it
        // throws -- and UMEXT's body is fully commented out. Declared for completeness.
        {"UMEXT", {SM::cauchy}},
        {"UMABA", {SM::cauchy}},
    };
    return conventions;
}

}  // namespace

umat_convention output_convention_of(const std::string &umat_name)
{
    const auto &conventions = umat_conventions();
    auto it = conventions.find(umat_name);
    if (it == conventions.end()) {
        throw std::invalid_argument(
            "output_convention_of: the umat '" + umat_name + "' has not declared what its raw "
            "outputs are expressed in. Add it to the table in umat_smart.cpp: a missing stress "
            "measure is an error of exactly J, a missing tangent rate is a wrong rate, and both "
            "are silent.");
    }
    return it->second;
}

bool stress_output_is_kirchhoff(const std::string &umat_name)
{
    // Total where output_convention_of throws, and it shares that function's table rather than
    // calling it: the python wrapper asks this for EVERY name it serves, including
    // small-strain-only ones (SMADI, ZENER, the micro family) that the finite dispatch never
    // sees and that therefore declare no finite convention. For those the answer is "leave the
    // stress alone", not "refuse the call" -- and it must not cost a thrown exception per call.
    // Inside the finite dispatch, where a missing declaration IS a bug, use
    // output_convention_of and let it throw.
    const auto &conventions = umat_conventions();
    auto it = conventions.find(umat_name);
    return it != conventions.end() && it->second.stress == StressMeasure::kirchhoff;
}

void select_umat_T(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    UNUSED(solver_type);
    static const std::map<string, int> list_umat = {{"ELISO",1},{"ELIST",2},{"ELORT",3},{"EPICP",4},{"EPKCP",5},{"ZENER",6},{"ZENNK",7},{"PRONK",8},{"SMADI",9},{"SMADC",9},{"SMAAI",9},{"SMAAC",9}};

    // Same frame handling as select_umat_M_finite: the caller's DR is global;
    // rotate it with the other state variables so the local UMAT receives the
    // material-frame increment.
    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_local);
    const mat &DR = umat_T->DR;
    
    auto it_umat = list_umat.find(rve.sptr_matprops->umat_name);
    switch (it_umat != list_umat.end() ? it_umat->second : -1) {
        case 0: {
//            umat_external_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 1: {
            umat_elasticity_iso_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 2: {
            umat_elasticity_trans_iso_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 3: {
            umat_elasticity_ortho_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 4: {
            umat_plasticity_iso_CCP_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 5: {
            umat_plasticity_kin_iso_CCP_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
           break;
        }
        case 6: {
            umat_zener_fast_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 7: {
            umat_zener_Nfast_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 8: {
            umat_prony_Nfast_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        case 9: {
            umat_sma_unified_T_T(rve.sptr_matprops->umat_name, umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt, umat_T->tangent_mode);
            break;
        }
        default: {
            throw std::invalid_argument("Unknown umat name in the thermomechanical dispatch: " + rve.sptr_matprops->umat_name);
        }
    }
    // Thermomechanical route is small strain (F = I, J = 1): the UMATs return Cauchy sigma;
    // mirror the mechanical stress-route convention (see select_umat_M/_finite) so the
    // Kirchhoff-route consumers (output tau/J, set_start, sinks) stay correct: tau = J sigma = sigma,
    // PKII = sigma. Without this the thermomechanical outputs read an empty tau (zero stress).
    umat_T->tau = umat_T->sigma;
    umat_T->PKII = umat_T->sigma;
    rve.local2global();
    
}

namespace {

// One tensorial internal variable carried by the rotation-free increment M of corates 4 and 5,
// with its variance: Truesdell (4) is tensor2's own push-forward (stress F X F^T, strain
// F^-T X F^-1, no Piola factor); log_F (5) the similarity M X M^-1 for both, symmetrised. Same
// rules as the solver's transport of tau and etot.
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

}  // namespace

void select_umat_M_finite(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const int &corate_type, double &tnew_dt)
{
    const std::map<string, int> &list_umat = finite_umat_names();

    // guarded lookup: operator[] would default-insert 0 (=UMEXT, a no-op) for an
    // unknown name and silently return zero stress; -1 falls to the default case
    auto it_umat = list_umat.find(rve.sptr_matprops->umat_name);
    const int id_umat = (it_umat != list_umat.end()) ? it_umat->second : -1;

    // The caller provides DR in global coordinates; rotate it with the other
    // state variables so that the local UMAT receives a consistent increment.
    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    const mat &DR = umat_M->DR;

    // Transport the start state with THIS increment's DR before the kernel sees it, so the
    // total strain, the start stress and the internal variables (rotated by DR inside the
    // kernels) all live in the same configuration. The stored etot is left untransported:
    // set_start transports it at commit, together with the increment it belongs to.
    const vec etot_stored = umat_M->etot;
    if (corate_type != 4 && corate_type != 5) {
        umat_M->etot = rotate_strain(etot_stored, DR);
        umat_M->sigma = rotate_stress(umat_M->sigma, DR);
    }
    else if (corate_type == 4) {   // Truesdell, DR = DF: strain lower-, Kirchhoff stress upper-convected
        mat DR_inv;
        if (!inv(DR_inv, DR))
            throw simcoon::exception_inv("select_umat_M_finite: DF is not invertible");
        umat_M->etot = t2v_strain(DR_inv.t()*v2t_strain(etot_stored)*DR_inv);
        umat_M->sigma = t2v_stress(DR*v2t_stress(umat_M->sigma)*DR.t());
    }
    else {   // log_F, DR = DF: similarity transport, as in set_start
        mat DR_inv;
        if (!inv(DR_inv, DR))
            throw simcoon::exception_inv("select_umat_M_finite: DF is not invertible");
        umat_M->etot = t2v_strain(DR*v2t_strain(etot_stored)*DR_inv);
        umat_M->sigma = t2v_stress(DR*v2t_stress(umat_M->sigma)*DR_inv);
    }

    const vec tau_start_tr = umat_M->sigma;   // transported start stress, for the work correction

    // ONE declaration per kernel (output_convention_of) answers the questions below. It THROWS
    // for a kernel that has not declared, so a new one is caught instead of inheriting a
    // default; unknown names fall to the switch's own error.
    umat_convention conv{StressMeasure::kirchhoff};
    if (id_umat >= 0)
        conv = output_convention_of(rve.sptr_matprops->umat_name);

    if (id_umat == 200) {   // MODUL
        // tau = d psi / d eps_el is a stored-energy law only when the accumulated strain is
        // ln V, i.e. corate 3 (log_R); any other corate makes it a non-integrable rate law.
        if (corate_type != 3) {
            throw simcoon::exception_solver(
                "MODUL under finite strain requires corate_type = 3 "
                "(log_R): the modular composition is hyperelastic in "
                "the logarithmic strain; got corate_type = "
                + std::to_string(corate_type));
        }
    }

    if (id_umat == 7 && corate_type == 4) {   // EPCHA
        // X_i = 2/3 C_i a_i is stored: Truesdell convects a (strain) and X (stress) differently,
        // so the pair would drift apart. log_F transports both by the same similarity.
        throw simcoon::exception_solver(
            "EPCHA stores its back stresses X_i = 2/3 C_i a_i, which the Truesdell rate "
            "(corate_type 4) cannot keep consistent; use corate_type 3 (or 5)");
    }

    // A box kernel runs in the frame that follows the material, R_hat (see umat_convention), so
    // its anisotropy axes rotate with the body and its statev lives in that frame. For corates
    // 0-3, R_hat_{n+1} = DR R_hat_n and the frame-relative increment is I; for 4 and 5 it is the
    // rotation-free stretch M, applied here to the declared statev tensors.
    const vec Detot_stored = umat_M->Detot;
    if (corate_type == 4) {
        // Truesdell: the kernel's increment is the closed-form Almansi increment of DR = DF,
        // whatever the solver's control variable is (logarithmic control increments ln V).
        mat bDF_inv;
        if (!inv_sympd(bDF_inv, DR*DR.t()))
            throw simcoon::exception_inv("select_umat_M_finite: DF DF^T is not invertible");
        umat_M->Detot = t2v_strain(0.5*(eye(3,3) - bDF_inv));
    }
    mat R_hat = eye(3,3);
    mat DR_kernel = DR;
    if (conv.material_frame) {
        const bool convected = (corate_type == 4 || corate_type == 5);
        mat R_hat_n;
        if (!convected) {
            R_hat_n = umat_M->R;
            R_hat = DR*umat_M->R;
        }
        else {
            mat U;
            RU_decomposition(R_hat_n, U, umat_M->F0);
            RU_decomposition(R_hat, U, umat_M->F1);
        }
        umat_M->etot = rotate_strain(umat_M->etot, R_hat.t());
        umat_M->Detot = rotate_strain(Detot_stored, R_hat.t());
        umat_M->sigma = rotate_stress(umat_M->sigma, R_hat.t());
        DR_kernel = eye(3,3);
        if (convected) {
            if (!conv.layout_declared) {
                throw simcoon::exception_solver(
                    "corate_type " + std::to_string(corate_type) + " requires the tensorial "
                    "internal variables of '" + rve.sptr_matprops->umat_name + "' to be declared "
                    "(umat_conventions, umat_smart.cpp); use corate_type 3");
            }
            const mat M = R_hat.t()*DR*R_hat_n;
            for (const StatevTensor &t : conv.statev_tensors) {
                if (t.offset + 6 > static_cast<int>(umat_M->statev.n_elem)) {
                    throw simcoon::exception_solver(
                        "'" + rve.sptr_matprops->umat_name + "' declares a statev tensor at offset "
                        + std::to_string(t.offset) + ", past nstatev = " + std::to_string(umat_M->statev.n_elem));
                }
                const vec x = umat_M->statev.subvec(t.offset, t.offset + 5);
                umat_M->statev.subvec(t.offset, t.offset + 5) = transport_convected(x, M, corate_type, t.type);
            }
        }
    }

    switch (id_umat) {

            case 0: {
/*                umat_external_M(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
*/
                break;
            }
            case 5: {
                umat_hypoelasticity_ortho(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_type, umat_M->tangent_mode);
                break;
            }
             case 6: {
                umat_plasticity_iso_CCP(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                 break;
             }
             case 7: {
                // Chaboche on the log-strain/Kirchhoff box route, same as EPICP
                umat_plasticity_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                 break;
             }
            case 8: {
                umat_saint_venant(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_type, umat_M->tangent_mode);
                break;
            }
            case 201: {
                // Legacy names served by the modular engine on the log-strain box, like EPICP:
                // run in the material frame, so their anisotropy axes follow the body.
                umat_legacy_modular(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                break;
            }
            case 200: {
                // MODUL under NLGEOM is a Hencky hyperelastic composition:
                umat_modular(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                break;
            }
            case 9: {
                umat_neo_hookean_incomp(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_type, umat_M->tangent_mode);
                break;
            }                         
            case 10: case 11: case 12: case 13: case 14: case 15: case 16: {
                umat_generic_hyper_invariants(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_type, umat_M->tangent_mode);
                break;
            }
            case 22: {
                umat_generic_hyper_pstretch(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, corate_type, umat_M->tangent_mode);
                break;
            }
            case 300: {
                // PYEXT: process-wide callback UMAT (umat_callback.hpp; registered by the Python
                // bindings). Small-strain convention on the log-strain / Kirchhoff box, as EPICP.
                umat_callback_M(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR_kernel, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
                break;
            }
            default: {
                throw std::invalid_argument("Unknown umat name in the finite-strain dispatch: " + rve.sptr_matprops->umat_name);
            }
        }
    
        // Back from the material frame: the stress and the box tangent (both in the kernel's
        // frame) to the lab; statev stays in the material frame.
        if (conv.material_frame) {
            umat_M->sigma = rotate_stress(umat_M->sigma, R_hat);
            umat_M->Lt = rotate_stiffness(umat_M->Lt, R_hat);
            umat_M->L = rotate_stiffness(umat_M->L, R_hat);
        }

        // tau is the canonical Kirchhoff route stress (Wm stays on the Kirchhoff route). Every
        // NATIVE kernel returns it directly, so this is a pass-through for all but the
        // foreign-convention ones (the UMEXT/UMABA plugin adapters, whose contract belongs to
        // the host code). Cauchy is a derived OUTPUT (tau/J), never on the route.
        if (conv.stress == StressMeasure::kirchhoff)
            umat_M->tau = umat_M->sigma;                                                      // the kernel's output IS the Kirchhoff stress
        else
            umat_M->tau = t2v_stress(Cauchy2Kirchoff(v2t_stress(umat_M->sigma), umat_M->F1));  // foreign Cauchy -> Kirchhoff
        umat_M->PKII = t2v_stress(Kirchoff2PKII(v2t_stress(umat_M->tau), umat_M->F1));

        // Log corates: the kernel's work is on the corate strain, not D; restore the true work.
        if (corate_type == 2 || corate_type == 3 || corate_type == 5) {
            const double dW = Delta_work_conjugacy(tau_start_tr, umat_M->tau, Detot_stored, umat_M->F0, umat_M->F1);
            umat_M->Wm(0) += dW;
            umat_M->Wm(1) += dW;
        }

        // No tangent conversion: every kernel emits Lt in corate_type already.
        umat_M->etot = etot_stored;
        umat_M->Detot = Detot_stored;
        rve.local2global();
}
    
void select_umat_M(phase_characteristics &rve, const mat &DR_global,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{

    static const std::map<string, int> list_umat = {{"UMEXT",0},{"UMABA",1},{"ELISO",201},{"ELIST",201},{"ELORT",201},{"EPICP",5},{"EPKCP",201},{"EPCHA",7},{"SMADI",8},{"SMADC",8},{"SMAAI",8},{"SMAAC",8},{"SMRDI",9},{"SMRDC",9},{"SMRAI",9},{"SMRAC",9},{"LLDM0",10},{"ZENER",11},{"ZENNK",12},{"PRONK",13},{"EPHIL",201},{"EPTRI",201},{"EPHAC",201},{"EPANI",201},{"EPDFA",201},{"EPCHG",201},{"EPHIN",201},{"SMAMO",23},{"SMAMC",24},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104},{"MODUL",200},{"PYEXT",300}};

    // Same frame handling as select_umat_M_finite: the caller's DR is global;
    // rotate it with the other state variables so the local UMAT receives the
    // material-frame increment.
    rve.sptr_sv_global->DR = DR_global;
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    const mat &DR = umat_M->DR;

    auto it_umat = list_umat.find(rve.sptr_matprops->umat_name);
    switch (it_umat != list_umat.end() ? it_umat->second : -1) {

        case 0: {
            //umat_external(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);

            static dylib::library ext_lib("external/umat_plugin_ext", dylib::decorations::os_default());

            using create_fn = umat_plugin_ext_api*();
            using destroy_fn = void(umat_plugin_ext_api*);

            static create_fn* ext_create = ext_lib.get_function<create_fn>("create_api");
            static destroy_fn* ext_destroy = ext_lib.get_function<destroy_fn>("destroy_api");

            static std::unique_ptr<umat_plugin_ext_api, destroy_fn*> external_umat(
                ext_create(),
                ext_destroy
            );

            // ABI-preserving dispatch: external user plugins (umat_plugin_ext_api
            // in umat_plugin_api.hpp) are loaded as dylibs at runtime and were
            // compiled against the pre-tangent_mode signature. Do NOT forward
            // umat_M->tangent_mode here — it would mismatch the plugin's
            // umat_external_M symbol and crash on dispatch. To enable
            // tangent_mode for external plugins, bump umat_plugin_ext_api's
            // pure-virtual signature and recompile every plugin DSO.
            external_umat->umat_external_M(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);

            break;
        }
        case 1: {
            //
            static dylib::library aba_lib("external/umat_plugin_aba", dylib::decorations::os_default());

            using create_fn = umat_plugin_aba_api*();
            using destroy_fn = void(umat_plugin_aba_api*);

            static create_fn* aba_create = aba_lib.get_function<create_fn>("create_api");
            static destroy_fn* aba_destroy = aba_lib.get_function<destroy_fn>("destroy_api");

            static std::unique_ptr<umat_plugin_aba_api, destroy_fn*> abaqus_umat(
                aba_create(),
                aba_destroy
            );

            abaqus_umat->umat_abaqus(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 5: {
            umat_plasticity_iso_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 7: {
            umat_plasticity_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 8: {
            umat_sma_unified_T(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 9: {
            umat_sma_unified_TR(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 10: {
            umat_damage_LLD_0(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 11: {
            umat_zener_fast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 12: {
            umat_zener_Nfast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 13: {
            umat_prony_Nfast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 23: case 24: {
            // SMAMO (isotropic) and SMAMC (cubic) both use unified umat_sma_mono
            umat_sma_mono(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 200: {
            // Modular UMAT - composable constitutive model
            umat_modular(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 201: {
            // Legacy names served by the modular engine: props are translated
            // by the per-name adapter (legacy_adapters.cpp), state/outputs
            // flow through umat_modular untouched.
            umat_legacy_modular(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        case 100: case 101: case 103: case 104: {
            umat_multi(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt, it_umat->second);
            break;
        }
        case 300: {
            // PYEXT: process-wide callback UMAT (umat_callback.hpp; registered by the Python
            // bindings). Small-strain convention on the small-strain box, as EPICP.
            umat_callback_M(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt, umat_M->tangent_mode);
            break;
        }
        default: {
            throw std::invalid_argument("Unknown umat name in the small-strain dispatch: " + rve.sptr_matprops->umat_name);
        }
    }
    // Small-strain control (control_type 1): F = I, J = 1, so the box stress is
    // simultaneously the Cauchy and the Kirchhoff stress. Mirror it onto tau/PKII
    // so the route (the Newton iterates on tau) and the outputs stay consistent
    // with the finite dispatcher.
    umat_M->tau = umat_M->sigma;
    umat_M->PKII = t2v_stress(Kirchoff2PKII(v2t_stress(umat_M->tau), umat_M->F1));
    rve.local2global();
}
    
void run_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, double &tnew_dt)
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
//        case 2: case 3: case 4: case 5: {
//        select_umat_T_finite(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
//        break;
//        }
        default: {
            throw simcoon::exception_solver("run_umat: control type " + std::to_string(control_type)
                                            + " is not supported by this block type");
        }
    }
}

void run_umat_M(phase_characteristics &rve, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, const int &corate_type, double &tnew_dt)
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
            select_umat_M_finite(rve, DR, Time, DTime, ndi, nshr, start, solver_type, corate_type, tnew_dt);
            break;
        }
        default: {
            throw simcoon::exception_solver("run_umat: control type " + std::to_string(control_type)
                                            + " is not supported by this block type");
        }
    }
}
    
} //namespace simcoon
