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
///@brief Selection of constitutive laws and transfer to between Abaqus and simcoon formats
///@brief Thin wrapper around UmatDispatch singleton for phase_characteristics interface
///@version 2.0

#include <iostream>
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
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>

// UmatDispatch is the single source of truth for UMAT dispatch
#include <simcoon/Simulation/Solver_optimized/umat_dispatch.hpp>

// Finite strain UMATs (needed for direct calls with F0/F1)
#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>

// Hill_isoh_Nfast not in UmatDispatch (case 22 only used here)
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh_Nfast.hpp>

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

void select_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    UNUSED(solver_type);

    rve.global2local();
    auto umat_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_local);

    // Delegate to UmatDispatch singleton (single source of truth)
    bool success = UmatDispatch::instance().call_umat_T(
        rve.sptr_matprops->umat_name,
        umat_T->Etot, umat_T->DEtot,
        umat_T->sigma, umat_T->r,
        umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT,
        DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props,
        umat_T->nstatev, umat_T->statev,
        umat_T->T, umat_T->DT, Time, DTime,
        umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3),
        umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2),
        ndi, nshr, start, tnew_dt);

    if (!success) {
        cout << "Error: The choice of Thermomechanical Umat could not be found in the umat library: " << rve.sptr_matprops->umat_name << "\n";
        exit(0);
    }

    rve.local2global();
}

void select_umat_M_finite(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);

    // Delegate to UmatDispatch singleton (single source of truth)
    bool success = UmatDispatch::instance().call_umat_M_finite(
        rve.sptr_matprops->umat_name,
        umat_M->etot, umat_M->Detot,
        umat_M->F0, umat_M->F1,
        umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in,
        DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props,
        umat_M->nstatev, umat_M->statev,
        umat_M->T, umat_M->DT, Time, DTime,
        umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3),
        ndi, nshr, start, solver_type, tnew_dt);

    if (!success) {
        cout << "Error: The choice of finite strain Umat could not be found in the umat library: " << rve.sptr_matprops->umat_name << "\n";
        exit(0);
    }

    // Post-process: compute PKII and Kirchhoff stress from Cauchy
    umat_M->PKII = t2v_stress(Cauchy2PKII(v2t_stress(umat_M->sigma), umat_M->F1));
    umat_M->tau = t2v_stress(Cauchy2Kirchoff(v2t_stress(umat_M->sigma), umat_M->F1));

    rve.local2global();
}
    
void select_umat_M(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);

    // Get UMAT ID from UmatDispatch singleton (single source of truth)
    int umat_id = UmatDispatch::instance().get_umat_M_id(rve.sptr_matprops->umat_name);

    // Special cases that need rve context or plugin loading
    switch (umat_id) {
        case 0: {
            // UMEXT - External UMAT plugin
            static dylib::library ext_lib("external/umat_plugin_ext", dylib::decorations::os_default());
            using create_fn = umat_plugin_ext_api*();
            using destroy_fn = void(umat_plugin_ext_api*);
            static create_fn* ext_create = ext_lib.get_function<create_fn>("create_api");
            static destroy_fn* ext_destroy = ext_lib.get_function<destroy_fn>("destroy_api");
            static std::unique_ptr<umat_plugin_ext_api, destroy_fn*> external_umat(ext_create(), ext_destroy);

            external_umat->umat_external_M(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 1: {
            // UMABA - Abaqus UMAT plugin
            static dylib::library aba_lib("external/umat_plugin_aba", dylib::decorations::os_default());
            using create_fn = umat_plugin_aba_api*();
            using destroy_fn = void(umat_plugin_aba_api*);
            static create_fn* aba_create = aba_lib.get_function<create_fn>("create_api");
            static destroy_fn* aba_destroy = aba_lib.get_function<destroy_fn>("destroy_api");
            static std::unique_ptr<umat_plugin_aba_api, destroy_fn*> abaqus_umat(aba_create(), aba_destroy);

            abaqus_umat->umat_abaqus(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 22: {
            // EPHIN - Hill_isoh_Nfast (not in UmatDispatch)
            umat_plasticity_hill_isoh_CCP_N(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 100: case 101: case 103: case 104: {
            // Multiphase UMATs - need rve context
            umat_multi(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt, umat_id);
            break;
        }
        default: {
            // Delegate to UmatDispatch for all other UMATs
            bool success = UmatDispatch::instance().call_umat_M(
                rve.sptr_matprops->umat_name,
                umat_M->Etot, umat_M->DEtot,
                umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in,
                DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props,
                umat_M->nstatev, umat_M->statev,
                umat_M->T, umat_M->DT, Time, DTime,
                umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3),
                ndi, nshr, start, solver_type, tnew_dt);

            if (!success) {
                cout << "Error: The choice of Umat could not be found in the umat library: " << rve.sptr_matprops->umat_name << "\n";
                exit(0);
            }
            break;
        }
    }
    rve.local2global();
}
    
void run_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, double &tnew_dt)
{
    
    tnew_dt = 1.;
    
    if (Time > sim_limit) {
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
            cout << "Error: The control type of the block does not correspond" << endl;
            exit(0);
        }
    }
}

void run_umat_M(phase_characteristics &rve, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, double &tnew_dt)
{
    
    tnew_dt = 1.;
    if (Time > sim_limit) {
        start = false;
    }
    switch (control_type) {
        case 1: {
            select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 2: case 3: case 4: case 5: case 6: {
            select_umat_M_finite(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        default: {
            cout << "Error: The control type of the block does not correspond" << endl;
            exit(0);
        }
    }
}
    
} //namespace simcoon
