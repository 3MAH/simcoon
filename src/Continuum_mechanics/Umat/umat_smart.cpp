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
///@brief Implemented in 1D-2D-3D
///@version 1.0

#include <iostream>
#include <map>
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

#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_comp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Ani_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/DFA_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Generic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/aniso_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>

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

void select_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    UNUSED(solver_type);
    std::map<string, int> list_umat;
    list_umat = {{"UMEXT",0},{"ELISO",1},{"ELIST",2},{"ELORT",3},{"EPICP",4},{"EPKCP",5},{"ZENER",6},{"ZENNK",7},{"PRONK",8},{"SMAUT",9}};

    rve.global2local();
    auto umat_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_local);
    
    switch (list_umat[rve.sptr_matprops->umat_name]) {
        case 0: {
//            umat_external_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 1: {
            umat_elasticity_iso_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 2: {
            umat_elasticity_trans_iso_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 3: {
            umat_elasticity_ortho_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 4: {
            umat_plasticity_iso_CCP_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 5: {
            umat_plasticity_kin_iso_CCP_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
           break;
        }
        case 6: {
            umat_zener_fast_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 7: {
            umat_zener_Nfast_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 8: {
            umat_prony_Nfast_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        case 9: {
            umat_sma_unified_T_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
            break;
        }
        default: {
            cout << "Error: The choice of Thermomechanical Umat could not be found in the umat library :" << rve.sptr_matprops->umat_name << "\n";
                exit(0);
        }
    }
    rve.local2global();
    
}

void select_umat_M_finite(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    std::map<string, int> list_umat;
    
    list_umat = {{"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"HYPOO",5},{"EPICP",6},{"EPKCP",7},{"SNTVE",8},{"NEOHI",9},{"NEOHC",10},{"MOORI",11},{"YEOHH",12},{"ISHAH",13},{"GETHH",14},{"SWANH",15}};
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    
    switch (list_umat[rve.sptr_matprops->umat_name]) {

            case 0: {
/*                umat_external_M(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
*/
                break;
            }
            case 2: {
                umat_elasticity_iso(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                 break;
             }
            case 3: {
                umat_elasticity_trans_iso(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                 break;
             }
            case 4: {
                umat_elasticity_ortho(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                 break;
             }
            case 5: {
                umat_hypoelasticity_ortho(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
             case 6: {
                umat_plasticity_iso_CCP(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                 break;
             }
             case 7: {
                umat_plasticity_kin_iso_CCP(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                 break;
             }
            case 8: {
                umat_saint_venant(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 9: {
                umat_neo_hookean_incomp(umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }                         
            case 10: case 11: case 12: case 13: case 14: case 15: {
                umat_generic_hyper_invariants(rve.sptr_matprops->umat_name, umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }     
            default: {
                cout << "Error: The choice of Umat could not be found in the umat library :" << rve.sptr_matprops->umat_name << "\n";
                exit(0);
            }                   
        }
    
        umat_M->PKII = t2v_stress(Cauchy2PKII(v2t_stress(umat_M->sigma), umat_M->F1));
        umat_M->tau = t2v_stress(Cauchy2Kirchoff(v2t_stress(umat_M->sigma), umat_M->F1));
        rve.local2global();
}
    
void select_umat_M(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{

    std::map<string, int> list_umat;
    
    list_umat = {{"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"SMAUT",8},{"SMANI",9},{"LLDM0",10},{"ZENER",11},{"ZENNK",12},{"PRONK",13},{"EPHIL",17},{"EPHAC",18},{"EPANI",19},{"EPDFA",20},{"EPCHG",21},{"EPHIN",22},{"SMAMO",23},{"SMAMC",24},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104}};
    
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);

    switch (list_umat[rve.sptr_matprops->umat_name]) {

        case 0: {
            //umat_external(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);

            static dylib::library ext_lib("external/umat_plugin_ext", dylib::decorations::os_default());

            using create_fn = umat_plugin_ext_api*();
            using destroy_fn = void(umat_plugin_ext_api*);

            static create_fn* ext_create = ext_lib.get_function<create_fn>("create_api");
            static destroy_fn* ext_destroy = ext_lib.get_function<destroy_fn>("destroy_api");

            static std::unique_ptr<umat_plugin_ext_api, destroy_fn*> external_umat(
                ext_create(),
                ext_destroy
            );

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
        case 2: {
            umat_elasticity_iso(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
        case 3: {
            umat_elasticity_trans_iso(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 4: {
            umat_elasticity_ortho(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 5: {
            umat_plasticity_iso_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 6: {
            umat_plasticity_kin_iso_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 7: {
            umat_plasticity_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 8: {
            umat_sma_unified_T(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 9: {
            umat_sma_aniso_T(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 10: {
            umat_damage_LLD_0(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 11: {
            umat_zener_fast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 12: {
            umat_zener_Nfast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 13: {
            umat_prony_Nfast(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 17: {
            umat_plasticity_hill_isoh_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 18: {
            umat_hill_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
        }
        case 19: {
            umat_ani_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
        }
        case 20: {
            umat_dfa_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
        }
        case 21: {
            umat_generic_chaboche_CCP(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
        }
        case 22: {
            umat_plasticity_hill_isoh_CCP_N(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }            
        case 23: case 24: {
            // SMAMO (isotropic) and SMAMC (cubic) both use unified umat_sma_mono
            umat_sma_mono(rve.sptr_matprops->umat_name, umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 100: case 101: case 103: case 104: {
            umat_multi(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt, list_umat[rve.sptr_matprops->umat_name]);
            break;
        }
        default: {
            cout << "Error: The choice of Umat could not be found in the umat library :" << rve.sptr_matprops->umat_name << "\n";
            exit(0);
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
