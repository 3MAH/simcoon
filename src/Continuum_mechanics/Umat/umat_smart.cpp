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
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono_cubic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>

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

namespace simcoon{

///@param stress array containing the components of the stress tensor (dimension ntens)
///@param stran array containing total strain component (dimension ntens) at the beginning of increment
///@param dstran array containing the component of total strain increment (dimension ntens)
///@param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
///@param dtime time increment
///@param temperature temperature avlue at the beginning of increment
///@param Dtemperature temperature increment
///@param ndi number of direct stress components
///@param nshr number of shear stress components
///@param drot rotation increment matrix (dimension 3*3)
///@param tnewdt ratio of suggested new time increment
///@param celent characteristic element length
///@param dfgrd0 array containing the deformation gradient at the beginning of increment (dimension 3*3)
///@param dfgrd1 array containing the deformation gradient at the end of increment (dimension 3*3)
///@param noel element number
///@param npt integration point number
///@param layer layer number - not used
///@param kspt section point number within the current layer - not used
///@param kstep step number
///@param kinc increment number


void size_statev(phase_characteristics &rve, unsigned int &size) {
    
    for (auto r:rve.sub_phases) {
        size+=r.sptr_sv_local->nstatev;
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

        shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
        unsigned int nstatev = umat_phase_M->statev.n_elem;

        vec vide = zeros(6);
        umat_phase_M->Etot = statev.subvec(pos,size(vide));
        umat_phase_M->DEtot = statev.subvec(pos+6,size(vide));
        umat_phase_M->sigma = statev.subvec(pos+12,size(vide));
        umat_phase_M->sigma_start = statev.subvec(pos+18,size(vide));
        umat_phase_M->T = statev(pos+19);
        umat_phase_M->DT = statev(pos+20);
        for (int i=0; i<6; i++) {
            umat_phase_M->L.col(i) = statev.subvec(pos+21+i*6,pos+21+i*6+6);
        }
        umat_phase_M->statev = statev.subvec(pos+62,pos+62+nstatev);
    
        pos+=62+nstatev;
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
        
        shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
        unsigned int nstatev = umat_phase_M->statev.n_elem;
        
        vec vide = zeros(6);
        statev.subvec(pos,size(vide)) = umat_phase_M->Etot;
        statev.subvec(pos+6,size(vide)) = umat_phase_M->DEtot;
        statev.subvec(pos+12,size(vide)) = umat_phase_M->sigma;
        statev(pos+19) = umat_phase_M->T;
        statev(pos+20) = umat_phase_M->DT;
        for (int i=0; i<6; i++) {
            statev.subvec(pos+21+i*6,pos+21+i*6+6) = umat_phase_M->L.col(i);
        }
        statev.subvec(pos+62,pos+62+nstatev) = umat_phase_M->statev;
        
        pos+=62+nstatev;
        phases_2_statev(statev,pos,r);
    }
    
}
    

void abaqus2smart_M(double *stress, double *ddsdde, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops,const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &Lt, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &statev_smart, mat &DR, bool &start)
{
    
    if(ndi == 1){						// 1D
        sigma(0) = stress[0];
        Etot(0) = stran[0];
        DEtot(0) = dstran[0];
        Lt(0,0) = ddsdde[0];
        
    }
    else if(ndi == 2){					// 2D Plane Stress
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(3) = stress[2];
        
        Etot(0) = stran[0];
        Etot(1) = stran[1];
        Etot(3) = stran[2];
        
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(3) = dstran[2];
        
        for(int i=0 ; i<3 ; i++)
        {
            for(int j=0 ; j<3 ; j++)
                Lt(j,i) = ddsdde[i*3+j];
        }
    }
    else if(ndi == 3){
        if(nshr == 1) {
            sigma(0) = stress[0];		// 2D Generalized Plane Strain (Plane Strain, Axisymetric)
            sigma(1) = stress[1];
            sigma(2) = stress[2];
            sigma(3) = stress[3];
            
            Etot(0) = stran[0];
            Etot(1) = stran[1];
            Etot(2) = stran[2];
            Etot(3) = stran[3];
            
            DEtot(0) = dstran[0];
            DEtot(1) = dstran[1];
            DEtot(2) = dstran[2];
            DEtot(3) = dstran[3];
            
            for(int i=0 ; i<4 ; i++)
            {
                for(int j=0 ; j<4 ; j++)
                    Lt(j,i) = ddsdde[i*4+j];
            }
            
        }
        else {							// 3D
            sigma(0) = stress[0];
            sigma(1) = stress[1];
            sigma(2) = stress[2];
            sigma(3) = stress[3];
            sigma(4) = stress[4];
            sigma(5) = stress[5];
            
            Etot(0) = stran[0];
            Etot(1) = stran[1];
            Etot(2) = stran[2];
            Etot(3) = stran[3];
            Etot(4) = stran[4];
            Etot(5) = stran[5];
            
            DEtot(0) = dstran[0];
            DEtot(1) = dstran[1];
            DEtot(2) = dstran[2];
            DEtot(3) = dstran[3];
            DEtot(4) = dstran[4];
            DEtot(5) = dstran[5];
            
            for(int i=0 ; i<6 ; i++)
            {
                for(int j=0 ; j<6 ; j++)
                    Lt(j,i) = ddsdde[i*6+j];
            }
            
        }
    }
    
    ///@brief rotation matrix
    DR(0,0) = drot[0];
    DR(0,1) = drot[3];
    DR(0,2) = drot[6];
    DR(1,0) = drot[1];
    DR(1,1) = drot[4];
    DR(1,2) = drot[7];
    DR(2,0) = drot[2];
    DR(2,1) = drot[5];
    DR(2,2) = drot[8];
    
    ///@brief Temperature and temperature increment creation
    T = temperature;
    DT = Dtemperature;
    
    ///@brief Time
    Time = time[1];
    DTime = dtime;
    
    ///@brief Initialization
    if(Time < 1E-12)
    {
        start = true;
    }
    else
    {
        start = false;
    }
    
    ///@brief : Pass the material properties and the internal variables
    for (int i=0; i<nprops; i++) {
        props_smart(i) = props[i];
    }
    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    for (int i=0; i<nstatev-4; i++) {
        statev_smart(i) = statev[i+4];
    }
    
}

void abaqus2smart_T(double *stress, double *ddsdde, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops, const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &Wt, vec &statev_smart, mat &DR, bool &start) {
    
    if(ndi == 1){						// 1D
        sigma(0) = stress[0];
        Etot(0) = stran[0];
        DEtot(0) = dstran[0];
        dSdE(0,0) = ddsdde[0];
        dSdT(0,0) = ddsddt[0];
        drpldE(0,0) = drplde[0];
    }
    else if(ndi == 2){					// 2D Plane Stress
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(3) = stress[2];
        
        Etot(0) = stran[0];
        Etot(1) = stran[1];
        Etot(3) = stran[2];
        
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(3) = dstran[2];
        
        for(int i=0 ; i<3 ; i++)
        {
            dSdT(0,i) = ddsddt[i];
            drpldE(0,i) = drplde[i];
            
            for(int j=0 ; j<3 ; j++)
                dSdE(j,i) = ddsdde[i*3+j];
        }
    }
    else if(ndi == 3){
        if(nshr == 1) {
            sigma(0) = stress[0];		// 2D Generalized Plane Strain (Plane Strain, Axisymetric)
            sigma(1) = stress[1];
            sigma(2) = stress[2];
            sigma(3) = stress[3];
            
            Etot(0) = stran[0];
            Etot(1) = stran[1];
            Etot(2) = stran[2];
            Etot(3) = stran[3];
            
            DEtot(0) = dstran[0];
            DEtot(1) = dstran[1];
            DEtot(2) = dstran[2];
            DEtot(3) = dstran[3];
            
            for(int i=0 ; i<4 ; i++)
            {
                dSdT(0,i) = ddsddt[i];
                drpldE(0,i) = drplde[i];
                
                for(int j=0 ; j<4 ; j++)
                    dSdE(j,i) = ddsdde[i*4+j];
            }
            
        }
        else {							// 3D
            sigma(0) = stress[0];
            sigma(1) = stress[1];
            sigma(2) = stress[2];
            sigma(3) = stress[3];
            sigma(4) = stress[4];
            sigma(5) = stress[5];
            
            Etot(0) = stran[0];
            Etot(1) = stran[1];
            Etot(2) = stran[2];
            Etot(3) = stran[3];
            Etot(4) = stran[4];
            Etot(5) = stran[5];
            
            DEtot(0) = dstran[0];
            DEtot(1) = dstran[1];
            DEtot(2) = dstran[2];
            DEtot(3) = dstran[3];
            DEtot(4) = dstran[4];
            DEtot(5) = dstran[5];
            
            for(int i=0 ; i<6 ; i++)
            {
                dSdT(0,i) = ddsddt[i];
                drpldE(0,i) = drplde[i];                
                
                for(int j=0 ; j<6 ; j++)
                    dSdE(j,i) = ddsdde[i*6+j];
            }
        }
    }
    
    drpldT(0) = drpldt;
    
    ///@brief rotation matrix
    DR(0,0) = drot[0];
    DR(0,1) = drot[3];
    DR(0,2) = drot[6];
    DR(1,0) = drot[1];
    DR(1,1) = drot[4];
    DR(1,2) = drot[7];
    DR(2,0) = drot[2];
    DR(2,1) = drot[5];
    DR(2,2) = drot[8];
    
    ///@brief Temperature and temperature increment creation
    T = temperature;
    DT = Dtemperature;
    
    ///@brief Time
    Time = time[1];
    DTime = dtime;
    
    ///@brief Initialization
    if(Time < 1E-12)
    {
        start = true;
    }
    else
    {
        start = false;
    }
    
    ///@brief : Pass the material properties and the internal variables
    for (int i=0; i<nprops; i++) {
        props_smart(i) = props[i];
    }
    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    for (int i=0; i<3; i++) {
        Wt(i) = statev[i+4];
    }
    for (int i=0; i<nstatev-4; i++) {
        statev_smart(i) = statev[i+7];
    }
    
}
    
void select_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
{
    UNUSED(solver_type);
    std::map<string, int> list_umat;
    list_umat = {{"ELISO",1},{"ELIST",2},{"ELORT",3},{"EPICP",4},{"EPKCP",5},{"ZENER",6},{"ZENNK",7},{"PRONK",8},{"SMAUT",9}};

    rve.global2local();
    auto umat_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_local);
    
    switch (list_umat[rve.sptr_matprops->umat_name]) {
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
    
void select_umat_M(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
{
	
    std::map<string, int> list_umat;
    
    list_umat = {{"ELISO",1},{"ELIST",2},{"ELORT",3},{"EPICP",4},{"EPKCP",5},{"SMAUT",6},{"LLDM0",8},{"ZENER",10},{"ZENNK",11},{"PRONK",12},{"EPHIC",17},{"EPHIN",18},{"SMAMO",19},{"SMAMC",20},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104}};
    
        rve.global2local();
        auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    
        switch (list_umat[rve.sptr_matprops->umat_name]) {
                
            case 1: {
                umat_elasticity_iso(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
             case 2: {
                umat_elasticity_trans_iso(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
             case 3: {
                umat_elasticity_ortho(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
             case 4: {
                umat_plasticity_iso_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
             case 5: {
                umat_plasticity_kin_iso_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
            case 6: {
                umat_sma_unified_T(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 8: {
                umat_damage_LLD_0(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
            }
            case 10: {
                umat_zener_fast(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 11: {
                umat_zener_Nfast(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 12: {
                umat_prony_Nfast(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 17: {
                umat_plasticity_hill_isoh_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 18: {
                umat_plasticity_hill_isoh_CCP_N(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 19: {
                umat_sma_mono(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 20: {
                umat_sma_mono_cubic(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
                break;
            }
            case 100: case 101: case 103: case 104: {
                umat_multi(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt, list_umat[rve.sptr_matprops->umat_name]);
                break;
            }
            default: {
                cout << "Error: The choice of Umat could not be found in the umat library :" << rve.sptr_matprops->umat_name << "\n";
                exit(0);
            }
        }
        rve.local2global();

}
    
void run_umat_T(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    
    tnew_dt = 1.;
    
    select_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
    
    if (Time + DTime > limit) {
        start = false;
    }
}

void run_umat_M(phase_characteristics &rve, const mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt)
{
    
    tnew_dt = 1.;
    
    select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
    
    if (Time + DTime > limit) {
        start = false;
    }
}
    
void smart2abaqus_M(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const vec &Wm, const mat &Lt)
{
    
    if(ndi == 1) {							// 1D
        stress[0] = sigma(0);
        ddsdde[0] = Lt(0,0);
    }
    else if(ndi == 2) {						// 2D Plane Stress
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(3);
        
        assert(Lt(2,2) > 0.);
        ddsdde[0] = Lt(0,0)-Lt(0,2)*Lt(2,0)/Lt(2,2);
        ddsdde[4] = Lt(1,1)-Lt(1,2)*Lt(2,1)/Lt(2,2);
        ddsdde[3] = Lt(0,1)-Lt(0,2)*Lt(2,1)/Lt(2,2);
        ddsdde[1] = Lt(1,0)-Lt(1,2)*Lt(2,0)/Lt(2,2);
        ddsdde[6] = Lt(0,3)-Lt(0,2)*Lt(2,3)/Lt(2,2);
        ddsdde[7] = Lt(1,3)-Lt(1,2)*Lt(2,3)/Lt(2,2);
        ddsdde[2] = Lt(3,0)-Lt(3,2)*Lt(2,0)/Lt(2,2);
        ddsdde[5] = Lt(3,1)-Lt(3,2)*Lt(2,1)/Lt(2,2);
        ddsdde[8] = Lt(3,3)-Lt(3,2)*Lt(2,3)/Lt(2,2);
    }
    else if(ndi == 3){
        if (nshr == 1) {					// 2D Generalized Plane Strain (Plane Strain, Axisymetric)
            stress[0] = sigma(0);
            stress[1] = sigma(1);
            stress[2] = sigma(2);
            stress[3] = sigma(3);
            
            ddsdde[0] = Lt(0,0);
            ddsdde[4] = Lt(0,1);
            ddsdde[8] = Lt(0,2);
            ddsdde[3] = Lt(0,3);
            ddsdde[1] = Lt(1,0);
            ddsdde[5] = Lt(1,1);
            ddsdde[9] = Lt(1,2);
            ddsdde[7] = Lt(1,3);
            ddsdde[2] = Lt(2,0);
            ddsdde[6] = Lt(2,1);
            ddsdde[10] = Lt(2,2);
            ddsdde[11] = Lt(2,3);
            ddsdde[12] = Lt(3,0);
            ddsdde[13] = Lt(3,1);
            ddsdde[14] = Lt(3,2);
            ddsdde[15] = Lt(3,3);
        }
        else {								// 3D
            stress[0] = sigma(0);
            stress[1] = sigma(1);
            stress[2] = sigma(2);
            stress[3] = sigma(3);
            stress[4] = sigma(4);
            stress[5] = sigma(5);
            
            ddsdde[0] = Lt(0,0);
            ddsdde[6] = Lt(0,1);
            ddsdde[12] = Lt(0,2);
            ddsdde[18] = Lt(0,3);
            ddsdde[24] = Lt(0,4);
            ddsdde[30] = Lt(0,5);
            ddsdde[1] = Lt(1,0);
            ddsdde[7] = Lt(1,1);
            ddsdde[13] = Lt(1,2);
            ddsdde[19] = Lt(1,3);
            ddsdde[25] = Lt(1,4);
            ddsdde[31] = Lt(1,5);
            ddsdde[2] = Lt(2,0);
            ddsdde[8] = Lt(2,1);
            ddsdde[14] = Lt(2,2);
            ddsdde[20] = Lt(2,3);
            ddsdde[26] = Lt(2,4);
            ddsdde[32] = Lt(2,5);
            ddsdde[3] = Lt(3,0);
            ddsdde[9] = Lt(3,1);
            ddsdde[15] = Lt(3,2);
            ddsdde[21] = Lt(3,3);
            ddsdde[27] = Lt(3,4);
            ddsdde[33] = Lt(3,5);
            ddsdde[4] = Lt(4,0);
            ddsdde[10] = Lt(4,1);
            ddsdde[16] = Lt(4,2);
            ddsdde[22] = Lt(4,3);
            ddsdde[28] = Lt(4,4);
            ddsdde[34] = Lt(4,5);
            ddsdde[5] = Lt(5,0);
            ddsdde[11] = Lt(5,1);
            ddsdde[17] = Lt(5,2);
            ddsdde[23] = Lt(5,3);
            ddsdde[29] = Lt(5,4);
            ddsdde[35] = Lt(5,5);
        }
    }
    
    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+4] = statev_smart(i);
    }
}
    
void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, double *drplde, double &drpldt, double &rpl, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const double &r, const vec &Wm, const vec &Wt, const mat &dSdE, const mat &dSdT, const mat &drpldE, const mat &drpldT) {
    
    if(ndi == 1) {							// 1D
        stress[0] = sigma(0);
        ddsdde[0] = dSdE(0,0);
        ddsddt[0] = dSdT(0,0);
        drplde[0] = drpldE(0,0);
    }
    else if(ndi == 2) {						// 2D Plane Stress
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(3);
        
        ddsddt[0] = dSdT(0,0);
        ddsddt[1] = dSdT(0,1);
        ddsddt[2] = dSdT(0,2);
        
        drplde[0] = drpldE(0,0);
        drplde[1] = drpldE(0,1);
        drplde[2] = drpldE(0,2);
        
        assert(dSdE(2,2) > 0.);
        ddsdde[0] = dSdE(0,0)-dSdE(0,2)*dSdE(2,0)/dSdE(2,2);
        ddsdde[4] = dSdE(1,1)-dSdE(1,2)*dSdE(2,1)/dSdE(2,2);
        ddsdde[3] = dSdE(0,1)-dSdE(0,2)*dSdE(2,1)/dSdE(2,2);
        ddsdde[1] = dSdE(1,0)-dSdE(1,2)*dSdE(2,0)/dSdE(2,2);
        ddsdde[6] = dSdE(0,3)-dSdE(0,2)*dSdE(2,3)/dSdE(2,2);
        ddsdde[7] = dSdE(1,3)-dSdE(1,2)*dSdE(2,3)/dSdE(2,2);
        ddsdde[2] = dSdE(3,0)-dSdE(3,2)*dSdE(2,0)/dSdE(2,2);
        ddsdde[5] = dSdE(3,1)-dSdE(3,2)*dSdE(2,1)/dSdE(2,2);
        ddsdde[8] = dSdE(3,3)-dSdE(3,2)*dSdE(2,3)/dSdE(2,2);
    }
    else if(ndi == 3){
        if (nshr == 1) {					// 2D Generalized Plane Strain (Plane Strain, Axisymetric)
            stress[0] = sigma(0);
            stress[1] = sigma(1);
            stress[2] = sigma(2);
            stress[3] = sigma(3);
            
            ddsddt[0] = dSdT(0,0);
            ddsddt[1] = dSdT(0,1);
            ddsddt[2] = dSdT(0,2);
            ddsddt[3] = dSdT(0,3);
            
            drplde[0] = drpldE(0,0);
            drplde[1] = drpldE(0,1);
            drplde[2] = drpldE(0,2);
            drplde[3] = drpldE(0,3);
            
            ddsdde[0] = dSdE(0,0);
            ddsdde[4] = dSdE(0,1);
            ddsdde[8] = dSdE(0,2);
            ddsdde[3] = dSdE(0,3);
            ddsdde[1] = dSdE(1,0);
            ddsdde[5] = dSdE(1,1);
            ddsdde[9] = dSdE(1,2);
            ddsdde[7] = dSdE(1,3);
            ddsdde[2] = dSdE(2,0);
            ddsdde[6] = dSdE(2,1);
            ddsdde[10] = dSdE(2,2);
            ddsdde[11] = dSdE(2,3);
            ddsdde[12] = dSdE(3,0);
            ddsdde[13] = dSdE(3,1);
            ddsdde[14] = dSdE(3,2);
            ddsdde[15] = dSdE(3,3);
        }
        else {								// 3D
            stress[0] = sigma(0);
            stress[1] = sigma(1);
            stress[2] = sigma(2);
            stress[3] = sigma(3);
            stress[4] = sigma(4);
            stress[5] = sigma(5);
            
            ddsddt[0] = dSdT(0,0);
            ddsddt[1] = dSdT(0,1);
            ddsddt[2] = dSdT(0,2);
            ddsddt[3] = dSdT(0,3);
            ddsddt[4] = dSdT(0,4);
            ddsddt[5] = dSdT(0,5);
            
            drplde[0] = drpldE(0,0);
            drplde[1] = drpldE(0,1);
            drplde[2] = drpldE(0,2);
            drplde[3] = drpldE(0,3);
            drplde[4] = drpldE(0,4);
            drplde[5] = drpldE(0,5);
            
            ddsdde[0] = dSdE(0,0);
            ddsdde[6] = dSdE(0,1);
            ddsdde[12] = dSdE(0,2);
            ddsdde[18] = dSdE(0,3);
            ddsdde[24] = dSdE(0,4);
            ddsdde[30] = dSdE(0,5);
            ddsdde[1] = dSdE(1,0);
            ddsdde[7] = dSdE(1,1);
            ddsdde[13] = dSdE(1,2);
            ddsdde[19] = dSdE(1,3);
            ddsdde[25] = dSdE(1,4);
            ddsdde[31] = dSdE(1,5);
            ddsdde[2] = dSdE(2,0);
            ddsdde[8] = dSdE(2,1);
            ddsdde[14] = dSdE(2,2);
            ddsdde[20] = dSdE(2,3);
            ddsdde[26] = dSdE(2,4);
            ddsdde[32] = dSdE(2,5);
            ddsdde[3] = dSdE(3,0);
            ddsdde[9] = dSdE(3,1);
            ddsdde[15] = dSdE(3,2);
            ddsdde[21] = dSdE(3,3);
            ddsdde[27] = dSdE(3,4);
            ddsdde[33] = dSdE(3,5);
            ddsdde[4] = dSdE(4,0);
            ddsdde[10] = dSdE(4,1);
            ddsdde[16] = dSdE(4,2);
            ddsdde[22] = dSdE(4,3);
            ddsdde[28] = dSdE(4,4);
            ddsdde[34] = dSdE(4,5);
            ddsdde[5] = dSdE(5,0);
            ddsdde[11] = dSdE(5,1);
            ddsdde[17] = dSdE(5,2);
            ddsdde[23] = dSdE(5,3);
            ddsdde[29] = dSdE(5,4);
            ddsdde[35] = dSdE(5,5);
        }
    }
    
    drpldt = drpldT(0,0);
    rpl = r;
    
    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    for (int i=0; i<3; i++) {
        statev[i+4] = Wt(i);
    }
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+7] = statev_smart(i);
    }
}
    
} //namespace simcoon
