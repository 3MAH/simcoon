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
#define BOOST_DLL_USE_STD_FS  // Forces Boost.DLL to use std::filesystem
#include <boost/dll.hpp>
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
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono_cubic.hpp>
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

void abaqus2smart_M_light(const double *stress, const double *ddsdde, const int &nstatev, double *statev, const int &ndi, const int &nshr, vec &sigma, mat &Lt, vec &Wm, vec &statev_smart)
{
    
    if(ndi == 1){						// 1D
        sigma(0) = stress[0];
        Lt(0,0) = ddsdde[0];
    }
    else if(ndi == 2){					// 2D Plane Stress
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(3) = stress[2];
        
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
            
            Lt(0,0) = ddsdde[0];
            Lt(0,1) = ddsdde[4];
            Lt(0,2) = ddsdde[8];
            Lt(0,3) = ddsdde[3];
            Lt(1,0) = ddsdde[1];
            Lt(1,1) = ddsdde[5];
            Lt(1,2) = ddsdde[9];
            Lt(1,3) = ddsdde[7];
            Lt(2,0) = ddsdde[2];
            Lt(2,1) = ddsdde[6];
            Lt(2,2) = ddsdde[10];
            Lt(2,3) = ddsdde[11];
            Lt(3,0) = ddsdde[12];
            Lt(3,1) = ddsdde[13];
            Lt(3,2) = ddsdde[14];
            Lt(3,3) = ddsdde[15];
        }
        else {							// 3D
            sigma(0) = stress[0];
            sigma(1) = stress[1];
            sigma(2) = stress[2];
            sigma(3) = stress[3];
            sigma(4) = stress[4];
            sigma(5) = stress[5];
            
            Lt(0,0) = ddsdde[0];
            Lt(0,1) = ddsdde[6];
            Lt(0,2) = ddsdde[12];
            Lt(0,3) = ddsdde[18];
            Lt(0,4) = ddsdde[24];
            Lt(0,5) = ddsdde[30];
            Lt(1,0) = ddsdde[1];
            Lt(1,1) = ddsdde[7];
            Lt(1,2) = ddsdde[13];
            Lt(1,3) = ddsdde[19];
            Lt(1,4) = ddsdde[25];
            Lt(1,5) = ddsdde[31];
            Lt(2,0) = ddsdde[2];
            Lt(2,1) = ddsdde[8];
            Lt(2,2) = ddsdde[14];
            Lt(2,3) = ddsdde[20];
            Lt(2,4) = ddsdde[26];
            Lt(2,5) = ddsdde[32];
            Lt(3,0) = ddsdde[3];
            Lt(3,1) = ddsdde[9];
            Lt(3,2) = ddsdde[15];
            Lt(3,3) = ddsdde[21];
            Lt(3,4) = ddsdde[27];
            Lt(3,5) = ddsdde[33];
            Lt(4,0) = ddsdde[4];
            Lt(4,1) = ddsdde[10];
            Lt(4,2) = ddsdde[16];
            Lt(4,3) = ddsdde[22];
            Lt(4,4) = ddsdde[28];
            Lt(4,5) = ddsdde[34];
            Lt(5,0) = ddsdde[5];
            Lt(5,1) = ddsdde[11];
            Lt(5,2) = ddsdde[17];
            Lt(5,3) = ddsdde[23];
            Lt(5,4) = ddsdde[29];
            Lt(5,5) = ddsdde[35];
        }
    }
    
    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    for (int i=0; i<nstatev-4; i++) {
        statev_smart(i) = statev[i+4];
    }
    
}


void abaqus2smart_M(const double *stress, const double *ddsdde, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops,const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &Lt, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &statev_smart, mat &DR, bool &start)
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
            
            Lt(0,0) = ddsdde[0];
            Lt(0,1) = ddsdde[4];
            Lt(0,2) = ddsdde[8];
            Lt(0,3) = ddsdde[3];
            Lt(1,0) = ddsdde[1];
            Lt(1,1) = ddsdde[5];
            Lt(1,2) = ddsdde[9];
            Lt(1,3) = ddsdde[7];
            Lt(2,0) = ddsdde[2];
            Lt(2,1) = ddsdde[6];
            Lt(2,2) = ddsdde[10];
            Lt(2,3) = ddsdde[11];
            Lt(3,0) = ddsdde[12];
            Lt(3,1) = ddsdde[13];
            Lt(3,2) = ddsdde[14];
            Lt(3,3) = ddsdde[15];
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
            
            Lt(0,0) = ddsdde[0];
            Lt(0,1) = ddsdde[6];
            Lt(0,2) = ddsdde[12];
            Lt(0,3) = ddsdde[18];
            Lt(0,4) = ddsdde[24];
            Lt(0,5) = ddsdde[30];
            Lt(1,0) = ddsdde[1];
            Lt(1,1) = ddsdde[7];
            Lt(1,2) = ddsdde[13];
            Lt(1,3) = ddsdde[19];
            Lt(1,4) = ddsdde[25];
            Lt(1,5) = ddsdde[31];
            Lt(2,0) = ddsdde[2];
            Lt(2,1) = ddsdde[8];
            Lt(2,2) = ddsdde[14];
            Lt(2,3) = ddsdde[20];
            Lt(2,4) = ddsdde[26];
            Lt(2,5) = ddsdde[32];
            Lt(3,0) = ddsdde[3];
            Lt(3,1) = ddsdde[9];
            Lt(3,2) = ddsdde[15];
            Lt(3,3) = ddsdde[21];
            Lt(3,4) = ddsdde[27];
            Lt(3,5) = ddsdde[33];
            Lt(4,0) = ddsdde[4];
            Lt(4,1) = ddsdde[10];
            Lt(4,2) = ddsdde[16];
            Lt(4,3) = ddsdde[22];
            Lt(4,4) = ddsdde[28];
            Lt(4,5) = ddsdde[34];
            Lt(5,0) = ddsdde[5];
            Lt(5,1) = ddsdde[11];
            Lt(5,2) = ddsdde[17];
            Lt(5,3) = ddsdde[23];
            Lt(5,4) = ddsdde[29];
            Lt(5,5) = ddsdde[35];
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

void abaqus2smart_T(const double *stress, const double *ddsdde, const double *ddsddt, const double *drplde, const double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops, const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &Wt, vec &statev_smart, mat &DR, bool &start) {
    
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
            dSdT(i,0) = ddsddt[i];
            drpldE(i,0) = drplde[i];
            
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
                dSdT(i,0) = ddsddt[i];
                drpldE(i,0) = drplde[i];
            }
            
            dSdE(0,0) = ddsdde[0];
            dSdE(0,1) = ddsdde[4];
            dSdE(0,2) = ddsdde[8];
            dSdE(0,3) = ddsdde[3];
            dSdE(1,0) = ddsdde[1];
            dSdE(1,1) = ddsdde[5];
            dSdE(1,2) = ddsdde[9];
            dSdE(1,3) = ddsdde[7];
            dSdE(2,0) = ddsdde[2];
            dSdE(2,1) = ddsdde[6];
            dSdE(2,2) = ddsdde[10];
            dSdE(2,3) = ddsdde[11];
            dSdE(3,0) = ddsdde[12];
            dSdE(3,1) = ddsdde[13];
            dSdE(3,2) = ddsdde[14];
            dSdE(3,3) = ddsdde[15];
            
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
                dSdT(i,0) = ddsddt[i];
                drpldE(i,0) = drplde[i];
            }
            
            dSdE(0,0) = ddsdde[0];
            dSdE(0,1) = ddsdde[6];
            dSdE(0,2) = ddsdde[12];
            dSdE(0,3) = ddsdde[18];
            dSdE(0,4) = ddsdde[24];
            dSdE(0,5) = ddsdde[30];
            dSdE(1,0) = ddsdde[1];
            dSdE(1,1) = ddsdde[7];
            dSdE(1,2) = ddsdde[13];
            dSdE(1,3) = ddsdde[19];
            dSdE(1,4) = ddsdde[25];
            dSdE(1,5) = ddsdde[31];
            dSdE(2,0) = ddsdde[2];
            dSdE(2,1) = ddsdde[8];
            dSdE(2,2) = ddsdde[14];
            dSdE(2,3) = ddsdde[20];
            dSdE(2,4) = ddsdde[26];
            dSdE(2,5) = ddsdde[32];
            dSdE(3,0) = ddsdde[3];
            dSdE(3,1) = ddsdde[9];
            dSdE(3,2) = ddsdde[15];
            dSdE(3,3) = ddsdde[21];
            dSdE(3,4) = ddsdde[27];
            dSdE(3,5) = ddsdde[33];
            dSdE(4,0) = ddsdde[4];
            dSdE(4,1) = ddsdde[10];
            dSdE(4,2) = ddsdde[16];
            dSdE(4,3) = ddsdde[22];
            dSdE(4,4) = ddsdde[28];
            dSdE(4,5) = ddsdde[34];
            dSdE(5,0) = ddsdde[5];
            dSdE(5,1) = ddsdde[11];
            dSdE(5,2) = ddsdde[17];
            dSdE(5,3) = ddsdde[23];
            dSdE(5,4) = ddsdde[29];
            dSdE(5,5) = ddsdde[35];
            
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
    for (int i=0; i<nstatev-7; i++) {
        statev_smart(i) = statev[i+7];
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
                umat_elasticity_iso(umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
            case 3: {
                umat_elasticity_trans_iso(umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
            case 4: {
                umat_elasticity_ortho(umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
            case 5: {
                umat_hypoelasticity_ortho(umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
            }                                    
             case 6: {
                umat_plasticity_iso_CCP(umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
             case 7: {
                umat_plasticity_kin_iso_CCP(umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                 break;
             }
            case 8: {
                umat_saint_venant(umat_M->etot, umat_M->Detot, umat_M->F0, umat_M->F1, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
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
    
    list_umat = {{"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"SMAUT",8},{"LLDM0",9},{"ZENER",10},{"ZENNK",11},{"PRONK",12},{"EPHIL",17},{"EPHAC",18},{"EPANI",19},{"EPDFA",20},{"EPCHG",21},{"EPHIN",22},{"SMAMO",23},{"SMAMC",24},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104}};
    
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);

    switch (list_umat[rve.sptr_matprops->umat_name]) {

        case 0: {
            //umat_external(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);

            fs::path lib_path("external");  // Path to the directory with our plugin library
            fs::path boost_lib_path = lib_path / "umat_plugin_ext";            

            boost::dll::shared_library lib(boost_lib_path, boost::dll::load_mode::append_decorations);
            auto& external_umat = lib.get<umat_plugin_ext_api>("external_umat");     

            external_umat.umat_external_M(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);

            break;
        }
        case 1: {
            //
            fs::path lib_path("external");  // Path to the directory with our plugin library
            fs::path boost_lib_path = lib_path / "umat_plugin_aba";            

            boost::dll::shared_library lib2(boost_lib_path, boost::dll::load_mode::append_decorations);
            auto& abaqus_umat = lib2.get<umat_plugin_aba_api>("abaqus_umat");
            abaqus_umat.umat_abaqus(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 2: {
            umat_elasticity_iso(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
            }
        case 3: {
            umat_elasticity_trans_iso(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 4: {
            umat_elasticity_ortho(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 5: {
            umat_plasticity_iso_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 6: {
            umat_plasticity_kin_iso_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 7: {
            umat_plasticity_chaboche_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
            break;
        }
        case 8: {
            umat_sma_unified_T(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 9: {
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
            umat_hill_chaboche_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
        }
        case 19: {
            umat_ani_chaboche_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
        }
        case 20: {
            umat_dfa_chaboche_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
        }
        case 21: {
            umat_generic_chaboche_CCP(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);
                break;
        }                                                
        case 22: {
            umat_plasticity_hill_isoh_CCP_N(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }            
        case 23: {
            umat_sma_mono(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
            break;
        }
        case 24: {
            umat_sma_mono_cubic(umat_M->Etot, umat_M->DEtot, umat_M->sigma, umat_M->Lt, umat_M->L, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, tnew_dt);
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

void smart2abaqus_M_full(double *stress, double *ddsdde, double *stran, double *dstran, double *time, double &dtime, double &temperature, double &Dtemperature, int &nprops, double *props,  int &nstatev, double *statev, const int &ndi, const int &nshr, double *drot, const vec &sigma, const mat &Lt, const vec &Etot, const vec &DEtot, const double &T, const double &DT, const double &Time, const double &DTime, const vec &props_smart, const vec &Wm, const vec &statev_smart, const mat &DR, bool &start)
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
            
            stran[0] = Etot(0);
            stran[1] = Etot(1);
            stran[2] = Etot(2);
            stran[3] = Etot(3);
            stran[4] = Etot(4);
            stran[5] = Etot(5);
            
            dstran[0] = DEtot(0);
            dstran[1] = DEtot(1);
            dstran[2] = DEtot(2);
            dstran[3] = DEtot(3);
            dstran[4] = DEtot(4);
            dstran[5] = DEtot(5);
            
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

    
    ///@brief rotation matrix
    drot[0] = DR(0,0);
    drot[3] = DR(0,1);
    drot[6] = DR(0,2);
    drot[1] = DR(1,0);
    drot[4] = DR(1,1);
    drot[7] = DR(1,2);
    drot[2] = DR(2,0);
    drot[5] = DR(2,1);
    drot[8] = DR(2,2);
    
    ///@brief Temperature and temperature increment creation
    temperature = T;
    Dtemperature = DT;
    
    ///@brief Time
    time[1] = Time;
    dtime = DTime;
    
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
        props[i] = props_smart(i);
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
        ddsddt[1] = dSdT(1,0);
        ddsddt[2] = dSdT(2,0);
        
        drplde[0] = drpldE(0,0);
        drplde[1] = drpldE(1,0);
        drplde[2] = drpldE(2,0);
        
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
            ddsddt[1] = dSdT(1,0);
            ddsddt[2] = dSdT(2,0);
            ddsddt[3] = dSdT(3,0);
            
            drplde[0] = drpldE(0,0);
            drplde[1] = drpldE(1,0);
            drplde[2] = drpldE(2,0);
            drplde[3] = drpldE(3,0);
            
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
            ddsddt[1] = dSdT(1,0);
            ddsddt[2] = dSdT(2,0);
            ddsddt[3] = dSdT(3,0);
            ddsddt[4] = dSdT(4,0);
            ddsddt[5] = dSdT(5,0);
            
            drplde[0] = drpldE(0,0);
            drplde[1] = drpldE(1,0);
            drplde[2] = drpldE(2,0);
            drplde[3] = drpldE(3,0);
            drplde[4] = drpldE(4,0);
            drplde[5] = drpldE(5,0);
            
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
