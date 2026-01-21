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

///@file fea_transfer.cpp
///@brief Functions to transfer data between FEA software formats and simcoon internal format
///@version 1.0

#include <iostream>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

//=============================================================================
// ABAQUS UMAT TRANSFER FUNCTIONS
//=============================================================================

void abaqus2smart_M_light(const double *stress, const double *ddsdde, const int &nstatev, double *statev, const int &ndi, const int &nshr, vec &sigma, mat &Lt, vec &Wm, vec &statev_smart)
{
    
    if(ndi == 1){                       // 1D
        sigma(0) = stress[0];
        Lt(0,0) = ddsdde[0];
    }
    else if(ndi == 2){                  // 2D Plane Stress
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
            sigma(0) = stress[0];       // 2D Generalized Plane Strain (Plane Strain, Axisymetric)
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
        else {                          // 3D
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
    
    if(ndi == 1){                       // 1D
        sigma(0) = stress[0];
        Etot(0) = stran[0];
        DEtot(0) = dstran[0];
        Lt(0,0) = ddsdde[0];
        
    }
    else if(ndi == 2){                  // 2D Plane Stress
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
            sigma(0) = stress[0];       // 2D Generalized Plane Strain (Plane Strain, Axisymetric)
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
        else {                          // 3D
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
    
    if(ndi == 1){                       // 1D
        sigma(0) = stress[0];
        Etot(0) = stran[0];
        DEtot(0) = dstran[0];
        dSdE(0,0) = ddsdde[0];
        dSdT(0,0) = ddsddt[0];
        drpldE(0,0) = drplde[0];
    }
    else if(ndi == 2){                  // 2D Plane Stress
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
            sigma(0) = stress[0];       // 2D Generalized Plane Strain (Plane Strain, Axisymetric)
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
        else {                          // 3D
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
    
void smart2abaqus_M(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const vec &Wm, const mat &Lt)
{
    
    if(ndi == 1) {                          // 1D
        stress[0] = sigma(0);
        ddsdde[0] = Lt(0,0);
    }
    else if(ndi == 2) {                     // 2D Plane Stress
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
        if (nshr == 1) {                    // 2D Generalized Plane Strain (Plane Strain, Axisymetric)
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
        else {                              // 3D
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
    
    if(ndi == 1) {                          // 1D
        stress[0] = sigma(0);
        ddsdde[0] = Lt(0,0);
    }
    else if(ndi == 2) {                     // 2D Plane Stress
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
        if (nshr == 1) {                    // 2D Generalized Plane Strain (Plane Strain, Axisymetric)
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
        else {                              // 3D
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
    
    if(ndi == 1) {                          // 1D
        stress[0] = sigma(0);
        ddsdde[0] = dSdE(0,0);
        ddsddt[0] = dSdT(0,0);
        drplde[0] = drpldE(0,0);
    }
    else if(ndi == 2) {                     // 2D Plane Stress
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
        if (nshr == 1) {                    // 2D Generalized Plane Strain (Plane Strain, Axisymetric)
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
        else {                              // 3D
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

//=============================================================================
// ANSYS USERMAT TRANSFER FUNCTIONS
//=============================================================================

void ansys2smart_M(const double *stress, const double *dstran,
    const double &sedEl, const double &sedPl, const double &epseq,
    const double *statev, const double *props,
    const double &Time, const double &DTime,
    const double &temperature, const double &Dtemperature,
    const int &ncomp, const int &nprops, const int &nstatev,
    vec &sigma, vec &DEtot, vec &Wm, 
    vec &statev_smart, vec &props_smart,
    double &T, double &DT)
{
    // Ansys uses Voigt notation: (11, 22, 33, 12, 23, 13)
    // simcoon uses:              (11, 22, 33, 12, 13, 23)
    // So components 4 and 5 are swapped
    
    if (ncomp == 6) {  // 3D
        // Stress: swap components 4 and 5
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(2) = stress[2];
        sigma(3) = stress[3];
        sigma(4) = stress[5];  // Ansys 13 → simcoon index 4
        sigma(5) = stress[4];  // Ansys 23 → simcoon index 5
        
        // Strain increment: swap components 4 and 5
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(2) = dstran[2];
        DEtot(3) = dstran[3];
        DEtot(4) = dstran[5];  // Ansys 13 → simcoon index 4
        DEtot(5) = dstran[4];  // Ansys 23 → simcoon index 5
    }
    else if (ncomp == 4) {  // 2D (plane strain/axisymmetric)
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(2) = stress[2];
        sigma(3) = stress[3];
        
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(2) = dstran[2];
        DEtot(3) = dstran[3];
    }
    else if (ncomp == 3) {  // 2D plane stress
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(3) = stress[2];
        
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(3) = dstran[2];
    }
    else if (ncomp == 1) {  // 1D
        sigma(0) = stress[0];
        DEtot(0) = dstran[0];
    }
    
    // Temperature
    T = temperature;
    DT = Dtemperature;
    
    // Material properties
    for (int i = 0; i < nprops; i++) {
        props_smart(i) = props[i];
    }
    
    // Work quantities from sedEl and sedPl
    Wm(0) = sedEl + sedPl;  // Total strain energy
    Wm(1) = sedEl;          // Elastic (reversible) strain energy
    Wm(2) = sedPl;          // Plastic (irreversible) strain energy
    Wm(3) = 0.;             // Damage energy (not directly available)
    
    // State variables (first 4 are Wm, then user statev)
    int nstatev_smart = nstatev - 4;
    for (int i = 0; i < nstatev_smart; i++) {
        statev_smart(i) = statev[i + 4];
    }
}

void smart2ansys_M(double *stress, double *ddsdde,
    double &sedEl, double &sedPl, double *statev,
    const int &ncomp,
    const vec &sigma, const mat &Lt,
    const vec &Wm, const vec &statev_smart)
{
    // Ansys uses Voigt notation: (11, 22, 33, 12, 23, 13)
    // simcoon uses:              (11, 22, 33, 12, 13, 23)
    // So components 4 and 5 are swapped
    
    if (ncomp == 6) {  // 3D
        // Stress: swap components 4 and 5
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(2);
        stress[3] = sigma(3);
        stress[4] = sigma(5);  // simcoon 23 → Ansys index 4
        stress[5] = sigma(4);  // simcoon 13 → Ansys index 5
        
        // Tangent matrix: need to swap rows 4↔5 and columns 4↔5
        // Create permutation: simcoon indices [0,1,2,3,4,5] → Ansys [0,1,2,3,5,4]
        int perm[6] = {0, 1, 2, 3, 5, 4};
        
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                // Ansys uses row-major storage
                ddsdde[i * 6 + j] = Lt(perm[i], perm[j]);
            }
        }
    }
    else if (ncomp == 4) {  // 2D (plane strain/axisymmetric)
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(2);
        stress[3] = sigma(3);
        
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                ddsdde[i * 4 + j] = Lt(i, j);
            }
        }
    }
    else if (ncomp == 3) {  // 2D plane stress
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(3);
        
        // Condense tangent for plane stress
        assert(Lt(2,2) > 0.);
        ddsdde[0] = Lt(0,0) - Lt(0,2)*Lt(2,0)/Lt(2,2);
        ddsdde[4] = Lt(1,1) - Lt(1,2)*Lt(2,1)/Lt(2,2);
        ddsdde[3] = Lt(0,1) - Lt(0,2)*Lt(2,1)/Lt(2,2);
        ddsdde[1] = Lt(1,0) - Lt(1,2)*Lt(2,0)/Lt(2,2);
        ddsdde[6] = Lt(0,3) - Lt(0,2)*Lt(2,3)/Lt(2,2);
        ddsdde[7] = Lt(1,3) - Lt(1,2)*Lt(2,3)/Lt(2,2);
        ddsdde[2] = Lt(3,0) - Lt(3,2)*Lt(2,0)/Lt(2,2);
        ddsdde[5] = Lt(3,1) - Lt(3,2)*Lt(2,1)/Lt(2,2);
        ddsdde[8] = Lt(3,3) - Lt(3,2)*Lt(2,3)/Lt(2,2);
    }
    else if (ncomp == 1) {  // 1D
        stress[0] = sigma(0);
        ddsdde[0] = Lt(0,0);
    }
    
    // Work quantities
    sedEl = Wm(1);  // Elastic strain energy
    sedPl = Wm(2);  // Plastic strain energy
    
    // State variables
    for (int i = 0; i < 4; i++) {
        statev[i] = Wm(i);
    }
    for (unsigned int i = 0; i < statev_smart.n_elem; i++) {
        statev[i + 4] = statev_smart(i);
    }
}
    
} // namespace simcoon
