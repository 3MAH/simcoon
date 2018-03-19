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

///@file umat_singleM.cpp
///@brief umat template to run smart subroutines using Abaqus
///@brief Implemented in 1D-2D-3D
///@author Chemisky & Despringre
///@version 1.0
///@date 12/04/2013

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>

///@param stress array containing the components of the stress tensor (dimension ntens)
///@param statev array containing the evolution variables (dimension nstatev)
///@param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
///@param sse unused
///@param spd unused
///@param scd unused
///@param rpl unused
///@param ddsddt array containing the thermal tangent operator
///@param drple unused
///@param drpldt unused
///@param stran array containing total strain component (dimension ntens) at the beginning of increment
///@param dstran array containing the component of total strain increment (dimension ntens)
///@param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
///@param dtime time increment
///@param temperature temperature avlue at the beginning of increment
///@param Dtemperature temperature increment
///@param predef unused
///@param dpred unused
///@param cmname user-defined material name
///@param ndi number of direct stress components
///@param nshr number of shear stress components
///@param ntens number stress and strain components
///@param nstatev number of evolution variables
///@param props array containing material properties
///@param nprops number of material properties
///@param coords coordinates of the considered point
///@param drot rotation increment matrix (dimension 3*3)
///@param pnewdt ratio of suggested new time increment
///@param celent characteristic element length
///@param dfgrd0 array containing the deformation gradient at the beginning of increment (dimension 3*3)
///@param dfgrd1 array containing the deformation gradient at the end of increment (dimension 3*3)
///@param noel element number
///@param npt integration point number
///@param layer layer number - not used
///@param kspt section point number within the current layer - not used
///@param kstep step number
///@param kinc increment number

using namespace std;
using namespace arma;
using namespace simcoon;

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc)
{
    UNUSED(sse);
    UNUSED(spd);
	UNUSED(scd);
	UNUSED(rpl);
	UNUSED(ddsddt);
	UNUSED(drplde);
	UNUSED(drpldt);
	UNUSED(predef);
	UNUSED(dpred);
	UNUSED(ntens);
	UNUSED(coords);
	UNUSED(celent);
	UNUSED(dfgrd0);
	UNUSED(dfgrd1);
	UNUSED(noel);
	UNUSED(npt);
	UNUSED(layer);
	UNUSED(kspt);
	UNUSED(kstep);
	UNUSED(kinc);
	
	bool start = false;
	double Time = 0.;
	double DTime = 0.;
    int solver_type = 0;
    
	int nstatev_smart = nstatev-4;
    vec props_smart = zeros(nprops);
    vec statev_smart = zeros(nstatev_smart);
    
	vec sigma = zeros(6);
	vec Etot = zeros(6);
	vec DEtot = zeros(6);
	mat Lt = zeros(6,6);
	mat DR = zeros(3,3);
    vec Wm = zeros(4);
    
	string umat_name(cmname);
	umat_name = umat_name.substr(0, 5);
	   
    phase_characteristics rve;
    rve.construct(0,1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_smart);
    rve.sptr_matprops->resize(nprops);
    
	abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_smart, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
    select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
    
	smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
}
