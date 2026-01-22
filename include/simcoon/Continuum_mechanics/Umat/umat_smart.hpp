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
* @file umat_smart.hpp
* @author Yves Chemisky 
* @brief Implemented in 1D-2D-3D
* @section Selection of constitutive laws and transfer to between Abaqus and simcoon formats
*/

#pragma once
#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>

extern "C"{
/**
 * @brief run a Abaqus umat_ using standard Abaqus UMAT arguments
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param sse unused
 * @param spd unused
 * @param scd unused
 * @param rpl unused
 * @param ddsddt array containing the thermal tangent operator
 * @param drplde unused
 * @param drpldt unused
 * @param stran array containing total strain component (dimension ntens) at the beginning of increment
 * @param dstran array containing the component of total strain increment (dimension ntens)
 * @param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
 * @param dtime time increment
 * @param temperature temperature value at the beginning of increment
 * @param Dtemperature temperature increment
 * @param predef unused
 * @param dpred unused
 * @param cmname user-defined material name
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param ntens number stress and strain components
 * @param nstatev number of evolution variables
 * @param props array containing material properties
 * @param nprops number of material properties
 * @param coords coordinates of the considered point
 * @param drot rotation increment matrix (dimension 3*3)
 * @param pnewdt ratio of suggested new time increment
 * @param celent characteristic element length
 * @param dfgrd0 array containing the deformation gradient at the beginning of increment (dimension 3*3)
 * @param dfgrd1 array containing the deformation gradient at the end of increment (dimension 3*3)
 * @param noel element number
 * @param npt integration point number
 * @param layer layer number - not used
 * @param kspt section point number within the current layer - not used
 * @param kstep step number
 * @param kinc increment number
 * @details Example: 
 * @code         
        ///Macroscopic state variables and control increments
        double stress[6];
        double stran[6];
        double dstran[6];
        double temperature;
        double Dtemperature;
        
        // Umat variable list unused here
        double sse = 0.;
        double spd = 0.;
        double scd = 0.;
        double rpl = 0.;
        double drpldt = 0.;
        double predef = 0.;
        double dpred = 0.;
        int layer = 0;
        int kspt = 0;
        double celent = 0.;
        double dfgrd0[9];
        double dfgrd1[3];
        double drplde[6];
        double coords = 0;
        double drot[9];
            
        // Usefull UMAT variables
        int ntens = 6;
        int noel = 1;
        int npt = 1;
        int kstep = 0;
        int kinc = 0;
        double ddsdde[36];
        double ddsddt[6];
        char cmname[5];
        int nprops = 2
        double props[nprops];
        int nstatev = 1
        double statev[nstatev];
        double pnewdt = 1.
        double time[2];
        double dtime;
        
        umat_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temperature, Dtemperature, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc);        
    }
 * @endcode
*/
	void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc);
}

namespace simcoon{

/**
 * @brief Set the size of the statev vector, required to store information of an phase_characteristics object
 * @param rve Reference to the phase characteristics object.
 * @param size Reference to the size of the statev vector
 * @details Example: 
 * @code 
    phase_characteristics rve;
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve.construct(0,1);
    natural_basis nb;
    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(3,3), zeros(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
    unsigned int statev_abaqus = 0;
    size_statev(rve, statev_abaqus);
    cout << "The Umat has a number of statev equal to: " << statev_abaqus << endl;	
 * @endcode
*/
void size_statev(phase_characteristics &rve, unsigned int &size);
    
/**
 * @brief From mechanical variables stored in a statev vector, parse a rve object that contains all information about the phases. 
 * @param rve Reference to the phase characteristics object.
 * @param pos Reference to the position in the statev vector that contain subphase information
 * @param statev Reference to the statev vector that hold the information.
 * @details Example: 
 * @code 
    unsigned int nstatev_multi = nstatev-nstatev_macro-4;
    vec statev_macro = zeros(nstatev_macro);
    vec statev_multi(&statev[nstatev_macro+4], nstatev_multi, false, false);
    
    phase_characteristics rve;
    rve.construct(0,1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_macro);
    rve.sptr_matprops->resize(nprops_macro);
    unsigned int pos=0;    
    statev_2_phases(rve, pos, statev_multi);
 * @endcode
*/
void statev_2_phases(phase_characteristics &rve, unsigned int &pos, const arma::vec &statev);

/**
 * @brief From a rve object, write a statev vector that contains all information about the phases
 * @param statev Reference to the statev vector that hold the information.
 * @param pos Reference to the position in the statev vector that contain subphase information
 * @param rve Reference to the phase characteristics.
 * @details Example: 
 * @code 
    unsigned int nstatev_multi = nstatev-nstatev_macro-4;
    vec statev_macro = zeros(nstatev_macro);
    vec statev_multi(&statev[nstatev_macro+4], nstatev_multi, false, false);
    
    phase_characteristics rve;
    rve.construct(0,1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_macro);
    rve.sptr_matprops->resize(nprops_macro);
    
    unsigned int pos=0;
    
    statev_2_phases(rve, pos, statev_multi);
	abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev_macro, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_macro, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops_macro, props_macro);
    start = true;
    select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
    phases_2_statev(statev_multi, pos, rve);
 * @endcode
*/
void phases_2_statev(arma::vec &statev, unsigned int &pos, const phase_characteristics &rve);

/**
 * @brief From a abaqus umat_ function, get the essential simcoon variables for a mechanical simulation 
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param nstatev number of evolution variables
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param sigma arma::vec containing the components of the stress tensor
 * @param Lt arma::mat containing the mechanical tangent operator
 * @param Wm variables related to mechanical work
 * @param statev_smart evolution variables in the form of an arma::vec
 * @details Example: 
 * @code 
        // Macroscopic state variables and control increments
        double stress[6];
        double stran[6];
        double dstran[6];
        double temperature;
        double Dtemperature;
        
        // Umat variable list unused here
        double sse = 0.;
        double spd = 0.;
        double scd = 0.;
        double rpl = 0.;
        double drpldt = 0.;
        double predef = 0.;
        double dpred = 0.;
        int layer = 0;
        int kspt = 0;
        double celent = 0.;
        double dfgrd0[9];
        double dfgrd1[3];
        double drplde[6];
        double coords = 0;
        double drot[9];
            
        // Usefull UMAT variables
        int ntens = 6;
        int noel = 1;
        int npt = 1;
        int kstep = 0;
        int kinc = 0;
        double ddsdde[36];
        double ddsddt[6];
        char cmname[5];
        strcpy(cmname, rve.sptr_matprops->umat_name.c_str());
        int nprops = rve.sptr_matprops->nprops;
        double props[nprops];
        int nstatev = rve.sptr_sv_global->nstatev;
        double statev[nstatev+4];                   //+4 for a mechanical response to store Wm components
        double pnewdt = tnew_dt;
        double time[2];
        double dtime;
        
        auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
        
        smart2abaqus_M_full(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, rve.sptr_matprops->props, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
//        smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
        umat_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temperature, Dtemperature, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc);
        
        abaqus2smart_M_light(stress, ddsdde, nstatev, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Wm, rve_sv_M->statev);
 * @endcode
*/
void abaqus2smart_M_light(const double *stress, const double *ddsdde, const int &nstatev, double *statev, const int &ndi, const int &nshr, arma::vec &sigma, arma::mat &Lt, arma::vec &Wm, arma::vec &statev_smart);
    
/**
 * @brief From a abaqus umat_ function, get the simcoon variables for a mechanical simulation 
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param stran array containing total strain component (dimension ntens) at the beginning of increment
 * @param dstran array containing the component of total strain increment (dimension ntens)
 * @param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
 * @param dtime time increment
 * @param temperature temperature value at the beginning of increment
 * @param Dtemperature temperature increment
 * @param nstatev number of evolution variables
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param nprops number of material properties
 * @param props array containing material properties
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param drot rotation increment matrix (dimension 3*3)
 * @param sigma arma::vec containing the components of the stress tensor
 * @param Lt arma::mat containing the mechanical tangent operator
 * @param Etot arma::vec containing total strain component at the beginning of increment
 * @param DEtot arma::vec containing the component of total strain increment
 * @param T temperature value at the beginning of increment
 * @param DT temperature increment
 * @param Time value of step time at the beginning of the current increment 
 * @param DTime Increment of time
 * @param props_smart arma::vec containing material properties
 * @param Wm variables related to mechanical work
 * @param statev_smart evolution variables in the form of an arma::vec
 * @param DR arma::mat increment of rigid body rotation
 * @param start bolean that states if it is the beginning of the simulation (or not)
 * @details Example: 
 * @code 
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
 * @endcode
*/
void abaqus2smart_M(const double *stress, const double *ddsdde, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops,const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, arma::vec &sigma, arma::mat &Lt, arma::vec &Etot, arma::vec &DEtot, double &T, double &DT, double &Time, double &DTime, arma::vec &props_smart, arma::vec &Wm, arma::vec &statev_smart, arma::mat &DR, bool &start);

/**
 * @brief From a abaqus umat_ function, get the simcoon variables for a thermo-mechanical simulation 
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param ddsddt array containing the mechanical/thermal tangent operator
 * @param drplde array containing the heat source/strain tangent operator
 * @param drpldt array containing the heat source/temperature tangent operator
 * @param stran array containing total strain component (dimension ntens) at the beginning of increment
 * @param dstran array containing the component of total strain increment (dimension ntens)
 * @param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
 * @param dtime time increment
 * @param temperature temperature value at the beginning of increment
 * @param Dtemperature temperature increment
 * @param nstatev number of evolution variables
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param nprops number of material properties
 * @param props array containing material properties
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param drot rotation increment matrix (dimension 3*3)
 * @param sigma arma::vec containing the components of the stress tensor
 * @param dSdE arma::mat containing the mechanical tangent operator
 * @param dSdT arma::mat containing the mechanical tangent operator
 * @param drpldE arma::mat containing the heat source/strain tangent operator
 * @param drpldT arma::mat containing the heat source/temperature tangent operator
 * @param Etot arma::vec containing total strain component at the beginning of increment
 * @param DEtot arma::vec containing the component of total strain increment
 * @param T temperature value at the beginning of increment
 * @param DT temperature increment
 * @param Time value of step time at the beginning of the current increment 
 * @param DTime Increment of time
 * @param props_smart arma::vec containing material properties
 * @param Wm variables related to mechanical work
 * @param Wt variables related to thermal work
 * @param statev_smart evolution variables in the form of an arma::vec
 * @param DR arma::mat increment of rigid body rotation
 * @param start bolean that states if it is the beginning of the simulation (or not)
 * @details Example: 
 * @code 
	bool start = false;
	double Time = 0.;
	double DTime = 0.;
    int solver_type = 0;
    
	int nstatev_smart = nstatev-7;
    vec props_smart = zeros(nprops);
    vec statev_smart = zeros(nstatev_smart);
    
    vec sigma = zeros(6);
    vec Etot = zeros(6);
    vec DEtot = zeros(6);    
    mat dSdE = zeros(6,6);
    mat dSdT = zeros(6,1);
    mat drpldE = zeros(6,1);
    mat drpldT = zeros(1,1);
    mat DR = zeros(3,3);
    vec Wm = zeros(4);
    vec Wt = zeros(3);
    double T = 0.;
    double DT = 0.;
    double r = 0.;
    
	abaqus2smart_T(stress, ddsdde, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, sigma, dSdE, dSdT, drpldE, drpldT, Etot, DEtot, T, DT, Time, DTime, props_smart, Wm, Wt, statev_smart, DR, start);
 * @endcode
*/
void abaqus2smart_T(const double *stress, const double *ddsdde, const double *ddsddt, const double *drplde, const double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops, const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, arma::vec &sigma, arma::mat &dSdE, arma::mat &dSdT, arma::mat &drpldE, arma::mat &drpldT, arma::vec &Etot, arma::vec &DEtot, double &T, double &DT, double &Time, double &DTime, arma::vec &props_smart, arma::vec &Wm, arma::vec &Wt, arma::vec &statev_smart, arma::mat &DR, bool &start);

/**
 * @brief From the name of the umat, select the appropriate function to determine the thermo-mechanical response
 * @param rve Reference to the phase characteristics.
 * @param DR arma::mat increment of rigid body rotation
 * @param Time value of step time at the beginning of the current increment 
 * @param DTime Increment of time
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param start bolean that states if it is the beginning of the simulation (or not)
 * @param solver_type type of solver to be used (0 = Newton-Raphson)
 * @param tnew_dt New increment of time if the max number of iteration has not converged
 * @details Example: 
 * @code 
    bool start = false;
    double Time = 0.;
    double DTime = 0.;
    int solver_type = 0;
    
    int nstatev_smart = nstatev-7;
    vec props_smart = zeros(nprops);
    vec statev_smart = zeros(nstatev_smart);
    
    vec sigma = zeros(6);
    vec Etot = zeros(6);
    vec DEtot = zeros(6);
    mat dSdE = zeros(6,6);
    mat dSdT = zeros(1,6);
    mat drpldE = zeros(6,1);
    mat drpldT = zeros(1,1);
    mat DR = zeros(3,3);
    vec Wm = zeros(4);
    vec Wt = zeros(3);
    
    string umat_name(cmname);
    umat_name = umat_name.substr(0, 5);
	   
    phase_characteristics rve;
    rve.construct(0,2);
    auto rve_sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
    rve_sv_T->resize(nstatev_smart);
    rve.sptr_matprops->resize(nprops);

    
    abaqus2smart_T(stress, ddsdde, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_T->sigma, rve_sv_T->dSdE, rve_sv_T->dSdT, rve_sv_T->drdE, rve_sv_T->drdT, rve_sv_T->Etot, rve_sv_T->DEtot, rve_sv_T->T, rve_sv_T->DT, Time, DTime, props_smart, rve_sv_T->Wm, rve_sv_T->Wt, rve_sv_T->statev, DR, start);
    
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
    select_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
 * @endcode
*/
void select_umat_T(phase_characteristics &rve, const arma::mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt);

/**
 * @brief From the name of the umat, select the appropriate function to determine the mechanical response considering non-linear kinematics
 * @param rve Reference to the phase characteristics.
 * @param DR arma::mat increment of rigid body rotation
 * @param Time value of step time at the beginning of the current increment 
 * @param DTime Increment of time
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param start bolean that states if it is the beginning of the simulation (or not)
 * @param solver_type type of solver to be used (0 = Newton-Raphson)
 * @param tnew_dt New increment of time if the max number of iteration has not converged
 * @details Example: 
 * @code 
    bool start = false;
    double Time = 0.;
    double DTime = 0.;
    int solver_type = 0;
    
    int nstatev_smart = nstatev-7;
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

    abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev_macro, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_macro, rve_sv_M->Wm, rve_sv_M->statev, DR, start);    
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
    select_umat_M_finite(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
 * @endcode
*/
void select_umat_M_finite(phase_characteristics &rve, const arma::mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt);    

/**
 * @brief From the name of the umat, select the appropriate function to determine the mechanical response considering small strain assumption
 * @param rve Reference to the phase characteristics.
 * @param DR arma::mat increment of rigid body rotation
 * @param Time value of step time at the beginning of the current increment 
 * @param DTime Increment of time
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param start bolean that states if it is the beginning of the simulation (or not)
 * @param solver_type type of solver to be used (0 = Newton-Raphson)
 * @param tnew_dt New increment of time if the max number of iteration has not converged
 * @details Example: 
 * @code 
    bool start = false;
    double Time = 0.;
    double DTime = 0.;
    int solver_type = 0;
    
    int nstatev_smart = nstatev-7;
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

    abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev_macro, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_macro, rve_sv_M->Wm, rve_sv_M->statev, DR, start);    
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
    select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
 * @endcode
*/    
void select_umat_M(phase_characteristics &rve, const arma::mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
void run_umat_T(phase_characteristics &rve, const arma::mat &DR,const double &Time,const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, const unsigned int &control_type, double &tnew_dt);

void run_umat_M(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, const unsigned int &, double &);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */    
    
/**
 * @brief Transfer variables from simcoon to Abaqus format, considering a mechanical constitutive law, returoning only updated variables
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param sigma arma::vec containing the components of the stress tensor
 * @param statev_smart evolution variables in the form of an arma::vec
 * @param Wm variables related to mechanical work
 * @param Lt arma::mat containing the mechanical tangent operator
 * @details Example: 
 * @code 
	bool start = false;
	double Time = 0.;
	double DTime = 0.;
    int solver_type = 0;
    
	vec sigma = zeros(6);
	vec Etot = zeros(6);
	vec DEtot = zeros(6);
	mat Lt = zeros(6,6);
	mat DR = zeros(3,3);
    vec Wm = zeros(4);
    
	string umat_name(cmname);
	umat_name = umat_name.substr(0, 5);
	   
    unsigned int nstatev_macro = 0;
    unsigned int nprops_macro = 0;
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props_macro = zeros(nprops);
    read_matprops(umat_name, nprops_macro, props_macro, nstatev_macro, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    
    unsigned int nstatev_multi = nstatev-nstatev_macro-4;
    vec statev_macro = zeros(nstatev_macro);
    vec statev_multi(&statev[nstatev_macro+4], nstatev_multi, false, false);
    
    phase_characteristics rve;
    rve.construct(0,1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_macro);
    rve.sptr_matprops->resize(nprops_macro);
    
    unsigned int pos=0;
    
    statev_2_phases(rve, pos, statev_multi);
	abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev_macro, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_macro, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops_macro, props_macro);
    start = true;
    select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
    phases_2_statev(statev_multi, pos, rve);
    
	smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
 * @endcode
*/	
void smart2abaqus_M(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const arma::vec &sigma, const arma::vec &statev_smart, const arma::vec &Wm, const arma::mat &Lt);

/**
 * @brief Transfer variables from simcoon to Abaqus format, considering a mechanical constitutive law, returning all Abaqus variables
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param stran array containing total strain component (dimension ntens) at the beginning of increment
 * @param dstran array containing the component of total strain increment (dimension ntens)
 * @param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
 * @param dtime time increment
 * @param temperature temperature value at the beginning of increment
 * @param Dtemperature temperature increment
 * @param nstatev number of evolution variables
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param nprops number of material properties
 * @param props array containing material properties
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param drot rotation increment matrix (dimension 3*3)
 * @param sigma arma::vec containing the components of the stress tensor
 * @param Lt arma::mat containing the mechanical tangent operator
 * @param Etot arma::vec containing total strain component at the beginning of increment
 * @param DEtot arma::vec containing the component of total strain increment
 * @param T temperature value at the beginning of increment
 * @param DT temperature increment
 * @param Time value of step time at the beginning of the current increment 
 * @param DTime Increment of time
 * @param props_smart arma::vec containing material properties
 * @param Wm variables related to mechanical work
 * @param statev_smart evolution variables in the form of an arma::vec
 * @param DR arma::mat increment of rigid body rotation
 * @param start bolean that states if it is the beginning of the simulation (or not)
 * @details Example: 
 * @code 
	bool start = false;
	double Time = 0.;
	double DTime = 0.;
    int solver_type = 0;
    
	vec sigma = zeros(6);
	vec Etot = zeros(6);
	vec DEtot = zeros(6);
	mat Lt = zeros(6,6);
	mat DR = zeros(3,3);
    vec Wm = zeros(4);
    
	string umat_name(cmname);
	umat_name = umat_name.substr(0, 5);
	   
    unsigned int nstatev_macro = 0;
    unsigned int nprops_macro = 0;
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props_macro = zeros(nprops);
    read_matprops(umat_name, nprops_macro, props_macro, nstatev_macro, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    
    unsigned int nstatev_multi = nstatev-nstatev_macro-4;
    vec statev_macro = zeros(nstatev_macro);
    vec statev_multi(&statev[nstatev_macro+4], nstatev_multi, false, false);
    
    phase_characteristics rve;
    rve.construct(0,1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_macro);
    rve.sptr_matprops->resize(nprops_macro);
    
    unsigned int pos=0;
    
    statev_2_phases(rve, pos, statev_multi);
	abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev_macro, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_macro, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops_macro, props_macro);
    start = true;
    select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
    phases_2_statev(statev_multi, pos, rve);
    
	smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
 * @endcode
*/ 
void smart2abaqus_M_full(double *stress, double *ddsdde, double *stran, double *dstran, double *time, double &dtime, double &temperature, double &Dtemperature, int &nprops, double *props,  int &nstatev, double *statev, const int &ndi, const int &nshr, double *drot, const arma::vec &sigma, const arma::mat &Lt, const arma::vec &Etot, const arma::vec &DEtot, const double &T, const double &DT, const double &Time, const double &DTime, const arma::vec &props_smart, const arma::vec &Wm, const arma::vec &statev_smart, const arma::mat &DR, bool &start);

/**
 * @brief Transfer variables from simcoon to Abaqus format, considering a mechanical constitutive law, returning all Abaqus variables
 * @param stress array containing the components of the stress tensor (dimension ntens)
 * @param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
 * @param ddsddt array containing the mechanical/thermal tangent operator
 * @param drplde array containing the heat source/strain tangent operator
 * @param drpldt array containing the heat source/temperature tangent operator
 * @param rpl Volumetric heat generation per unit time at the end of the increment caused by mechanical work
 * @param statev array containing the evolution variables (dimension nstatev)
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param sigma arma::vec containing the components of the stress tensor
 * @param statev_smart evolution variables in the form of an arma::vec
 * @param r Volumetric heat generation per unit time (total)
 * @param Wm variables related to mechanical work
 * @param Wt variables related to thermal work 
 * @param dSdE arma::mat containing the mechanical tangent operator
 * @param dSdT arma::mat containing the mechanical tangent operator
 * @param drpldE arma::mat containing the heat source/strain tangent operator
 * @param drpldT arma::mat containing the heat source/temperature tangent operator
 * @details Example: 
 * @code 
    bool start = false;
    double Time = 0.;
    double DTime = 0.;
    int solver_type = 0;
    
    int nstatev_smart = nstatev-7;
    vec props_smart = zeros(nprops);
    vec statev_smart = zeros(nstatev_smart);
    
    vec sigma = zeros(6);
    vec Etot = zeros(6);
    vec DEtot = zeros(6);
    mat dSdE = zeros(6,6);
    mat dSdT = zeros(1,6);
    mat drpldE = zeros(6,1);
    mat drpldT = zeros(1,1);
    mat DR = zeros(3,3);
    vec Wm = zeros(4);
    vec Wt = zeros(3);
    
    string umat_name(cmname);
    umat_name = umat_name.substr(0, 5);
	   
    phase_characteristics rve;
    rve.construct(0,2);
    auto rve_sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
    rve_sv_T->resize(nstatev_smart);
    rve.sptr_matprops->resize(nprops);

    
    abaqus2smart_T(stress, ddsdde, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_T->sigma, rve_sv_T->dSdE, rve_sv_T->dSdT, rve_sv_T->drdE, rve_sv_T->drdT, rve_sv_T->Etot, rve_sv_T->DEtot, rve_sv_T->T, rve_sv_T->DT, Time, DTime, props_smart, rve_sv_T->Wm, rve_sv_T->Wt, rve_sv_T->statev, DR, start);
    
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
    select_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);

    smart2abaqus_T(stress, ddsdde, ddsddt, drplde, drpldt, rpl, statev, ndi, nshr, rve_sv_T->sigma, rve_sv_T->statev, rve_sv_T->r, rve_sv_T->Wm, rve_sv_T->Wt, rve_sv_T->dSdE, rve_sv_T->dSdT, rve_sv_T->drdE, rve_sv_T->drdT);
 * @endcode
*/ 
void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, double *drplde, double &drpldt, double &rpl, double *statev, const int &ndi, const int &nshr, const arma::vec &sigma, const arma::vec &statev_smart, const double &r, const arma::vec &Wm, const arma::vec &Wt, const arma::mat &dSdE, const arma::mat &dSdT, const arma::mat &drpldE, const arma::mat &drpldT);
            
} //namespace simcoon
