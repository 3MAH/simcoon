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

///@file multiphase.cpp
///@brief User subroutine for non-linear N-phases heterogeneous materials
///@version 1.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <stdexcept>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/schemes.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    

int sub_phase_shape(const std::string &umat_name) {
    if (umat_name == "MIHEN" || umat_name == "MIMTN" || umat_name == "MISCN") {
        return 2;
    }
    if (umat_name == "MIPLN") {
        return 1;
    }
    return 0;
}

int self_consistent_start(const phase_characteristics &phase) {
    const vec &props = phase.sptr_matprops->props;
    return (props.n_elem > 3) ? static_cast<int>(props(3)) : 1;
}

void check_sub_phases(const phase_characteristics &phase) {
    const std::string &name = phase.sptr_matprops->umat_name;
    const vec &props = phase.sptr_matprops->props;
    const int nphases = static_cast<int>(phase.sub_phases.size());
    if (nphases == 0) {
        throw std::invalid_argument(name + " needs its sub-phases: pass them (phases=) to the solver or to L_eff.");
    }

    //The props are the scheme's settings only. An exact length turns the pre-2.0 layout
    //[nphases, file number, mp, np, n_matrix] into an error instead of a misread n_matrix.
    const bool self_consistent = (name == "MISCN");
    const uword n_min = (name == "MIPLN") ? 0 : ((name == "MIHEN") ? 2 : 3);
    const uword n_max = self_consistent ? 4 : n_min;
    if (props.n_elem < n_min || props.n_elem > n_max) {
        throw std::invalid_argument(name + " takes props = " + std::string(name == "MIPLN" ? "[]" : (name == "MIHEN" ? "[mp, np]"
                                    : (self_consistent ? "[mp, np, n_matrix] or [mp, np, n_matrix, start]" : "[mp, np, n_matrix]")))
                                    + ", got " + std::to_string(props.n_elem) + " values (the leading [nphases, file number] slots no longer exist)");
    }
    //comparisons written so that NaN fails them; the casts to int come after
    if (n_min >= 2 && !(props(0) >= 1. && props(0) <= 10000. && props(1) >= 1. && props(1) <= 10000.)) {
        throw std::invalid_argument(name + ": mp and np, the integration points of the Eshelby integrals, must be in [1, 10000]");
    }
    if (n_min >= 3) {
        if (!(props(2) >= -1. && props(2) < static_cast<double>(nphases))) {
            throw std::invalid_argument(name + ": n_matrix = " + std::to_string(props(2)) + " is not one of the "
                                        + std::to_string(nphases) + " phases given");
        }
        const int n_matrix = static_cast<int>(props(2));
        const int start = self_consistent ? self_consistent_start(phase) : 1;
        if (start != 0 && start != 1) {
            throw std::invalid_argument(name + ": start = " + std::to_string(start) + " is neither 0 (homogeneous strain) nor 1 (Mori-Tanaka)");
        }
        if (start == 1 && (n_matrix < 0 || n_matrix >= nphases)) {
            throw std::invalid_argument(name + ": n_matrix = " + std::to_string(n_matrix) + " is not one of the "
                                        + std::to_string(nphases) + " phases given");
        }
        if (start == 0 && n_matrix >= 0) {
            throw std::invalid_argument(name + ": the homogeneous-strain start (start = 0) takes n_matrix < 0");
        }
    }
}

void umat_multi(phase_characteristics &phase, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const unsigned int &solver_type, double &tnew_dt, const int &method)
{

    check_sub_phases(phase);
    const int nphases = static_cast<int>(phase.sub_phases.size());

    shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local); //shared_ptr on state variables of the rve
    shared_ptr<state_variables_M> umat_sub_phases_M; //shared_ptr on state variables
    
    //1 - Quadrature points of the Eshelby integrals (ellipsoidal schemes only)
    if (start && (sub_phase_shape(phase.sptr_matprops->umat_name) == 2)) {
        ellipsoid_multi::mp = phase.sptr_matprops->props(0);
        ellipsoid_multi::np = phase.sptr_matprops->props(1);
        ellipsoid_multi::x.set_size(ellipsoid_multi::mp);
        ellipsoid_multi::wx.set_size(ellipsoid_multi::mp);
        ellipsoid_multi::y.set_size(ellipsoid_multi::np);
        ellipsoid_multi::wy.set_size(ellipsoid_multi::np);
        points(ellipsoid_multi::x, ellipsoid_multi::wx, ellipsoid_multi::y, ellipsoid_multi::wy,ellipsoid_multi::mp, ellipsoid_multi::np);
    }
    
	//Initialization
	if (start) {

        for (int i=0; i<nphases; i++) {
            // The localization schemes (Hill interaction tensors, concentration
            // A, Lt_eff assembly) are built on the CONTINUUM phase tangents —
            // the incremental Mori-Tanaka / self-consistent / periodic-layer
            // formulations. Pin the sub-phase mode so the caller's tangent_mode
            // (algorithmic by default since 2.0) never leaks into the
            // localization. Persistent member: set once at start.
            phase.sub_phases[i].sptr_sv_global->tangent_mode = simcoon::tangent_continuum;
            //Run the appropriate constitutive model
            select_umat_M(phase.sub_phases[i], DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
        }
    }
    
    // Preliminaries of the convergence loop
    int nbiter = 0;
    double error = 1.;
    std::vector<vec> DE_N(nphases); //Table that stores all the previous increments of strain
    
	//Convergence loop, localization
	while ((error > simcoon::precision_micro)&&(nbiter <= simcoon::maxiter_micro)) {
	  
        for(int i=0; i<nphases; i++) {
            auto sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
            DE_N[i] = sv_r->DEtot;
        }

        //Compute the strain concentration tensor for each phase:
        //Also update of all the local strain increment
        switch (method) {
                
            case 100: {
                DE_Homogeneous_E(phase);
                break;
            }
            case 101: {
                int n_matrix = phase.sptr_matprops->props(2);
                DE_Mori_Tanaka(phase, n_matrix);
                break;
            }
            case 102: {
                int n_matrix = phase.sptr_matprops->props(2);
                DE_Mori_Tanaka_iso(phase, n_matrix);
                break;
            }
            case 103: {
                int n_matrix = phase.sptr_matprops->props(2);
                DE_Self_Consistent(phase, n_matrix, start, self_consistent_start(phase));
                break;
            }
            case 104: {
                dE_Periodic_Layer(phase, nbiter);
                break;
            }
        
        }
    
        for (unsigned int i=0; i<phase.sub_phases.size(); i++) {
            phase.sub_phases[i].sptr_sv_global->to_start();
            
            //Theta method for the tangent modulus
            //mat Lt_start = umat_sub_phases_M->Lt
            select_umat_M(phase.sub_phases[i], DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);

            //Theta method for the tangent modulus
            //umat_sub_phases_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
            //Lt* = (1 - (2./3.))*Lt_start + 2./3.*Lt;
        }
        
        error = 0.;
        for(int i=0; i<nphases; i++) {
            auto sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
            error += norm(DE_N[i] - sv_r->DEtot,2);
        }
        error*=(1./nphases);
        nbiter++;
	}
    
    //Now we can calculate the concentration tensors only for the tangent modulus
    switch (method) {
            
            case 100: {
                Lt_Homogeneous_E(phase);
                break;
            }
            case 101: {
                int n_matrix = phase.sptr_matprops->props(2);
                Lt_Mori_Tanaka(phase, n_matrix);
                break;
            }
            case 102: {
                int n_matrix = phase.sptr_matprops->props(2);
                Lt_Mori_Tanaka_iso(phase, n_matrix);
                break;
            }
            case 103: {
                int n_matrix = phase.sptr_matprops->props(2);
                Lt_Self_Consistent(phase, n_matrix, start, self_consistent_start(phase));
                break;
            }
            case 104: {
                Lt_Periodic_Layer(phase);
                break;
            }
            
    }
    
    //	Homogenization
	//Compute the effective stress
	umat_phase_M->sigma = zeros(6);
    for (auto r : phase.sub_phases) {
		umat_phase_M->sigma += r.sptr_shape->concentration*r.sptr_sv_global->sigma;
	}
    
    umat_phase_M->Lt = zeros(6,6);
	// Compute the effective tangent modulus, and the effective stress
    for (auto r : phase.sub_phases) {
        umat_sub_phases_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
		umat_phase_M->Lt += r.sptr_shape->concentration*(umat_sub_phases_M->Lt*r.sptr_multi->A);
	}
    
}

} //namespace simcoon
