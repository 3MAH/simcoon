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
#include <cmath>
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

mean_field_options mean_field_options_of(const phase_characteristics &phase) {
    const std::string &name = phase.sptr_matprops->umat_name;
    const vec &props = phase.sptr_matprops->props;
    mean_field_options opt;
    //first option slot: after n_matrix (MIMTN), after start (MISCN)
    const uword first = (name == "MIMTN") ? 3 : ((name == "MISCN") ? 4 : props.n_elem);
    if (props.n_elem > first) opt.closed_form = static_cast<int>(props(first));
    if (props.n_elem > first + 1) opt.rule = static_cast<int>(props(first + 1));
    if (name == "MISCN") {
        if (props.n_elem > first + 2) opt.anderson = static_cast<int>(props(first + 2));
        if (props.n_elem > first + 3) opt.rtol = props(first + 3);
    }
    return opt;
}

void check_sub_phases(const phase_characteristics &phase) {
    const std::string &name = phase.sptr_matprops->umat_name;
    const vec &props = phase.sptr_matprops->props;
    const int nphases = static_cast<int>(phase.sub_phases.size());
    if (nphases == 0) {
        throw std::invalid_argument(name + " needs its sub-phases: pass them (phases=) to the solver or to L_eff.");
    }

    //The props are the scheme's settings only, [mp, np, n_matrix, ...] then the options of
    //mean_field_options. The length and the option values turn the pre-2.0 layout
    //[nphases, file number, mp, np, n_matrix] into an error instead of a misread n_matrix.
    const bool self_consistent = (name == "MISCN");
    const uword n_min = (name == "MIPLN") ? 0 : ((name == "MIHEN") ? 2 : 3);
    const uword n_max = self_consistent ? 8 : ((name == "MIMTN") ? 5 : n_min);
    if (props.n_elem < n_min || props.n_elem > n_max) {
        throw std::invalid_argument(name + " takes props = " + std::string(name == "MIPLN" ? "[]" : (name == "MIHEN" ? "[mp, np]"
                                    : (self_consistent ? "[mp, np, n_matrix, (start, closed_form, rule, anderson, rtol)]" : "[mp, np, n_matrix, (closed_form, rule)]")))
                                    + " (trailing slots optional), got " + std::to_string(props.n_elem) + " values (the leading [nphases, file number] slots no longer exist)");
    }
    //options: comparisons written so that NaN fails them
    const uword first_option = self_consistent ? 4 : 3;
    for (uword k = first_option; k < props.n_elem; k++) {
        const uword slot = k - first_option;
        const double v = props(k);
        if (slot <= 1 && !(v == 0. || v == 1.)) {
            throw std::invalid_argument(name + ": " + (slot == 0 ? "closed_form" : "rule") + " (props[" + std::to_string(k)
                                        + "]) must be 0 or 1, got " + std::to_string(v)
                                        + " (the leading [nphases, file number] slots of the pre-2.0 layout no longer exist)");
        }
        if (slot == 2 && !(v >= 0. && v <= 10. && v == std::floor(v))) {
            throw std::invalid_argument(name + ": anderson (props[" + std::to_string(k) + "]), the depth of the Anderson "
                                        "acceleration, must be an integer in [0, 10], got " + std::to_string(v));
        }
        if (slot == 3 && !(v >= 0. && v < 1.)) {
            throw std::invalid_argument(name + ": rtol (props[" + std::to_string(k) + "]) must be in [0, 1), got " + std::to_string(v));
        }
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
    
    //The self-consistent scheme takes the effective tangent of the previous increment as its
    //reference medium. The solver clears that tangent (Lt = 0) before the priming call of each
    //block: a cleared tangent is seeded again, as at start, instead of being inverted.
    const bool sc_seed = start || umat_phase_M->Lt.is_zero();

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
                //Isotropized Mori-Tanaka (Doghri & Ouaar): no UMAT name maps to code 102 yet
                //(umat_smart), so this scheme is not reachable from the solver nor from L_eff
                int n_matrix = phase.sptr_matprops->props(2);
                DE_Mori_Tanaka_iso(phase, n_matrix);
                break;
            }
            case 103: {
                int n_matrix = phase.sptr_matprops->props(2);
                DE_Self_Consistent(phase, n_matrix, sc_seed, self_consistent_start(phase));
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
                Lt_Self_Consistent(phase, n_matrix, sc_seed, self_consistent_start(phase));
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
