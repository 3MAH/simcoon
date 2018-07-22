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
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/schemes.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
///@brief The elastic isotropic UMAT requires 4 constants, plus the number of constants for each materials:
///@brief props[0] : Number of phases
///@brief props[1] : File # that stores the microstructure properties
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

///@brief The table Nphases.dat will store the necessary informations about the geometry of the phases and the material properties

void umat_multi(phase_characteristics &phase, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, const bool &start, const unsigned int &solver_type, double &tnew_dt, const int &method)
{

    int nphases = phase.sptr_matprops->props(0); // Number of phases
    string path_data = "data";
    string inputfile; //file # that stores the microstructure properties
    
    shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local); //shared_ptr on state variables of the rve
    shared_ptr<state_variables_M> umat_sub_phases_M; //shared_ptr on state variables
    
    //1 - We need to figure out the type of geometry and read the phase
    if(start) {
        switch (method) {
                
            case 100: case 101: case 102: case 103: {
                //Definition of the static vectors x,wx,y,wy
                ellipsoid_multi::mp = phase.sptr_matprops->props(2);
                ellipsoid_multi::np = phase.sptr_matprops->props(3);
                ellipsoid_multi::x.set_size(ellipsoid_multi::mp);
                ellipsoid_multi::wx.set_size(ellipsoid_multi::mp);
                ellipsoid_multi::y.set_size(ellipsoid_multi::np);
                ellipsoid_multi::wy.set_size(ellipsoid_multi::np);
                points(ellipsoid_multi::x, ellipsoid_multi::wx, ellipsoid_multi::y, ellipsoid_multi::wy,ellipsoid_multi::mp, ellipsoid_multi::np);
                
                inputfile = "Nellipsoids" + to_string(int(phase.sptr_matprops->props(1))) + ".dat";
                read_ellipsoid(phase, path_data, inputfile);
                break;
            }
            case 104: {
                inputfile = "Nlayers" + to_string(int(phase.sptr_matprops->props(1))) + ".dat";
                read_layer(phase, path_data, inputfile);
                break;
            }
        }
    }
    
	//Initialization
	if (start) {
        
        for (int i=0; i<nphases; i++) {
            //Run the appropriate constitutive model
            select_umat_M(phase.sub_phases[i], DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
        }
    }

    // Preliminaries of the convergence loop
    int nbiter = 0;
    double error = 1.;
    std::vector<vec> DE_N(nphases); //Table that stores all the previous increments of strain
    
	//Convergence loop, localization
	while ((error > precision_micro)&&(nbiter <= maxiter_micro)) {
	  
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
                int n_matrix = phase.sptr_matprops->props(4);
                DE_Mori_Tanaka(phase, n_matrix);
                break;
            }
            case 102: {
                int n_matrix = phase.sptr_matprops->props(4);
                DE_Mori_Tanaka_iso(phase, n_matrix);
                break;
            }
            case 103: {
                int n_matrix = phase.sptr_matprops->props(4);
                DE_Self_Consistent(phase, n_matrix, start, phase.sptr_matprops->props(5));
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
                int n_matrix = phase.sptr_matprops->props(4);
                Lt_Mori_Tanaka(phase, n_matrix);
                break;
            }
            case 102: {
                int n_matrix = phase.sptr_matprops->props(4);
                Lt_Mori_Tanaka_iso(phase, n_matrix);
                break;
            }
            case 103: {
                int n_matrix = phase.sptr_matprops->props(4);
                Lt_Self_Consistent(phase, n_matrix, start, phase.sptr_matprops->props(5));
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
