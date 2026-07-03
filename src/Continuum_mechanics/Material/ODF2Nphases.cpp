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

///@file ODF2Nphases.cpp
///@brief ODF2Nphases discretization of ODFs
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Simulation/Maths/random.hpp>
#include <simcoon/Simulation/Maths/stats.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF2Nphases.hpp>
#include <simcoon/Continuum_mechanics/Material/read.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
vec get_densities_ODF(const vec &x, const string &path_data, const string &input_peaks, const bool &radian) {
    
    vec y = zeros(x.n_elem);
    vec x_rad;
    if (!radian) {
        x_rad = simcoon::deg2rad(x);
        if(x_rad.min() < 0.) {
            cout << "Error : x.min() < 0. Please provide an angle vector with all angles >=0.";
            return y;
        }
        if(x_rad.max() > simcoon::pi) {
            cout << "Error : x.max() > pi Please provide an angle vector with all angles <=pi/180";
            return y;
        }
    }
    else {
        if(x.min() < 0.) {
            cout << "Error : x.min() < 0. Please provide an angle vector with all angles >=0.";
            return y;
        }
        if(x.max() > simcoon::pi) {
            cout << "Error : x.max() > pi Please provide an angle vector with all angles <=pi";
            return y;
        }
    }
    
    ODF odf_rve(0, radian, x.min(), x.max());
    read_peak(odf_rve, path_data, input_peaks);
    
    for(unsigned int i=0; i<x.n_elem; i++) {
        if (radian)
            y(i) = odf_rve.density(x(i));
        else
            y(i) = odf_rve.density(x_rad(i));
    }
    return y;
}
    
void fill_angles(const double &alpha, phase_characteristics &phase, const ODF &odf_rve, const int &angles_mat) {
    
    double temp_psi_geom = 0.;
    double temp_theta_geom = 0.;
    double temp_phi_geom = 0.;
    
    switch (odf_rve.Angle) {
        case 0: {
            if (angles_mat)
                phase.sptr_matprops->psi_mat = alpha;
            temp_psi_geom = alpha;
            break;
        }
        case 1: {
            if (angles_mat)
                phase.sptr_matprops->theta_mat = alpha;
            temp_theta_geom = alpha;
            break;
        }
        case 2: {
            if (angles_mat)
                phase.sptr_matprops->phi_mat = alpha;
            temp_phi_geom = alpha;
            break;
        }
        default: {
            break;
        }
    }
    
    //Switch case for the geometry of the phase
    switch (phase.shape_type) {
        case 0: {
            break;
        }
        case 1: {
            std::shared_ptr<layer> lay = std::dynamic_pointer_cast<layer>(phase.sptr_shape);
            lay->psi_geom = temp_psi_geom;
            lay->theta_geom = temp_theta_geom;
            lay->phi_geom = temp_phi_geom;
            break;
        }
        case 2: {
            std::shared_ptr<ellipsoid> elli = std::dynamic_pointer_cast<ellipsoid>(phase.sptr_shape);
            elli->psi_geom = temp_psi_geom;
            elli->theta_geom = temp_theta_geom;
            elli->phi_geom = temp_phi_geom;
            break;
        }
        case 3: {
            std::shared_ptr<cylinder> cyl = std::dynamic_pointer_cast<cylinder>(phase.sptr_shape);
            cyl->psi_geom = temp_psi_geom;
            cyl->theta_geom = temp_theta_geom;
            cyl->phi_geom = temp_phi_geom;
            break;
        }
    }
}
    
phase_characteristics discretize_ODF(const phase_characteristics &rve_init, ODF &odf_rve, const int &num_phase_disc, const int &nb_phases_disc, const int &angles_mat) {
    
    phase_characteristics rve;
    rve.copy(rve_init);
    
    int number = 0;
    odf_rve.norm = 0.;
    
    double angle_range = odf_rve.limits(1) - odf_rve.limits(0);
    assert(angle_range > 0.);
    
    double dalpha = angle_range/double(nb_phases_disc);
    double alpha = odf_rve.limits(0);
    
    rve.sub_phases.erase(rve.sub_phases.begin()+num_phase_disc);
    for (int i=0; i<nb_phases_disc; i++) {
        
        phase_characteristics temp;
        temp.copy(rve_init.sub_phases[num_phase_disc]);
        
        fill_angles(alpha, temp, odf_rve, angles_mat);
        
        if(alpha - odf_rve.limits(0) < dalpha)
            temp.sptr_shape->concentration = dalpha/6. * (odf_rve.density(simcoon::pi+alpha-dalpha/2) + 4.*odf_rve.density(alpha) + odf_rve.density(alpha+dalpha/2));
        else
            temp.sptr_shape->concentration = dalpha/6. * (odf_rve.density(alpha-dalpha/2) + 4.*odf_rve.density(alpha) + odf_rve.density(alpha+dalpha/2));

        odf_rve.norm += temp.sptr_shape->concentration;
        
        rve.sub_phases.insert(rve.sub_phases.begin()+i+num_phase_disc, temp);
        alpha += dalpha;
    }
    
    ///Normalization
    for(int i=0; i<nb_phases_disc; i++) {
        rve.sub_phases[i+num_phase_disc].sptr_shape->concentration *= (rve_init.sub_phases[num_phase_disc].sptr_shape->concentration / odf_rve.norm);
    }
    
    for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
        rve.sub_phases[i].sptr_matprops->number = number;
        number++;
    }
    
//    cout << "rve = " << rve;
    
    return rve;
    
}

} //namespace simcoon
