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
#include <simcoon/Continuum_Mechanics/Material/PDF.hpp>
#include <simcoon/Continuum_Mechanics/Material/PDF2Nphases.hpp>
#include <simcoon/Continuum_Mechanics/Material/read.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
vec get_densities_PDF(const vec &x, const string &path_data, const string &input_peaks) {
    
    vec y = zeros(x.n_elem);
    
    PDF pdf_rve(0, x.min(), x.max());
    read_peak(pdf_rve, path_data, input_peaks);
    
    for(unsigned int i=0; i<x.n_elem; i++) {
		y(i) = pdf_rve.density(x(i));
	}
    
    return y;
}
    
void fill_parameters(const double &alpha, phase_characteristics &phase, const PDF &pdf_rve) {
    
    phase.sptr_matprops->props(pdf_rve.Parameter) = alpha;
    
}
    
phase_characteristics discretize_PDF(const phase_characteristics &rve_init, PDF &pdf_rve, const int &num_phase_disc, const int &nb_phases_disc) {
    
    phase_characteristics rve;
    rve.copy(rve_init);
    
    int number = 0;
    pdf_rve.norm = 0.;
    
    double parameter_range = pdf_rve.limits(1) - pdf_rve.limits(0);
    assert(parameter_range > 0.);
    
    double dalpha = parameter_range/double(nb_phases_disc);
    double alpha = pdf_rve.limits(0);
    
    rve.sub_phases.erase(rve.sub_phases.begin()+num_phase_disc);
    for (int i=0; i<nb_phases_disc; i++) {
        
        phase_characteristics temp;
        temp.copy(rve_init.sub_phases[num_phase_disc]);
        
        fill_parameters(alpha, temp, pdf_rve);
        
        // if(alpha - pdf_rve.limits(0) < dalpha)
            // temp.sptr_shape->concentration = dalpha/6. * (pdf_rve.density(pi+alpha-dalpha/2) + 4.*pdf_rve.density(alpha) + pdf_rve.density(alpha+dalpha/2));
        // else
            // temp.sptr_shape->concentration = dalpha/6. * (pdf_rve.density(alpha-dalpha/2) + 4.*pdf_rve.density(alpha) + pdf_rve.density(alpha+dalpha/2));
            
        temp.sptr_shape->concentration = dalpha/6. * (pdf_rve.density(alpha-dalpha/2) + 4.*pdf_rve.density(alpha) + pdf_rve.density(alpha+dalpha/2));

        pdf_rve.norm += temp.sptr_shape->concentration;
        
        rve.sub_phases.insert(rve.sub_phases.begin()+i+num_phase_disc, temp);
        alpha += dalpha;
    }
    
    ///Normalization
    for(int i=0; i<nb_phases_disc; i++) {
        rve.sub_phases[i+num_phase_disc].sptr_shape->concentration *= (rve_init.sub_phases[num_phase_disc].sptr_shape->concentration / pdf_rve.norm);
    }
    
    for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
        rve.sub_phases[i].sptr_matprops->number = number;
        number++;
    }
    
//    cout << "rve = " << rve;
    
    return rve;
    
}

} //namespace simcoon
