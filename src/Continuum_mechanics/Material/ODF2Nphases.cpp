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
///@brief Densities of an ODF at given angles (its discretisation into phases is done in Python)
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF2Nphases.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
vec get_densities_ODF(const vec &x, const std::vector<peak> &peaks, const bool &radian) {
    
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
    odf_rve.peaks = peaks;
    if (!radian) {
        for (auto &p : odf_rve.peaks) {
            p.mean = simcoon::deg2rad(p.mean);
            p.s_dev = simcoon::deg2rad(p.s_dev);
            p.width = simcoon::deg2rad(p.width);
        }
    }
    
    for(unsigned int i=0; i<x.n_elem; i++) {
        if (radian)
            y(i) = odf_rve.density(x(i));
        else
            y(i) = odf_rve.density(x_rad(i));
    }
    return y;
}
    
} //namespace simcoon
