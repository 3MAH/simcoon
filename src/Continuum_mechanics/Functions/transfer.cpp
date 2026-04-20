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

///@file objective_rate.cpp
///@brief A set of function that help to define different quantities, depending on a selected objective rate
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//This function transforms the strain Voigt vector into a 3*3 strain matrix
mat v2t_strain(const vec &v) {
    assert(v.size()==6);
    mat strain(3,3);
    
    for (int i=0; i<3; i++)
    { strain (i,i) = v(i);
        for (int j=i+1; j<3; j++) {
            strain(i,j) = 0.5 * v(i+j+2);
            strain(j,i) = 0.5 * v(i+j+2);
        }
    }
    
    return strain;
}

//This function transforms a 3*3 strain matrix into a strain Voigt vector
vec t2v_strain (const mat &strain) {
    assert((strain.n_cols==3)&&(strain.n_rows==3));
    vec v(6);
    
    for (int i=0; i<3; i++)
    { v(i) = strain (i,i);
        for (int j=i+1; j<3; j++)
            v(i+j+2) =  strain(i,j) + strain(j,i);
    }
    
    return v;
}

//This function transforms the stress Voigt vector into a 3*3 stress matrix
mat v2t_stress(const vec &v) {
    assert(v.size()==6);
    mat stress(3,3);
    
    for (int i=0; i<3; i++)
    { stress (i,i) = v(i);
        for (int j=i+1; j<3; j++) {
            stress(i,j) = v(i+j+2);
            stress(j,i) = v(i+j+2);
        }
    }
    
    return stress;
}

//This function transforms a 3*3 stress matrix into a stress Voigt vector
vec t2v_stress (const mat &stress) {
    assert((stress.n_cols==3)&&(stress.n_rows==3));	
    vec v(6);
    
    for (int i=0; i<3; i++)
    { v(i) = stress (i,i);
        for (int j=i+1; j<3; j++)
            v(i+j+2) =  0.5*(stress(i,j) + stress(j,i));
    }
    
    return v;
}

//This function transforms a 3x3 symmetric matrix into a vector (6 components 11,22,33,12,13,23)
vec t2v_sym (const mat &m) {
    return t2v_stress(m);
}

//This function transforms a vector (6 components 11,22,33,12,13,23) into a symmetric 3x3 stress matrix
mat v2t_sym (const vec &v) {
    return v2t_stress(v);
}

//This function transforms a vector (6 components 11,22,33,12,13,23) into a skew-symmetric 3x3 stress matrix
mat v2t_skewsym (const vec &v) {
    assert(v.size()==6);
    mat w(3,3);
    
    for (int i=0; i<3; i++)
    { w (i,i) = v(i);
        for (int j=i+1; j<3; j++) {
            w(i,j) = v(i+j+2);
            w(j,i) = -1.*v(i+j+2);
        }
    }
    
    return w;
}

mat v2t (const vec &v) {
    assert(v.n_elem == 9);
    mat t = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            t(i,j) = v(i*3+j);
        }
    }
    return t;
}

} //namespace simcoon
