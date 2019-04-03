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
#include <FTensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;
using namespace FTensor;

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

vec t2v_sym (const mat &m) {
    return t2v_stress(m);
}

//This function transforms an armadillo 3x3 matrix to a FTensor Tensor of the 2nd rank
Tensor2<double,3,3> mat_FTensor2(const mat &m) {
    
    Tensor2<double,3,3> T;
    for(unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            T(i,j) = m(i,j);
        }
    }
    return T;
}
    
//This function transforms a strain Voigt vector to a FTensor Tensor of the 2nd rank
Tensor2<double,3,3> v_FTensor2_strain(const vec &v) {
    
    assert(v.size()==6);
    Tensor2<double,3,3> T;
    
    for (unsigned int i=0; i<3; i++) {
        T (i,i) = v(i);
        for (unsigned int j=i+1; j<3; j++) {
            T(i,j) = 0.5 * v(i+j+2);
            T(j,i) = 0.5 * v(i+j+2);
        }
    }	
    
    return T;
}

//This function transforms a stress Voigt vector to a FTensor Tensor of the 2nd rank
Tensor2<double,3,3> v_FTensor2_stress(const vec &v) {
    
    assert(v.size()==6);
    Tensor2<double,3,3> T;
    
    for (unsigned int i=0; i<3; i++) {
        T(i,i) = v(i);
        for (unsigned int j=i+1; j<3; j++) {
            T(i,j) = v(i+j+2);
            T(j,i) = v(i+j+2);
        }
    }	
    
    return T;
}
   
//This function transforms an armadillo 3x3 matrix to a FTensor Tensor of the 2nd rank
mat FTensor_mat(const Tensor2<double,3,3> &T) {
    
    mat m(3,3);
    for(unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            m(i,j) = T(i,j);
        }
    }
    return m;
}
    
//This function transforms a FTensor Tensor of the 2nd rank to a strain Voigt vector
vec FTensor2_v_strain(const Tensor2<double,3,3> &T) {
    
    vec v = zeros(6);
    
    for (unsigned int i=0; i<3; i++) {
        v(i) = T(i,i);
        for (unsigned int j=i+1; j<3; j++) {
            v(i+j+2) = T(i,j) + T(j,i);
        }
    }
    return v;
}

//This function transforms a FTensor Tensor of the 2nd rank to a stress Voigt vector
vec FTensor2_v_stress(const Tensor2<double,3,3> &T) {
    
    vec v = zeros(6);
    
    for (int i=0; i<3; i++) {
        v(i) = T(i,i);
        for (int j=i+1; j<3; j++) {
            v(i+j+2)=0.5*(T(i,j) + T(j,i));
        }
    }
    return v;
}

//This function transforms a FTensor Tensor of the 4nd rank to a stiffness symmetric matrix 6x6
mat FTensor4_mat(const Tensor4<double,3,3,3,3> &C) {
    
    mat L = zeros(6,6);
    L(0,0) = C(0,0,0,0); //C_1111
    L(0,1) = C(0,0,1,1); //C_1122
    L(0,2) = C(0,0,2,2); //C_1133
    L(0,3) = 0.5*C(0,0,0,1)+0.5*C(0,0,1,0); //.5*C_1112+.5*C_1121
    L(0,4) = 0.5*C(0,0,0,2)+0.5*C(0,0,2,0); //.5*C_1113+.5*C_1131
    L(0,5) = 0.5*C(0,0,1,2)+0.5*C(0,0,2,1); //.5*C_1123+.5*C_1123

    L(1,0) = C(1,1,0,0); //C_2211
    L(1,1) = C(1,1,1,1); //C_2222
    L(1,2) = C(1,1,2,2); //C_2233
    L(1,3) = 0.5*C(1,1,0,1)+0.5*C(1,1,1,0); //.5*C_2212+.5*C_2221
    L(1,4) = 0.5*C(1,1,0,2)+0.5*C(1,1,2,0); //.5*C_2213+.5*C_2231
    L(1,5) = 0.5*C(1,1,1,2)+0.5*C(1,1,2,1); //.5*C_2223+.5*C_2232

    L(2,0) = C(2,2,0,0); //C_3311
    L(2,1) = C(2,2,1,1); //C_3322
    L(2,2) = C(2,2,2,2); //C_3333
    L(2,3) = 0.5*C(2,2,0,1)+0.5*C(2,2,1,0); //.5*C_3312+.5*C_3321
    L(2,4) = 0.5*C(2,2,0,2)+0.5*C(2,2,2,0); //.5*C_3313+.5*C_3331
    L(2,5) = 0.5*C(2,2,1,2)+0.5*C(2,2,2,1); //.5*C_3323+.5*C_3332
    
    L(3,0) = 0.5*C(0,1,0,0)+0.5*C(1,0,0,0); //.5*C_1211+.5*C_2111
    L(3,1) = 0.5*C(0,1,1,1)+0.5*C(1,0,1,1); //.5*C_1222+.5*C_2122
    L(3,2) = 0.5*C(0,1,2,2)+0.5*C(1,0,2,2); //.5*C_1233+.5*C_2133
    L(3,3) = 0.5*C(0,1,0,1)+0.5*C(1,0,1,0); //.5*C_1212+.5*C_2112
    L(3,4) = 0.5*C(0,1,0,2)+0.5*C(1,0,2,0); //.5*C_1213+.5*C_2131
    L(3,5) = 0.5*C(0,1,1,2)+0.5*C(1,0,2,1); //.5*C_1223+.5*C_2132

    L(4,0) = 0.5*C(0,2,0,0)+0.5*C(2,0,0,0); //.5*C_1311+.5*C_3111
    L(4,1) = 0.5*C(0,2,1,1)+0.5*C(2,0,1,1); //.5*C_1322+.5*C_3122
    L(4,2) = 0.5*C(0,2,2,2)+0.5*C(2,0,2,2); //.5*C_1333+.5*C_3133
    L(4,3) = 0.5*C(0,2,0,1)+0.5*C(2,0,1,0); //.5*C_1312+.5*C_3112
    L(4,4) = 0.5*C(0,2,0,2)+0.5*C(2,0,2,0); //.5*C_1313+.5*C_3131
    L(4,5) = 0.5*C(0,2,1,2)+0.5*C(2,0,2,1); //.5*C_1323+.5*C_3132

    L(5,0) = 0.5*C(1,2,0,0)+0.5*C(2,1,0,0); //.5*C_2311+.5*C_3211
    L(5,1) = 0.5*C(1,2,1,1)+0.5*C(2,1,1,1); //.5*C_2322+.5*C_3222
    L(5,2) = 0.5*C(1,2,2,2)+0.5*C(2,1,2,2); //.5*C_2333+.5*C_3233
    L(5,3) = 0.5*C(1,2,0,1)+0.5*C(2,1,1,0); //.5*C_2312+.5*C_3212
    L(5,4) = 0.5*C(1,2,0,2)+0.5*C(2,1,2,0); //.5*C_2313+.5*C_3231
    L(5,5) = 0.5*C(1,2,1,2)+0.5*C(2,1,2,1); //.5*C_2323+.5*C_3232
    return L;
}

//This function transforms a stiffness symmetric matrix 6x6 to a FTensor Tensor of the 4nd rank    
Tensor4<double,3,3,3,3> mat_FTensor4(const mat &L) {
    Tensor4<double,3,3,3,3> C;
    
    int ij=0;
    int kl=0;
    
    umat Id(3,3);
    Id(0,0) = 0;
    Id(0,1) = 3;
    Id(0,2) = 4;
    Id(1,0) = 3;
    Id(1,1) = 1;
    Id(1,2) = 5;
    Id(2,0) = 4;
    Id(2,1) = 5;
    Id(2,2) = 2;
    
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=i; j<3; j++) {
            ij = Id(i,j);
            for (unsigned int k=0; k<3; k++) {
                for (unsigned int l=k; l<3; l++) {
                    kl = Id(k,l);
                    C(i,j,k,l) = L(ij,kl);
                    C(j,i,l,k) = C(i,j,k,l);
                }
            }
        }
    }
    return C;
}
    
    
} //namespace simcoon
