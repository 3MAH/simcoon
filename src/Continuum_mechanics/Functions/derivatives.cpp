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

///@file derivatives.cpp
///@brief A set of functions that performs the derivative of specific tensors attributes (Invariants, determinant, inverse...)
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
using namespace std;
using namespace arma;

namespace simcoon{

mat dI1DS(const mat &S) {
    UNUSED(S);
    return eye(3,3);
}

mat dI2DS(const mat &S) {
    return S;
}

mat dI3DS(const mat &S) {
    return (S*S).t();
}

mat dJ2DS(const mat &S) {

    mat S_dev = dev(S);  
    return S_dev;
}

mat dJ3DS(const mat &S) {

    mat S_dev = dev(S);

    mat S2 = (S_dev * S_dev).t();
    double trS2 = trace(S2);
    return S2 - (trS2 / 3.0) * eye<mat>(3,3);
}


mat dtrSdS(const mat &S) {
    UNUSED(S);
    return eye(3,3);
}

mat ddetSdS(const mat &S) {

    try {
        return det(S)*(inv(S)).t();
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv (combined expression), error inv is thrown: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside ddetSdS.");
    }    
}

mat dinvSdSsym(const mat &S) {

    mat invS;

    try {
        invS = inv(S);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv : " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside dinvSdSsym.");
    }    

    // dinvSdS_ijkl = -0.5*(invS_ik*invS_jl + invS_il*invS_jk): the symmetric dyadic operator of invS, negated
    return -auto_sym_dyadic_operator(mat::fixed<3,3>(invS));
}
    
} //namespace simcoon
