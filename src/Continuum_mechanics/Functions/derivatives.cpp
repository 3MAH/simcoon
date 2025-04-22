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
#include <simcoon/FTensor.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>

using namespace std;
using namespace arma;
using namespace FTensor;

namespace simcoon{

//This function returns the derivative of the first invariant (trace) of a tensor
mat dI1DS(const mat &S) {
    UNUSED(S);
    return eye(3,3);
}

//This function returns the derivative of the second invariant of a tensor : I_2 = 1/2 S_ij S_ij
mat dI2DS(const mat &S) {
    return S;
}

//This function returns the derivative of the third invariant of a tensor : I_3 = 1/3 S_ij S_jk S_ki
mat dI3DS(const mat &S) {
    return (S*S).t();
}

//This function returns the derivative of the trace of a tensor
mat dtrSdS(const mat &S) {
    UNUSED(S);
    return eye(3,3);
}

//This function returns the derivative of the determinant of a tensor
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
        throw simcoon::exception_inv("Error in eig_sym function inside dinvSdSsym.");
    }    

    Tensor2<double,3,3> invS_ = mat_FTensor2(invS);
    Tensor4<double,3,3,3,3> dinvSdSsym_;
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
        
    dinvSdSsym_(i,j,k,l) = invS_(i,k)*invS_(j,l)+invS_(i,l)*invS_(j,k);
    return 0.5*FTensor4_mat(dinvSdSsym_);
}
    
} //namespace simcoon
