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
#include <FTensor.hpp>
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
    return det(S)*(inv(S)).t();
}

//This function returns the derivative of the inverse of a symmetric tensor
/*mat dinvSdS(const mat &S) {

    Tensor2<double,3,3> invS_ = mat_FTensor2(inv(S));
    Tensor4<double,3,3,3,3> I_;

    Index<'I', 3> I;
    Index<'J', 3> J;
    Index<'K', 3> K;
    Index<'L', 3> L;

    I_(I,J,K,L) = invS_(I,K)*invS_(J,L) + invS_(I,L)*invS_(J,K);

    return 0.5*FTensor4_mat(I_);
}*/

mat dinvSdS(const mat &S) {

	mat invS = inv(S);
	mat F = zeros(6,6);

	int ij=0;
	int kl=0;
	int ik=0;
	int jl=0;
	int il=0;
	int jk=0;
    
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
	
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			ij = Id(i,j);
			for (int k=0; k<3; k++) {
				for (int l=k; l<3; l++) {
					kl = Id(k,l);
					ik = Id(i,k);
					jl = Id(j,l);
					il = Id(i,l);
					jk = Id(j,k);
					
					F(ij,kl) += -0.5*(S(ik,jl)+S(il,jk));
				}
			}
		}
	}
	
	return F;
}

mat sym_dyadic(const mat &A, const mat &B) {

    vec A_v = t2v_sym(A);
    vec B_v = t2v_sym(B);
    mat C = zeros(6,6);

    for (unsigned int i=0; i<6; i++) {
        for (unsigned int j=0; j<6; i++) {
            C(i,j) = A_v(i)*A_v(j);
        }
    }
    return C;
}
    
} //namespace simcoon
