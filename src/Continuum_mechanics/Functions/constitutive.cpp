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

///@file constitutive.cpp
///@brief Constitutive tensors in Voigt notation
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <string>
#include <stdexcept>
#include <armadillo>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

mat Ireal(){
	
	mat Ireal = eye(6,6);
	for (int i=3; i<6; i++)
		Ireal(i,i) = 0.5;
		
	return Ireal;
}

mat Ivol(){
	
	mat Ivol = zeros(6,6);
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			Ivol(i,j) = (1./3.);
		}
	}
	return Ivol;
}

// a*Ivol + b*X entry by entry: one allocation, same rounding as the expression
static mat iso_combination(const double &a, const double &b, const double &X_shear) {
	const double v = a*(1./3.);             // a Ivol_ij on the normal block
	const double d_ii = b*(1.-1./3.);       // b X_ii, normal block
	const double d_ij = b*(0.-1./3.);       // b X_ij, normal block, i != j
	const double d_ss = b*X_shear;          // b X_ii, shear block
	const double n_ii = v + d_ii;
	const double n_ij = v + d_ij;
	mat L(6, 6, fill::zeros);
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			L(i,j) = (i == j) ? n_ii : n_ij;
		}
		L(i+3,i+3) = d_ss;
	}
	return L;
}

mat Idev(){
	
	return	iso_combination(0., 1., 0.5);
}

mat Ireal2(){
	
	mat Ireal2 = eye(6,6);
	for (int i=3; i<6; i++)
		Ireal2(i,i) = 2.;
		
	return Ireal2;
}

mat Idev2(){
	
	return	iso_combination(0., 1., 2.);
}

vec Ith(){
	
	vec Ith = zeros(6);	
	for (int i=0; i<3; i++)
		Ith(i) = 1.;
		
	return Ith;
}

vec Ir2(){
	
	vec Ir2 = ones(6);	
	for (int i=3; i<6; i++)
		Ir2(i) = 2.;
	
	return Ir2;
}

vec Ir05(){
	
	vec Ir05 = ones(6);	
	for (int i=3; i<6; i++)
		Ir05(i) = 0.5;
	
	return Ir05;
}

mat L_iso(const double &C1, const double &C2, const std::string &conv) {
	
	double K = 0.;
	double mu = 0.;	
	
	if(conv == "Enu") {
		assert(C1 > 0.);
		assert(C2 < 0.5);
		K = C1/(3.*(1.-2*C2));
		mu = C1/(2.*(1.+C2));			
	}
	else if(conv == "nuE") {	
		assert(C2 > 0.);
		assert(C1 < 0.5);
		K = C2/(3.*(1.-2*C1));
		mu = C2/(2.*(1.+C1));
	}
	else if((conv == "Kmu")||(conv == "KG")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C1;
		mu = C2;
	}
	else if((conv == "muK")||(conv == "GK")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C2;
		mu = C1;
	}
	else if((conv == "lambdamu")||(conv == "lambdaG")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C1+(2./3.)*C2;
		mu = C2;
	}
	else if((conv == "mulambda")||(conv == "Glambda")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C2+(2./3.)*C1;
		mu = C1;
	}	
	else {
		throw std::invalid_argument("L_iso: unknown convention '" + conv + "'; expected Enu, nuE, Kmu/KG, muK/GK, lambdamu/lambdaG, mulambda/Glambda");
	}
	
	return iso_combination(3.*K, 2.*mu, 0.5);   // 3K Ivol + 2mu Idev
}

mat M_iso(const double &C1, const double &C2, const string &conv) {
	
	double K = 0.;
	double mu = 0.;	
	
	if(conv == "Enu") {
		assert(C1 > 0.);
		assert(C2 < 0.5);
		K = C1/(3.*(1.-2*C2));
		mu = C1/(2.*(1.+C2));		
	}
	else if(conv == "nuE") {	
		assert(C2 > 0.);
		assert(C1 < 0.5);
		K = C2/(3.*(1.-2*C1));
		mu = C2/(2.*(1.+C1));
	}
	else if((conv == "Kmu")||(conv == "KG")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C1;
		mu = C2;
	}
	else if((conv == "muK")||(conv == "GK")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C2;
		mu = C1;
	}
	else if((conv == "lambdamu")||(conv == "lambdaG")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C1+(2./3.)*C2;
		mu = C2;
	}
	else if((conv == "mulambda")||(conv == "Glambda")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C2+(2./3.)*C1;
		mu = C1;
	}	
	else {
		throw std::invalid_argument("M_iso: unknown convention '" + conv + "'; expected Enu, nuE, Kmu/KG, muK/GK, lambdamu/lambdaG, mulambda/Glambda");
	}
	
	return iso_combination(1/(3.*K), 1/(2.*mu), 2.);   // 1/(3K) Ivol + 1/(2mu) Idev2
}

mat L_cubic(const double &C1, const double &C2, const double &C4, const string &conv){
	
    mat L = zeros(6,6);
	if (conv == "Cii") {
        L(0,0) = C1;
        L(0,1) = C2;
        L(0,2) = C2;
        L(1,0) = C2;
        L(1,1) = C1;
        L(1,2) = C2;
        L(2,0) = C2;
        L(2,1) = C2;
        L(2,2) = C1;
        L(3,3) = C4;
        L(4,4) = C4;
        L(5,5) = C4;
    }
    else if(conv == "EnuG") {
        double E = C1;
        double nu = C2;
        double G = C4;
    
        double C11 = E*( 1-pow(nu,2))/(1 - 3*pow(nu,2) - 2*pow(nu,3) );
        double C12 = C11/(1./nu - 1.);
        double C44 = G;

        L(0,0) = C11;
        L(0,1) = C12;
        L(0,2) = C12;
        L(1,0) = C12;
        L(1,1) = C11;
        L(1,2) = C12;
        L(2,0) = C12;
        L(2,1) = C12;
        L(2,2) = C11;
        L(3,3) = C44;
        L(4,4) = C44;
        L(5,5) = C44;
    }
    else {
        throw std::invalid_argument("L_cubic: unknown convention '" + conv + "'; expected Cii, EnuG");
    }
    
	return L;
}

mat M_cubic(const double &C1, const double &C2, const double &C4, const string &conv){
	
    mat M = zeros(6,6);
	if (conv == "Cii") {
        double delta = (C1-C2)*(C1+2.*C2);
        M(0,0) = (C1+C2)/delta;
        M(0,1) = -C2/delta;
        M(0,2) = -C2/delta;
        M(1,0) = -C2/delta;
        M(1,1) = (C1+C2)/delta;
        M(1,2) = -C2/delta;
        M(2,0) = -C2/delta;
        M(2,1) = -C2/delta;
        M(2,2) = (C1+C2)/delta;
        M(3,3) = 1./C4;
        M(4,4) = 1./C4;
        M(5,5) = 1./C4;
    }
    else if (conv == "EnuG") {
        double E = C1;
        double nu = C2;
        double G = C4;
        
        M(0,0) = 1./E;
        M(0,1) = -nu/E;
        M(0,2) = -nu/E;
        M(1,0) = -nu/E;
        M(1,1) = 1./E;
        M(1,2) = -nu/E;
        M(2,0) = -nu/E;
        M(2,1) = -nu/E;
        M(2,2) = 1./E;
        M(3,3) = 1./G;
        M(4,4) = 1./G;
        M(5,5) = 1./G;
    }
    else {
        throw std::invalid_argument("M_cubic: unknown convention '" + conv + "'; expected Cii, EnuG");
    }
	
	return M;
}

bool cubic_coefficients(const mat &A, double &a, double &b, double &c, const double &tol) {
    if (A.n_rows != 6 || A.n_cols != 6) {
        return false;
    }
    a = A(0,0);
    b = A(0,1);
    c = A(3,3);
    const double thr = (tol > 0.) ? tol*abs(A).max() : 0.;
    for (int j=0; j<6; j++) {
        for (int i=0; i<6; i++) {
            const double expected = (i < 3 && j < 3) ? ((i == j) ? a : b) : ((i == j) ? c : 0.);
            if (!(std::fabs(A(i,j) - expected) <= thr)) {   // a NaN fails too
                return false;
            }
        }
    }
    return true;
}

bool is_isotropic(const mat &L, const Tensor4Type &type, const double &tol) {
    double a, b, c, row, col;
    mandel_factors(type, row, col);   // isotropic: Mandel shear c row col = a - b
    return cubic_coefficients(L, a, b, c, tol) && std::fabs(c*row*col - (a - b)) <= tol*abs(L).max();
}

// closed form: eigenvalues a+2b, a-b of the normal block, shear 1/c
mat inv_cubic(const mat &A) {

    double a = 0., b = 0., c = 0.;
    if (!cubic_coefficients(A, a, b, c)) {
        mat A_inv;
        if (!inv(A_inv, A)) {
            throw simcoon::exception_inv("Error in inv function inside inv_cubic.");
        }
        return A_inv;
    }
    const double l_vol = a + 2.*b;
    const double l_dev = a - b;
    if (l_vol == 0. || l_dev == 0. || c == 0. || !std::isfinite(l_vol*l_dev*c)) {
        throw simcoon::exception_inv("Error in inv function inside inv_cubic (singular matrix).");
    }
    const double a_inv = (1./l_vol + 2./l_dev)/3.;
    const double b_inv = (1./l_vol - 1./l_dev)/3.;
    mat A_inv(6, 6, fill::zeros);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            A_inv(i,j) = (i == j) ? a_inv : b_inv;
        }
        A_inv(i+3,i+3) = 1./c;
    }
    return A_inv;
}

mat L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv){
	mat L = zeros(6,6);
	
	if (conv == "Cii") {
	    L(0,0) = C11;
	    L(0,1) = C12;
	    L(0,2) = C13;
	    L(1,0) = C12;
	    L(1,1) = C22;
	    L(1,2) = C23;
	    L(2,0) = C13;
	    L(2,1) = C23;
	    L(2,2) = C33;
	    L(3,3) = C44;
	    L(4,4) = C55;
	    L(5,5) = C66;
		
	}
	else if (conv == "EnuG") {
        double E1 = C11;
        double E2 = C12;
        double E3 = C13;
        double nu12 = C22;
        double nu13 = C23;
        double nu23 = C33;
        double G12 = C44;
        double G13 = C55;
        double G23 = C66;
        
        double nu21=nu12*(E2/E1);
        double nu32=nu23*(E3/E2);
        double nu31=nu13*(E3/E1);

        double invdelta =1./(1-nu12*nu21-nu23*nu32-nu31*nu13-2.*nu21*nu32*nu13);
         
        L(0,0) = E1*(1-nu23*nu32)*invdelta;
        L(0,1) = E1*(nu21+nu31*nu23)*invdelta;
        L(0,2) = E1*(nu31+nu21*nu32)*invdelta;
        L(1,0) = L(0,1);
        L(1,1) = E2*(1-nu13*nu31)*invdelta;
        L(1,2) = E2*(nu32+nu12*nu31)*invdelta;
        L(2,0) = L(0,2);
        L(2,1) = L(1,2);
        L(2,2) = E3*(1-nu12*nu21)*invdelta;
        L(3,3) = G12;
        L(4,4) = G13;
        L(5,5) = G23;
    }
    else {
        throw std::invalid_argument("L_ortho: unknown convention '" + conv + "'; expected Cii, EnuG");
    }
	
	return L;
}

mat M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv){
	
    mat M = zeros(6,6);
    mat L = zeros(6,6);
    
    if (conv == "Cii") {
	    L(0,0) = C11;
	    L(0,1) = C12;
	    L(0,2) = C13;
	    L(1,0) = C12;
	    L(1,1) = C22;
	    L(1,2) = C23;
	    L(2,0) = C13;
	    L(2,1) = C23;
	    L(2,2) = C33;
	    L(3,3) = C44;
	    L(4,4) = C55;
	    L(5,5) = C66;
		
		try {
			M = inv(L);
		} catch (const std::runtime_error &e) {
			cerr << "Error in inv: " << e.what() << endl;
			throw simcoon::exception_inv("Error in inv function inside M_ortho.");
		}     
	}
    
    else if (conv == "EnuG") {
		double E1 = C11;
		double E2 = C12;
		double E3 = C13;
		double nu12 = C22;
		double nu13 = C23;
		double nu23 = C33;
		double G12 = C44;
		double G13 = C55;
		double G23 = C66;
        
		M(0,0) = 1./E1;
		M(0,1) = -nu12/E1;
		M(0,2) = -nu13/E1;
		M(1,0) = -nu12/E1;
		M(1,1) = 1./E2;
		M(1,2) = -nu23/E2;
		M(2,0) = -nu13/E1;
		M(2,1) = -nu23/E2;
		M(2,2) = 1./E3;
		M(3,3) = 1./G12;
		M(4,4) = 1./G13;
		M(5,5) = 1./G23;
	}
    else {
        throw std::invalid_argument("M_ortho: unknown convention '" + conv + "'; expected Cii, EnuG");
    }
	
	return M;
}

mat L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis){
	
    mat L = zeros(6,6);
    double delta = (1.+nuTT)*(2*EL*nuTL*nuTL+ET*(nuTT-1.))/ET;

	switch(axis) {
		case 1: {	
		    L(0,0) = EL*(nuTT*nuTT-1.)/delta;
		    L(0,1) = -EL*nuTL*(1+nuTT)/delta;
		    L(0,2) = -EL*nuTL*(1+nuTT)/delta;
		    L(1,0) = -EL*nuTL*(1+nuTT)/delta;
		    L(1,1) = (EL*nuTL*nuTL-ET)/delta;
		    L(1,2) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(2,0) = -EL*nuTL*(1+nuTT)/delta;
		    L(2,1) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(2,2) = (EL*nuTL*nuTL-ET)/delta;
		    L(3,3) = GLT;
		    L(4,4) = GLT;
		    L(5,5) = ET/(2.*(1.+nuTT));
		    break;
		}
		case 2: {
		    L(0,0) = (EL*nuTL*nuTL-ET)/delta;
		    L(0,1) = -EL*nuTL*(1+nuTT)/delta;
		    L(0,2) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(1,0) = -EL*nuTL*(1+nuTT)/delta;
		    L(1,1) = EL*(nuTT*nuTT-1.)/delta;
		    L(1,2) = -EL*nuTL*(1+nuTT)/delta;
		    L(2,0) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(2,1) = -EL*nuTL*(1+nuTT)/delta;
		    L(2,2) = (EL*nuTL*nuTL-ET)/delta;
		    L(3,3) = GLT;
		    L(4,4) = ET/(2.*(1.+nuTT));
		    L(5,5) = GLT;
			break;
		}
		case 3: {
			L(0,0) = (EL*nuTL*nuTL-ET)/delta;
			L(0,1) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
			L(0,2) = -EL*nuTL*(1+nuTT)/delta;
			L(1,0) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
			L(1,1) = (EL*nuTL*nuTL-ET)/delta;
			L(1,2) = -EL*nuTL*(1+nuTT)/delta;
			L(2,0) = -EL*nuTL*(1+nuTT)/delta;
			L(2,1) = -EL*nuTL*(1+nuTT)/delta;
			L(2,2) = EL*(nuTT*nuTT-1.)/delta;
			L(3,3) = ET/(2.*(1.+nuTT));
			L(4,4) = GLT;
			L(5,5) = GLT;
			break;	
		}
        default : {
            throw std::invalid_argument("L_isotrans: the axis must be 1, 2 or 3; got " + std::to_string(axis));
        }
            
	}

    
	return L;
}

mat M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis){
	
	mat M = zeros(6,6);
	switch(axis) {
		case 1: {	
			M(0,0) = 1./EL;
			M(0,1) = -nuTL/ET;
			M(0,2) = -nuTL/ET;
			M(1,0) = -nuTL/ET;
			M(1,1) = 1./ET;
			M(1,2) = -nuTT/ET;
			M(2,0) = -nuTL/ET;
			M(2,1) = -nuTT/ET;
			M(2,2) = 1./ET;
			M(3,3) = 1./GLT;
			M(4,4) = 1./GLT;
			M(5,5) = (2.*(1.+nuTT))/ET;
			break;
		}
		case 2: {	
			M(0,0) = 1./ET;
			M(0,1) = -nuTL/ET;
			M(0,2) = -nuTT/ET;
			M(1,0) = -nuTL/ET;
			M(1,1) = 1./EL;
			M(1,2) = -nuTL/ET;
			M(2,0) = -nuTT/ET;
			M(2,1) = -nuTL/ET;
			M(2,2) = 1./ET;
			M(3,3) = 1./GLT;
			M(4,4) = (2.*(1.+nuTT))/ET;
			M(5,5) = 1./GLT;
			break;
		}
		case 3: {	
			M(0,0) = 1./ET;
			M(0,1) = -nuTT/ET;
			M(0,2) = -nuTL/ET;
			M(1,0) = -nuTT/ET;
			M(1,1) = 1./ET;
			M(1,2) = -nuTL/ET;
			M(2,0) = -nuTL/ET;
			M(2,1) = -nuTL/ET;
			M(2,2) = 1./EL;
			M(3,3) = (2.*(1.+nuTT))/ET;
			M(4,4) = 1./GLT;
			M(5,5) = 1./GLT;
			break;
		}
        default : {
            throw std::invalid_argument("M_isotrans: the axis must be 1, 2 or 3; got " + std::to_string(axis));
        }
	}			
	return M;
}

mat H_iso(const double &etaB, const double &etaS) {
		
	return iso_combination(3.*etaB, 2.*etaS, 0.5);   // 3 etaB Ivol + 2 etaS Idev
}

// ndi < 3 (Abaqus order): first ndi direct + first nshr shears kept, the rest stress-free;
// Q = L_KK - L_KC L_CC^-1 L_CK. nshr < 0: Abaqus default for ndi
static vec condensed_pred(const mat &L, const vec &E, const int &ndi, const int &nshr) {

    const int ns = (nshr < 0) ? ((ndi == 1) ? 0 : 1) : nshr;
    if ((ndi != 1 && ndi != 2) || ns > 3) {
        throw std::invalid_argument("el_pred: invalid condensed state (ndi = " + std::to_string(ndi) + ", nshr = " + std::to_string(nshr) + ")");
    }
    uvec K(ndi + ns), C(6 - ndi - ns);
    int nk = 0, nc = 0;
    for (int i = 0; i < 6; i++) {
        if ((i < 3) ? (i < ndi) : (i - 3 < ns)) K(nk++) = i;
        else C(nc++) = i;
    }
    const mat Q = L(K, K) - L(K, C) * solve(L(C, C), L(C, K));
    vec sigma = zeros(6);
    sigma(K) = Q * E(K);
    return sigma;
}

vec el_pred(const vec &sigma_start, const mat &L, const vec &DE_el, const int &ndi, const int &nshr) {

    if (ndi == 3)
        return sigma_start + L*DE_el;
    return sigma_start + condensed_pred(L, DE_el, ndi, nshr);
}

vec el_pred(const mat &L, const vec &E_el, const int &ndi, const int &nshr) {

    if (ndi == 3)
        return L*E_el;
    return condensed_pred(L, E_el, ndi, nshr);
}
    
mat Isotropize(const mat &Lt) {

    return accu(Ivol()%Lt)*Ivol() + (4./11.)*accu(Idev()%Lt)*Idev();

}
    
} //namespace simcoon
