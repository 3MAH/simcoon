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

///@file crystallo.cpp
///@brief Some definitions coming from crystallography
///@version 1.0

#include <iostream>
#include <math.h>
#include <string.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Material/crystallo.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

mat Schmid(const vec &n, const vec &m){
	
	return 0.5*(n*trans(m) + m*trans(n))*(1./(norm(n,2)*norm(m,2)));
}

vec Schmid_v(const vec &n, const vec &m){
		
	return t2v_strain(Schmid(n, m));
}

mat F_nm(const vec &N){

	vec b = zeros(6);	
	b = Ith() - t2v_stress(N*trans(N));
	
	return p_ikjl(b);
}

mat Q_nm(const vec &N, const double &mu, const double &lambda) {

	mat F = zeros(6,6);
	mat Q = zeros(6,6);	
	vec b = zeros(6);
	b = Ith() - t2v_stress(N*trans(N));
	
	F = p_ikjl(b);
	Q = 2.*mu*(F + (lambda/(lambda+2.*mu))*((b)*trans(b)));
	
	return Q;
}

} //namespace simcoon
