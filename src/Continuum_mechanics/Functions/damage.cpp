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

///@file damage.hpp
///@brief Functions that computes damage evolution laws
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/stats.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/damage.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//This function returns damage evolution (/dt) considering a Weibull damage law
double damage_weibull(const vec &stress, const double &damage, const double &alpha, const double &beta, const double &DTime, const string &criterion) {
		
	if ((criterion == "Mises")) {
		return ((1.-damage)/DTime)*(1. - exp( -1.* pow( Mises_stress(stress) / beta, alpha)));
	}
	else if ((criterion == "hydro")) {
		return ((1.-damage)/DTime)*(1. - exp( -1.* pow( tr(stress) /(3.* beta), alpha)));
	}
	else if ((criterion == "J3")) {
		return ((1.-damage)/DTime)*(1. - exp( -1.* pow( J3_stress(stress) / beta, alpha)));
	}
	else {
		cout << "ERROR : criterion name not found in damage.hpp";
		return 0.;
	}
}

//This function returns damage evolution (/dt) considering Kachanov's creep damage law
double damage_kachanov(const vec &stress, const vec &strain, const double &damage, const double &A0, const double &r, const string &criterion) {
		
	if ((criterion == "vonmises")) {
		return pow(Mises_stress(stress*(1.+strain))/(A0*(1.-damage)),r);
	}
	else if ((criterion == "hydro")) {
		return pow(tr(stress*(1.+strain))/(3.*A0*(1.-damage)),r);
	}
	else if ((criterion == "J3")) {
		return pow(J3_stress(stress*(1.+strain))/(A0*(1.-damage)),r);
	}
	else {
		cout << "ERROR : criterion name not found in damage.hpp";
		return 0.;
	}
}

//This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law
double damage_miner(const double &S_max, const double &S_mean, const double &S_ult, const double &b, const double &B0, const double &beta, const double &Sl_0) {
		
	return (S_max - (S_mean + Sl_0*(1.-b*S_mean)))/(S_ult-S_max)*pow((S_max-S_mean)/(B0*(1.-b*S_mean)),beta);
}

//This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law
double damage_manson(const double &S_amp, const double &C2, const double &gamma2) {
		
	return pow(S_amp/C2,gamma2);
}

} //namespace simcoon
