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

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//This function returns the first Piola-Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Cauchy2PKI(const mat &sigma, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(3,3);
    }
    else {
        return J*sigma*inv(F.t());
    }
}

//This function returns the Biot stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Cauchy2Biot(const mat &sigma, const mat &F, const mat &mR, const double &mJ) {

    mat PKI = Cauchy2PKI(sigma, F, mJ);
    mat R = mR;

    if (norm(R,1) < sim_iota) {
        mat U = zeros(3,3);
        RU_decomposition(R,U,F);        
    }
    return 0.5*(R.t()*PKI + PKI.t()*R);
}

//This function returns the second Piola-Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Cauchy2PKII(const mat &sigma, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(3,3);
    }
    else {
        return J*inv(F)*sigma*inv(F.t());
    }
}

//This function returns the Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Cauchy2Kirchoff(const mat &sigma, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(3,3);
    }
    else {
        return J*sigma;
    }

}

//This function returns the Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
vec Cauchy2Kirchoff(const vec &sigma, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(6);
    }
    else {
        return J*sigma;
    }

}

//This function returns the Cauchy stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Kirchoff2Cauchy(const mat &tau, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(3,3);
    }
    else {
        return (1./J)*tau;
    }
}

//This function returns the Cauchy stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
vec Kirchoff2Cauchy(const vec& tau, const mat& F, const double& mJ) {

	double J = mJ;
	if (fabs(mJ) < sim_iota) {
		J = det(F);
	}
	//If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
	if (fabs(J) < sim_iota) {
		return zeros(6);
	}
	else {
		return (1. / J) * tau;
	}
}


//This function returns the first Piola-Kirchoff stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Kirchoff2PKI(const mat &tau, const mat &F, const double &mJ) {

    UNUSED(mJ);
    return tau*inv(F.t());
}

//This function returns the second Piola-Kirchoff stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat Kirchoff2PKII(const mat &tau, const mat &F, const double &mJ) {

    UNUSED(mJ);
    return inv(F)*tau*inv(F.t());
}

vec Kirchoff2PKII(const vec &tau, const mat &F, const double &mJ) {

    UNUSED(mJ);
    return t2v_stress(inv(F)*v2t_stress(tau)*inv(F.t()));
}


//This function returns the Kirchoff stress tensor from the first Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat PKI2Kirchoff(const mat &Sigma, const mat &F, const double &mJ) {

    UNUSED(mJ);
    return Sigma*F.t();
}

//This function returns the Kirchoff stress tensor from the second Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat PKII2Kirchoff(const mat &S, const mat &F, const double &mJ) {

    UNUSED(mJ);
    return F*S*F.t();
}

//This function returns the Cauchy stress tensor from the second Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat PKI2Cauchy(const mat &Sigma, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(3,3);
    }
    else {
        return (1./J)*Sigma*F.t();
    }
}

//This function returns the Cauchy stress tensor from the second Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
mat PKII2Cauchy(const mat &S, const mat &F, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < sim_iota) {
        J = det(F);
    }
    //If J is still less than a small value, we assume that sigma=tau=PK1=PKII = 0
    if (fabs(J) < sim_iota) {
        return zeros(3,3);
    }
    else {
        return (1./J)*F*S*F.t();
    }
}
    
} //namespace simcoon
