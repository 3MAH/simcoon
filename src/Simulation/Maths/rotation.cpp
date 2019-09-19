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

///@file rotation.cpp
///@brief rotation of a Voigt tensor
///@version 1.0

#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

vec rotate_vec(const vec &v, const mat &DR) {
    return DR*v;
}
    
vec rotate_vec(const vec &v, const double &alpha, const int &axis) {
    mat DR = fillR(alpha, axis);
    return rotate_vec(v, DR);
}

mat rotate_mat(const mat &m, const mat &DR) {
    return trans(DR)*m*DR;
}

mat rotate_mat(const mat &m, const double &alpha, const int &axis) {
    mat DR = fillR(alpha, axis);
    return rotate_mat(m, DR);
}

mat fillR(const double &alpha, const int &axis, const bool &active) {
    
    double act = 0;
    if (active)
    act = 1.;
    else
    act = -1.;
    
    mat R = zeros(3,3);
    double c = cos(alpha);
    double s = act*sin(alpha);
    
    switch(axis) {
        case 1: {
            R = { {1,0,0}, {0,c,-s}, {0,s,c} };
            break;
        }
        case 2: {
            R = mat {{c, 0, s}, {0,1,0}, {-s,0,c}};
            break;
        }
        case 3: {
            R = {{c,-s,0}, {s,c, 0}, {0,0,1}};
            break;
        }
        default: {
            cout << "Please choose the axis 1,2 or 3\n";
        }
    }
    return R;
}

mat fillR(const double &psi, const double &theta, const double &phi, const bool &active, const std::string &conv) {
    
    mat R = zeros(3,3);
    
    if(conv == "zxz") {
        double c1 = cos(psi);
        double s1 = sin(psi);
        double c2 = cos(theta);
        double s2 = sin(theta);
        double c3 = cos(phi);
        double s3 = sin(phi);

        if (active)
            R = { {c3*c1-c2*s3*s1,-c3*s1-c2*c1*s3,s3*s2}, {c1*s3+c3*c2*s1,c3*c2*c1-s3*s1,-c3*s2}, {s2*s1,c1*s2,c2} };
        else
            R = { {c3*c1-c2*s3*s1,c1*s3+c3*c2*s1,s2*s1}, {-c3*s1-c2*c1*s3,c3*c2*c1-s3*s1,c1*s2}, {s3*s2,-c3*s2,c2} };
    }
    else if(conv == "zyz") {
        double c1 = cos(psi);
        double s1 = sin(psi);
        double c2 = cos(theta);
        double s2 = sin(theta);
        double c3 = cos(phi);
        double s3 = sin(phi);

        if (active)
            R = { {c3*c2*c1-s3*s1,-c1*s3-c3*c2*s1,c3*s2}, {c3*s1+c2*c1*s3,c3*c1-c2*s3*s1,s3*s2}, {-c1*s2,s2*s1,c2} };
        else
            R = { {c3*c2*c1-s3*s1,c3*s1+c2*c1*s3,-c1*s2}, {-c1*s3-c3*c2*s1,c3*c1-c2*s3*s1,s2*s1}, {c3*s2,s3*s2,c2} };
    }
    else if(conv == "") {
            
        mat R1 = fillR(psi, axis_psi, active);
        mat R2 = fillR(theta, axis_theta, active);
        mat R3 = fillR(phi, axis_phi, active);
        
        R = R3*R2*R1;
    }
    else {
        cout << "error in Simulation/Maths/rotation.cpp: please provide a consistent convention for Euler rotation, i.e. zxz or zyz. The by-default convention has been utilized (as defined in parameter.hpp)" << endl;
        mat R1 = fillR(psi, axis_psi, active);
        mat R2 = fillR(theta, axis_theta, active);
        mat R3 = fillR(phi, axis_phi, active);
        
        R = R3*R2*R1;
    }
    
    return R;
}
    
mat fillQS(const double &alpha, const int &axis, const bool &active) {

    double act = 0;
    if (active) {
        act = 1.;
    }
    else {
        act = -1.;
    }
    
	mat QS = zeros(6,6);
	double c = cos(alpha);
	double s = act*sin(alpha);
	
	switch(axis) {
		case 1: {
			QS(0,0) = 1.;
			QS(0,1) = 0.;
			QS(0,2) = 0.;
			QS(0,3) = 0.;
			QS(0,4) = 0.;
			QS(0,5) = 0.;
			QS(1,0) = 0.;
			QS(1,1) = c*c;
			QS(1,2) = s*s;
			QS(1,3) = 0.;
			QS(1,4) = 0.;
			QS(1,5) = -2.*s*c;
			QS(2,0) = 0.;
			QS(2,1) = s*s;
			QS(2,2) = c*c;
			QS(2,3) = 0.;
			QS(2,4) = 0.;
			QS(2,5) = 2.*s*c;
			QS(3,0) = 0.;
			QS(3,1) = 0.;
			QS(3,2) = 0.;
			QS(3,3) = c;
			QS(3,4) = -s;
			QS(3,5) = 0.;
			QS(4,0) = 0.;
			QS(4,1) = 0.;
			QS(4,2) = 0.;
			QS(4,3) = s;
			QS(4,4) = c;
			QS(4,5) = 0.;
			QS(5,0) = 0.;
			QS(5,1) = s*c;
			QS(5,2) = -s*c;
			QS(5,3) = 0.;
			QS(5,4) = 0.;
			QS(5,5) = c*c-s*s;
			
			break;
		}
		case 2: {
			QS(0,0) = c*c;
			QS(0,1) = 0.;
			QS(0,2) = s*s;
			QS(0,3) = 0.;
			QS(0,4) = 2*c*s;
			QS(0,5) = 0.;
			QS(1,0) = 0.;
			QS(1,1) = 1.;
			QS(1,2) = 0.;
			QS(1,3) = 0.;
			QS(1,4) = 0.;
			QS(1,5) = 0.;
			QS(2,0) = s*s;
			QS(2,1) = 0.;
			QS(2,2) = c*c;
			QS(2,3) = 0.;
			QS(2,4) = -2*c*s;
			QS(2,5) = 0.;
			QS(3,0) = 0.;
			QS(3,1) = 0.;
			QS(3,2) = 0.;
			QS(3,3) = c;
			QS(3,4) = 0.;
			QS(3,5) = s;
			QS(4,0) = -c*s;
			QS(4,1) = 0.;
			QS(4,2) = c*s;
			QS(4,3) = 0.;
			QS(4,4) = c*c-s*s;
			QS(4,5) = 0.;
			QS(5,0) = 0.;
			QS(5,1) = 0.;
			QS(5,2) = 0.;
			QS(5,3) = -s;
			QS(5,4) = 0.;
			QS(5,5) = c;
			
			break;
			
		}
		case 3: {
			QS(0,0) = c*c;
			QS(0,1) = s*s;
			QS(0,2) = 0.;
			QS(0,3) = -2*s*c;
			QS(0,4) = 0.;
			QS(0,5) = 0.;
			QS(1,0) = s*s;
			QS(1,1) = c*c;
			QS(1,2) = 0.;
			QS(1,3) = 2*s*c;
			QS(1,4) = 0.;
			QS(1,5) = 0.;
			QS(2,0) = 0.;
			QS(2,1) = 0.;
			QS(2,2) = 1.;
			QS(2,3) = 0.;
			QS(2,4) = 0.;
			QS(2,5) = 0.;
			QS(3,0) = s*c;
			QS(3,1) = -s*c;
			QS(3,2) = 0.;
			QS(3,3) = c*c-s*s;
			QS(3,4) = 0.;
			QS(3,5) = 0.;
			QS(4,0) = 0.;
			QS(4,1) = 0.;
			QS(4,2) = 0.;
			QS(4,3) = 0.;
			QS(4,4) = c;
			QS(4,5) = -s;
			QS(5,0) = 0.;
			QS(5,1) = 0.;
			QS(5,2) = 0.;
			QS(5,3) = 0.;
			QS(5,4) = s;
			QS(5,5) = c;
						
			break;
		}
		default: {
			cout << "Please choose the axis 1,2 or 3\n";
		}
	}
	return QS;
}

mat fillQS(const mat &DR, const bool &active) {

    double a = 0.;
    double d = 0.;
    double g = 0.;
    double b = 0.;
    double e = 0.;
    double h = 0.;
    double c = 0.;
    double f = 0.;
    double i = 0.;
    
    if(active) {
        a = DR(0,0);
        b = DR(0,1);
        c = DR(0,2);
        d = DR(1,0);
        e = DR(1,1);
        f = DR(1,2);
        g = DR(2,0);
        h = DR(2,1);
        i = DR(2,2);
    }
    else{
        a = DR(0,0);
        d = DR(0,1);
        g = DR(0,2);
        b = DR(1,0);
        e = DR(1,1);
        h = DR(1,2);
        c = DR(2,0);
        f = DR(2,1);
        i = DR(2,2);
    }
    
    mat QS= zeros(6,6);
    QS(0,0) = a*a;
    QS(0,1) = b*b;
    QS(0,2) = c*c;
    QS(0,3) = 2.*a*b;
    QS(0,4) = 2.*a*c;
    QS(0,5) = 2.*b*c;
    QS(1,0) = d*d;
    QS(1,1) = e*e;
    QS(1,2) = f*f;
    QS(1,3) = 2.*d*e;
    QS(1,4) = 2.*d*f;
    QS(1,5) = 2.*e*f;
    QS(2,0) = g*g;
    QS(2,1) = h*h;
    QS(2,2) = i*i;
    QS(2,3) = 2.*g*h;
    QS(2,4) = 2.*g*i;
    QS(2,5) = 2.*h*i;
    QS(3,0) = a*d;
    QS(3,1) = b*e;
    QS(3,2) = c*f;
    QS(3,3) = d*b+a*e;
    QS(3,4) = d*c+a*f;
    QS(3,5) = e*c+b*f;
    QS(4,0) = a*g;
    QS(4,1) = b*h;
    QS(4,2) = c*i;
    QS(4,3) = g*b+a*h;
    QS(4,4) = g*c+a*i;
    QS(4,5) = h*c+b*i;
    QS(5,0) = d*g;
    QS(5,1) = e*h;
    QS(5,2) = f*i;
    QS(5,3) = g*e+d*h;
    QS(5,4) = g*f+d*i;
    QS(5,5) = h*f+e*i;
    
    return QS;
}
    
mat fillQE(const double &alpha, const int &axis, const bool &active) {

    double act = 0;
    if (active) {
        act = 1.;
    }
    else {
        act = -1.;
    }
    
	mat QE = zeros(6,6);
    double c = cos(alpha);
    double s = act*sin(alpha);
	
	switch(axis) {
		case 1: {
			QE(0,0) = 1.;
			QE(0,1) = 0.;
			QE(0,2) = 0.;
			QE(0,3) = 0.;
			QE(0,4) = 0.;
			QE(0,5) = 0.;
			QE(1,0) = 0.;
			QE(1,1) = c*c;
			QE(1,2) = s*s;
			QE(1,3) = 0.;
			QE(1,4) = 0.;
			QE(1,5) = -s*c;
			QE(2,0) = 0.;
			QE(2,1) = s*s;
			QE(2,2) = c*c;
			QE(2,3) = 0.;
			QE(2,4) = 0.;
			QE(2,5) = s*c;
			QE(3,0) = 0.;
			QE(3,1) = 0.;
			QE(3,2) = 0.;
			QE(3,3) = c;
			QE(3,4) = -s;
			QE(3,5) = 0.;
			QE(4,0) = 0.;
			QE(4,1) = 0.;
			QE(4,2) = 0.;
			QE(4,3) = s;
			QE(4,4) = c;
			QE(4,5) = 0.;
			QE(5,0) = 0.;
			QE(5,1) = 2.*s*c;
			QE(5,2) = -2.*s*c;
			QE(5,3) = 0.;
			QE(5,4) = 0.;
			QE(5,5) = c*c-s*s;
			break;
		}
		case 2: {
			QE(0,0) = c*c;
			QE(0,1) = 0.;
			QE(0,2) = s*s;
			QE(0,3) = 0.;
			QE(0,4) = c*s;
			QE(0,5) = 0.;
			QE(1,0) = 0.;
			QE(1,1) = 1.;
			QE(1,2) = 0.;
			QE(1,3) = 0.;
			QE(1,4) = 0.;
			QE(1,5) = 0.;
			QE(2,0) = s*s;
			QE(2,1) = 0.;
			QE(2,2) = c*c;
			QE(2,3) = 0.;
			QE(2,4) = -c*s;
			QE(2,5) = 0.;
			QE(3,0) = 0.;
			QE(3,1) = 0.;
			QE(3,2) = 0.;
			QE(3,3) = c;
			QE(3,4) = 0.;
			QE(3,5) = s;
			QE(4,0) = -2.*c*s;
			QE(4,1) = 0.;
			QE(4,2) = 2.*c*s;
			QE(4,3) = 0.;
			QE(4,4) = c*c-s*s;
			QE(4,5) = 0.;
			QE(5,0) = 0.;
			QE(5,1) = 0.;
			QE(5,2) = 0.;
			QE(5,3) = -s;
			QE(5,4) = 0.;
			QE(5,5) = c;
			break;
		}
		case 3: {
			QE(0,0) = c*c;
			QE(0,1) = s*s;
			QE(0,2) = 0.;
			QE(0,3) = -s*c;
			QE(0,4) = 0.;
			QE(0,5) = 0.;
			QE(1,0) = s*s;
			QE(1,1) = c*c;
			QE(1,2) = 0.;
			QE(1,3) = s*c;
			QE(1,4) = 0.;
			QE(1,5) = 0.;
			QE(2,0) = 0.;
			QE(2,1) = 0.;
			QE(2,2) = 1.;
			QE(2,3) = 0.;
			QE(2,4) = 0.;
			QE(2,5) = 0.;
			QE(3,0) = 2.*s*c;
			QE(3,1) = -2.*s*c;
			QE(3,2) = 0.;
			QE(3,3) = c*c-s*s;
			QE(3,4) = 0.;
			QE(3,5) = 0.;
			QE(4,0) = 0.;
			QE(4,1) = 0.;
			QE(4,2) = 0.;
			QE(4,3) = 0.;
			QE(4,4) = c;
			QE(4,5) = -s;
			QE(5,0) = 0.;
			QE(5,1) = 0.;
			QE(5,2) = 0.;
			QE(5,3) = 0.;
			QE(5,4) = s;
			QE(5,5) = c;
			break;
		}
		default: {
			cout << "Please choose the axis 1,2 or 3\n";
		}
	}
	return QE;		
}
    
mat fillQE(const mat &DR, const bool &active) {
    
    double a = 0.;
    double d = 0.;
    double g = 0.;
    double b = 0.;
    double e = 0.;
    double h = 0.;
    double c = 0.;
    double f = 0.;
    double i = 0.;
    
    if(active) {
        a = DR(0,0);
        b = DR(0,1);
        c = DR(0,2);
        d = DR(1,0);
        e = DR(1,1);
        f = DR(1,2);
        g = DR(2,0);
        h = DR(2,1);
        i = DR(2,2);
    }
    else{
        a = DR(0,0);
        d = DR(0,1);
        g = DR(0,2);
        b = DR(1,0);
        e = DR(1,1);
        h = DR(1,2);
        c = DR(2,0);
        f = DR(2,1);
        i = DR(2,2);
    }

    mat QE= zeros(6,6);
    QE(0,0) = a*a;
    QE(0,1) = b*b;
    QE(0,2) = c*c;
    QE(0,3) = a*b;
    QE(0,4) = a*c;
    QE(0,5) = b*c;
    QE(1,0) = d*d;
    QE(1,1) = e*e;
    QE(1,2) = f*f;
    QE(1,3) = d*e;
    QE(1,4) = d*f;
    QE(1,5) = e*f;
    QE(2,0) = g*g;
    QE(2,1) = h*h;
    QE(2,2) = i*i;
    QE(2,3) = g*h;
    QE(2,4) = g*i;
    QE(2,5) = h*i;
    QE(3,0) = 2.*a*d;
    QE(3,1) = 2.*b*e;
    QE(3,2) = 2.*c*f;
    QE(3,3) = d*b+a*e;
    QE(3,4) = d*c+a*f;
    QE(3,5) = e*c+b*f;
    QE(4,0) = 2.*a*g;
    QE(4,1) = 2.*b*h;
    QE(4,2) = 2.*c*i;
    QE(4,3) = g*b+a*h;
    QE(4,4) = g*c+a*i;
    QE(4,5) = h*c+b*i;
    QE(5,0) = 2.*d*g;
    QE(5,1) = 2.*e*h;
    QE(5,2) = 2.*f*i;
    QE(5,3) = g*e+d*h;
    QE(5,4) = g*f+d*i;
    QE(5,5) = h*f+e*i;
    
    return QE;
}


//To rotate a stiffness matrix (6,6)
mat rotateL(const mat &L, const double &alpha, const int &axis, const bool &active) {

	mat QS = fillQS(alpha, axis, active);
	return QS*(L*trans(QS));
}

mat rotateL(const mat &L, const mat &DR, const bool &active) {
    
    mat QS = fillQS(DR, active);
    return QS*(L*trans(QS));
}
    
    
//To rotate a compliance matrix (6,6)
mat rotateM(const mat &M, const double &alpha, const int &axis, const bool &active) {

	mat QE = fillQE(alpha, axis, active);
	return QE*(M*trans(QE));
}

mat rotateM(const mat &M, const mat &DR, const bool &active) {
    
    mat QE = fillQE(DR, active);
    return QE*(M*trans(QE));
}
    
    
//To rotate an interaction matrix of type A (strain) (6,6)
mat rotateA(const mat &A, const double &alpha, const int &axis, const bool &active) {
        
    mat QE = fillQE(alpha, axis, active);
    mat QS = fillQS(alpha, axis, active);
    
    return QE*(A*trans(QS));
}

mat rotateA(const mat &A, const mat &DR, const bool &active) {
    
    mat QE = fillQE(DR, active);
    mat QS = fillQS(DR, active);
    
    return QE*(A*trans(QS));
}
    
//To rotate an interaction matrix of type B (stress) (6,6)
mat rotateB(const mat &B, const double &alpha, const int &axis, const bool &active) {
    
    mat QE = fillQE(alpha, axis, active);
    mat QS = fillQS(alpha, axis, active);
    
    return QS*(B*trans(QE));
}

mat rotateB(const mat &B, const mat &DR, const bool &active) {
    
    mat QE = fillQE(DR, active);
    mat QS = fillQS(DR, active);
    
    return QS*(B*trans(QE));
}
    
//To rotate a stress vector (6)
vec rotate_stress(const vec &V, const double &alpha, const int &axis, const bool &active) {

	mat QS = fillQS(alpha, axis, active);
	return QS*V;
}

//To rotate a stress vector (6)
vec rotate_stress(const vec &V, const mat &DR, const bool &active) {
    
    mat QS = fillQS(DR, active);
    return QS*V;
}
    
//To rotate a strain vector (6)
vec rotate_strain(const vec &V, const double &alpha, const int &axis, const bool &active) {

	mat QE = fillQE(alpha, axis, active);
	return QE*V;
}

vec rotate_strain(const vec &V, const mat &DR, const bool &active) {
    
    mat QE = fillQE(DR, active);
    return QE*V;
}
    

//To rotate from local to global a strain tensor (6) using Euler angles
mat rotate_l2g_strain(const vec &E, const double &psi, const double &theta, const double &phi) {
    
    mat E_temp = E;
    if(fabs(phi) > sim_iota) {
        E_temp = rotate_strain(E_temp, -phi, axis_phi, false);
    }
    if(fabs(theta) > sim_iota) {
        E_temp = rotate_strain(E_temp, -theta, axis_theta, false);
    }
    if(fabs(psi) > sim_iota) {
        E_temp = rotate_strain(E_temp, -psi, axis_psi, false);
    }
    
    return E_temp;
}
    
//To rotate from global to local a strain tensor (6) using Euler angles
mat rotate_g2l_strain(const vec &E, const double &psi, const double &theta, const double &phi) {
    
    mat E_temp = E;
    if(fabs(psi) > sim_iota) {
        E_temp = rotate_strain(E_temp, psi, axis_psi, false);
    }
    if(fabs(theta) > sim_iota) {
        E_temp = rotate_strain(E_temp, theta, axis_theta, false);
    }
    if(fabs(phi) > sim_iota) {
        E_temp = rotate_strain(E_temp, phi, axis_phi, false);
    }
    
    return E_temp;
}

//To rotate from local to global a stress tensor (6)
mat rotate_l2g_stress(const vec &S, const double &psi, const double &theta, const double &phi) {
    
    mat S_temp = S;
    if(fabs(phi) > sim_iota) {
        S_temp = rotate_stress(S_temp, -phi, axis_phi, false);
    }
    if(fabs(theta) > sim_iota) {
        S_temp = rotate_stress(S_temp, -theta, axis_theta, false);
    }
    if(fabs(psi) > sim_iota) {
        S_temp = rotate_stress(S_temp, -psi, axis_psi, false);
    }
    
    return S_temp;
}

//To rotate from global to local a stress tensor (6)
mat rotate_g2l_stress(const vec &S, const double &psi, const double &theta, const double &phi) {
    
    mat S_temp = S;
    if(fabs(psi) > sim_iota) {
        S_temp = rotate_stress(S_temp, psi, axis_psi, false);
    }
    if(fabs(theta) > sim_iota) {
        S_temp = rotate_stress(S_temp, theta, axis_theta, false);
    }
    if(fabs(phi) > sim_iota) {
        S_temp = rotate_stress(S_temp, phi, axis_phi, false);
    }
    
    return S_temp;
}
    
//To rotate from local to global a stiffness matrix (6,6)
mat rotate_l2g_L(const mat &Lt, const double &psi, const double &theta, const double &phi) {
    
    mat Lt_temp = Lt;
  	if(fabs(phi) > sim_iota) {
		Lt_temp = rotateL(Lt_temp, -phi, axis_phi, false);
	}
  	if(fabs(theta) > sim_iota) {
		Lt_temp = rotateL(Lt_temp, -theta, axis_theta, false);
	}
	if(fabs(psi) > sim_iota) {
		Lt_temp = rotateL(Lt_temp, -psi, axis_psi, false);
	}
    
	return Lt_temp;
}

//To rotate from global to local a stiffness matrix (6,6)
mat rotate_g2l_L(const mat &Lt, const double &psi, const double &theta, const double &phi) {
    
    mat Lt_temp = Lt;
  	if(fabs(psi) > sim_iota) {
		Lt_temp = rotateL(Lt_temp, psi, axis_psi, false);
	}
	if(fabs(theta) > sim_iota) {
		Lt_temp = rotateL(Lt_temp, theta, axis_theta, false);
	}
	if(fabs(phi) > sim_iota) {
		Lt_temp = rotateL(Lt_temp, phi, axis_phi, false);
    }
    
	return Lt_temp;
}

//To rotate from local to global a  localisation matrix (6,6)
mat rotate_l2g_M(const mat &M, const double &psi, const double &theta, const double &phi) {
    
    mat M_temp = M;
    if(fabs(phi) > sim_iota) {
        M_temp = rotateM(M_temp, -phi, axis_phi, false);
    }
    if(fabs(theta) > sim_iota) {
        M_temp = rotateM(M_temp, -theta, axis_theta, false);
    }
    if(fabs(psi) > sim_iota) {
        M_temp = rotateM(M_temp, -psi, axis_psi, false);
    }
    
    return M_temp;
}

//To rotate from global to local a localisation matrix (6,6)
mat rotate_g2l_M(const mat &M, const double &psi, const double &theta, const double &phi) {
    
    mat M_temp = M;
    if(fabs(psi) > sim_iota) {
        M_temp = rotateM(M_temp, psi, axis_psi, false);
    }
    if(fabs(theta) > sim_iota) {
        M_temp = rotateM(M_temp, theta, axis_theta, false);
    }
    if(fabs(phi) > sim_iota) {
        M_temp = rotateM(M_temp, phi, axis_phi, false);
    }
    
    return M_temp;
}

//To rotate from local to global a strain localisation matrix (6,6)
mat rotate_l2g_A(const mat &A, const double &psi, const double &theta, const double &phi) {
    
    mat A_temp = A;
  	if(fabs(phi) > sim_iota) {
		A_temp = rotateA(A_temp, -phi, axis_phi, false);
	}
  	if(fabs(theta) > sim_iota) {
		A_temp = rotateA(A_temp, -theta, axis_theta, false);
	}
	if(fabs(psi) > sim_iota) {
		A_temp = rotateA(A_temp, -psi, axis_psi, false);
	}
    
	return A_temp;
}

//To rotate from global to local a strain localisation matrix (6,6)
mat rotate_g2l_A(const mat &A, const double &psi, const double &theta, const double &phi) {
    
    mat A_temp = A;
  	if(fabs(psi) > sim_iota) {
		A_temp = rotateA(A_temp, psi, axis_psi, false);
	}
	if(fabs(theta) > sim_iota) {
		A_temp = rotateA(A_temp, theta, axis_theta, false);
	}
	if(fabs(phi) > sim_iota) {
		A_temp = rotateA(A_temp, phi, axis_phi, false);
    }
    
	return A_temp;
}

//To rotate from local to global a stress localisation matrix (6,6)
mat rotate_l2g_B(const mat &B, const double &psi, const double &theta, const double &phi) {
    
    mat B_temp = B;
  	if(fabs(phi) > sim_iota) {
		B_temp = rotateB(B_temp, -phi, axis_phi, false);
	}
  	if(fabs(theta) > sim_iota) {
		B_temp = rotateB(B_temp, -theta, axis_theta, false);
	}
	if(fabs(psi) > sim_iota) {
		B_temp = rotateB(B_temp, -psi, axis_psi, false);
	}
    
	return B_temp;
}

//To rotate from global to local a stress localisation matrix (6,6)
mat rotate_g2l_B(const mat &B, const double &psi, const double &theta, const double &phi) {
    
    mat B_temp = B;
  	if(fabs(psi) > sim_iota) {
		B_temp = rotateB(B_temp, psi, axis_psi, false);
	}
	if(fabs(theta) > sim_iota) {
		B_temp = rotateB(B_temp, theta, axis_theta, false);
	}
	if(fabs(phi) > sim_iota) {
		B_temp = rotateB(B_temp, phi, axis_phi, false);
    }
    
	return B_temp;
}
    
} //namespace simcoon
