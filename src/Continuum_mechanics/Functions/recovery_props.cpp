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

///@file recovery_props.cpp
///@brief A set of function that allow to get the Elastic properties from stiffness/compliance tensors
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/recovery_props.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void check_symetries(const mat &L, std::string &umat_type, int &axis, vec &props, int &maj_sym, const double &tol) {
 
    double max_tol_simcoon::limit = simcoon::limit;
    if(tol > simcoon::limit) {
        max_tol_simcoon::limit = tol;
    }

    axis = 0; //Indicate no preferential axis
    //First thing to do is to check how symtric is the tensor:
    mat L_sym = 0.5*(L + trans(L));
    double check_L_sym_diff = norm(L-L_sym,2);
    mat Q = zeros(3,3); //rotation/reflexion matrix
    mat L_test = zeros(6,6);
    
    
    mat check_sym = zeros(3,3);  //First col : reflections around x, y and z
                                //second col : 90deg rotation around x, y, and z
                                //Third col : Equality between constants
    Mat<int> type_sym(3,3);
    type_sym.zeros();
        
    //Everything is checked based on a "symmtrized" matrix
    //Check reflexions
    //reflexion around the z axis - symmetry plane xy
    Q = { {1,0,0}, {0,1,0}, {0,0,-1} };
    L_test = rotateL(L_sym, Q);
    check_sym(2,0) = norm(L_sym - L_test,2);
    //reflexion around the y axis - symmetry plane xz
    Q = { {1,0,0}, {0,-1,0}, {0,0,1} };
    L_test = rotateL(L_sym, Q);
    check_sym(1,0) = norm(L_sym - L_test,2);
    //reflexion around the x axis - symmetry plane yz
    Q = { {-1,0,0}, {0,1,0}, {0,0,1} };
    L_test = rotateL(L_sym, Q);
    check_sym(0,0) = norm(L_sym - L_test,2);
    
    //Check 90deg rotations
    //90 rotation around the z axis
    Q = { {0,1,0}, {-1,0,0}, {0,0,1} };
    L_test = rotateL(L_sym, Q);
    check_sym(2,1) = norm(L_sym - L_test,2);
    //90 rotation around the y axis
    Q = { {0,0,1}, {0,1,0}, {-1,0,0} };
    L_test = rotateL(L_sym, Q);
    check_sym(1,1) = norm(L_sym - L_test,2);
    //90 rotation around the x axis
    Q = { {1,0,0}, {0,0,1}, {0,-1,0} };
    L_test = rotateL(L_sym, Q);
    check_sym(0,1) = norm(L_sym - L_test,2);

    //Check equality between constants
    //All rotation around the z axis
    check_sym(2,2) = fabs(L_sym(3,3) - (0.5*(L_sym(0,0)-L_sym(0,1))));
    //All rotation around the y axis
    check_sym(1,2) = fabs(L_sym(4,4) - (0.5*(L_sym(0,0)-L_sym(0,2))));
    //All rotation around the x axis
    check_sym(0,2) = fabs(L_sym(5,5) - (0.5*(L_sym(1,1)-L_sym(1,2))));
    
    for (unsigned int i=0; i<check_sym.n_rows; i++) {
        for (unsigned int j=0; j<check_sym.n_cols; j++) {
            if(check_sym(i,j) < max_tol_simcoon::limit)
                type_sym(i,j) = 1;
            else
                type_sym(i,j) = 0;
        }
    }
        
    if (check_L_sym_diff < max_tol_simcoon::limit)
        maj_sym = 1;
    else
        maj_sym = 0;

    //Anisotropic / Triclinic: no symmetry planes, fully anisotropic.
    umat_type = "ELANI";
    
    //check Monoclinic
    //One symmetry plane (need to check all the three planes..)
    if(type_sym(0,0) + type_sym(1,0) + type_sym(2,0) == 1) {
        umat_type = "ELMON";
        if(type_sym(0,0) == 1)
            axis = 1;
        else if(type_sym(1,0) == 1)
            axis = 2;
        else if(type_sym(2,0) == 1)
            axis = 3;
        else {
            cout << "Error in the Monoclinic check";
            exit(0);
        }
    }
    
    //Check  orthotropic
    //Orthotropic: three mutually orthogonal planes of reflection. Symmetry. a!=b=c, alpha=beta=gamma=90 Number of independent coefficients: 9
    //Symmetry transformations: reflections about all three orthogonal planes
    if(type_sym(0,0) + type_sym(1,0) + type_sym(2,0) == 3) {
        umat_type = "ELORT";
    }
    
    //Check  cubic
    //three mutually orthogonal planes of reflection symmetry
    //Plus 90deg rotation symmetry with respect to those planes. a=b=c, alpha=beta=gamma=90 Number of independent coefficients: 3
    //Symmetry transformations: reflections and 90 rotations about all three orthogonal planes
    if((type_sym(0,0) + type_sym(1,0) + type_sym(2,0) == 3) && (type_sym(0,1) + type_sym(1,1) + type_sym(2,1) == 3)) {
        umat_type = "ELCUB";
    }

    //Check Transversely isotropic
    //Transversely isotropic: The physical properties are symmetric about an axis that is normal to a plane of isotropy (xy-plane in the figure).
    //Three mutually orthogonal planes of reflection symmetry and axial symmetry with respect to z-axis. Number of independent coefficients: 5
    //Symmetry transformations: reflections about all three orthogonal planes plus all rotations about z- axis, x- axis or y - axis
    if((type_sym(0,0) + type_sym(1,0) + type_sym(2,0) == 3) && (type_sym(0,1) + type_sym(1,1) + type_sym(2,1) == 1) && (type_sym(0,2) + type_sym(1,2) + type_sym(2,2) == 1)) {
        if(type_sym(0,1) + type_sym(0,2) == 2) {
            umat_type = "ELIST";
            axis = 1;
        }
        else if(type_sym(1,1) + type_sym(1,2) == 2) {
            umat_type = "ELIST";
            axis = 2;
        }
        else if(type_sym(2,1) + type_sym(2,2) == 2) {
            umat_type = "ELIST";
            axis = 3;
        }
        else {
            cout << "Error in the Transversely isotropic check";
            exit(0);
        }
    }

    //Check Isotropic
    //Isotropic: The physical properties are symmetric about all axis
    //Three mutually orthogonal planes of reflection symmetry and axial symmetry with respect to z-axis. Number of independent coefficients: 3
    //Symmetry transformations: reflections about all three orthogonal planes plus all rotations about z- axis, x- axis and y - axis
    if((type_sym(0,0) + type_sym(1,0) + type_sym(2,0) == 3) && (type_sym(0,1) + type_sym(1,1) + type_sym(2,1) == 3) && (type_sym(0,2) + type_sym(1,2) + type_sym(2,2) == 3)) {
        umat_type = "ELISO";
    }
        
    //Here is the full list of possibilities
    //Triclinic -- fully anisotropic - 21 independant parameters
    //Monoclinic - xy                - 15 independant parameters
    //Monoclinic - xz                - 15 independant parameters
    //Monoclinic - yz                - 15 independant parameters
    //Orthotropic                    - 9 independant parameters
    //Cubic                          - 3 independant parameters
    //Transversely isotropic x       - 5 independant parameters
    //Transversely isotropic y       - 5 independant parameters
    //Transversely isotropic z       - 5 independant parameters
    //Isotropic                      - 2 independant parameters
    
    if (umat_type == "ELISO") {
        props = L_iso_props(L_sym);
    }
    else if(umat_type == "ELIST") {
        props = L_isotrans_props(L_sym, axis);
    }
    else if(umat_type == "ELCUB") {
        props = L_cubic_props(L_sym);
    }
    else if(umat_type == "ELORT") {
        props = L_ortho_props(L_sym);
    }
    else
        props = zeros(1);
}
 
vec L_iso_props(const mat &Lt) {
    
    double mu = (1./3.)*(Lt(3,3) + Lt(4,4) + Lt(5,5));
    double lambda = (1./6)*(Lt(0,1)+Lt(0,2)+Lt(1,2)+Lt(1,0)+Lt(2,0)+Lt(2,1));
    
    double E = mu*(3.+2.*(mu/lambda))/(1+(mu/lambda));
    double nu = 0.5/(1+(mu/lambda));
    
    vec props = zeros(2);
    props = {E,nu};
    return props;
}

vec M_iso_props(const mat &Mt) {
    
    double E = (1./3.)*((1./Mt(0,0))+(1./Mt(1,1))+(1./Mt(2,2)));
    double nu = (-E/6.)*(Mt(0,1)+Mt(0,2)+Mt(1,2)+Mt(1,0)+Mt(2,0)+Mt(2,1));

    vec props = zeros(2);
    props = {E,nu};
    return props;
}

vec L_isotrans_props(const mat &Lt, const int &axis) {
    
    double yiiii;
    double yjjjj;
    double yiijj;
    double yjjkk;
    double yijij;
    
    double EL;
    double ET;
    double nuTL;
    double nuTT;
    double GLT;
    
    if (axis==1) {
        
        yiiii = Lt(0,0);
        yjjjj = ( Lt(1,1) + Lt(2,2) )/2.;
        
        yiijj = ( Lt(0,1) + Lt(0,2) + Lt(2,0) + Lt(1,0) )/4.;
        yjjkk = ( Lt(1,2) + Lt(2,1) )/2.;
        
        yijij = ( Lt(3,3) + Lt(4,4) )/2.;
    }
    else if (axis==2) {
        
        yiiii = Lt(1,1);
        yjjjj = ( Lt(0,0) + Lt(2,2) )/2.;
        
        yiijj = ( Lt(1,2) + Lt(1,0) + Lt(0,1) + Lt(2,1) )/4.;
        yjjkk = ( Lt(0,2) + Lt(2,0) )/2.;
        
        yijij = ( Lt(3,3) + Lt(5,5) )/2.;
    }
    else if (axis==3) {
        
        yiiii = Lt(2,2);
        yjjjj = ( Lt(0,0) + Lt(1,1) )/2.;
        
        yiijj = ( Lt(0,2) + Lt(1,2) + Lt(2,0) + Lt(2,1) )/4.;
        yjjkk = ( Lt(0,1) + Lt(1,0) )/2.;
        
        yijij = ( Lt(4,4) + Lt(5,5) )/2.;
    }
    else {
        cout << "Wrong axis !" << endl;
        exit(0);
    }
    
    EL = yiiii - 2*yiijj*yiijj/(yjjjj+yjjkk);
    ET = ( yjjjj-yjjkk )*( yjjjj*yiiii+yjjkk*yiiii-2*yiijj*yiijj )/( yjjjj*yiiii-yiijj*yiijj );
    
    double Delta = (yjjjj-yjjkk)*((yjjjj+yjjkk)*yiiii-2*pow(yiijj,2));
    nuTL = ET*yiijj*(yjjjj-yjjkk)/Delta;
    nuTT = ET*(yjjkk*yiiii-pow(yiijj,2))/Delta;
    
    GLT = yijij;
    
    vec props = zeros(5);
    
    props = {EL,ET,nuTL,nuTT,GLT};
    
    return props;
}

vec M_isotrans_props(const mat &Mt, const int &axis) {
    
    double EL;
    double ET;
    double nuTL;
    double nuTT;
    double GLT;
    
    if (axis==1) {
        EL = 1./Mt(0,0);
        ET = 1./Mt(1,1);
        
        nuTT = (( Mt(5,5)/2 - (Mt(1,2)+Mt(2,1)) )*ET - 1 )/3. ;
        nuTL = -ET*(Mt(1,0)+Mt(2,0)+Mt(0,1)+Mt(0,2))/4.;
        
        GLT = (1/Mt(3,3) + 1/Mt(4,4))/2;
    }
    else if (axis==2) {
        EL = 1./Mt(1,1);
        ET = (1/Mt(0,0)+1/Mt(2,2))/2;
        
        nuTT = (( Mt(4,4)/2 - (Mt(0,2)+Mt(2,0)) )*ET - 1 )/3. ;
        nuTL = -(Mt(0,1)+Mt(2,1)+Mt(1,0)+Mt(1,2))*ET/4.;
        
        GLT = (1/Mt(3,3) + 1/Mt(5,5))/2;
    }
    else if (axis==3) {
        EL = 1./Mt(2,2);
        ET = (1./Mt(0,0)+1./Mt(1,1))/2.;
        
        nuTT = (( Mt(3,3)/2 - (Mt(0,1)+Mt(1,0)) )*ET - 1 )/3. ;
        nuTL = -(Mt(0,2)+Mt(1,2)+Mt(2,0)+Mt(2,1))*ET/4.;
        
        GLT = (1/Mt(4,4) + 1/Mt(5,5))/2;
    }
    else {
        cout << "Error with recovery_MisoT_props : Wrong axis !" << endl;
        exit(0);
    }
    
    vec props = zeros(5);
    
    props = {EL,ET,nuTL,nuTT,GLT};
    
    return props;
}

vec L_cubic_props(const mat &Lt) {
    
    double yiiii = ( Lt(0,0) + Lt(1,1) + Lt(2,2) )/3.;
    double yiijj = ( Lt(0,1) + Lt(1,0) + Lt(0,2) + Lt(2,0) + Lt(1,2) + Lt(2,1) )/6.;
    double yijij = ( Lt(3,3) + Lt(4,4) + Lt(5,5))/3.;
    
    double nu = 1 / (1 + yiiii/yiijj);
    double E = yiiii*( 1 - 3*pow(nu,2) - 2*pow(nu,3) )/( 1-pow(nu,2) );
    double G = yijij;
    
    vec props = zeros(3);
    
    props = {E,nu,G};
    
    return props;
}

vec M_cubic_props(const mat &Mt) {
    
    double E = ( 1/Mt(0,0) + 1/Mt(1,1) + 1/Mt(2,2) )/3.;
    
    double nu = -E*( Mt(0,1) + Mt(1,0) + Mt(0,2) + Mt(2,0) + Mt(1,2) + Mt(2,1) )/6.;
    
    double G = ( 1/Mt(3,3) + 1/Mt(4,4) + 1/Mt(5,5))/3.;
    
    vec props = zeros(3);
    
    props = {E,nu,G};
    
    return props;
}

vec L_ortho_props(const mat &Lt) {
    
    mat Mt = inv(Lt);
    
    double E1 = 1/Mt(0,0);
    double E2 = 1/Mt(1,1);
    double E3 = 1/Mt(2,2);
    
    double nu12 = -E1*( Mt(0,1) + Mt(1,0) )/2.;
    double nu13 = -E1*( Mt(0,2) + Mt(2,0) )/2.;
    double nu23 = -E2*( Mt(1,2) + Mt(2,1) )/2.;
    
    double G23 = 1/Mt(5,5);
    double G13 = 1/Mt(4,4);
    double G12 = 1/Mt(3,3);
    
    vec props = zeros(9);
    props = {E1,E2,E3,nu12,nu13,nu23,G12,G13,G23};
    
    return props;
}

vec M_ortho_props(const mat &Mt) {
    
    double E1 = 1/Mt(0,0);
    double E2 = 1/Mt(1,1);
    double E3 = 1/Mt(2,2);
    
    double nu12 = -E1*( Mt(0,1) + Mt(1,0) )/2.;
    double nu13 = -E1*( Mt(0,2) + Mt(2,0) )/2.;
    double nu23 = -E2*( Mt(1,2) + Mt(2,1) )/2.;
    
    double G23 = 1/Mt(5,5);
    double G13 = 1/Mt(4,4);
    double G12 = 1/Mt(3,3);
    
    vec props = zeros(9);
    props = {E1,E2,E3,nu12,nu13,nu23,G12,G13,G23};
    
    return props;
}
    
vec M_aniso_props(const mat &Mt) {

	double E1 = 1/Mt(0,0);
	double E2 = 1/Mt(1,1);
	double E3 = 1/Mt(2,2);
	
    double nu12 = -E1*( Mt(0,1) + Mt(1,0) )/2.;
    double nu13 = -E1*( Mt(0,2) + Mt(2,0) )/2.;
    double nu23 = -E2*( Mt(1,2) + Mt(2,1) )/2.;
    
	double G23 = 1/Mt(5,5);
	double G13 = 1/Mt(4,4);
	double G12 = 1/Mt(3,3);
    
	double eta14 = E1*( Mt(0,3) + Mt(3,0) )/2.;
	double eta15 = E1*( Mt(0,4) + Mt(4,0) )/2.;
	double eta16 = E1*( Mt(0,5) + Mt(5,0) )/2.;
	
	double eta24 = E2*( Mt(1,3) + Mt(3,1) )/2.;
	double eta25 = E2*( Mt(1,4) + Mt(4,1) )/2.;
	double eta26 = E2*( Mt(1,5) + Mt(5,1) )/2.;
	
	double eta34 = E3*( Mt(2,3) + Mt(3,2) )/2.;
	double eta35 = E3*( Mt(2,4) + Mt(4,2) )/2.;
	double eta36 = E3*( Mt(2,5) + Mt(5,2) )/2.;
	
    double eta45 = G12*( Mt(3,4) + Mt(4,3) )/2.;
    double eta46 = G12*( Mt(3,5) + Mt(5,3) )/2.;
    double eta56 = G13*( Mt(4,5) + Mt(5,4) )/2.;
    
	vec props = zeros(21);
    props = {E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,eta14,eta15,eta16,eta24,eta25,eta26,eta34,eta35,eta36,eta45,eta46,eta56};
	
	return props;
}

} //end of namespace simcoon
