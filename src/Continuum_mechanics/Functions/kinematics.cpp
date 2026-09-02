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
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//This function returns F (in a vectorized), from E (Green-Lagrange strain) and R (Rotation matrix), according to a RU decomposition
mat ER_to_F(const mat &E, const mat &R) {

    assert(E.n_cols == 3);
    assert(E.n_rows == 3);
    assert(R.n_cols == 3);
    assert(R.n_rows == 3);

    //From E we compute C : E = 1/2 (C-I) --> C = U^2 = 2E+I
    //symmatu() removes any numerical asymmetry so a mathematically-SPD C is not
    //rejected by sqrtmat_sympd's symmetry check (the input E is symmetric by construction).
    mat C = symmatu(2.*E+eye(3,3));

    try {
        return (R*sqrtmat_sympd(C));
    } catch (const std::runtime_error &e) {
        cerr << "Error in sqrtmat_sympd : " << e.what() << endl;
        throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside ER_to_F.");
    }            
}

mat eR_to_F(const mat &e, const mat &R) {

    assert(e.n_cols == 3);
    assert(e.n_rows == 3);
    assert(R.n_cols == 3);
    assert(R.n_rows == 3);

    //From e we compute V : e = 1/2 (C-I) --> C = U^2 = 2E+I
    mat V;
    try {
        V = expmat_sym(e);
    } catch (const std::runtime_error &e) {
        cerr << "Error in expmat_sym : " << e.what() << endl;
        throw simcoon::exception_expmat_sym("Error in expmat_sym function inside eR_to_F.");
    }            

    return (V*R);
}
    
//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
mat G_UdX(const mat &F) {
    return F - eye(3,3);
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
mat G_Udx(const mat &F) {

    try {
        return eye(3,3) - inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv : " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside G_Udx.");
    }        
}
    
//This function computes the Right Cauchy-Green C
mat R_Cauchy_Green(const mat &F) {
    return F.t()*F;
}
    
//This function computes the Left Cauchy-Green B
mat L_Cauchy_Green(const mat &F) {
    return F*F.t();
}

//Provides the RU decomposition of the transformation gradient F
void RU_decomposition(mat &R, mat &U, const mat &F) {
    mat U2 = F.t()*F;

    try {
        U = sqrtmat_sympd(U2);
    } catch (const std::runtime_error &e) {
        cerr << "Error in sqrtmat_sympd : " << e.what() << endl;
        throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside RU_decomposition.");
    }                

    try {
        R = F*inv(U);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside RU_decomposition.");
    }     
}

//Provides the VR decomposition of the transformation gradient F
void VR_decomposition(mat &V, mat &R, const mat &F) {
    mat V2 = F*F.t();

    try {
        V = sqrtmat_sympd(V2);
    } catch (const std::runtime_error &e) {
        cerr << "Error in sqrtmat_sympd : " << e.what() << endl;
        throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside VR_decomposition.");
    }                

    try {
        R = inv(V)*F;
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside RU_decomposition.");
    }     
}
    
//This function computes the common Right (or Left) Cauchy-Green invariants
vec Inv_X(const mat &X) {
    
    vec I = zeros(3);
    I(0) = trace(X);
    I(1) = 0.5*(pow(trace(X),2.) - trace(X*X));
    I(2) = det(X);    
    return I;
}

//This function computes the Cauchy deformation tensor b
mat Cauchy(const mat &F) {
    try {
        return inv(L_Cauchy_Green(F));
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Cauchy.");
    }     
}
    
//This function computes the Green-Lagrange finite strain tensor E
mat Green_Lagrange(const mat &F) {
    return 0.5*(R_Cauchy_Green(F) - eye(3,3));
}

//This function computes the Euler-Almansi finite strain tensor e
mat Euler_Almansi(const mat &F) {
    return 0.5*(eye(3,3) - Cauchy(F));
}

//This function computes the Euler-Almansi finite strain tensor h
mat Log_strain(const mat &F) {
    try {
        return  0.5*logmat_sympd(L_Cauchy_Green(F));
    } catch (const std::runtime_error &e) {
        cerr << "Error in logmat_sympd: " << e.what() << endl;
        throw simcoon::exception_logmat_sympd("Error in logmat_sympd function inside Log_strain.");
    }     
}

//This function computes the velocity difference
mat finite_L(const mat &F0, const mat &F1, const double &DTime) {
    
    //Definition of L = dot(F)*F^-1
    try {
        return (1./DTime)*(F1-F0)*inv(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_L.");
    }   
}

//This function computes the deformation rate D
mat finite_D(const mat &F0, const mat &F1, const double &DTime) {
    
    //Definition of L = dot(F)*F^-1
    mat L;
    try {
        L = (1./DTime)*(F1-F0)*inv(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_D.");
    }   
    //Definition of the deformation rate D
    return 0.5*(L+L.t());
    
}

//This function computes the spin tensor W (correspond to Jaumann rate)
mat finite_W(const mat &F0, const mat &F1, const double &DTime) {

    //Definition of L = dot(F)*F^-1
    mat L;
    try {
        L = (1./DTime)*(F1-F0)*inv(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_W.");
    }   
    
    //Definition of the rotation matrix Q
    return 0.5*(L-L.t());

}
    
//This function computes the exact relative polar rotation and its midpoint spin
// Note : here R is the rigid body rotation in the polar decomposition of the deformation gradient F
void finite_rotation(const mat &F0, const mat &F1, const double &DTime, mat &DR, mat &Omega) {
    mat R0, U0, R1, U1;
    RU_decomposition(R0, U0, F0);
    RU_decomposition(R1, U1, F1);
    // Exact and DTime-independent: DR stays the true relative rotation even
    // when the time increment degenerates (the spin, a rate, does not).
    DR = R1*R0.t();

    if (DTime <= simcoon::iota) {
        Omega = zeros(3,3);
        return;
    }
    try {
        // Inverse Cayley: exactly skew for an orthogonal DR, and its
        // Hughes-Winget (Cayley) update recovers DR exactly. Singular at a
        // 180-degree increment (Cayley chart limit) and for improper inputs
        // (det F < 0) -> exception_inv, the solver's step-cut signal.
        Omega = (2./DTime)*(DR-eye(3,3))*inv(DR+eye(3,3));
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_rotation.");
    }
}

//This function computes the spin tensor Omega (correspond to Green-Naghdi rate)
mat finite_Omega(const mat &F0, const mat &F1, const double &DTime) {
    mat DR, Omega;
    finite_rotation(F0, F1, DTime, DR, Omega);
    return Omega;
}

//This function computes the increment of finite rotation
mat finite_DQ(const mat &Omega0, const mat &Omega1, const double &DTime) {
    
    try {
        return (eye(3,3)+0.5*DTime*Omega0)*(inv(eye(3,3)-0.5*DTime*Omega1));
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_DQ.");
    }   
    
}
    
    
} //namespace simcoon
