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
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

vec isochoric_invariants(const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside isochoric_invariants.");
        }            
    }
    mat b_bar = pow(J,-2./3.)*b;
    vec I = zeros(3);    
    I(0) = trace(b_bar);
    I(1) = 0.5*(pow(trace(b_bar),2.)-trace(powmat(b_bar,2)));
    I(2) = 1.;
    return I;
}

vec isochoric_invariants(const vec &lambda, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        J = prod(lambda);
    }
    vec lambda_bar = pow(J,-1./3.)*lambda;
    vec I = zeros(3);    
    I(0) = pow(lambda_bar(0),2.) + pow(lambda_bar(1),2.) + pow(lambda_bar(2),2.);
    I(1) = pow(lambda_bar(0),-2.) + pow(lambda_bar(1),-2.) + pow(lambda_bar(2),-2.);
    I(2) = 1.;
    return I;
}

vec isochoric_pstretch_from_V(const mat &V, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = det(V);
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside isochoric_pstretch_from_V.");
        }    
    }
    vec lambda;
    try {
        lambda = eig_sym(V);
    } catch (const std::runtime_error &e) {
        cerr << "Error in eig_sym: " << e.what() << endl;
        throw simcoon::exception_eig_sym("Failed to compute eigenvalues in isochoric_pstretch_from_V.");
    }
    vec lambda_bar = pow(J,-1./3.)*lambda;
    return lambda_bar;    
}

vec isochoric_pstretch_from_b(const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside isochoric_pstretch_from_b.");
        }    
    }
    vec lambda;
    try {
        lambda = eig_sym(b);
    } catch (const std::runtime_error &e) {
        cerr << "Error in eig_sym: " << e.what() << endl;
        throw simcoon::exception_eig_sym("Failed to compute eigenvalues in isochoric_pstretch_from_b.");
    }
    lambda.transform( [](double val) { return (sqrt(val)); } );
    vec lambda_bar = pow(J,-1./3.)*lambda;
    return lambda_bar;
}

vec isochoric_pstretch(const mat &input, const string &input_tensor, const double &mJ) {
    if (input_tensor == "b") {
        return isochoric_pstretch_from_b(input, mJ);
    }
    if (input_tensor == "v" || input_tensor == "V") {
        return isochoric_pstretch_from_V(input, mJ);
    }
    throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
}

void pstretch(vec &lambda, mat &n_pvectors, const mat &input, const string &input_tensor, const double &mJ) {

    double J=mJ;
    if (input_tensor == "b") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside pstretch.");
            }    
        }
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside pstretch.");
        }
        lambda.transform( [](double val) { return (sqrt(val)); } );
    }
    else if (input_tensor == "v" || input_tensor == "V") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside pstretch.");
            }   
        }
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);        
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside pstretch.");
        }
    }
    else {
        throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
    }
}

void pstretch(vec &lambda, mat &n_pvectors, std::vector<mat> &N_projectors, const mat &input, const string &input_tensor, const double &mJ) {

    pstretch(lambda, n_pvectors, input, input_tensor, mJ);
    N_projectors[0] = n_pvectors.col(0)*(n_pvectors.col(0)).t();
    N_projectors[1] = n_pvectors.col(1)*(n_pvectors.col(1)).t();
    N_projectors[2] = n_pvectors.col(2)*(n_pvectors.col(2)).t();        
}

void isochoric_pstretch(vec &lambda_bar, mat &n_pvectors, const mat &input, const string &input_tensor, const double &mJ) {

    double J=mJ;
    vec lambda = zeros(3);

    if (input_tensor == "b") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside isochoric_pstretch.");
            }   
        }        
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside isochoric_pstretch.");
        }        
        lambda.transform( [](double val) { return (sqrt(val)); } );
    }
    else if (input_tensor == "v" || input_tensor == "V") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside isochoric_pstretch.");
            }   
        }
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside isochoric_pstretch.");
        }        
    }
    else {
        throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
    }
    lambda_bar = pow(J,-1./3.)*lambda;
}

void isochoric_pstretch(vec &lambda_bar, mat &n_pvectors, std::vector<mat> &N_projectors, const mat &input, const string &input_tensor, const double &mJ) {

    isochoric_pstretch(lambda_bar, n_pvectors, input, input_tensor, mJ);
    N_projectors[0] = n_pvectors.col(0)*(n_pvectors.col(0)).t();
    N_projectors[1] = n_pvectors.col(1)*(n_pvectors.col(1)).t();
    N_projectors[2] = n_pvectors.col(2)*(n_pvectors.col(2)).t();        
}

vec beta_coefs(const vec &dWdlambda_bar, const vec &lambda_bar) {

    vec beta = lambda_bar%dWdlambda_bar - (1./3.)*sum(lambda_bar%dWdlambda_bar)*ones(3);
    return beta;
}

mat gamma_coefs(const vec &dWdlambda_bar, const mat &dW2dlambda_bar2, const vec &lambda_bar) {

    mat gamma = zeros(3,3);
    mat factor = dW2dlambda_bar2%(lambda_bar*lambda_bar.t()) + (dWdlambda_bar*lambda_bar.t())%ones(3,3);
    vec factor_sum_col = (-1./3.)*sum(factor,0);
    vec factor_sum_row = (-1./3.)*sum(factor,1);

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
                gamma(i,j) = factor_sum_col(i) + factor_sum_row(j);
        }            
    }    

    gamma += factor;
    gamma += (1./9.)*accu(factor)*ones(3,3);

    return gamma;
}

mat tau_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &b, const double &mJ) {
    vec lambda_bar = zeros(3);
    vec n_pvectors = zeros(3);
    std::vector<mat> N_projectors(3);
    isochoric_pstretch(lambda_bar, n_pvectors, N_projectors, b, "b", mJ);
    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);

    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return tau_iso;
}

mat tau_iso_hyper_pstretch(const vec &dWdlambda_bar, const vec &lambda_bar, const std::vector<mat> &N_projectors) {

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return tau_iso;
}

mat sigma_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside sigma_iso_hyper_pstretch.");
        }   
    }  
    vec lambda_bar = zeros(3);
    vec n_pvectors = zeros(3);
    std::vector<mat> N_projectors(3);
    isochoric_pstretch(lambda_bar, n_pvectors, N_projectors, b, "b", mJ);
    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);

    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return (1./J)*tau_iso;
}

mat sigma_iso_hyper_pstretch(const vec &dWdlambda_bar, const vec &lambda_bar, const std::vector<mat> &N_projectors, const double &J) {

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return (1./J)*tau_iso;
}

/*vec a_coefs(const double &dWdI_1_bar, const double &dWdI_2_bar, const vec &I_bar) {
    vec a = zeros(2);
    a(0) = 2.*(dWdI_1_bar + I_bar(0)*dWdI_2_bar);
    a(1) = 2.*(dWdI_2_bar);
    return a;
}

vec b_coefs(const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const vec &I_bar) {
    vec b = zeros(4);
    b(0) = 4.*(dW2dI_11_bar + dWdI_2_bar + pow(I_bar(0),2.)*dW2dI_22_bar + 2*I_bar(0)*dW2dI_12_bar);
    b(1) = 4.*(I_bar(0)*dW2dI_22_bar+dW2dI_12_bar);
    b(2) = 4.*dW2dI_22_bar;
    b(3) = 4.*dWdI_2_bar;
    return b;
}

vec delta_coefs(const vec &a_coefs, const vec &b_coefs, const mat &b) {

    vec I = Inv_X(b);    
    double trb2 = trace(powmat(b,2));
    vec delta = zeros(8);
    delta(0) = b_coefs(0);
    delta(1) = -b_coefs(1);
    delta(2) = (1./3.)*(2*a_coefs(0)-b_coefs(0)*I(0)+b_coefs(1)*trb2);
    delta(3) = b_coefs(2);
    delta(4) = 2.*a_coefs(1)+b_coefs(1)*I(0)-b_coefs(2)*trb2+b_coefs(3);
    delta(5) = (1./9.)*(2.*a_coefs(0)*I(0)-2.*a_coefs(1)*trb2+b_coefs(0)*pow(I(0),2.)) + 
               (1./9.)*(-2.*b_coefs(1)*I(0)*trb2+b_coefs(2)*pow(trb2,2.)-b_coefs(3)*trb2);
    delta(6) = (2./3.)*(a_coefs(0)*I(0)-a_coefs(1)*trb2);
    delta(7) = -b_coefs(3);
    return delta;
}*/

mat tau_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside tau_iso_hyper_invariants.");
        }   
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    return 2.*dWdI_1_bar*dev(b_bar) + 2.*dWdI_2_bar*(trace(b_bar)*dev(b_bar) - dev(b_bar2));
}

mat tau_vol_hyper(const double &dUdJ, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside tau_vol_hyper.");
        }   
    }    
    mat Id = eye(3,3);
    return J*dUdJ*eye(3,3);
}

mat sigma_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside sigma_iso_hyper_invariants.");
        }   
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    return (1./J)*(2.*dWdI_1_bar*dev(b_bar) + 2.*dWdI_2_bar*(trace(b_bar)*dev(b_bar) - dev(b_bar2)));
}

mat sigma_vol_hyper(const double &dUdJ, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside sigma_vol_hyper.");
        }   
    }    
    mat Id = eye(3,3);
    return dUdJ*eye(3,3);
}

/*mat L_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        J = sqrt(det(b));
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    vec I_bar = isochoric_invariants(b,J);
    mat Id = eye(3,3);

    mat H_1 = ;
    mat H_2 = ;

    mat gamma_1 = 2.*H_1 - (4./3.)*(I_bar(0)*Idev() - (sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar)));
    mat gamma_2 = (8./3.)*(I_bar(0)*(Ireal() - 2.*Ivol()) - I_bar(0)*(sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar))
                    + (sym_dyadic(dev_b_bar2,Id)+sym_dyadic(Id,dev_b_bar2))) + 4*(auto_sym_dyadic(dev_b_bar)-H_bar);
    mat gamma_11 = 4.*auto_sym_dyadic(dev_b_bar);
    mat gamma_22 = 4.*auto_sym_dyadic(devdevbb2);
    mat gamma_12 = 4.*(sym_dyadic(dev_b_bar,devdevbb2)+sym_dyadic(devdevbb2,dev_b_bar));

    mat L_iso = zeros (6,6);
    L_iso = (1./J)*(gamma_1*dWdI_1_bar+gamma_2*dWdI_2_bar+gamma_11*dW2dI_11_bar+gamma_12*dW2dI_12_bar+gamma_22*dW2dI_22_bar);

    return L_iso;
}*/

/*mat L_iso_hyper_invariants(const vec &delta_coefs, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        J = sqrt(det(b));
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    mat Id = eye(3,3);

    cout << "delta_coefs(0)*auto_sym_dyadic(b_bar)" << delta_coefs(0)*auto_sym_dyadic(b_bar) << endl;
    cout << "delta_coefs(1)*(sym_dyadic(b_bar,b_bar2)+sym_dyadic(b_bar2,b_bar))" << delta_coefs(1)*(sym_dyadic(b_bar,b_bar2)+sym_dyadic(b_bar2,b_bar)) << endl;
    cout << "delta_coefs(2)*(sym_dyadic(b_bar,Id) + sym_dyadic(Id,b_bar))" << delta_coefs(2)*(sym_dyadic(b_bar,Id) + sym_dyadic(Id,b_bar)) << endl;
    cout << "delta_coefs(3)*auto_sym_dyadic(b_bar2)" << delta_coefs(3)*auto_sym_dyadic(b_bar2) << endl;
    cout << "delta_coefs(4)*(sym_dyadic(b_bar2,Id) + sym_dyadic(Id,b_bar2))" << delta_coefs(4)*(sym_dyadic(b_bar2,Id) + sym_dyadic(Id,b_bar2)) << endl;
    cout << "delta_coefs(5)*auto_sym_dyadic(Id)" << delta_coefs(5)*auto_sym_dyadic(Id) << endl;    
    cout << "delta_coefs(6)*auto_sym_dyadic_operator(Id)" << delta_coefs(6)*auto_sym_dyadic_operator(Id) << endl;
    cout << "delta_coefs(7)*auto_sym_dyadic_operator(b_bar)" << delta_coefs(7)*auto_sym_dyadic_operator(b_bar) << endl;                

    mat L_iso = zeros (6,6);
    L_iso = delta_coefs(0)*auto_sym_dyadic(b_bar)+delta_coefs(1)*(sym_dyadic(b_bar,b_bar2)+sym_dyadic(b_bar2,b_bar))+ delta_coefs(2)*(sym_dyadic(b_bar,Id) + sym_dyadic(Id,b_bar)) + delta_coefs(3)*auto_sym_dyadic(b_bar2)
                + delta_coefs(4)*(sym_dyadic(b_bar2,Id) + sym_dyadic(Id,b_bar2)) + delta_coefs(5)*auto_sym_dyadic(Id)
                + delta_coefs(6)*auto_sym_dyadic_operator(Id) + delta_coefs(7)*auto_sym_dyadic_operator(b_bar);
    return L_iso;
}*/

mat L_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &dW2dlambda_bar2, const vec &lambda_bar, const mat &n_pvectors, const double &J) {

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat gamma = gamma_coefs(dWdlambda_bar, dW2dlambda_bar2, lambda_bar);
    mat c = zeros(6,6);
    mat delta = eye(3,3);

    double factor = 0.;

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            c += (gamma(i,j) - 2*delta(i,j)*beta(i))*dyadic_4vectors_sym(n_pvectors.col(i), n_pvectors.col(j), "aabb") ;
        }            

        for(unsigned int j=i; j<3; j++) {
            if(pow(lambda_bar(j),2.)-pow(lambda_bar(i),2.) < 1.E-6) {
                factor = pow(lambda_bar(i), 2.) * pow(lambda_bar(j), -2.) * (0.5*gamma(j,j) - beta(j)) - 0.5*gamma(i,j);
            }
            else {
                factor = (beta(j)*pow(lambda_bar(i),2.) - beta(i)*pow(lambda_bar(j),2.))/(pow(lambda_bar(j),2.)-pow(lambda_bar(i),2.));
            }

            c +=  factor*dyadic_4vectors_sym(n_pvectors.col(i), n_pvectors.col(j), "abab"); 
        }            
    }   
    return (1./J)*c;    
}

mat L_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &dW2dlambda_bar2, const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside L_iso_hyper_pstretch.");
        }   
    }    

    vec lambda;
    try {
        lambda = eig_sym(b);
    } catch (const std::runtime_error &e) {
        cerr << "Error in eig_sym: " << e.what() << endl;
        throw simcoon::exception_eig_sym("Failed to compute eigenvalues in L_iso_hyper_pstretch.");
    }    

    vec lambda_bar = zeros(3);
    vec n_pvectors = zeros(3);
    isochoric_pstretch(lambda_bar, n_pvectors, b, "b", J);

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat gamma = gamma_coefs(dWdlambda_bar, dW2dlambda_bar2, lambda_bar);
    mat c = zeros(6,6);
    mat delta = eye(3,3);

    double factor = 0.;

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            c += (gamma(i,j) - 2*delta(i,j)*beta(i))*dyadic_4vectors_sym(n_pvectors.col(i), n_pvectors.col(j), "aabb") ;
        }            

        for(unsigned int j=i; j<3; j++) {
            if(pow(lambda(j),2.)-pow(lambda(i),2.) < 1.E-6) {
                factor = pow(lambda(i), 2.) * pow(lambda(j), -2.) * (0.5*gamma(j,j) - beta(j)) - 0.5*gamma(i,j);
            }
            else {
                factor = (beta(j)*pow(lambda(i),2.) - beta(i)*pow(lambda(j),2.))/(pow(lambda(j),2.)-pow(lambda(i),2.));
            }

            c +=  factor*dyadic_4vectors_sym(n_pvectors.col(i), n_pvectors.col(j), "abab"); 
        }            
    }   
    return (1./J)*c;

}

mat L_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside L_iso_hyper_invariants.");
        }   
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    mat dev_b_bar = dev(b_bar);    
    mat dev_b_bar2 = dev(b_bar2);        
    vec I_bar = isochoric_invariants(b,J);
    mat Id = eye(3,3);
    mat H_bar = auto_sym_dyadic_operator(b_bar);

    mat devdevbb2 = I_bar(0)*dev(b_bar) - dev(b_bar2);

    mat gamma_1 = (4./3.)*(I_bar(0)*Idev() - (sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar)));
    mat gamma_2 = (8./3.)*(I_bar(1)*(Ireal() - 2.*Ivol()) - I_bar(0)*(sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar))
                    + (sym_dyadic(dev_b_bar2,Id)+sym_dyadic(Id,dev_b_bar2))) + 4*(auto_sym_dyadic(b_bar)-H_bar);
    mat gamma_11 = 4.*auto_sym_dyadic(dev_b_bar);
    mat gamma_22 = 4.*auto_sym_dyadic(devdevbb2);
    mat gamma_12 = 4.*(sym_dyadic(dev_b_bar,devdevbb2)+sym_dyadic(devdevbb2,dev_b_bar));

    return (1./J)*(gamma_1*dWdI_1_bar+gamma_2*dWdI_2_bar+gamma_11*dW2dI_11_bar+gamma_12*dW2dI_12_bar+gamma_22*dW2dI_22_bar);
}

mat L_vol_hyper(const double &dUdJ, const double &dU2dJ2, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside L_vol_hyper.");
        } 
    }
    return (dUdJ+dU2dJ2*J)*3.*Ivol() - 2.*dUdJ*Ireal();
}

} //namespace simcoon
