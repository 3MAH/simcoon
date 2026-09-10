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
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <stdexcept>

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
    try {
        I(1) = 0.5*(pow(trace(b_bar),2.)-trace(powmat(b_bar,2)));
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_powmat("Error in powmat function inside isochoric_invariants.");
    }        
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
    // g_ab = lambda_a lambda_b W_ab + delta_ab lambda_a W_a (the W_a term is DIAGONAL,
    // not an outer product): gamma = P.g.P with P the deviatoric projector, i.e.
    // gamma_ab = d beta_a / d ln(lambda_b)
    mat factor = dW2dlambda_bar2%(lambda_bar*lambda_bar.t()) + diagmat(dWdlambda_bar%lambda_bar);
    vec factor_sum_col = ((-1./3.)*sum(factor,0)).t();  // sum along dim 0 returns row vec, transpose to col
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
    mat n_pvectors = zeros(3,3);
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
    mat n_pvectors = zeros(3,3);
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

    mat b_bar2;
    try {
        b_bar2 = powmat(b_bar,2);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_powmat("Error in powmat function inside tau_iso_hyper_invariants.");
    }        

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
    mat b_bar2;
    try {
        b_bar2 = powmat(b_bar,2);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_powmat("Error in powmat function inside sigma_iso_hyper_invariants.");
    }            

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

        // shear term, summed over ordered pairs a!=b: the full 1/2-symmetrized dyadic
        // 2*(abab(n_i,n_j) + abab(n_j,n_i)) collapses to the rank-1 outer product s*s^T
        // with the symmetric dyadic of (n_i n_j^T + n_j n_i^T)
        for(unsigned int j=i+1; j<3; j++) {
            if(pow(lambda_bar(j),2.)-pow(lambda_bar(i),2.) < 1.E-6) {
                factor = pow(lambda_bar(i), 2.) * pow(lambda_bar(j), -2.) * (0.5*gamma(j,j) - beta(j)) - 0.5*gamma(i,j);
            }
            else {
                factor = (beta(j)*pow(lambda_bar(i),2.) - beta(i)*pow(lambda_bar(j),2.))/(pow(lambda_bar(j),2.)-pow(lambda_bar(i),2.));
            }

            c += factor*auto_sym_dyadic(n_pvectors.col(i)*(n_pvectors.col(j)).t() + n_pvectors.col(j)*(n_pvectors.col(i)).t());
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

    vec lambda_bar = zeros(3);
    mat n_pvectors = zeros(3,3);
    isochoric_pstretch(lambda_bar, n_pvectors, b, "b", J);
    // the shear coefficient is invariant under the J^(-1/3) bar-scaling, so the
    // lambda_bar-based assembly is exact for b as well
    return L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, lambda_bar, n_pvectors, J);
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

    mat b_bar2;
    try {
        b_bar2 = powmat(b_bar,2);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_powmat("Error in powmat function inside L_iso_hyper_invariants.");
    }            

    mat dev_b_bar = dev(b_bar);    
    mat dev_b_bar2 = dev(b_bar2);        
    vec I_bar = isochoric_invariants(b,J);
    mat Id = eye(3,3);
    mat H_bar = auto_sym_dyadic_operator(b_bar);

    mat devdevbb2 = I_bar(0)*dev(b_bar) - dev(b_bar2);

    // Plain locals on purpose: a function-local static of an armadillo type registers a
    // destructor that runs at DLL unload, and on Windows the unload order of the extension
    // module, the BLAS/LAPACK DLLs and the CRT is not defined — freeing the matrix after its
    // allocator is gone aborts the process once the tests are over. These are 6x6 identities,
    // negligible next to the six dyadic products this function already builds per call.
    mat I_real = Ireal();
    mat I_vol = Ivol();
    mat I_dev = Idev();

    mat gamma_1 = (4./3.)*(I_bar(0)*I_dev - (sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar)));
    mat gamma_2 = (8./3.)*(I_bar(1)*(I_real - 2.*I_vol) - I_bar(0)*(sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar))
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
    mat I_real = Ireal();                                // never static: see L_iso_hyper above
    mat I_vol = Ivol();
    return (dUdJ+dU2dJ2*J)*3.*I_vol - 2.*dUdJ*I_real;
}


namespace {

// props(i) is unchecked in release builds: a short props vector would be read
// out of bounds silently.
void require_props(const vec &props, const uword n, const char *name) {
    if (props.n_elem < n) {
        throw std::invalid_argument(std::string("hyper_potential_derivatives: ") + name + " needs "
                                    + std::to_string(n) + " parameters, got "
                                    + std::to_string(props.n_elem));
    }
}

}  // namespace

hyper_invariants_dW hyper_potential_derivatives(const HyperPotential &potential, const vec &props, const vec &I_bar, const double &J) {

    hyper_invariants_dW dW;
    double kappa = 0.;  // every potential shares U(J) = kappa (J ln J - J + 1)

    switch (potential) {
        case HyperPotential::NEOHC: {
            // \f$ W = \frac{\mu}{2}*\left(\bar{I}_1 -3 \right) + \kappa \left( J \]textrm{ln} J - J +1 \right) \f$
            require_props(props, 2, "NEOHC");
            double mu = props(0);
            kappa = props(1);
            dW.dWdI_1_bar = 0.5*mu;
            break;
        }
        case HyperPotential::MOORI: {
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{01} left(\bar{I}_2 -3\right) + \kappa \left( J textrm{ln} J - J +1 \right) \f$
            require_props(props, 3, "MOORI");
            double C_10 = props(0);
            double C_01 = props(1);
            kappa = props(2);
            dW.dWdI_1_bar = C_10;
            dW.dWdI_2_bar = C_01;
            break;
        }
        case HyperPotential::YEOHH: {
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{20} left(\bar{I}_1 -3\right)^2 + C_{30} left(\bar{I}_1 -3\right)^3 + \kappa \left( J textrm{ln} J - J +1 \right) \f$
            require_props(props, 4, "YEOHH");
            double C_10 = props(0);
            double C_20 = props(1);
            double C_30 = props(2);
            kappa = props(3);
            dW.dWdI_1_bar = C_10 + 2.*C_20*(I_bar(0)-3.) + 3.*C_30*pow((I_bar(0)-3.),2.);
            dW.dW2dI_11_bar = 2.*C_20 + 6.*C_30*(I_bar(0)-3.);
            break;
        }
        case HyperPotential::ISHAH: {
            // Isihara model (1951)
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{20} left(\bar{I}_1 -3\right)^2 + C_{01} left(\bar{I}_2 -3\right) + \kappa \left( J textrm{ln} J - J +1 \right) \f$
            require_props(props, 4, "ISHAH");
            double C_10 = props(0);
            double C_20 = props(1);
            double C_01 = props(2);
            kappa = props(3);
            dW.dWdI_1_bar = C_10 + 2.*C_20*(I_bar(0)-3.);
            dW.dW2dI_11_bar = 2.*C_20;
            dW.dWdI_2_bar = C_01;
            break;
        }
        case HyperPotential::GETHH: {
            // Gent-Thomas model (1958)
            // \f$ W = c_1 left(\bar{I}_1 -3\right) + c_2 \textrm{ln} left( \frac{\bar{I}_2}{3}\right) + \kappa \left( J textrm{ln} J - J +1 \right) \f$
            require_props(props, 3, "GETHH");
            double c_1 = props(0);
            double c_2 = props(1);
            kappa = props(2);
            dW.dWdI_1_bar = c_1;
            if(fabs(I_bar(1)) > simcoon::iota) {
                dW.dWdI_2_bar = c_2/I_bar(1);
                dW.dW2dI_22_bar = -1.*c_2/pow(I_bar(1),2.);
            }
            break;
        }
        case HyperPotential::SWANH: {
            // Swanson model (1985)
            // \f$ W = \frac{3}{2} \sum_{i=1}^n \frac{A_i}{1+\alpha_i} left(\frac{\bar{I}_1}{3}\right)^{1+\alpha_i} + \frac{3}{2} \sum_{i=1}^n \frac{B_i}{1+\beta_i} left(\frac{\bar{I}_2}{3}\right)^{1+\beta_i} + \kappa \left( J textrm{ln} J - J +1 \right) \f$
            require_props(props, 2, "SWANH");
            int N_Swanson = int(props(0));
            require_props(props, 2 + 4*std::max(N_Swanson, 0), "SWANH");
            kappa = props(1);
            for (int i=0; i<N_Swanson; i++) {
                const double A = props(2+i*4);
                const double B = props(2+i*4+1);
                const double alpha = props(2+i*4+2);
                const double beta = props(2+i*4+3);
                dW.dWdI_1_bar += 1./2.*A*pow((I_bar(0)/3.),alpha);
                dW.dW2dI_11_bar += A*alpha/6.*pow((I_bar(0)/3.),alpha-1.);
                dW.dWdI_2_bar += 1./2.*B*pow((I_bar(1)/3.),beta);
                dW.dW2dI_22_bar += B*beta/6.*pow((I_bar(1)/3.),beta-1.);
            }
            break;
        }
        default:
            throw std::invalid_argument("hyper_potential_derivatives: unknown potential "
                                        + std::to_string(static_cast<int>(potential)));
    }

    dW.dUdJ = kappa*log(J);
    dW.dU2dJ2 = kappa/J;
    return dW;
}

void hyper_invariants_response(const hyper_invariants_dW &dW, const mat &b, const double &J, const mat &F, vec &sigma, mat &Lt_box) {

    mat m_sigma_iso = sigma_iso_hyper_invariants(dW.dWdI_1_bar, dW.dWdI_2_bar, b, J);
    mat m_sigma_vol = sigma_vol_hyper(dW.dUdJ, b, J);
    mat m_sigma = m_sigma_iso + m_sigma_vol;
    sigma = t2v_stress(m_sigma);

    mat Lt_iso = L_iso_hyper_invariants(dW.dWdI_1_bar, dW.dWdI_2_bar, dW.dW2dI_11_bar, dW.dW2dI_12_bar, dW.dW2dI_22_bar, b, J);
    mat Lt_vol = L_vol_hyper(dW.dUdJ, dW.dU2dJ2, b, J);
    mat Lt_spatial = Lt_iso + Lt_vol;   // native hyperelastic tangent = Cauchy (Oldroyd/Lie) spatial elasticity, dsigma/dD

    // Standardize to the canonical box convention Lt = d(tau_hat)/d(De) (Kirchhoff, no-J,
    // XBM rate) -- identical object to the small-strain boxes and saint_venant.
    Lt_box = box_DtauDe_from_spatial(Lt_spatial, F, sigma);
}

} //namespace simcoon
