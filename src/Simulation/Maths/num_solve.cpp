/* This file is part of simcoon private.
 
 Only part of simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file num_solve.cpp
///@brief random number generators
///@author Chemisky & Anagnostou
///@version 1.0
///@date 10/23/2014

#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>

using namespace std;
using namespace arma;

/* Notice*/
// Recall that these methods are intended to solve an problem by means of the Kuhn-Tucker conditions
//This means : Phi <= 0 ; Dp >= 0; Phi*Dp = 0, with Dp = pdot * Dt, a non-zero time state

namespace simcoon{

void Newton_Raphon(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    int n=Phi.n_elem;
    
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(i)) > 0.);
    }
    
    try {
        if (fabs(det(denom)) > simcoon::limit) {
            dp = -1.*solve(denom,Phi);
        }
        else {
            dp = zeros(n);
        }
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Newton_Raphon.");
    }   
    
    //New update of the vector of the variable increment throughout the step
    Dp = Dp + dp;
    
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(Phi(i))/fabs(Y_crit(i));
    }
    
}
    
void Fischer_Burmeister(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    
    int n=Phi.n_elem;
    
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(n)) > 0);
    }
    
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Normalized Fischer-Burmeister set of equations
    for (int i=0; i<n; i++) {
        
        if ((fabs(Phi(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
            FB(i) = sqrt(pow(Phi(i),2.) + pow(Dp(i),2.)) + Phi(i) - Dp(i);
        }
        else if(fabs(Phi(i)) > 0.) {
            FB(i) = sqrt(pow(Phi(i),2.)) + Phi(i);
        }
        else if(fabs(Dp(i)) > 0.) {
            FB(i) = sqrt(pow(Dp(i),2.)) - Dp(i);
        }
        else {
            FB(i) = 0.;
        }
        
        for (int j=0; j<n; j++) {
            if ((fabs(Phi(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.) + pow(Dp(i),2.)))+1.)*denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Phi(i),2.) + pow(Dp(i),2.))) - 1.);
            }
            else if(fabs(Phi(i)) > 0.) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.)))+1.)*denom(i,j) - delta(i,j);
            }
            else if(fabs(Dp(i)) > 0.) {
                denomFB(i,j) = denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Dp(i),2.))) - 1.);
            }
            else {
                denomFB(i,j) = 1.E12;
            }
        }
        
    }
    
    try {
        if (fabs(det(denomFB)) > simcoon::limit) {
            dp = -1.*solve(denomFB,FB);
        }
        else {
            dp = zeros(n);
        }
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Fischer_Burmeister.");
    }   
    
    //New update of the vector of the variable increment throughout the step
    Dp = Dp + dp;
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
}

void Fischer_Burmeister_limits(const vec &Phi_p, const vec &Phi_l, const vec &Y_crit, const mat &denom, const mat &denom_l, vec &Dp, vec &dp, double &error)
{
    int n=Phi_p.n_elem;
    
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(n)) > 0);
    }
    
    vec Phi_1 = zeros(n);
    vec Phi_2 = zeros(n);
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    mat dPhi1dp = zeros(n,n);
    mat dPhi2dp = zeros(n,n);
    
    // Fischer-Burmeister set of equations
    //Phi_l is for now written as : Phi_l = p - 1 <= 0
    
    for (int i=0; i<n; i++) {
        
        if (fabs(Phi_l(i)) > 0.) {
            
            if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                Phi_1(i) = sqrt(pow(Phi_p(i),2.) + pow(Dp(i),2.)) + Phi_p(i) - Dp(i);
            }
            else if(fabs(Phi_p(i)) > 0.) {
                Phi_1(i) = sqrt(pow(Phi_p(i),2.)) + Phi_p(i);
            }
            else if(fabs(Dp(i)) > 0.) {
                Phi_1(i) = sqrt(pow(Dp(i),2.)) - Dp(i);
            }
            else {
                Phi_1(i) = 0.;
            }
            Phi_2(i) = sqrt(pow(Phi_l(i),2.)) + Phi_l(i);
            FB(i) = pow(Phi_1(i),2.) + pow(Phi_2(i),2.);
         
            for (int j=0; j<n; j++) {
                if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                    dPhi1dp(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dp(i),2.)))+1.)*denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dp(i),2.))) - 1.);
                }
                else if(fabs(Phi_p(i)) > 0.) {
                    dPhi1dp(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.)))+1.)*denom(i,j) - delta(i,j);
                }
                else if(fabs(Dp(i)) > 0.) {
                    dPhi1dp(i,j) = denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Dp(i),2.))) - 1.);
                }
                else {
                    if (i==j)
                        dPhi1dp(i,j) = 1.E12;
                    else
                        dPhi1dp(i,j) = 0.;
                }
            
                if(fabs(Phi_p(i)) > 0.) {
                    dPhi2dp(i,j) = (Phi_l(i)/(sqrt(pow(Phi_l(i),2.)))+1.)*denom_l(i,j);
                }
                else {
                    if (i==j)
                        dPhi2dp(i,j) = (Phi_l(i)/(sqrt(pow(Phi_l(i),2.)))+1.)*denom_l(i,j);
                    else
                        dPhi2dp(i,j) = 0.;
                }
                
                denomFB(i,j) = 2.*dPhi1dp(i,j)*Phi_1(i) + 2.*dPhi2dp(i,j)*Phi_2(i);
            }
        }
        else {
            if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                FB(i) = sqrt(pow(Phi_p(i),2.) + pow(Dp(i),2.)) + Phi_p(i) - Dp(i);
            }
            else if(fabs(Phi_p(i)) > 0.) {
                FB(i) = sqrt(pow(Phi_p(i),2.)) + Phi_p(i);
            }
            else if(fabs(Dp(i)) > 0.) {
                FB(i) = sqrt(pow(Dp(i),2.)) - Dp(i);
            }
            else {
                FB(i) = 0.;
            }
            for (int j=0; j<n; j++) {
                if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                    denomFB(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dp(i),2.)))+1.)*denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dp(i),2.))) - 1.);
                }
                else if(fabs(Phi_p(i)) > 0.) {
                    denomFB(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.)))+1.)*denom(i,j) - delta(i,j);
                }
                else if(fabs(Dp(i)) > 0.) {
                    denomFB(i,j) = denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Dp(i),2.))) - 1.);
                }
                else {
                    if (i==j)
                        denomFB(i,j) = 1.E12;
                    else
                        denomFB(i,j) = 0.;
                }
            }
        }
    }
    
    try {
        if (fabs(det(denomFB)) > simcoon::limit) {
            dp = -1.*solve(denomFB,FB);
        }
        else {
            dp = zeros(n);
        }
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Fischer_Burmeister_limits.");
    }   

    //New update of the vector of the variable increment throughout the step
    Dp = Dp + dp;
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
}
    
double Fischer_Burmeister_residual(const vec &Phi, const vec &Y_crit, const mat &denom, const vec &Dp)
{
    // Same branches and expressions as the FB vector of Fischer_Burmeister_m below.
    double error = 0.;
    for (uword i = 0; i < Phi.n_elem; i++) {
        const double Dpstar = Dp(i) * fabs(denom(i, i));
        double FB;
        if ((fabs(Phi(i)) > 0.) && (fabs(Dpstar) > 0.)) {
            FB = sqrt(pow(Phi(i), 2.) + pow(Dpstar, 2.)) + Phi(i) - Dpstar;
        }
        else if (fabs(Phi(i)) > 0.) {
            FB = sqrt(pow(Phi(i), 2.)) + Phi(i);
        }
        else if (fabs(Dpstar) > 0.) {
            FB = sqrt(pow(Dpstar, 2.)) - Dpstar;
        }
        else {
            FB = 0.;
        }
        error += fabs(FB) / fabs(Y_crit(i));
    }
    return error;
}

void Fischer_Burmeister_m(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    int n=Phi.n_elem;
        
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(i)) > 0);
    }
    
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Determine the eigenvalues of denom
    vec factor_denom = zeros(n);
    for (int i=0; i<n; i++) {
        factor_denom(i) = fabs(denom(i,i));
    }
    
    vec Dpstar = Dp%factor_denom;
    vec dDpstar = factor_denom;
    
    //Normalized Fischer-Burmeister set of equations
    for (int i=0; i<n; i++) {
        
        if ((fabs(Phi(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
            FB(i) = sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.)) + Phi(i) - Dpstar(i);
        }
        else if(fabs(Phi(i)) > 0.) {
            FB(i) = sqrt(pow(Phi(i),2.)) + Phi(i);
        }
        else if(fabs(Dpstar(i)) > 0.) {
            FB(i) = sqrt(pow(Dpstar(i),2.)) - Dpstar(i);
        }
        else {
            FB(i) = 0.;
        }
        
        for (int j=0; j<n; j++) {
            if ((fabs(Phi(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.)))+1.)*denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.))) - 1.);
            }
            else if(fabs(Phi(i)) > 0.) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.)))+1.)*denom(i,j) - delta(i,j)*dDpstar(i);
            }
            else if(fabs(Dpstar(i)) > 0.) {
                denomFB(i,j) = denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Dpstar(i),2.))) - 1.);
            }
            else {
                if (i==j)
                    denomFB(i,j) = 1.E12;
                else
                    denomFB(i,j) = 0.;
            }
        }
        
    }
    
    // One or two mechanisms (the common cases): closed forms, no LAPACK call. A det + solve pair
    // on a 1x1 or 2x2 costs ~0.3 us, the bulk of a local iteration.
    if (n == 1) {
        const double a = denomFB(0,0);
        if (fabs(a) > simcoon::limit) {
            dp.set_size(1);
            dp(0) = -FB(0)/a;
        }
        else {
            dp = zeros(n);
        }
    }
    else if (n == 2) {
        const double det_FB = denomFB(0,0)*denomFB(1,1) - denomFB(0,1)*denomFB(1,0);
        if (fabs(det_FB) > simcoon::limit) {
            dp.set_size(2);
            dp(0) = -(denomFB(1,1)*FB(0) - denomFB(0,1)*FB(1))/det_FB;   // Cramer
            dp(1) = -(denomFB(0,0)*FB(1) - denomFB(1,0)*FB(0))/det_FB;
        }
        else {
            dp = zeros(n);
        }
    }
    else {
        try {
            if (fabs(det(denomFB)) > simcoon::limit) {
                dp = -1.*solve(denomFB, FB);
            }
            else {
                dp = zeros(n);
            }
        } catch (const std::runtime_error &e) {
            cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside Fischer_Burmeister_m.");
        }
    }
    
    //New update of the transformation/orientation multipliers
    Dp = Dp + dp;

    // the entering iterate's residual: Fischer_Burmeister_residual, kept in step
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
}
                     
NCPMethod ncp_method_of(const int &code, const string &where)
{
    if (code == 0) return NCPMethod::fischer_burmeister;
    if (code == 1) return NCPMethod::min_fb_fallback;
    throw std::invalid_argument(where + ": ncp must be 0 (Fischer-Burmeister) or 1 "
                                "(min-NCP with a Fischer-Burmeister fallback); got " + std::to_string(code));
}

void ncp_step::operator()(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    // Several surfaces stay on Fischer-Burmeister: the active set of coupled surfaces makes the
    // plain semi-smooth Newton cycle (SMA with three surfaces)
    if (method_ == NCPMethod::fischer_burmeister || fallback_ || Phi.n_elem != 1) {
        Fischer_Burmeister_m(Phi, Y_crit, denom, Dp, dp, error);
        return;
    }
    // min(Dp*, -Phi) = 0, Dp* = Dp |B|: an active row solves Phi + B dp = 0 (classic cutting
    // plane), an inactive one resets its multiplier
    const double B = denom(0, 0);
    const double f = fabs(B);
    const bool active = (Phi(0) > 0.) || (-Phi(0) < Dp(0)*f);
    const double r = active ? -Phi(0) : Dp(0)*f;
    const double J = active ? -B : ((f > 0.) ? f : 1.);
    const double err = 2.*fabs(r)/fabs(Y_crit(0));
    // a growing error: the semi-smooth Newton overshoots or cycles, Fischer-Burmeister from here
    if ((prev_error_ >= 0. && err > prev_error_) || !(fabs(J) > simcoon::iota)) {
        fallback_ = true;
        Fischer_Burmeister_m(Phi, Y_crit, denom, Dp, dp, error);
        return;
    }
    prev_error_ = err;
    dp.set_size(1);
    dp(0) = -r/J;
    Dp(0) += dp(0);
    error = err;
}

void Fischer_Burmeister_m_limits(const vec &Phi_p, const vec &Phi_l, const vec &Y_crit, const mat &denom, const mat &denom_l, vec &Dp, vec &dp, double &error)
{
    int n=Phi_p.n_elem;
    
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(i)) > 0);
    }
    
    vec Phi_1 = zeros(n);
    vec Phi_2 = zeros(n);
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Determine the eigenvalues of denom
    vec factor_denom = zeros(n);
    for (int i=0; i<n; i++) {
        factor_denom(i) = fabs(denom(i,i));
    }
    
    vec Dpstar = Dp%factor_denom;
    vec dDpstar = factor_denom;
    
    mat dPhi1dp = zeros(n,n);
    mat dPhi2dp = zeros(n,n);
    
    // Fischer-Burmeister set of equations
    //Phi_l is for now written as : Phi_l = p - 1 <= 0
    
    for (int i=0; i<n; i++) {
        
        if (fabs(Phi_l(i)) > 0.) {
            
            if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                Phi_1(i) = sqrt(pow(Phi_p(i),2.) + pow(Dpstar(i),2.)) + Phi_p(i) - Dpstar(i);
            }
            else if(fabs(Phi_p(i)) > 0.) {
                Phi_1(i) = sqrt(pow(Phi_p(i),2.)) + Phi_p(i);
            }
            else if(fabs(Dp(i)) > 0.) {
                Phi_1(i) = sqrt(pow(Dpstar(i),2.)) - Dpstar(i);
            }
            else {
                Phi_1(i) = 0.;
            }
            Phi_2(i) = sqrt(pow(Phi_l(i),2.)) + Phi_l(i);
            FB(i) = pow(Phi_1(i),2.) + pow(Phi_2(i),2.);
            
            for (int j=0; j<n; j++) {
                if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                    dPhi1dp(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dpstar(i),2.)))+1.)*denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dpstar(i),2.))) - 1.);
                }
                else if(fabs(Phi_p(i)) > 0.) {
                    dPhi1dp(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.)))+1.)*denom(i,j) - delta(i,j)*dDpstar(i);
                }
                else if(fabs(Dp(i)) > 0.) {
                    dPhi1dp(i,j) = denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Dpstar(i),2.))) - 1.);
                }
                else {
                    if (i==j)
                        dPhi1dp(i,j) = 1.E12;
                    else
                        dPhi1dp(i,j) = 0.;
                }
                
                if(fabs(Phi_p(i)) > 0.) {
                    dPhi2dp(i,j) = (Phi_l(i)/(sqrt(pow(Phi_l(i),2.)))+1.)*denom_l(i,j);
                }
                else {
                    if (i==j)
                        dPhi2dp(i,j) = (Phi_l(i)/(sqrt(pow(Phi_l(i),2.)))+1.)*denom_l(i,j);
                    else
                        dPhi2dp(i,j) = 0.;
                }
                
                denomFB(i,j) = 2.*dPhi1dp(i,j)*Phi_1(i) + 2.*dPhi2dp(i,j)*Phi_2(i);
            }
        }
        else {
            
            if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
                FB(i) = sqrt(pow(Phi_p(i),2.) + pow(Dpstar(i),2.)) + Phi_p(i) - Dpstar(i);
            }
            else if(fabs(Phi_p(i)) > 0.) {
                FB(i) = sqrt(pow(Phi_p(i),2.)) + Phi_p(i);
            }
            else if(fabs(Dpstar(i)) > 0.) {
                FB(i) = sqrt(pow(Dpstar(i),2.)) - Dpstar(i);
            }
            else {
                FB(i) = 0.;
            }
            
            for (int j=0; j<n; j++) {
                if ((fabs(Phi_p(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
                    denomFB(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dpstar(i),2.)))+1.)*denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Phi_p(i),2.) + pow(Dpstar(i),2.))) - 1.);
                }
                else if(fabs(Phi_p(i)) > 0.) {
                    denomFB(i,j) = (Phi_p(i)/(sqrt(pow(Phi_p(i),2.)))+1.)*denom(i,j) - delta(i,j)*dDpstar(i);
                }
                else if(fabs(Dpstar(i)) > 0.) {
                    denomFB(i,j) = denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Dpstar(i),2.))) - 1.);
                }
                else {
                    if (i==j)
                        denomFB(i,j) = 1.E12;
                    else
                        denomFB(i,j) = 0.;
                }
            }
        }
    }
    
    try {
        if (fabs(det(denomFB)) > simcoon::limit) {
            dp = -1.*solve(denomFB, FB);
        }
        else {
            dp = zeros(n);
        }
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Fischer_Burmeister_m_limits.");
    }  
    
    //New update of the transformation/orientation multipliers
    Dp = Dp + dp;
    
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
        
}
    
mat denom_FB_m(const vec &Phi, const mat &denom, const vec &Dp)
{
    
    int n=Phi.n_elem;
    
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Determine the eigenvalues of denom
    vec factor_denom = zeros(n);
    for (int i=0; i<n; i++) {
        factor_denom(i) = fabs(denom(i,i));
    }
    
    vec Dpstar = Dp%factor_denom;
    vec dDpstar = factor_denom;
    
    //Normalized Fischer-Burmeister set of equations
    for (int i=0; i<n; i++) {
        
        for (int j=0; j<n; j++) {
            if ((fabs(Phi(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.)))+1.)*denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.))) - 1.);
            }
            else if(fabs(Phi(i)) > 0.) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.)))+1.)*denom(i,j) - delta(i,j)*dDpstar(i);
            }
            else if(fabs(Dpstar(i)) > 0.) {
                denomFB(i,j) = denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Dpstar(i),2.))) - 1.);
            }
            else {
                if (i==j)
                    denomFB(i,j) = 1.E12;
                else
                    denomFB(i,j) = 0.;
            }
        }
        
    }
    return denomFB;
}
    
} //namespace simcoon
