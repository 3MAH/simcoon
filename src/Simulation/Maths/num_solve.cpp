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
#include <armadillo>
#include <simcoon/parameter.hpp>

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
    
    if (fabs(det(denom)) > sim_limit) {
        dp = -1.*solve(denom,Phi);
    }
    else {
        dp = zeros(n);
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
    
    if (fabs(det(denomFB)) > sim_limit) {
        dp = -1.*solve(denomFB,FB);
    }
    else {
        dp = zeros(n);
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
    
    if (fabs(det(denomFB)) > sim_limit) {
        dp = -1.*solve(denomFB,FB);
    }
    else {
        dp = zeros(n);
    }
    
    //New update of the vector of the variable increment throughout the step
    Dp = Dp + dp;
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
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
    
    if (fabs(det(denomFB)) > sim_limit) {
        dp = -1.*solve(denomFB, FB);
    }
    else {
        dp = zeros(n);
    }
    
    //New update of the transformation/orientation multipliers
    Dp = Dp + dp;

    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
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
    
    if (fabs(det(denomFB)) > sim_limit) {
        dp = -1.*solve(denomFB, FB);
    }
    else {
        dp = zeros(n);
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
