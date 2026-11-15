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

///@file unified_T.cpp
///@brief Unified thermomechanical SMA model supporting multiple elastic symmetries and transformation criteria

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <assert.h>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/recovery_props.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/sma_flow.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

void sma_thermal_face(const SMAIncrement &inc, const vec &sigma, const vec &DEtot, const double &rho, const double &c_pA, const double &c_pM, const double &T, const double &DT, const int &tangent_mode, double &r, mat &dSdE, vec &dSdT, vec &drdE, double &drdT, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir) {
    // none clamped to continuum: P_epsilon, invBhat feed r (see header)
    const int mode = inc.closest_point ? (tangent_mode == tangent_closest_point_none ? tangent_closest_point_continuum : tangent_mode)
                                       : tangent_continuum;
    const ContinuumTangent ct = inc.tangent(mode);
    dSdE = ct.Lt;
    const uword N = inc.Ds.n_elem;
    // closest point without a solve (elastic trial, step cut): the elastic operator, zero sensitivities
    const bool solved = !inc.closest_point || (inc.rm.converged && !inc.rm.kappa_j.is_empty());
    const mat P_epsilon = solved ? ct.P_epsilon : mat(6, N, fill::zeros);
    const mat invBhat = solved ? ct.invBhat : mat(N, N, fill::zeros);
    const vec &Ds = inc.Ds;
    const vec &sigma_start = inc.sigma_start;
    const vec &alpha = inc.alpha;
    const double xi = inc.xi;
    const double xi_start = inc.xi_start;
    const double T_init = inc.T_init;
    const double rhoDs0 = inc.rhoDs0;

    // trial stress rate of a temperature move: the closest point freezes alpha at xi_n and folds
    // Dalpha into the flows, so its thermal strain is alpha(xi_n+1); the algorithmic operator acts through M^-1
    const bool algo = solved && !ct.kappa_tilde.is_empty();
    const vec alpha_q = inc.closest_point ? vec(alpha + (xi - xi_start)*inc.Dalpha) : alpha;
    const vec L_alpha = algo ? vec(ct.L_tilde*alpha_q) : vec(inc.L_tan*alpha_q);
    const mat kappa_q = !solved ? mat(6, N, fill::zeros) : (algo ? ct.kappa_tilde : inc.flux());
    const mat dPhi = solved ? inc.gradients() : mat(6, N, fill::zeros);

    // Delta s = invBhat * rhs, Bhat(l, j) = criterion l, mechanism j
    // the transformation rows (0, 1) in the expressions of the two-mechanism kernel, the reorientation row added after
    vec rhs(N), P_theta(N);
    for (uword l = 0; l < N; l++) rhs(l) = inc.dPhi_dtheta(l) - sum(dPhi.col(l)%L_alpha);
    for (uword l = 0; l < N; l++) {
        P_theta(l) = invBhat(l, 0)*rhs(0) + invBhat(l, 1)*rhs(1);
        if (N > 2) P_theta(l) += invBhat(l, 2)*rhs(2);
    }
    vec kP = kappa_q.col(0)*P_theta(0);
    for (uword j = 1; j < N; j++) kP += kappa_q.col(j)*P_theta(j);
    dSdT = -1.*L_alpha - kP;

    //Preliminaries for the computation of mechanical and thermal work
    const double c_0A = rho*c_pA;
    const double c_0M = rho*c_pM;
    const double c_0 = c_0A*(1-xi) + c_0M*xi;
    const double Dc_0 = c_0M - c_0A;
    const double T1 = T + DT;

    const double eta_r = (c_0A + Dc_0*xi)*log(T1/T_init) + sum(alpha%sigma) + rhoDs0*xi;
    const double eta_r_start = (c_0A + Dc_0*xi_start)*log(T/T_init) + sum(alpha%sigma_start) + rhoDs0*xi_start;
    const double Deta = eta_r - eta_r_start;

    // entropy forces of the mechanisms: the transformation rows, none for the reorientation
    vec dA_theta_ds(N, fill::zeros);
    dA_theta_ds(0) = Dc_0*log(T1/T_init) + sum(inc.Dalpha%sigma) + rhoDs0;
    dA_theta_ds(1) = -1.*Dc_0*log(T1/T_init) - sum(inc.Dalpha%sigma) - rhoDs0;

    // A_Reo = 0: the reorientation work is in DET
    const double Dgamma_loc = 0.5*sum((sigma_start+sigma)%inc.DET) + 0.5*(inc.A_start(0) + inc.A_end(0))*Ds(0) + 0.5*(inc.A_start(1) + inc.A_end(1))*Ds(1);

    // dissipation derivatives, flow directions held fixed over the increment
    vec dD = dSdE.t()*inc.DET;
    vec g_T(N);
    for (uword j = 0; j < N; j++) {
        vec g = dSdE.t()*inc.dA_dsigma.col(j);
        for (uword k = 0; k < N; k++) g += inc.dA_ds(j, k)*P_epsilon.col(k);
        dD += g*Ds(j);
        g_T(j) = inc.dA_dtheta(j) + sum(inc.dA_dsigma.col(j)%dSdT) + inc.dA_ds(j, 0)*P_theta(0) + inc.dA_ds(j, 1)*P_theta(1);
        if (N > 2) g_T(j) += inc.dA_ds(j, 2)*P_theta(2);
    }
    double dD_T_turn = 0.;
    vec f(N);
    for (uword j = 0; j < N; j++) f(j) = 0.5*(sum((sigma_start+sigma)%inc.flows.col(j)) + inc.A_start(j) + inc.A_end(j));
    // closest point: the turning of the flows and their multiplier chain, d DET = sum_k (Lambda_k dDs_k + Ds_k dLambda_k)
    if (!inc.dflow_dsigma.is_empty()) {
        vec turn = zeros(6);
        for (uword k = 0; k < N; k++) turn += Ds(k)*(inc.dflow_dsigma.cols(6*k, 6*k + 5).t()*(sigma_start+sigma));
        dD += dSdE.t()*turn;
        for (uword j = 0; j < N; j++) f(j) += 0.5*sum((sigma_start+sigma)%inc.flow_chain.col(j));
        dD_T_turn = sum(dSdT%turn);
    }
    vec dDgamma_dE = 0.5*dD;
    for (uword j = 0; j < N; j++) dDgamma_dE += f(j)*P_epsilon.col(j);
    double dD_T = sum(dSdT%inc.DET) + g_T(0)*Ds(0) + g_T(1)*Ds(1);
    if (N > 2) dD_T += g_T(2)*Ds(2);
    dD_T += dD_T_turn;
    double dDgamma_dT = 0.5*dD_T + f(0)*P_theta(0) + f(1)*P_theta(1);
    if (N > 2) dDgamma_dT += f(2)*P_theta(2);

    // Entropy derivatives at the end of the increment
    const vec deta_dE = dSdE.t()*alpha + dA_theta_ds(0)*P_epsilon.col(0) + dA_theta_ds(1)*P_epsilon.col(1);
    const double deta_dT = c_0/T1 + sum(dSdT%alpha) + dA_theta_ds(0)*P_theta(0) + dA_theta_ds(1)*P_theta(1);

    // heat of the increment and its derivatives (heat_source_increment)
    heat_source_increment(Dgamma_loc, dDgamma_dE, dDgamma_dT, Deta, deta_dE, deta_dT, T1, r, drdE, drdT);

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-inc.DET)) - 0.5*(inc.A_start(0) + inc.A_end(0))*Ds(0) - 0.5*(inc.A_start(1) + inc.A_end(1))*Ds(1);
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;

    Wt += T1*Deta;
    Wt_r += T1*Deta;
    Wt_ir += 0.;
}

void umat_sma_unified_T_T(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, vec &dSdT, vec &drdE, double &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode, const int &ncp) {
    const NCPMethod ncp_m = sma_ncp_method(ncp, tangent_mode, "unified_T_T");

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(tnew_dt);

    ///@brief Determine elasticity type and criteria type from umat_name
    bool cubic_elasticity = false;
    bool aniso_criteria = false;

    if (umat_name == "SMADI") {
        cubic_elasticity = false;
        aniso_criteria = false;
    }
    else if (umat_name == "SMADC") {
        cubic_elasticity = true;
        aniso_criteria = false;
    }
    else if (umat_name == "SMAAI") {
        cubic_elasticity = false;
        aniso_criteria = true;
    }
    else if (umat_name == "SMAAC") {
        cubic_elasticity = true;
        aniso_criteria = true;
    }
    else {
        throw std::invalid_argument("Unknown umat_name in umat_sma_unified_T_T: " + umat_name
                                    + " (valid: SMADI, SMADC, SMAAI, SMAAC)");
    }

    ///@brief Property offset depends on elastic symmetry type
    int offset = 0;

    // Taking the values from the material created
    double rho = props(0);
    double c_pA = props(1);
    double c_pM = props(2);
    int flagT = props(3);

    // closest point: the mechanical increment on the mechanical props, the thermal face on its result
    if (integrates_closest_point(tangent_mode)) {
        if (start) {
            Wm = 0.;
            Wm_r = 0.;
            Wm_ir = 0.;
            Wm_d = 0.;
        }
        mat L_out;
        const SMAIncrement inc = sma_unified_T_increment(umat_name, Etot, DEtot, sigma, L_out, DR, props.subvec(3, props.n_elem - 1), statev, T, DT, start, tnew_dt, tangent_mode, ncp_m);
        sma_thermal_face(inc, sigma, DEtot, rho, c_pA, c_pM, T, DT, tangent_mode, r, dSdE, dSdT, drdE, drdT, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir);
        return;
    }

    mat M_A;
    mat M_M;

    if (!cubic_elasticity) {
        // Isotropic elasticity: EA, EM, nuA, nuM
        M_A = M_iso(props(4), props(6), "Enu");
        M_M = M_iso(props(5), props(7), "Enu");
        offset = 8;
    }
    else {
        // Cubic elasticity: EA, EM, nuA, nuM, GA, GM
        M_A = M_cubic(props(4), props(6), props(8), "EnuG");
        M_M = M_cubic(props(5), props(7), props(9), "EnuG");
        offset = 10;
    }

    // Common parameters after elasticity
    double alphaA_iso = props(offset);
    double alphaM_iso = props(offset + 1);
    //parameters for Hcur
    double Hmin = props(offset + 2);
    double Hmax = props(offset + 3);
    double k1 = props(offset + 4);
    double sigmacrit = props(offset + 5);
    //parameters for phase transformation
    double C_A = props(offset + 6);
    double C_M = props(offset + 7);
    double Ms0 = props(offset + 8);
    double Mf0 = props(offset + 9);
    double As0 = props(offset + 10);
    double Af0 = props(offset + 11);
    //Additional parameters
    double n1 = props(offset + 12);
    double n2 = props(offset + 13);
    double n3 = props(offset + 14);
    double n4 = props(offset + 15);
    double sigmacaliber = props(offset + 16);
    double prager_b = props(offset + 17);
    double prager_n = props(offset + 18);
    double sigmastar = 0.;

    //Set the Lagrange multipliers coefficients
    double c_lambda = props(offset + 19);
    double p0_lambda = props(offset + 20);
    double n_lambda = props(offset + 21);
    double alpha_lambda = props(offset + 22);

    // DFA parameters for anisotropic criteria (only used if aniso_criteria == true)
    vec DFA_params = zeros(7);
    if (aniso_criteria) {
        double F_dfa = props(offset + 23);
        double G_dfa = props(offset + 24);
        double H_dfa = props(offset + 25);
        double L_dfa = props(offset + 26);
        double M_dfa = props(offset + 27);
        double N_dfa = props(offset + 28);
        double K_dfa = props(offset + 29);
        DFA_params = {F_dfa, G_dfa, H_dfa, L_dfa, M_dfa, N_dfa, K_dfa};
    }

    ///@brief Temperature initialization
    double T_init = statev(0);
    ///@brief Martensite volume fraction initialization
    double xi = statev(1);
    ///@brief mean strain tensor creation
    vec ET = zeros(6);
    ET(0) = statev(2);
    ET(1) = statev(3);
    ET(2) = statev(4);
    ET(3) = statev(5);
    ET(4) = statev(6);
    ET(5) = statev(7);

    ///@brief ETMax allow the definition of the lambdaTR
    vec DETF = zeros(6);
    vec DETR = zeros(6);
    vec ETMean = zeros(6);

    double xiF = statev(8);
    double xiR = statev(9);

    double rhoDs0 = statev(10);
    double rhoDE0 = statev(11);
    double D = statev(12);
    double a1 = statev(13);
    double a2 = statev(14);
    double a3 = statev(15);
    double Y0t = statev(16);

    // ######################  Elastic compliance and stiffness #################################
    // Reuss mixing on compliance tensor
    mat DM = M_M - M_A;
    mat L;
    mat M_eff = xi*M_M + (1. - xi)*M_A;
    L = inv_cubic(M_eff);

    //definition of the CTE tensor
    vec alpha = (alphaM_iso*xi + alphaA_iso*(1.-xi))*Ith();
    vec Dalpha = (alphaM_iso - alphaA_iso)*Ith();

    ///@brief Initialization
    if(start) {

        T_init = T;
        vec vide = zeros(6);
        sigma = zeros(6);
        ET = zeros(6);
        xiF = simcoon::limit;
        xiR = 0.;
        xi = xiF;

        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;

        //Get Hcurstar and its derivative to determine rhoDs0
        if (sigmacaliber > sigmacrit)
            sigmastar = sigmacaliber - sigmacrit;
        else
            sigmastar = 0.;

        double Hcurstar = Hmin + (Hmax - Hmin)*(1. - exp(-k1*sigmastar));
        assert(Hcurstar > 1E-12);

        double dHcurstar = (Hmax - Hmin)*k1*exp(-k1*sigmastar);
        rhoDs0 = -2.*C_M*C_A*(Hcurstar + sigmacaliber*(dHcurstar + (M_M(0,0)-M_A(0,0))))/(C_M+C_A);

        //Determination of the Smooth temperatures
        double MsSmooth = 0.;
        double MfSmooth = 0.;
        double AsSmooth = 0.;
        double AfSmooth = 0.;

        if(flagT == 0) {
            MsSmooth = 0.5*Ms0*(1.+ (n1+1)*pow(2.,-n1) + (n2-1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2)) + 0.5*Mf0*(-1.+ (n1-1)*pow(2.,-n1) + (n2+1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2));
            MfSmooth = 0.5*Ms0*(-1.+ (n1+1)*pow(2.,-n1) + (n2-1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2)) + 0.5*Mf0*(1.+ (n1-1)*pow(2.,-n1) + (n2+1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2));
            AsSmooth = 0.5*As0*(1.+ (n3-1)*pow(2.,-n3) + (n4+1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4)) + 0.5*Af0*(-1.+ (n3+1)*pow(2.,-n3) + (n4-1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4));
            AfSmooth = 0.5*As0*(-1.+ (n3-1)*pow(2.,-n3) + (n4+1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4)) + 0.5*Af0*(1.+ (n3+1)*pow(2.,-n3) + (n4-1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4));
        }
        else {
            MsSmooth = Ms0;
            MfSmooth = Mf0;
            AsSmooth = As0;
            AfSmooth = Af0;
        }

        rhoDE0 = 0.5*rhoDs0*(MsSmooth + AfSmooth);
        D = (C_M-C_A)*(Hcurstar + sigmacaliber*(dHcurstar + (M_M(0,0)-M_A(0,0))))/((C_A+C_M)*(Hcurstar + sigmacaliber*dHcurstar));
        a1 = rhoDs0*(MfSmooth - MsSmooth);
        a2 = rhoDs0*(AsSmooth - AfSmooth);
        a3 = -0.25*a1*(1 + 1./(n1 + 1.) - 1./(n2 + 1.))+0.25*a2*(1. + 1./(n3 + 1.) - 1./(n4 + 1.));
        Y0t = 0.5*rhoDs0*(MsSmooth - AfSmooth) - a3;
    }

    //Rotation of internal variables (tensors)
    ET = rotate_strain(ET, DR);

    //Variables values at the start of the increment
    vec	sigma_start = sigma;
    vec ET_start = ET;
    double xi_start = xi;
    vec alpha_start = alpha;

    // Find Hcur explicit
    if (Mises_stress(sigma) > sigmacrit)
        sigmastar = Mises_stress(sigma) - sigmacrit;
    else
        sigmastar = 0.;

    double Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

    //definition of Lambdas associated to transformation
    vec lambdaTF;
    if (aniso_criteria) {
        lambdaTF = Hcur*dDrucker_ani_stress(sigma, DFA_params, prager_b, prager_n);
    }
    else {
        lambdaTF = Hcur*dDrucker_stress(sigma, prager_b, prager_n);
    }

    if (Mises_strain(ET) > 1E-6)
        ETMean = dev(ET) / (xi);
    else if (Mises_stress(sigma) < 1.E-6)
        ETMean = lambdaTF;
    else
        ETMean = 0.*Ith();

    vec lambdaTR = -1.*ETMean;

    //Modified Y functions, Lagrange multipliers and criteria: set by the transformation loop before any use
    double YtF = 0.;
    double YtR = 0.;

    double HfF = 0.;
    double HfR = 0.;

    //Hardening function definition (Smooth hardening functions)
    if ((xi > 0.)&&((1. - xi) > 0.)) {
        HfF = 0.5*a1*(1. + pow(xi,n1) - pow(1. - xi,n2)) + a3;
        HfR = 0.5*a2*(1. + pow(xi,n3) - pow((1. - xi),n4)) - a3;
    }
    else if ((xi <= 0.)&&((1. - xi) > 0.)) {
        HfF = 0.5*a1*(1. - pow(1. - xi,n2)) + a3;
        HfR = 0.5*a2*(1. - pow((1. - xi),n4)) - a3;
    }
    else if ((xi > 0.)&&((1. - xi) <= 0.)) {
        HfF = 0.5*a1*(1. + pow(xi,n1)) + a3;
        HfR = 0.5*a2*(1. + pow(xi,n3)) - a3;
    }
    else {
        HfF =  a3;
        HfR = - a3;
    }

    //Lagrange multiplliers due to Physical limitations
    double lambda0 = 0.;
    double lambda1 = 0.;

    //Define the value of DM_sig
    vec DM_sig = (DM*sigma_start);
    //Define the value of Dalpha_T
    vec Dalpha_T = Dalpha*(T+DT-T_init);

    //Set the thermo forces
    double A_xiF = rhoDs0*(T+DT) - rhoDE0 + 0.5*sum(sigma%DM_sig) + sum(sigma%Dalpha)*(T+DT-T_init) - HfF;
    double A_xiF_start = rhoDs0*(T) - rhoDE0 + 0.5*sum(sigma_start%DM_sig) + sum(sigma_start%Dalpha)*(T-T_init) - HfF;
    double A_xiR = -1.*rhoDs0*(T+DT) + rhoDE0 - 0.5*sum(sigma%DM_sig) - sum(sigma%Dalpha)*(T+DT-T_init) + HfR;
    double A_xiR_start = -1.*rhoDs0*(T) + rhoDE0 - 0.5*sum(sigma_start%DM_sig) - sum(sigma_start%Dalpha)*(T-T_init) + HfR;

    //Transformation criteria
    double PhihatF = 0.;
    double PhihatR = 0.;

    //Variables required for the loop
    vec s_j = zeros(2);
    s_j(0) = xiF;
    s_j(1) = xiR;
    vec Ds_j = zeros(2);
    vec ds_j = zeros(2);

    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - ET;
    sigma = el_pred(L, Eel);

    //Define the functions for the system to solve
    vec Phi = zeros(2);
    mat B = zeros(2,2);
    vec Y_crit = zeros(2);

    //Define the function for the system to solve
    double dHfF = 0.;
    double dHfR = 0.;
    vec dHcurdsigma = zeros(6);
    //Relative to forward transformation
    vec dPhihatFdsigma = zeros(6);
    double dPhihatFdxiF = 0.;
    double dPhihatFdxiR = 0.;

    vec dA_xiFdsigma = zeros(6);
    double dA_xiFdxiF = 0.;
    double dA_xiFdxiR = 0.;

    vec dlambda1dsigma = zeros(6);
    double dlambda1dxiF = 0.;
    double dlambda1dxiR = 0.;

    vec dYtFdsigma = zeros(6);
    double dYtFdxiF = 0.;
    double dYtFdxiR = 0.;

    vec dPhiFdsigma = zeros(6);
    double dPhiFdxiF = 0.;
    double dPhiFdxiR = 0.;

    //Relative to reverse transformation: only the finite sigma-gradients (K(1,.) goes through ETMean)
    vec dPhihatRdsigma = zeros(6);
    vec dA_xiRdsigma = zeros(6);
    double dA_xiRdxiF = 0.;
    double dA_xiRdxiR = 0.;
    vec dlambda0dsigma = zeros(6);
    double dlambda0dxiF = 0.;
    double dlambda0dxiR = 0.;
    vec dYtRdsigma = zeros(6);
    vec dPhiRdsigma = zeros(6);

    //Compute the explicit flow direction
    std::vector<vec> kappa_j(2);
    mat K = zeros(2,2);

    //Loop parameters
    int compteur = 0;
    double error = 1.;

    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        // Reuss mixing on compliance tensor
        M_eff = xi*M_M + (1. - xi)*M_A;
        L = inv_cubic(M_eff);

        DM_sig = DM*sigma;

        // equivalent stresses and gradients at the current stress, shared below (Dalpha_T is loop-invariant)
        const double sig_eq = Mises_stress(sigma);
        const vec eta_sig = eta_stress(sigma);
        double Dr_sig;
        vec dDr_sig;
        if (aniso_criteria) {
            Dr_sig = Drucker_ani_stress(sigma, DFA_params, prager_b, prager_n);
            dDr_sig = dDrucker_ani_stress(sigma, DFA_params, prager_b, prager_n);
        }
        else {
            Dr_sig = Drucker_stress(sigma, prager_b, prager_n);
            dDr_sig = dDrucker_stress(sigma, prager_b, prager_n);
        }

        lambdaTF = Hcur * dDr_sig;
        lambdaTR = -1. * ETMean;

        kappa_j[0] = L*(lambdaTF + DM_sig + Dalpha_T);
        kappa_j[1] = L*(lambdaTR - DM_sig - Dalpha_T); //derivative w/ xiR is minus derivative w/ xi

        //Hardening function definition (Smooth hardening functions)
        if ((xi > 0.) && ((1. - xi) > 0.)) {
            HfF = 0.5*a1*(1. + pow(xi, n1) - pow(1. - xi, n2)) + a3;
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            HfF = 0.5*a1*(1. - pow(1. - xi, n2)) + a3;
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            HfF = 0.5*a1*(1. + pow(xi, n1)) + a3;
        }
        else  {
            HfF = 0.5*a1 + a3;
        }

        if ((xi > 0.) && ((1. - xi) > 0.)) {
            HfR = 0.5*a2*(1. + pow(xi, n3) - pow((1. - xi), n4)) - a3;
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            HfR = 0.5*a2*(1. - pow((1. - xi), n4)) - a3;
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            HfR = 0.5*a2*(1. + pow(xi, n3)) - a3;
        }
        else
            HfR = - a3;

        // Find Hcur explicit
        if (sig_eq > sigmacrit)
            sigmastar = sig_eq - sigmacrit;
        else
            sigmastar = 0.;

        Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

        //Forward transformation thermodynamic force
        PhihatF = Hcur*Dr_sig;
        A_xiF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sigma%DM_sig) + sum(sigma%Dalpha)*(T + DT - T_init) - HfF;
        lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        const double DHM = D*Hcur*sig_eq;
        YtF = Y0t + DHM;
        Phi(0) = PhihatF + A_xiF - lambda1 - YtF;

        //Reverse transformation thermodynamic force
        PhihatR = sum(sigma%ETMean);
        A_xiR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sigma%DM_sig) - sum(sigma%Dalpha)*(T + DT - T_init) + HfR;
        lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        YtR = Y0t + D*PhihatR;
        Phi(1) = -1.*PhihatR + A_xiR + lambda0 - YtR;  // PhiR < 0.

        if ((xi > 0.) && ((1. - xi) > 0.)) {
            dHfF = 0.5*a1*(n1*pow(xi, n1 - 1.) + n2*pow(1. - xi, n2 - 1.));
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            dHfF = 0.5*a1*(n2*pow(1. - xi, n2 - 1.));
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            dHfF = 0.5*a1*(n1*pow(xi, n1 - 1.));
        }
        else
        {
            dHfF = 0.;
        }

        if ((xi > 0.) && ((1. - xi) > 0.)) {
            dHfR = 0.5*a2*(n3*pow(xi, n3 - 1.) + n4*pow(1. - xi, n4 - 1.));
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            dHfR = 0.5*a2*(n4*pow(1. - xi, n4 - 1.));
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            dHfR = 0.5*a2*(n3*pow(xi, n3 - 1.));
        }
        else
        {
            dHfR = 0.;
        }

        dHcurdsigma = k1*(Hmax - Hmin)*exp(-1.*k1*sigmastar)*eta_sig;

        //Related to forward transformation
        dPhihatFdsigma = dHcurdsigma * Dr_sig + Hcur * dDr_sig;
        dPhihatFdxiF = 0.;
        dPhihatFdxiR = 0.;

        dA_xiFdsigma = DM_sig + Dalpha_T;
        dA_xiFdxiF = -dHfF;
        dA_xiFdxiR = dHfF;

        dlambda1dsigma = zeros(6);
        dlambda1dxiF = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda1dxiR = -1.*dlambda1dxiF;

        dYtFdsigma = D*(dHcurdsigma * sig_eq + Hcur * eta_sig);
        dYtFdxiF = 0.;
        dYtFdxiR = 0.;

        dPhiFdsigma = dPhihatFdsigma + dA_xiFdsigma - dlambda1dsigma - dYtFdsigma;
        dPhiFdxiF = dPhihatFdxiF + dA_xiFdxiF - dlambda1dxiF - dYtFdxiF;
        dPhiFdxiR = dPhihatFdxiR + dA_xiFdxiR - dlambda1dxiR - dYtFdxiR;

        //Relative to reverse transformation
        dPhihatRdsigma = ETMean;
        dA_xiRdsigma = -1.*DM_sig - Dalpha_T;
        dA_xiRdxiF =  dHfR;
        dA_xiRdxiR = -dHfR;

        dlambda0dsigma = zeros(6);
        dlambda0dxiR = dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda0dxiF = -1.*dlambda0dxiR;

        dYtRdsigma = D*ETMean;
        dPhiRdsigma = -1.*dPhihatRdsigma + dA_xiRdsigma + dlambda0dsigma - dYtRdsigma;

        // K(1,.) through ETMean: avoids a ~12-digit cancellation at small xi (Mechanical/SMA/unified_T.hpp)
        const vec dPhiRdETMean = -(1. + D) * sigma;
        vec Lambda_ETMean_F = zeros(6);
        if (Mises_strain(lambdaTF - ETMean) > simcoon::iota) {
            Lambda_ETMean_F = (lambdaTF - ETMean) / xi;
        }
        const vec Lambda_ETMean_R = zeros(6);
        const double dPhiRdxiF_finite = dA_xiRdxiF + dlambda0dxiF;   // = +dHfR - d\lambda_0/d\xi 
        const double dPhiRdxiR_finite = dA_xiRdxiR + dlambda0dxiR;   // = -dHfR + d\lambda_0/d\xi 

        K(0,0) = dPhiFdxiF;
        K(0,1) = dPhiFdxiR;
        K(1,0) = dPhiRdxiF_finite + sum(dPhiRdETMean % Lambda_ETMean_F);
        K(1,1) = dPhiRdxiR_finite + sum(dPhiRdETMean % Lambda_ETMean_R);

        B(0,0) = -1.*sum(dPhiFdsigma%kappa_j[0]) + K(0,0);
        B(0,1) = -1.*sum(dPhiFdsigma%kappa_j[1]) + K(0,1);
        B(1,0) = -1.*sum(dPhiRdsigma%kappa_j[0]) + K(1,0);
        B(1,1) = -1.*sum(dPhiRdsigma%kappa_j[1]) + K(1,1);

        // magnitude scales, not the signed YtF/YtR (they cross zero when D < 0)
        Y_crit(0) = std::max(fabs(Y0t) + fabs(DHM), simcoon::iota);
        Y_crit(1) = std::max(fabs(Y0t) + fabs(D*PhihatR), simcoon::iota);

        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);

        s_j(0) += ds_j(0);
        s_j(1) += ds_j(1);

        ET = ET + ds_j(0)*lambdaTF + ds_j(1)*lambdaTR;

        xiF = s_j(0);
        xiR = s_j(1);
        xi = xi + ds_j(0) - 1.*ds_j(1);

        DETF += ds_j(0)*lambdaTF;
        DETR += -1.*ds_j(1)*lambdaTR;
        
        if((Mises_strain(ET) > simcoon::precision_umat)&&(xi > simcoon::precision_umat))
        {
            ETMean = dev(ET) / (xi);
        }
        else {
            ETMean = lambdaTF;
        }

        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - ET;
        sigma = el_pred(L, Eel);
    }

    SMAIncrement inc;
    inc.L_tan = L;
    inc.Ds = Ds_j;
    inc.kappa = join_rows(kappa_j[0], kappa_j[1]);
    inc.dPhi = join_rows(dPhiFdsigma, dPhiRdsigma);
    inc.Bhat = zeros(2, 2);
    inc.Bhat(0,0) = sum(dPhiFdsigma%kappa_j[0]) - K(0,0);
    inc.Bhat(0,1) = sum(dPhiFdsigma%kappa_j[1]) - K(0,1);
    inc.Bhat(1,0) = sum(dPhiRdsigma%kappa_j[0]) - K(1,0);
    inc.Bhat(1,1) = sum(dPhiRdsigma%kappa_j[1]) - K(1,1);
    inc.flows = join_rows(lambdaTF, lambdaTR);
    inc.A_start = {A_xiF_start, A_xiR_start};
    inc.A_end = {A_xiF, A_xiR};
    inc.dA_dsigma = join_rows(dA_xiFdsigma, dA_xiRdsigma);
    inc.dA_ds = {{dA_xiFdxiF, dA_xiFdxiR}, {dA_xiRdxiF, dA_xiRdxiR}};
    const double dA_xiFdtheta = rhoDs0 + sum(sigma%Dalpha);
    inc.dA_dtheta = {dA_xiFdtheta, -1.*rhoDs0 - sum(sigma%Dalpha)};
    inc.dPhi_dtheta = inc.dA_dtheta;
    inc.sigma_start = sigma_start;
    inc.DET = ET - ET_start;
    inc.alpha = alpha_start;
    inc.Dalpha = Dalpha;
    inc.xi_start = xi_start;
    inc.xi = xi;
    inc.rhoDs0 = rhoDs0;
    inc.T_init = T_init;
    sma_thermal_face(inc, sigma, DEtot, rho, c_pA, c_pM, T, DT, tangent_mode, r, dSdE, dSdT, drdE, drdT, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir);

    ///@brief statev evolving variables
    statev(0) = T_init;
    statev(1) = xi;

    statev(2) = ET(0);
    statev(3) = ET(1);
    statev(4) = ET(2);
    statev(5) = ET(3);
    statev(6) = ET(4);
    statev(7) = ET(5);

    statev(8) = xiF;
    statev(9) = xiR;

    if(start) {
        statev(10) = rhoDs0;
        statev(11) = rhoDE0;
        statev(12) = D;

        statev(13) = a1;
        statev(14) = a2;
        statev(15) = a3;

        statev(16) = Y0t;
    }

}

} //namespace simcoon
