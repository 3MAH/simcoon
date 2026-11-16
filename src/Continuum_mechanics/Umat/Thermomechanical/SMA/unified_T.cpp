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
#include <simcoon/Continuum_mechanics/Functions/sma.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

void umat_sma_unified_T_T(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, vec &dSdT, vec &drdE, double &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode, const int &ncp) {
    const NCPMethod ncp_m = sma_ncp_method(ncp, tangent_mode, "unified_T_T");

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(ndi);
    UNUSED(nshr);

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
    auto drucker = [&](const vec &s, const int level) {
        return aniso_criteria ? Drucker_ani_stress_response(s, DFA_params, prager_b, prager_n, level)
                              : Drucker_stress_response(s, prager_b, prager_n, level);
    };

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

    // Find Hcur explicit
    if (Mises_stress(sigma) > sigmacrit)
        sigmastar = Mises_stress(sigma) - sigmacrit;
    else
        sigmastar = 0.;

    double Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

    //definition of Lambdas associated to transformation
    vec lambdaTF = Hcur*drucker(sigma, 1).df;

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

    //Hardening functions at xi clamped to [0,1]
    SMAHardening hf = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
    HfF = hf.F;
    HfR = hf.R;

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

    // closest point: L(xi), alpha(xi) frozen at xi_n, exact since both are linear in xi (see header)
    ReturnMappingResult rm;
    mat dflow_dsigma;   // closest point: total dLambda^j/dsigma of the transformation flows (6x12)
    vec chain_F;        // closest point: Ds_R dLambda^R/dDs_F through ETMean (6)
    const bool use_cpp = integrates_closest_point(tangent_mode);
    const mat L_cpp = use_cpp ? L : mat();   // frozen stiffness used consistently by the helper and the tangent
    if (use_cpp) {
        const double xi_n = xi;
        const double thp = (T + DT - T_init);
        vec Dl_cur = zeros(2);        // multipliers of the refreshed state
        bool ETMean_is_flow = false;  // the guard branch ETMean = lambdaTF is active
        bool lin_stale = true;        // the linearisation below is not the one of the refreshed state

        // Backward-Euler state at (sig, Dl): xi, Hcur, lambdaTF, ETMean, ET, hardening, Lagrange.
        auto refresh_state = [&](const vec &sig, const vec &Dl) {
            Dl_cur = Dl;
            lin_stale = true;
            xi = xi_n + Dl(0) - Dl(1);
            Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*std::max(Mises_stress(sig) - sigmacrit, 0.)));
            lambdaTF = sma_transformation_flow(sig, sigmacrit, Hmin, Hmax, k1, aniso_criteria, DFA_params, prager_b, prager_n);
            // closed form ETMean (xi + Dl_R) = dev(ET_n) + Dl_F dev(lambdaTF): the sweep diverges when Dl_R > xi
            const double xi_before = xi + Dl(1);
            ETMean_is_flow = !(xi_before > simcoon::precision_umat);
            ETMean = ETMean_is_flow ? lambdaTF : vec((dev(ET_start) + Dl(0)*dev(lambdaTF))/xi_before);
            ET = ET_start + Dl(0)*lambdaTF - Dl(1)*ETMean;
            if (!((Mises_strain(ET) > simcoon::precision_umat) && (xi > simcoon::precision_umat))) {
                ETMean_is_flow = true;
                ETMean = lambdaTF;
                ET = ET_start + Dl(0)*lambdaTF - Dl(1)*ETMean;
            }
            hf = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
            HfF = hf.F;
            HfR = hf.R;
            lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
            lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
            return true;
        };
        // Criterion values at (sig, current state).
        auto PhiF_val = [&](const vec &sig) {
            const double PhF = Hcur*drucker(sig, 0).f;
            const vec DMsig = DM*sig;
            const double AxF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sig%DMsig) + sum(sig%Dalpha)*thp - HfF;
            return PhF + AxF - lambda1 - (Y0t + D*Hcur*Mises_stress(sig));
        };
        auto PhiR_val = [&](const vec &sig) {
            const vec DMsig = DM*sig;
            const double AxR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sig%DMsig) - sum(sig%Dalpha)*thp + HfR;
            return -1.*sum(sig%ETMean) + AxR + lambda0 - (Y0t + D*sum(sig%ETMean));
        };
        // effective flows: stiffness/CTE mixture variation folded in
        auto LamF_val = [&](const vec &sig) { return vec(lambdaTF + DM*sig + Dalpha*thp); };
        auto LamR_val = [&](const vec &sig) { return vec(-1.*ETMean - DM*sig - Dalpha*thp); };
        // analytic total derivatives at (sig, Dl_cur), once per iterate (derivation in header)
        mat dev_op = eye(6, 6);   // dev() of a Voigt vector
        dev_op.submat(0, 0, 2, 2) -= ones(3, 3)/3.;
        vec dPhiF_ds = zeros(6), dPhiR_ds = zeros(6), c_F = zeros(6);
        mat dLamF_ds = zeros(6, 6), dLamR_ds = zeros(6, 6), K_cpp = zeros(2, 2);
        auto linearise = [&](const vec &sig) {
            if (!lin_stale) return;
            lin_stale = false;
            const double sig_eq = Mises_stress(sig);
            const vec eta = eta_stress(sig);
            vec dHcur = zeros(6);
            if (sig_eq > sigmacrit) dHcur = (Hmax - Hmin)*k1*exp(-1.*k1*(sig_eq - sigmacrit))*eta;
            const CriterionResponse dr = drucker(sig, 2);
            const double Dr = dr.f;
            const vec &g = dr.df;
            const mat &ddDr = dr.ddf;
            const mat dlTF_ds = g*dHcur.t() + Hcur*ddDr;

            mat dETM_ds = dlTF_ds;
            vec dETM_dF = zeros(6);
            if (!ETMean_is_flow) {
                const double xi_before = xi_n + Dl_cur(0);
                dETM_ds = (Dl_cur(0)/xi_before)*(dev_op*dlTF_ds);
                dETM_dF = (dev(lambdaTF) - ETMean)/xi_before;
            }

            // slopes of the clamped hardening functions and of the penalties
            const SMAHardening hs = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
            const double dHfF_dxi = hs.dF, dHfR_dxi = hs.dR;
            const double dl1_dxi = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
            const double dl0_dxi = -1.*dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);

            const vec DMsig = DM*sig;
            dPhiF_ds = Dr*dHcur + Hcur*g + DMsig + Dalpha*thp - D*(sig_eq*dHcur + Hcur*eta);
            dPhiR_ds = -1.*(1. + D)*(ETMean + dETM_ds.t()*sig) - DMsig - Dalpha*thp;
            dLamF_ds = dlTF_ds + DM;
            dLamR_ds = -1.*dETM_ds - DM;
            K_cpp(0, 0) = -1.*(dHfF_dxi + dl1_dxi);
            K_cpp(0, 1) = dHfF_dxi + dl1_dxi;
            K_cpp(1, 0) = -1.*(1. + D)*sum(sig%dETM_dF) + dHfR_dxi + dl0_dxi;
            K_cpp(1, 1) = -1.*(dHfR_dxi + dl0_dxi);
            c_F = -1.*Dl_cur(1)*dETM_dF;
        };
        std::vector<ReturnMechanism> mechs(2);
        mechs[0].Phi = PhiF_val;
        mechs[1].Phi = PhiR_val;
        mechs[0].Lambda = [&](const vec &sig, vec &Lam) { Lam = LamF_val(sig); };
        mechs[1].Lambda = [&](const vec &sig, vec &Lam) { Lam = LamR_val(sig); };
        mechs[0].dPhi_dsigma = [&](const vec &sig, vec &n) { linearise(sig); n = dPhiF_ds; };
        mechs[1].dPhi_dsigma = [&](const vec &sig, vec &n) { linearise(sig); n = dPhiR_ds; };
        mechs[0].dLambda_dsigma = [&](const vec &sig, mat &D) { linearise(sig); D = dLamF_ds; };
        mechs[1].dLambda_dsigma = [&](const vec &sig, mat &D) { linearise(sig); D = dLamR_ds; };
        ReturnStateHooks hooks;
        hooks.update_state = refresh_state;
        hooks.K = [&](const vec &sig, const vec &, mat &K) { linearise(sig); K = K_cpp; };
        hooks.flow_state_coupling = [&](const vec &sig, const vec &, mat &dLdDs) { linearise(sig); dLdDs.col(0) = c_F; dLdDs.col(1).zeros(); };
        vec Ycrit_cpp(2);
        // FB scales floored at an absolute material scale: both thresholds vanish in virgin states
        const double Yfloor = 0.1*std::max(Y0t + D*Hmax*sigmacaliber, 1.e-2);
        Ycrit_cpp(0) = std::max(std::fabs(Y0t) + std::fabs(D*Hcur*Mises_stress(sigma)), Yfloor);
        Ycrit_cpp(1) = std::max(std::fabs(Y0t) + std::fabs(D*sum(sigma%ETMean)), Yfloor);

        ReturnMappingControl control;
        control.ncp = ncp_m;
        rm = closest_point_return_mapping(sigma, L_cpp, mechs, hooks, Ycrit_cpp, control);
        if (rm.converged) {
            sigma = rm.sigma;
            Ds_j = rm.Dlambda;
            s_j += Ds_j;
            xiF = s_j(0);
            xiR = s_j(1);
            DETF = Ds_j(0)*lambdaTF;
            DETR = Ds_j(1)*ETMean;
            // output L at the converged xi; the tangent keeps the frozen L_cpp
            M_eff = xi*M_M + (1. - xi)*M_A;
            L = inv_cubic(M_eff);
            A_xiF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sigma%(DM*sigma)) + sum(sigma%Dalpha)*thp - HfF;
            A_xiR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sigma%(DM*sigma)) - sum(sigma%Dalpha)*thp + HfR;
            hf = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
            dHfF = hf.dF;
            dHfR = hf.dR;
            dA_xiFdsigma = DM*sigma + Dalpha*thp;
            dA_xiRdsigma = -1.*dA_xiFdsigma;
            dA_xiFdxiF = -dHfF;
            dA_xiFdxiR = dHfF;
            dA_xiRdxiF = dHfR;
            dA_xiRdxiR = -dHfR;
            lambdaTR = -1.*ETMean;
            // flows of the transformation strain (mixture terms of the effective flows taken out), for the heat
            dflow_dsigma = join_rows(dLamF_ds - DM, dLamR_ds + DM);
            chain_F = c_F;
        }
        else {
            // step cut, never a silent cutting-plane fallback
            tnew_dt = 0.5;
            sigma = sigma_start;
            refresh_state(sigma, zeros(2));
        }
    }
    else {
    //Loop
    // backtracking of the cutting plane: the state before the last step
    double bt_error = -1., bt_xi = xi;
    int n_bt = 0;
    vec bt_s = s_j, bt_ET = ET, bt_DETF = DETF, bt_DETR = DETR, bt_Ds = Ds_j, bt_lF = zeros(6), bt_lR = zeros(6), bt_ds = zeros(2);
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        // Reuss mixing on compliance tensor
        M_eff = xi*M_M + (1. - xi)*M_A;
        L = inv_cubic(M_eff);

        DM_sig = DM*sigma;

        // equivalent stresses and gradients at the current stress, shared below (Dalpha_T is loop-invariant)
        const double sig_eq = Mises_stress(sigma);
        const vec eta_sig = eta_stress(sigma);
        const CriterionResponse dr_sig = drucker(sigma, 1);
        const double Dr_sig = dr_sig.f;
        const vec &dDr_sig = dr_sig.df;

        // Hcur at the current stress: the flow and the criterion of an iterate share it
        if (sig_eq > sigmacrit)
            sigmastar = sig_eq - sigmacrit;
        else
            sigmastar = 0.;

        Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

        lambdaTF = Hcur * dDr_sig;
        lambdaTR = -1. * ETMean;

        kappa_j[0] = L*(lambdaTF + DM_sig + Dalpha_T);
        kappa_j[1] = L*(lambdaTR - DM_sig - Dalpha_T); //derivative w/ xiR is minus derivative w/ xi

        //Hardening functions at xi clamped to [0,1]
        hf = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
        HfF = hf.F;
        HfR = hf.R;

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

        dHfF = hf.dF;
        dHfR = hf.dR;

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
        if ((bt_error >= 0.) && (error > 2.*bt_error) && (n_bt < 5)) {
            // the merit grew: back to the state before the last step, half of it along its flows
            s_j = bt_s; ET = bt_ET; xi = bt_xi; DETF = bt_DETF; DETR = bt_DETR; Ds_j = bt_Ds;
            lambdaTF = bt_lF; lambdaTR = bt_lR;
            ds_j = 0.5*bt_ds;
            bt_ds = ds_j;
            Ds_j += ds_j;
            L = inv_cubic(xi*M_M + (1. - xi)*M_A);
            n_bt++;
        }
        else {
            bt_error = error; n_bt = 0;
            bt_s = s_j; bt_ET = ET; bt_xi = xi; bt_DETF = DETF; bt_DETR = DETR; bt_Ds = Ds_j - ds_j;
            bt_lF = lambdaTF; bt_lR = lambdaTR; bt_ds = ds_j;
        }

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
    }

    // bound penalties of xi in the transformation forces, with their slopes: stored like the hardening
    A_xiF_start -= lagrange_pow_1(xi_start, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    A_xiR_start -= lagrange_pow_0(xi_start, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    A_xiF -= lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    A_xiR -= lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    {
        const double dl1 = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        const double dl0 = -1.*dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dA_xiFdxiF -= dl1;
        dA_xiFdxiR += dl1;
        dA_xiRdxiF += dl0;
        dA_xiRdxiR -= dl0;
    }

    // tangent: closest point on its converged solve (frozen L_cpp), cutting plane on its last iterate;
    // none clamped to continuum: P_epsilon, invBhat feed r
    const mat &L_tan = use_cpp ? L_cpp : L;
    mat kappa, dPhi;
    ContinuumTangent ct;
    if (use_cpp) {
        ct = cpp_tangent(rm, L_cpp, tangent_mode == tangent_closest_point_none ? tangent_closest_point_continuum : tangent_mode);
        kappa = rm.kappa_j;
        dPhi = rm.dPhidsigma_l;
    }
    else {
        kappa = join_rows(kappa_j[0], kappa_j[1]);
        dPhi = join_rows(dPhiFdsigma, dPhiRdsigma);
        mat Bhat = zeros(2, 2);
        Bhat(0,0) = sum(dPhiFdsigma%kappa_j[0]) - K(0,0);
        Bhat(0,1) = sum(dPhiFdsigma%kappa_j[1]) - K(0,1);
        Bhat(1,0) = sum(dPhiRdsigma%kappa_j[0]) - K(1,0);
        Bhat(1,1) = sum(dPhiRdsigma%kappa_j[1]) - K(1,1);
        ct = compute_tangent_operator(tangent_continuum, Bhat, kappa, dPhi, Ds_j, L);
    }
    dSdE = ct.Lt;
    // closest point without a solve (elastic trial, step cut): the elastic operator, zero sensitivities
    const bool solved = !use_cpp || (rm.converged && !rm.kappa_j.is_empty());
    const mat P_epsilon = solved ? ct.P_epsilon : mat(6, 2, fill::zeros);
    const mat invBhat = solved ? ct.invBhat : mat(2, 2, fill::zeros);

    // trial stress rate of a temperature move: the closest point freezes alpha at xi_n and folds
    // Dalpha into the flows, so its thermal strain is alpha(xi_n+1); the algorithmic operator acts through M^-1
    const bool algo = solved && !ct.kappa_tilde.is_empty();
    const vec alpha_q = use_cpp ? vec(alpha + (xi - xi_start)*Dalpha) : alpha;
    const vec L_alpha = algo ? vec(ct.L_tilde*alpha_q) : vec(L_tan*alpha_q);
    const mat kappa_q = !solved ? mat(6, 2, fill::zeros) : (algo ? ct.kappa_tilde : kappa);
    const mat dPhi_q = solved ? dPhi : mat(6, 2, fill::zeros);

    // P_theta = invBhat (dPhi/dtheta - dPhi/dsigma : L alpha)
    const double dA_xiFdtheta = rhoDs0 + sum(sigma%Dalpha);
    const double dA_xiRdtheta = -1.*rhoDs0 - sum(sigma%Dalpha);
    const double rhs0 = dA_xiFdtheta - sum(dPhi_q.col(0)%L_alpha);
    const double rhs1 = dA_xiRdtheta - sum(dPhi_q.col(1)%L_alpha);
    vec P_theta(2);
    P_theta(0) = invBhat(0, 0)*rhs0 + invBhat(0, 1)*rhs1;
    P_theta(1) = invBhat(1, 0)*rhs0 + invBhat(1, 1)*rhs1;
    vec kP = kappa_q.col(0)*P_theta(0);
    kP += kappa_q.col(1)*P_theta(1);
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

    // entropy forces of the transformation rows
    const double dA_theta_dsF = Dc_0*log(T1/T_init) + sum(Dalpha%sigma) + rhoDs0;
    const double dA_theta_dsR = -1.*Dc_0*log(T1/T_init) - sum(Dalpha%sigma) - rhoDs0;

    const vec DET = ET - ET_start;
    const double DxiF = Ds_j(0);
    const double DxiR = Ds_j(1);
    const double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DET) + 0.5*(A_xiF_start + A_xiF)*DxiF + 0.5*(A_xiR_start + A_xiR)*DxiR;

    // dissipation derivatives, flow directions held fixed over the increment
    vec dD = dSdE.t()*DET;
    vec gF = dSdE.t()*dA_xiFdsigma;
    gF += dA_xiFdxiF*P_epsilon.col(0);
    gF += dA_xiFdxiR*P_epsilon.col(1);
    dD += gF*DxiF;
    vec gR = dSdE.t()*dA_xiRdsigma;
    gR += dA_xiRdxiF*P_epsilon.col(0);
    gR += dA_xiRdxiR*P_epsilon.col(1);
    dD += gR*DxiR;
    const double g_TF = dA_xiFdtheta + sum(dA_xiFdsigma%dSdT) + dA_xiFdxiF*P_theta(0) + dA_xiFdxiR*P_theta(1);
    const double g_TR = dA_xiRdtheta + sum(dA_xiRdsigma%dSdT) + dA_xiRdxiF*P_theta(0) + dA_xiRdxiR*P_theta(1);
    double f_F = 0.5*(sum((sigma_start+sigma)%lambdaTF) + A_xiF_start + A_xiF);
    const double f_R = 0.5*(sum((sigma_start+sigma)%lambdaTR) + A_xiR_start + A_xiR);
    double dD_T_turn = 0.;
    // closest point: the turning of the flows and their multiplier chain, d DET = sum_k (Lambda_k dDs_k + Ds_k dLambda_k)
    if (!dflow_dsigma.is_empty()) {
        vec turn = zeros(6);
        turn += DxiF*(dflow_dsigma.cols(0, 5).t()*(sigma_start+sigma));
        turn += DxiR*(dflow_dsigma.cols(6, 11).t()*(sigma_start+sigma));
        dD += dSdE.t()*turn;
        f_F += 0.5*sum((sigma_start+sigma)%chain_F);
        dD_T_turn = sum(dSdT%turn);
    }
    vec dDgamma_dE = 0.5*dD;
    dDgamma_dE += f_F*P_epsilon.col(0);
    dDgamma_dE += f_R*P_epsilon.col(1);
    double dD_T = sum(dSdT%DET) + g_TF*DxiF + g_TR*DxiR;
    dD_T += dD_T_turn;
    const double dDgamma_dT = 0.5*dD_T + f_F*P_theta(0) + f_R*P_theta(1);

    // Entropy derivatives at the end of the increment
    const vec deta_dE = dSdE.t()*alpha + dA_theta_dsF*P_epsilon.col(0) + dA_theta_dsR*P_epsilon.col(1);
    const double deta_dT = c_0/T1 + sum(dSdT%alpha) + dA_theta_dsF*P_theta(0) + dA_theta_dsR*P_theta(1);

    // heat of the increment and its derivatives (heat_source_increment)
    heat_source_increment(Dgamma_loc, dDgamma_dE, dDgamma_dT, Deta, deta_dE, deta_dT, T1, r, drdE, drdT);

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DET)) - 0.5*(A_xiF_start + A_xiF)*DxiF - 0.5*(A_xiR_start + A_xiR)*DxiR;
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;

    Wt += T1*Deta;
    Wt_r += T1*Deta;
    Wt_ir += 0.;

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
