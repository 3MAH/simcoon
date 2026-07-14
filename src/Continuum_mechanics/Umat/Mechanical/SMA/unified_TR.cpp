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

///@file unified_TR.cpp
///@brief Unified SMA model with transformation AND martensite reorientation
///@brief Based on Chatziathanasiou Ph.D Thesis (2016)
///@brief Implemented by Y. Chemisky and D. Chatziathanasiou
///@brief
///@brief The variant is determined by umat_name:
///@brief   SMRDI: isotropic elasticity, Drucker criteria
///@brief   SMRDC: cubic elasticity, Drucker criteria
///@brief   SMRAI: isotropic elasticity, anisotropic Drucker (DFA) criteria
///@brief   SMRAC: cubic elasticity, anisotropic Drucker (DFA) criteria
///@brief
///@brief Implemented in 1D-2D-3D.

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/recovery_props.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

void umat_sma_unified_TR(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode) {

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);

    ///@brief Determine elasticity type and criteria type from umat_name
    bool cubic_elasticity = false;
    bool aniso_criteria = false;

    if (umat_name == "SMRDI") {
        cubic_elasticity = false;
        aniso_criteria = false;
    }
    else if (umat_name == "SMRDC") {
        cubic_elasticity = true;
        aniso_criteria = false;
    }
    else if (umat_name == "SMRAI") {
        cubic_elasticity = false;
        aniso_criteria = true;
    }
    else if (umat_name == "SMRAC") {
        cubic_elasticity = true;
        aniso_criteria = true;
    }
    else {
        throw std::invalid_argument("Unknown umat_name in umat_sma_unified_TR: '" + umat_name
                                    + "'. Valid options: SMRDI, SMRDC, SMRAI, SMRAC.");
    }

    ///@brief Property offset depends on elastic symmetry type
    int offset = 0;
    int flagT = props(0);

    mat M_A;
    mat M_M;

    if (!cubic_elasticity) {
        M_A = inv(L_iso(props(1), props(3), "Enu"));
        M_M = inv(L_iso(props(2), props(4), "Enu"));
        offset = 5;
    }
    else {
        M_A = inv(L_cubic(props(1), props(3), props(5), "EnuG"));
        M_M = inv(L_cubic(props(2), props(4), props(6), "EnuG"));
        offset = 7;
    }

    // Common transformation parameters (mirror umat_sma_unified_T)
    double alphaA_iso = props(offset);
    double alphaM_iso = props(offset + 1);
    double Hmin = props(offset + 2);
    double Hmax = props(offset + 3);
    double k1 = props(offset + 4);
    double sigmacrit = props(offset + 5);
    double C_A = props(offset + 6);
    double C_M = props(offset + 7);
    double Ms0 = props(offset + 8);
    double Mf0 = props(offset + 9);
    double As0 = props(offset + 10);
    double Af0 = props(offset + 11);
    double n1 = props(offset + 12);
    double n2 = props(offset + 13);
    double n3 = props(offset + 14);
    double n4 = props(offset + 15);
    double sigmacaliber = props(offset + 16);
    double prager_b = props(offset + 17);
    double prager_n = props(offset + 18);
    double sigmastar = 0.;

    double c_lambda = props(offset + 19);
    double p0_lambda = props(offset + 20);
    double n_lambda = props(offset + 21);
    double alpha_lambda = props(offset + 22);

    int aniso_offset = offset + 23;
    vec DFA_params = zeros(7);
    if (aniso_criteria) {
        DFA_params = {props(aniso_offset), props(aniso_offset + 1), props(aniso_offset + 2),
                      props(aniso_offset + 3), props(aniso_offset + 4), props(aniso_offset + 5),
                      props(aniso_offset + 6)};
        aniso_offset += 7;
    }

    // Reorientation parameters (Form B: saturation-coupled Phi_Reo)
    double YReo = props(aniso_offset);
    double HReo = props(aniso_offset + 1);
    double ETRmax = props(aniso_offset + 2);
    double c_lambdaReo = props(aniso_offset + 3);
    double p0_lambdaReo = props(aniso_offset + 4);
    double n_lambdaReo = props(aniso_offset + 5);
    double alpha_lambdaReo = props(aniso_offset + 6);

    ///@brief Statev unpack (mirror unified_T then append reorientation block)
    double T_init = statev(0);
    double xi = statev(1);
    vec ET(6);
    for (int i = 0; i < 6; ++i) ET(i) = statev(2 + i);

    double xiF = statev(8);
    double xiR = statev(9);

    double rhoDs0 = statev(10);
    double rhoDE0 = statev(11);
    double D = statev(12);
    double a1 = statev(13);
    double a2 = statev(14);
    double a3 = statev(15);
    double Y0t = statev(16);

    double pTR = statev(17);
    vec areo(6), EReo(6);
    for (int i = 0; i < 6; ++i) {
        areo(i) = statev(18 + i);
        EReo(i) = statev(24 + i);
    }

    // ######################  Elastic compliance and stiffness #################################
    mat DM = M_M - M_A;
    mat M_eff = xi*M_M + (1. - xi)*M_A;
    L = inv(M_eff);

    vec alpha = (alphaM_iso*xi + alphaA_iso*(1.-xi))*Ith();
    vec Dalpha = (alphaM_iso - alphaA_iso)*Ith();

    ///@brief Initialisation block (start=true)
    if(start) {
        T_init = T;
        stress = zeros(6);
        ET = zeros(6);
        EReo = zeros(6);
        areo = zeros(6);
        pTR = simcoon::limit;
        xiF = simcoon::limit;
        xiR = 0.;
        xi = xiF;

        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;

        if (sigmacaliber > sigmacrit) sigmastar = sigmacaliber - sigmacrit;
        else                          sigmastar = 0.;

        double Hcurstar = Hmin + (Hmax - Hmin)*(1. - exp(-k1*sigmastar));
        if (Hcurstar <= 1E-12) {
            throw std::runtime_error("umat_sma_unified_TR: degenerate Hcurstar <= 0 at "
                                     "initialisation; check Hmin/Hmax/k1/sigmacaliber parameters.");
        }

        double dHcurstar = (Hmax - Hmin)*k1*exp(-k1*sigmastar);
        rhoDs0 = -2.*C_M*C_A*(Hcurstar + sigmacaliber*(dHcurstar + (M_M(0,0)-M_A(0,0))))/(C_M+C_A);

        double MsSmooth = 0., MfSmooth = 0., AsSmooth = 0., AfSmooth = 0.;
        if(flagT == 0) {
            MsSmooth = 0.5*Ms0*(1.+ (n1+1)*pow(2.,-n1) + (n2-1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2)) + 0.5*Mf0*(-1.+ (n1-1)*pow(2.,-n1) + (n2+1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2));
            MfSmooth = 0.5*Ms0*(-1.+ (n1+1)*pow(2.,-n1) + (n2-1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2)) + 0.5*Mf0*(1.+ (n1-1)*pow(2.,-n1) + (n2+1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2));
            AsSmooth = 0.5*As0*(1.+ (n3-1)*pow(2.,-n3) + (n4+1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4)) + 0.5*Af0*(-1.+ (n3+1)*pow(2.,-n3) + (n4-1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4));
            AfSmooth = 0.5*As0*(-1.+ (n3-1)*pow(2.,-n3) + (n4+1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4)) + 0.5*Af0*(1.+ (n3+1)*pow(2.,-n3) + (n4-1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4));
        }
        else {
            MsSmooth = Ms0; MfSmooth = Mf0; AsSmooth = As0; AfSmooth = Af0;
        }

        rhoDE0 = 0.5*rhoDs0*(MsSmooth + AfSmooth);
        D = (C_M-C_A)*(Hcurstar + sigmacaliber*(dHcurstar + (M_M(0,0)-M_A(0,0))))/((C_A+C_M)*(Hcurstar + sigmacaliber*dHcurstar));
        a1 = rhoDs0*(MfSmooth - MsSmooth);
        a2 = rhoDs0*(AsSmooth - AfSmooth);
        a3 = -0.25*a1*(1 + 1./(n1 + 1.) - 1./(n2 + 1.))+0.25*a2*(1. + 1./(n3 + 1.) - 1./(n4 + 1.));
        Y0t = 0.5*rhoDs0*(MsSmooth - AfSmooth) - a3;
    }

    // Rotate strain-like internal variables. The back-strain accumulator areo
    // store the same flow direction as the transformation/reorientation flow \Lambda_\bullet 
    // (gradient of Prager/Drucker stress potential), which has engineering shear
    // convention — they MUST rotate as strain (legacy SmartPlus used rotate_stress
    // by mistake, producing factor-2 errors on shears after any finite rotation).
    rotate_strain(ET, DR);
    rotate_strain(EReo, DR);
    rotate_strain(areo, DR);

    // Variables values at the start of the increment
    vec stress_start = stress;
    vec ET_start = ET;
    double xiF_start = xiF, xiR_start = xiR;
    // xi_start is the martensite volume fraction available for reorientation at the
    // start of this increment. It is the "existing" martensite that can be reoriented;
    // any new martensite formed by forward transformation in this increment forms
    // already aligned with \sigma and does not need reorienting. xi_start is constant
    // during the local Newton loop -> no chain-rule contributions to K(2,F) or K(2,R).
    const double xi_start = xi;

    // xi-floor for divisions
    auto xi_safe = [](double x) { return std::max(x, simcoon::limit); };

    // Hcur explicit
    if (Mises_stress(stress) > sigmacrit) sigmastar = Mises_stress(stress) - sigmacrit;
    else                                  sigmastar = 0.;
    double Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

    // Forward flow direction
    auto build_lambdaTF = [&](const vec& s) -> vec {
        return aniso_criteria
            ? Hcur * dDrucker_ani_stress(s, DFA_params, prager_b, prager_n)
            : Hcur * dDrucker_stress(s, prager_b, prager_n);
    };
    auto build_PhihatF = [&](const vec& s) -> double {
        return aniso_criteria
            ? Hcur * Drucker_ani_stress(s, DFA_params, prager_b, prager_n)
            : Hcur * Drucker_stress(s, prager_b, prager_n);
    };

    vec lambdaTF = build_lambdaTF(stress);

    vec ETMean = zeros(6);
    if (Mises_strain(ET) > simcoon::precision_umat) ETMean = dev(ET) / xi_safe(xi);
    else if (Mises_stress(stress) < 1.E-6)          ETMean = lambdaTF;
    else                                            ETMean = 0.*Ith();

    vec lambdaTR = -1.*ETMean;

    // Back-stress + reorientation flow.
    // Per Chatziathanasiou §4.3, the only kinematic back-strain in the model
    // is v^re (= `areo`). aF/aR were SmartPlus legacy artifacts and have been
    // removed: forward and reverse transformation use isotropic hardening only.
    vec X = HReo * (areo % Ir05());

    double a_eq = Mises_strain(areo);
    double lambda1Reo = lagrange_pow_1(a_eq/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
    vec sigma_eff = stress - (1. + lambda1Reo) * X;

    vec lambdaReo;
    vec etaReo;
    if (Mises_stress(sigma_eff) > simcoon::iota) {
        if (aniso_criteria) {
            lambdaReo = dDrucker_ani_stress(sigma_eff, DFA_params, prager_b, prager_n);
        } else {
            lambdaReo = dDrucker_stress(sigma_eff, prager_b, prager_n);
        }
        etaReo = lambdaReo;
    } else {
        lambdaReo = zeros(6);
        etaReo = zeros(6);
    }

    // Back-strain channel flows: each shadows its own leading-mechanism flow
    vec lambda_areo = lambdaReo;

    // Hardening function definition (Smooth hardening functions) — same as unified_T
    auto compute_Hf = [&](double xi_, double& HfF_, double& HfR_) {
        if ((n1==1.)&&(n2==1.)) {
            HfF_ = a1*xi_ + a3;
        } else {
            if ((xi_ > 0.) && ((1. - xi_) > 0.))   HfF_ = 0.5*a1*(1. + pow(xi_, n1) - pow(1. - xi_, n2)) + a3;
            else if (xi_ <= 0.)                     HfF_ = 0.5*a1-1.+a3;
            else                                     HfF_ = a1+a3;
        }
        if ((n3==1.)&&(n4==1.)) {
            HfR_ = a2*xi_ - a3;
        } else {
            if ((xi_ > 0.) && ((1. - xi_) > 0.))   HfR_ = 0.5*a2*(1. + pow(xi_, n3) - pow((1. - xi_), n4)) - a3;
            else if (xi_ <= 0.)                     HfR_ = 0.5*a2-1.-a3;
            else                                     HfR_ = a2-a3;
        }
    };
    double HfF = 0., HfR = 0.;
    compute_Hf(xi, HfF, HfR);

    double lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);

    vec DM_sig = DM*stress_start;
    vec Dalpha_T = Dalpha*(T+DT);

    double A_xiF       = rhoDs0*(T+DT) - rhoDE0 + 0.5*sum(stress%DM_sig)        + sum(stress%Dalpha)*(T+DT-T_init) - HfF;
    double A_xiF_start = rhoDs0*(T)    - rhoDE0 + 0.5*sum(stress_start%DM_sig)  + sum(stress_start%Dalpha)*(T-T_init) - HfF;
    double A_xiR       = -1.*rhoDs0*(T+DT) + rhoDE0 - 0.5*sum(stress%DM_sig)       - sum(stress%Dalpha)*(T+DT-T_init) + HfR;
    double A_xiR_start = -1.*rhoDs0*(T)    + rhoDE0 - 0.5*sum(stress_start%DM_sig) - sum(stress_start%Dalpha)*(T-T_init) + HfR;

    double PhihatF = build_PhihatF(stress);
    double PhihatR = sum(stress%ETMean);

    // Variables required for the local loop — 3 leading mechanisms
    vec s_j = zeros(3);
    s_j(0) = xiF; s_j(1) = xiR; s_j(2) = pTR;
    vec Ds_j = zeros(3);
    vec ds_j = zeros(3);

    // Increment accumulators (signed to match the energy bookkeeping convention)
    vec DETF = zeros(6);
    vec DETR = zeros(6);  // stored as positive: DETR = +\sum ds_R \cdot ETMean
    vec DEReo = zeros(6);

    // Elastic prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - ET;
    stress = el_pred(L, Eel, ndi);

    // Residual / Jacobian containers
    vec Phi = zeros(3);
    mat B = zeros(3,3);
    vec Y_crit = zeros(3);

    double dHfF = 0., dHfR = 0.;
    vec dHcurdsigma = zeros(6);

    // Forward
    vec dPhihatFdsigma = zeros(6);
    vec dA_xiFdsigma   = zeros(6);
    vec dlambda1dsigma = zeros(6);
    vec dYtFdsigma     = zeros(6);
    vec dPhiFdsigma    = zeros(6);
    double dPhiFdxiF = 0., dPhiFdxiR = 0., dPhiFdpTR = 0.;

    // Reverse — only finite \sigma -gradients survive; K(1,\cdot ) uses the ETMean chain rule.
    vec dPhihatRdsigma = zeros(6);
    vec dA_xiRdsigma = zeros(6);
    double dA_xiRdxiF = 0., dA_xiRdxiR = 0.;
    vec dlambda0dsigma = zeros(6);
    double dlambda0dxiF = 0., dlambda0dxiR = 0.;
    vec dYtRdsigma = zeros(6);
    vec dPhiRdsigma = zeros(6);
    double dPhiRdpTR = 0.;

    // Reorientation — aF/aR are not part of the Chatziathanasiou framework
    // (only v^re = areo). dPhiReoda{F,R} and dPhiReodxi{F,R} are zero by
    // construction with the thesis-aligned X = HReo\cdot areo\cdot Ir05.
    double dlambda1Reo_arg = 0.;     // d\lambda_1Reo/d(a_eq/ETRmax)
    vec dPhiReodsigma = zeros(6);
    vec dPhiReodareo = zeros(6);
    double dPhiReodpTR = 0.;

    std::vector<vec> kappa_j(3);
    mat K = zeros(3,3);

    int compteur = 0;
    double error = 1.;

    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        // Reuss mixing
        M_eff = xi*M_M + (1. - xi)*M_A;
        L = inv(M_eff);

        DM_sig   = DM*stress;
        Dalpha_T = Dalpha*(T+DT);

        // Hcur explicit (depends on \sigma )
        if (Mises_stress(stress) > sigmacrit) sigmastar = Mises_stress(stress) - sigmacrit;
        else                                  sigmastar = 0.;
        Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

        // Re-evaluate flows at current \sigma , areo, \lambda_1Reo
        lambdaTF = build_lambdaTF(stress);
        lambdaTR = -1. * ETMean;

        // Back-stress for the reorientation surface — see the framework
        // rationale at the first computation of X above.
        X = HReo * (areo % Ir05());
        a_eq = Mises_strain(areo);
        // Saturation argument is the *intrinsic* (per-unit-martensite) ratio
        // ||v^re|| / ETRmax, NOT ||v^re|| / (\xi \cdot ETRmax). v^re is the per-unit back-strain
        // (its evolution \dot{v}^re = \dot{p}^re \cdot \Lambda_V^re is intrinsic per the thesis); the
        // macroscopic \varepsilon^re reads \xi \cdot \bar{\varepsilon}^re. Putting \xi in the denominator of the
        // saturation argument (legacy SmartPlus form) makes the Lagrange penalty
        // blow up when \xi -> 0 even if ||v^re|| is tiny, producing astronomical
        // \sigma_eff = \sigma - (1+1e20)\cdot X garbage.
        lambda1Reo = lagrange_pow_1(a_eq/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
        sigma_eff = stress - (1. + lambda1Reo) * X;

        // Raw (per-unit-martensite) reorientation flow: \bar{\Lambda}^re = dDrucker(\sigma_eff).
        vec lambdaReo_raw;
        if (Mises_stress(sigma_eff) > simcoon::iota) {
            if (aniso_criteria) {
                lambdaReo_raw = dDrucker_ani_stress(sigma_eff, DFA_params, prager_b, prager_n);
            } else {
                lambdaReo_raw = dDrucker_stress(sigma_eff, prager_b, prager_n);
            }
        } else {
            lambdaReo_raw = zeros(6);
        }
        // Thesis §4 (eq 503): the macroscopic flow contributing to \varepsilon^T from
        // reorientation is \xi -scaled — \Lambda_\varepsilon^re = \xi \cdot \bar{\Lambda}^re — because it only
        // moves the variant orientation in the pre-existing martensite of
        // volume \xi . The back-strain flow \Lambda_V^re = \bar{\Lambda}^re is NOT \xi -scaled: the
        // back-stress tracks the intrinsic reorientation history within the
        // martensite domain, evolving per unit of \dot{p}^re regardless of how much
        // martensite there is. Criterion gradient \partial Phi^re/\partial \sigma also stays
        // un-scaled (it's just the stress gradient of the yield function).
        // This asymmetric structure keeps B(2,2) finite at \xi =0 through the
        // K(2,2) entry from the back-strain channel.
        lambdaReo   = xi_start * lambdaReo_raw;   // macroscopic flow into \varepsilon^T
        etaReo      = lambdaReo_raw;              // \partial Phi^re/\partial \sigma — un-scaled
        lambda_areo = lambdaReo_raw;              // back-strain v^re evolution — un-scaled

        kappa_j[0] = L*(lambdaTF  + DM_sig + Dalpha_T);
        kappa_j[1] = L*(lambdaTR  - DM_sig - Dalpha_T);
        kappa_j[2] = L*(lambdaReo);   // \xi -scaled — vanishes at \xi =0 as it should

        // Hardening
        compute_Hf(xi, HfF, HfR);

        // Forward thermodynamic force
        PhihatF = build_PhihatF(stress);
        A_xiF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(stress%DM_sig) + sum(stress%Dalpha)*(T + DT - T_init) - HfF;
        lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        double YtF = Y0t + D*Hcur*Mises_stress(stress);
        Phi(0) = PhihatF + A_xiF - lambda1 - YtF;

        // Reverse thermodynamic force
        PhihatR = sum(stress%ETMean);
        A_xiR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(stress%DM_sig) - sum(stress%Dalpha)*(T + DT - T_init) + HfR;
        lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        double YtR = Y0t + D*sum(stress%ETMean);
        Phi(1) = -1.*PhihatR + A_xiR + lambda0 - YtR;

        // Reorientation criterion (Chatziathanasiou thesis eq 664, Form B):
        //   Phi^re = Drucker(\sigma - (1+\lambda_1Reo)\cdot X) - Y^re        (UN-scaled in \xi )
        // No \xi multiplier on the criterion — the yield condition is per unit
        // martensite volume. The \xi -scaling lives in the macroscopic flow
        // (\Lambda_\varepsilon^re = \xi \cdot \bar{\Lambda}^re, applied above to lambdaReo and kappa_j[2]).
        // At \xi_start = 0, \kappa_2 = 0 so Phi^re cannot be driven down through \sigma ,
        // but it can still be driven down by growing v^re (the un-scaled
        // back-strain channel) — which is the correct physical statement:
        // the kinematic back-strain tracks the intrinsic reorientation
        // history within the martensite domain even when no macroscopic strain
        // is yet produced. B(2,2) stays finite via the back-strain entry in K(2,2).
        double Prager_sigma_eff = aniso_criteria
            ? Drucker_ani_stress(sigma_eff, DFA_params, prager_b, prager_n)
            : Drucker_stress(sigma_eff, prager_b, prager_n);
        Phi(2) = Prager_sigma_eff - YReo;

        // Hardening derivatives w.r.t. xi
        if ((n1==1.)&&(n2==1.))                       dHfF = a1;
        else if ((xi > 0.) && ((1. - xi) > 0.))        dHfF = 0.5*a1*(n1*pow(xi, n1 - 1.) + n2*pow(1. - xi, n2 - 1.));
        else                                            dHfF = 1.E12;
        if ((n3==1.)&&(n4==1.))                       dHfR = a2;
        else if ((xi > 0.) && ((1. - xi) > 0.))        dHfR = 0.5*a2*(n3*pow(xi, n3 - 1.) + n4*pow(1. - xi, n4 - 1.));
        else                                            dHfR = 1.E12;

        // d(Hcur)/d\sigma 
        dHcurdsigma = k1*(Hmax - Hmin)*exp(-1.*k1*sigmastar)*eta_stress(stress);

        // ---- Forward derivatives ----
        if (aniso_criteria) {
            dPhihatFdsigma = dHcurdsigma * Drucker_ani_stress(stress, DFA_params, prager_b, prager_n) + Hcur * dDrucker_ani_stress(stress, DFA_params, prager_b, prager_n);
        } else {
            dPhihatFdsigma = dHcurdsigma * Drucker_stress(stress, prager_b, prager_n) + Hcur * dDrucker_stress(stress, prager_b, prager_n);
        }
        dA_xiFdsigma = DM_sig + Dalpha*(T+DT);
        dlambda1dsigma = zeros(6);
        dYtFdsigma = D*(dHcurdsigma * Mises_stress(stress) + Hcur * eta_stress(stress));

        dPhiFdsigma = dPhihatFdsigma + dA_xiFdsigma - dlambda1dsigma - dYtFdsigma;
        double dlambda1dxi = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dPhiFdxiF = -dHfF - dlambda1dxi;
        dPhiFdxiR =  dHfF + dlambda1dxi;
        dPhiFdpTR = 0.;

        // ---- Reverse derivatives ----
        const double xi_div = xi_safe(xi);
        dPhihatRdsigma = ETMean;
        dA_xiRdsigma   = -1.*DM_sig - 1.*Dalpha*(T+DT-T_init);
        double dlambda0dxi = dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda0dsigma = zeros(6);
        dlambda0dxiF = -1.*dlambda0dxi;
        dlambda0dxiR =     dlambda0dxi;
        dYtRdsigma = D*ETMean;
        dPhiRdsigma = -1.*dPhihatRdsigma + dA_xiRdsigma + dlambda0dsigma - dYtRdsigma;
        dA_xiRdxiF =  dHfR;
        dA_xiRdxiR = -dHfR;
        dPhiRdpTR = 0.;

        // Phi_R partial w.r.t. ETMean is well-defined and finite:
        // Phi_R = -\sigma :ETMean + A_xiR + \lambda_0 - YtR, with YtR = Y0t + D\cdot \sigma :ETMean.
        // => \partial Phi_R/\partial ETMean = -(1+D)\cdot \sigma 
        const vec dPhiRdETMean = -(1. + D) * stress;
        // "Finite" Phi_R partials w.r.t. \xi (only the hardening + Lagrange-penalty
        // pieces; the \sigma :dev(\varepsilon^T)/\xi = \sigma :ETMean parts that would otherwise blow up
        // as \xi ->0 are absorbed into the \Lambda_ETMean chain-rule contributions below).
        const double dPhiRdxi_finite_F =  dA_xiRdxiF + dlambda0dxiF;   //  +dHfR  - d\lambda_0/d\xi 
        const double dPhiRdxi_finite_R =  dA_xiRdxiR + dlambda0dxiR;   //  -dHfR  + d\lambda_0/d\xi 
        // The legacy dPhiRdET and dPhihatRdxiF symbols are no longer assembled
        // here — K(1,\cdot ) is built directly from dPhiRdETMean \cdot \Lambda_ETMean^j below
        // to avoid the +(1+D)\cdot \sigma :ETMean/\xi - (1+D)\cdot \sigma :\Lambda_\varepsilon^j/\xi = (finite) cancellation
        // that loses ~12 digits at \xi \approx 1e-12.

        // ---- Reorientation derivatives (Form B, consistent linearisation) ----
        // Saturation argument: ||v^re||/ETRmax (intrinsic, no \xi ).
        dlambda1Reo_arg = dlagrange_pow_1(a_eq/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
        const double eta_dot_X = sum(etaReo % X);                  // \eta_Reo : X (tensor double contraction)
        const vec   eta_strain_areo = (a_eq > simcoon::iota) ? eta_strain(areo) : zeros(6);

        // d\lambda_1Reo/dareo = d\lambda_1Reo/d(arg) \cdot (1/ETRmax) \cdot \eta_strain(areo)
        const vec dlambda1Reo_dareo = dlambda1Reo_arg * (1./ETRmax) * eta_strain_areo;

        dPhiReodsigma = etaReo;

        // \partial \Phi^re/\partial areo = -(1+\lambda_1)\cdot H^Reo\cdot (\eta_Reo \odot Ir05)  -  (\eta_Reo:X)\cdot d\lambda_1Reo/dareo
        dPhiReodareo = -(1. + lambda1Reo) * HReo * (etaReo % Ir05())
                       - eta_dot_X * dlambda1Reo_dareo;

        // Saturation argument no longer depends on \xi => d\lambda_1Reo/d\xi_{F,R} = 0.
        dPhiReodpTR = 0.;

        // ---- Assemble K (state-coupling block, excluding the \sigma -flow term) ----
        // mech 0 (forward) evolves \xi^F (+1 in \xi ), \varepsilon^F via lambdaTF
        // mech 1 (reverse) evolves \xi^R (-1 in \xi ), \varepsilon^R via lambdaTR
        // mech 2 (reorient.) evolves pTR (+1), \varepsilon^re_macro = \xi \cdot \bar{\Lambda}^re into \varepsilon^T,
        //                    EReo via lambdaReo (= \xi \cdot \bar{\Lambda}^re), areo via \bar{\Lambda}^re
        // \varepsilon^T is shared by all mechanisms; the chain rule for row 1 routes
        // through ETMean as the natural intermediate (see K(1,\cdot ) below).
        K(0,0) = dPhiFdxiF;
        K(0,1) = dPhiFdxiR;
        K(0,2) = dPhiFdpTR;

        // ---- K(1,\cdot ) via \Lambda_ETMean^j chain rule (no \sigma /\xi blow-up) ----
        // For a mechanism j with kinetic flow \Lambda_\varepsilon^j and \xi -growth rate \Lambda_\xi^j,
        // the rate of change of ETMean = dev(\varepsilon^T)/\xi is
        //   \Lambda_ETMean^j = (dev(\Lambda_\varepsilon^j) - \Lambda_\xi^j \cdot ETMean) / \xi 
        // Closed forms:
        //   \Lambda_ETMean^F  = (\Lambda_\varepsilon^F - ETMean)/\xi (\Lambda_\varepsilon^F deviatoric, \Lambda_\xi^F = +1)
        //   \Lambda_ETMean^R  = 0                          (\Lambda_\varepsilon^R = -ETMean, \Lambda_\xi^R = -1: cancels)
        //   \Lambda_ETMean^re = \bar{\Lambda}^re                      (\Lambda_\varepsilon^re_macro = \xi \cdot \bar{\Lambda}^re, \Lambda_\xi^re = 0)
        // The ETMean fallback (= lambdaTF when \varepsilon^T and \xi are both small) makes
        // (\Lambda_\varepsilon^F - ETMean) = 0 identically there, so \Lambda_ETMean^F = 0/0 -> 0 in
        // floating point as long as we don't actually divide by \xi when both
        // numerator and \xi are below their respective floors.
        vec Lambda_ETMean_F = zeros(6);
        if (Mises_strain(lambdaTF - ETMean) > simcoon::iota) {
            Lambda_ETMean_F = (lambdaTF - ETMean) / xi_div;
        }
        const vec Lambda_ETMean_R = zeros(6);
        const vec Lambda_ETMean_re = lambdaReo_raw;
        K(1,0) = dPhiRdxi_finite_F + sum(dPhiRdETMean % Lambda_ETMean_F);
        K(1,1) = dPhiRdxi_finite_R + sum(dPhiRdETMean % Lambda_ETMean_R);
        K(1,2) = dPhiRdpTR        + sum(dPhiRdETMean % Lambda_ETMean_re);

        // K(2,F) = K(2,R) = 0 by construction: \Phi^re depends on internal state
        // only through areo (via X) and \sigma ; aF and aR are not in the model, and
        // \xi no longer appears in the saturation argument after the v^re fix.
        K(2,0) = 0.;
        K(2,1) = 0.;
        K(2,2) = dPhiReodpTR + sum(dPhiReodareo % lambda_areo);

        // ---- Build B = -\partial \Phi :\kappa + K  (Bhat in the helper convention is +\partial \Phi :\kappa - K) ----
        for (int l = 0; l < 3; ++l) {
            const vec& dPhi = (l == 0 ? dPhiFdsigma : (l == 1 ? dPhiRdsigma : dPhiReodsigma));
            for (int j = 0; j < 3; ++j) {
                B(l, j) = -1.*sum(dPhi % kappa_j[j]) + K(l, j);
            }
        }

        Y_crit(0) = Y0t + D*Hcur*Mises_stress(stress);
        Y_crit(1) = Y0t + D*sum(stress%ETMean);
        Y_crit(2) = std::max(YReo, simcoon::iota);   // guard the FB /Y_crit normalization against a (degenerate) YReo = 0

        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);

        s_j(0) += ds_j(0);
        s_j(1) += ds_j(1);
        s_j(2) += ds_j(2);

        // State updates
        ET   = ET   + ds_j(0)*lambdaTF + ds_j(1)*lambdaTR + ds_j(2)*lambdaReo;
        EReo = EReo + ds_j(2)*lambdaReo;
        // aF and aR are not updated (legacy accumulators kept at zero).
        areo = areo + ds_j(2)*lambda_areo;

        xiF = s_j(0);
        xiR = s_j(1);
        pTR = s_j(2);
        xi  = xi + ds_j(0) - ds_j(1);

        DETF  += ds_j(0)*lambdaTF;
        DETR  += -1.*ds_j(1)*lambdaTR;   // DETR stored positive (DETR = +\sum ds_R \cdot ETMean)
        DEReo += ds_j(2)*lambdaReo;

        if ((Mises_strain(ET) > simcoon::precision_umat) && (xi > simcoon::precision_umat)) {
            ETMean = dev(ET) / xi_safe(xi);
        } else {
            ETMean = lambdaTF;
        }

        Eel = Etot + DEtot - alpha*(T + DT - T_init) - ET;
        stress = el_pred(L, Eel, ndi);
    }

    // If the local Newton did not converge within maxiter_umat, request a sub-step
    // (tnew_dt < 1) from the global solver rather than committing a non-converged increment.
    if (error > simcoon::precision_umat) {
        tnew_dt = 0.5;
    }

    // ---- Tangent assembly via shared 3-mechanism helper ----
    mat Bhat = zeros(3, 3);
    for (int l = 0; l < 3; ++l) {
        const vec& dPhi = (l == 0 ? dPhiFdsigma : (l == 1 ? dPhiRdsigma : dPhiReodsigma));
        for (int j = 0; j < 3; ++j) {
            Bhat(l, j) = sum(dPhi % kappa_j[j]) - K(l, j);
        }
    }

    const std::vector<vec> dPhidsigma_l = { dPhiFdsigma, dPhiRdsigma, dPhiReodsigma };
    const ContinuumTangent ct = compute_tangent_operator(
        tangent_mode, Bhat, kappa_j, dPhidsigma_l, Ds_j, L,
        [&]() -> std::vector<mat> {  // lazy: evaluated only in algorithmic mode
            // Simo-Hughes algorithmic tangent (closest-point), 3-mechanism. dLambda/dsigma by central
            // FD of the transformation flows: forward Hcur(sigma)*dDrucker(sigma); reorientation
            // xi_start*dDrucker(sigma_eff), sigma_eff = sigma-(1+lambda1Reo)*X; plus the analytic linear
            // DM term. ETMean, X, lambda1Reo held fixed -> transformation/back-strain state-coupling is
            // the deferred term (closest-point/CPP rework, future release).
            auto lambdaTF_at = [&](const vec &s) -> vec {
                double sstar = Mises_stress(s) - sigmacrit;
                if (sstar < 0.) sstar = 0.;
                double Hc = Hmin + (Hmax - Hmin) * (1. - exp(-1. * k1 * sstar));
                return aniso_criteria ? Hc * dDrucker_ani_stress(s, DFA_params, prager_b, prager_n)
                                      : Hc * dDrucker_stress(s, prager_b, prager_n);
            };
            auto lambdaReo_raw_at = [&](const vec &s) -> vec {
                vec se = s - (1. + lambda1Reo) * X;
                if (Mises_stress(se) > simcoon::iota) {
                    return aniso_criteria ? dDrucker_ani_stress(se, DFA_params, prager_b, prager_n)
                                          : dDrucker_stress(se, prager_b, prager_n);
                }
                return zeros(6);
            };
            const double hfd = 1.e-5 * (norm(stress, 2) + 1.);
            mat dLambdaF = zeros(6, 6), dLambdaReo = zeros(6, 6);
            for (int c = 0; c < 6; c++) {
                vec sp = stress, sm = stress;
                sp(c) += hfd;
                sm(c) -= hfd;
                dLambdaF.col(c)   = (lambdaTF_at(sp) - lambdaTF_at(sm)) / (2. * hfd);
                dLambdaReo.col(c) = xi_start * (lambdaReo_raw_at(sp) - lambdaReo_raw_at(sm)) / (2. * hfd);
            }
            dLambdaF += DM;
            const std::vector<mat> dLambda_dsigma_l = { dLambdaF, -1. * DM, dLambdaReo };
            return dLambda_dsigma_l;
        });
    Lt = ct.Lt;

    // ---- Energy partition (SMA: Wm_ir is always 0; everything is either recoverable or dissipated) ----
    // Back-strain accumulators (aF, aR, areo) are kinematic *history* variables that only
    // shift the reorientation criterion — they do NOT store energy in the free energy \psi , so
    // there are no X:\Delta a contributions. Partition mirrors umat_sma_unified_T extended with
    // the reorientation strain increment.
    double DxiF = Ds_j[0];
    double DxiR = Ds_j[1];

    double Dgamma_loc = 0.5*sum((stress_start+stress) % (DETF - DETR + DEReo))
                      + 0.5*(A_xiF_start + A_xiF)*DxiF
                      + 0.5*(A_xiR_start + A_xiR)*DxiR;

    Wm    += 0.5*sum((stress_start+stress) % DEtot);
    Wm_r  += 0.5*sum((stress_start+stress) % (DEtot - DETF + DETR - DEReo))
           - 0.5*(A_xiF_start + A_xiF)*DxiF
           - 0.5*(A_xiR_start + A_xiR)*DxiR;
    Wm_ir += 0.;   // SMA has no irrecoverable channel by construction
    Wm_d  += Dgamma_loc;

    // ---- Statev write-back ----
    statev(0) = T_init;
    statev(1) = xi;
    for (int i = 0; i < 6; ++i) statev(2 + i) = ET(i);
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

    statev(17) = pTR;
    for (int i = 0; i < 6; ++i) {
        statev(18 + i) = areo(i);
        statev(24 + i) = EReo(i);
    }
}

} //namespace simcoon
