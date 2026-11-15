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
///@brief Thermomechanical face of the unified SMA model with martensite reorientation

#include <iostream>
#include <stdexcept>
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
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

void umat_sma_unified_TR_T(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, vec &dSdT, vec &drdE, double &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode, const int &ncp) {
    const NCPMethod ncp_m = sma_ncp_method(ncp, tangent_mode, "unified_TR_T");

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(ndi);
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
        throw std::invalid_argument("Unknown umat_name in umat_sma_unified_TR_T: '" + umat_name
                                    + "'. Valid options: SMRDI, SMRDC, SMRAI, SMRAC.");
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
        M_A = M_iso(props(4), props(6), "Enu");
        M_M = M_iso(props(5), props(7), "Enu");
        offset = 8;
    }
    else {
        M_A = M_cubic(props(4), props(6), props(8), "EnuG");
        M_M = M_cubic(props(5), props(7), props(9), "EnuG");
        offset = 10;
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
    mat L = inv_cubic(M_eff);

    vec alpha = (alphaM_iso*xi + alphaA_iso*(1.-xi))*Ith();
    vec Dalpha = (alphaM_iso - alphaA_iso)*Ith();

    ///@brief Initialisation block (start=true)
    if(start) {
        T_init = T;
        sigma = zeros(6);
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
            throw std::runtime_error("umat_sma_unified_TR_T: degenerate Hcurstar <= 0 at "
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

    // areo carries engineering shears like the flows: rotate as strain, not sigma
    ET = rotate_strain(ET, DR);
    EReo = rotate_strain(EReo, DR);
    areo = rotate_strain(areo, DR);

    // Variables values at the start of the increment
    vec sigma_start = sigma;
    vec ET_start = ET;
    double xiF_start = xiF, xiR_start = xiR;
    // martensite available for reorientation, constant over the local loop: K(2,F) = K(2,R) = 0
    const double xi_start = xi;

    // xi-floor for divisions
    auto xi_safe = [](double x) { return std::max(x, simcoon::limit); };

    // Hcur explicit
    if (Mises_stress(sigma) > sigmacrit) sigmastar = Mises_stress(sigma) - sigmacrit;
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

    vec lambdaTF = build_lambdaTF(sigma);

    vec ETMean = zeros(6);
    if (Mises_strain(ET) > simcoon::precision_umat) ETMean = dev(ET) / xi_safe(xi);
    else if (Mises_stress(sigma) < 1.E-6)          ETMean = lambdaTF;
    else                                            ETMean = 0.*Ith();

    vec lambdaTR = -1.*ETMean;

    // back-sigma from v^re (areo) only, Chatziathanasiou §4.3
    vec X = HReo * (areo % Ir05());

    double a_eq = Mises_strain(areo);
    double lambda1Reo = lagrange_pow_1(a_eq/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
    vec sigma_eff = sigma - (1. + lambda1Reo) * X;
    const vec S_start = (1. + lambda1Reo) * X;   // back sigma of the reorientation at the start

    // reorientation flow and criterion gradient: set by the local loop before any use
    vec lambdaReo = zeros(6);
    vec etaReo = zeros(6);

    // Back-strain channel flows: each shadows its own leading-mechanism flow
    vec lambda_areo = lambdaReo;

    // Hardening function definition (Smooth hardening functions) — same as unified_T
    // hardening functions at xi clamped to [0,1] (sma_hardening)
    auto compute_Hf = [&](double xi_, double& HfF_, double& HfR_) {
        const SMAHardening h = sma_hardening(xi_, a1, a2, a3, n1, n2, n3, n4);
        HfF_ = h.F;
        HfR_ = h.R;
    };
    double HfF = 0., HfR = 0.;
    compute_Hf(xi, HfF, HfR);

    // Lagrange multipliers and criteria: set by the local loop before any use
    double lambda0 = 0.;
    double lambda1 = 0.;

    vec DM_sig = DM*sigma_start;
    vec Dalpha_T = Dalpha*(T+DT-T_init);

    double A_xiF       = rhoDs0*(T+DT) - rhoDE0 + 0.5*sum(sigma%DM_sig)        + sum(sigma%Dalpha)*(T+DT-T_init) - HfF;
    double A_xiF_start = rhoDs0*(T)    - rhoDE0 + 0.5*sum(sigma_start%DM_sig)  + sum(sigma_start%Dalpha)*(T-T_init) - HfF;
    double A_xiR       = -1.*rhoDs0*(T+DT) + rhoDE0 - 0.5*sum(sigma%DM_sig)       - sum(sigma%Dalpha)*(T+DT-T_init) + HfR;
    double A_xiR_start = -1.*rhoDs0*(T)    + rhoDE0 - 0.5*sum(sigma_start%DM_sig) - sum(sigma_start%Dalpha)*(T-T_init) + HfR;

    double PhihatF = 0.;
    double PhihatR = 0.;

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
    sigma = el_pred(L, Eel);

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

    // Reorientation
    double dlambda1Reo_arg = 0.;     // d\lambda_1Reo/d(a_eq/ETRmax)
    vec dPhiReodsigma = zeros(6);
    vec dPhiReodareo = zeros(6);
    double dPhiReodpTR = 0.;

    std::vector<vec> kappa_j(3);
    mat K = zeros(3,3);

    int compteur = 0;
    double error = 1.;

    // closest point: frozen-mixture construction of unified_T (see header)
    ReturnMappingResult rm;
    // dS/dareo of the back sigma S = (1 + lambda1Reo) X at the current areo
    auto backstress_slope = [&]() -> mat {
        const vec T_sh = Ir05();
        const double a_eq_cur = Mises_strain(areo);
        const vec dl1Reo = (dlagrange_pow_1(a_eq_cur/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo)/ETRmax)
                           *(eta_strain(areo) % T_sh);
        mat dS = HReo*(areo % T_sh)*dl1Reo.t();
        dS.diag() += HReo*(1. + lagrange_pow_1(a_eq_cur/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo))*T_sh;
        return dS;
    };
    // closest point, for the heat: total dLambda^j/dsigma of the flows of the transformation strain (6x18) and their
    // multiplier chain (6x3); sensitivities of the work of the back-strain force at fixed multipliers
    // (closest point: exact; cutting plane: flows held fixed)
    mat dflow_dsigma, flow_chain;
    vec dWReo_dsigma, dWReo_ds;

    const bool use_cpp = integrates_closest_point(tangent_mode);
    const mat L_cpp = use_cpp ? L : mat();   // frozen stiffness used consistently by the helper and the tangent
    if (use_cpp) {
        const double xi_n = xi;
        const vec EReo_n = EReo, areo_n = areo;
        const double thp = (T + DT - T_init);
        vec lambdaReo_raw = zeros(6);
        vec Dl_cur = zeros(3);        // multipliers of the refreshed state
        bool ETMean_is_flow = false;  // the guard branch ETMean = lambdaTF is active
        bool lin_stale = true;        // the linearisation below is not the one of the refreshed state
        const vec T_shear = Ir05();   // engineering -> tensorial shear

        // r = dDrucker(sigma_eff), its Hessian, and dsigma_eff/dareo
        auto reo_direction = [&](const vec &seff) -> vec {
            if (!(Mises_stress(seff) > simcoon::iota)) return zeros(6);
            return aniso_criteria ? dDrucker_ani_stress(seff, DFA_params, prager_b, prager_n)
                                  : dDrucker_stress(seff, prager_b, prager_n);
        };
        auto reo_hessian = [&](const vec &seff) -> mat {
            if (!(Mises_stress(seff) > simcoon::iota)) return zeros(6, 6);
            return aniso_criteria ? ddDrucker_ani_stress(seff, DFA_params, prager_b, prager_n)
                                  : ddDrucker_stress(seff, prager_b, prager_n);
        };
        auto dsigma_eff_dareo = [&]() -> mat { return -1.*backstress_slope(); };
        auto set_sigma_eff = [&](const vec &sig) {
            X = HReo * (areo % T_shear);
            lambda1Reo = lagrange_pow_1(Mises_strain(areo)/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
            sigma_eff = sig - (1. + lambda1Reo) * X;
        };

        auto refresh_state = [&](const vec &sig, const vec &Dl) {
            Dl_cur = Dl;
            lin_stale = true;
            xi = xi_n + Dl(0) - Dl(1);
            // Hcur, forward flow (explicit in sigma)
            Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*std::max(Mises_stress(sig) - sigmacrit, 0.)));
            lambdaTF = build_lambdaTF(sig);
            // back-strain by Newton (a sweep does not contract); no root rejects the trial
            areo = areo_n;
            set_sigma_eff(sig);
            lambdaReo_raw = reo_direction(sigma_eff);
            if (Dl(2) > 0.) {
                bool areo_converged = false;
                for (int it_nr = 0; it_nr < 25; it_nr++) {
                    const vec R_a = areo - areo_n - Dl(2)*lambdaReo_raw;
                    if (norm(R_a, 2) <= 1.e-13*(norm(areo, 2) + Dl(2))) { areo_converged = true; break; }
                    vec d_a;
                    if (!solve(d_a, eye(6, 6) - Dl(2)*reo_hessian(sigma_eff)*dsigma_eff_dareo(), -1.*R_a, solve_opts::fast)) return false;
                    areo += d_a;
                    set_sigma_eff(sig);
                    lambdaReo_raw = reo_direction(sigma_eff);
                }
                if (!areo_converged || !areo.is_finite()) return false;
            }
            lambdaReo = xi_start * lambdaReo_raw;
            EReo = EReo_n + Dl(2)*lambdaReo;
            // ETMean in closed form: the sweep diverges when Dl_R > xi
            const double xi_before = xi + Dl(1);
            ETMean_is_flow = !(xi_before > simcoon::precision_umat);
            ETMean = ETMean_is_flow ? lambdaTF
                : vec((dev(ET_start) + Dl(0)*dev(lambdaTF) + Dl(2)*dev(lambdaReo))/xi_before);
            ET = ET_start + Dl(0)*lambdaTF - Dl(1)*ETMean + Dl(2)*lambdaReo;
            if (!((Mises_strain(ET) > simcoon::precision_umat) && (xi > simcoon::precision_umat))) {
                ETMean_is_flow = true;
                ETMean = lambdaTF;
                ET = ET_start + Dl(0)*lambdaTF - Dl(1)*ETMean + Dl(2)*lambdaReo;
            }
            compute_Hf(xi, HfF, HfR);
            lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
            lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
            return true;
        };
        auto PhiF_val = [&](const vec &sig) {
            const vec DMsig = DM*sig;
            const double AxF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sig%DMsig) + sum(sig%Dalpha)*thp - HfF;
            return build_PhihatF(sig) + AxF - lambda1 - (Y0t + D*Hcur*Mises_stress(sig));
        };
        auto PhiR_val = [&](const vec &sig) {
            const vec DMsig = DM*sig;
            const double AxR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sig%DMsig) - sum(sig%Dalpha)*thp + HfR;
            return -1.*sum(sig%ETMean) + AxR + lambda0 - (Y0t + D*sum(sig%ETMean));
        };
        auto PhiReo_val = [&](const vec &) {
            const double Pse = aniso_criteria ? Drucker_ani_stress(sigma_eff, DFA_params, prager_b, prager_n)
                                              : Drucker_stress(sigma_eff, prager_b, prager_n);
            return Pse - YReo;
        };
        auto LamF_val   = [&](const vec &sig) { return vec(lambdaTF + DM*sig + Dalpha*thp); };
        auto LamR_val   = [&](const vec &sig) { return vec(-1.*ETMean - DM*sig - Dalpha*thp); };
        auto LamReo_val = [&](const vec &) { return lambdaReo; };
        // analytic total derivatives at (sig, Dl_cur), once per iterate (derivation in header)
        mat dev_op = eye(6, 6);   // dev() of a Voigt vector
        dev_op.submat(0, 0, 2, 2) -= ones(3, 3)/3.;
        vec dPhiF_ds = zeros(6), dPhiR_ds = zeros(6), dPhiReo_ds = zeros(6), c_F = zeros(6), c_Reo = zeros(6);
        mat dLamF_ds = zeros(6, 6), dLamR_ds = zeros(6, 6), dLamReo_ds = zeros(6, 6), K_cpp = zeros(3, 3);
        mat dseff_ds_l = eye(6, 6);
        vec dseff_dre_l = zeros(6), dr_dre_l = zeros(6);
        auto linearise = [&](const vec &sig) {
            if (!lin_stale) return;
            lin_stale = false;
            const double sig_eq = Mises_stress(sig);
            const vec eta = eta_stress(sig);
            vec dHcur = zeros(6);
            if (sig_eq > sigmacrit) dHcur = (Hmax - Hmin)*k1*exp(-1.*k1*(sig_eq - sigmacrit))*eta;
            const double Dr = aniso_criteria ? Drucker_ani_stress(sig, DFA_params, prager_b, prager_n)
                                             : Drucker_stress(sig, prager_b, prager_n);
            const vec g = aniso_criteria ? dDrucker_ani_stress(sig, DFA_params, prager_b, prager_n)
                                         : dDrucker_stress(sig, prager_b, prager_n);
            const mat ddDr = aniso_criteria ? ddDrucker_ani_stress(sig, DFA_params, prager_b, prager_n)
                                            : ddDrucker_stress(sig, prager_b, prager_n);
            const mat dlTF_ds = g*dHcur.t() + Hcur*ddDr;

            // reorientation chain
            const vec &r = lambdaReo_raw;
            const mat G = reo_hessian(sigma_eff);
            const mat A = dsigma_eff_dareo();
            const mat Jinv = inv(eye(6, 6) - Dl_cur(2)*(G*A));
            const mat dseff_ds = eye(6, 6) + A*(Jinv*(Dl_cur(2)*G));
            const vec dseff_dre = A*(Jinv*r);
            const mat dr_ds = G*dseff_ds;
            const vec dr_dre = G*dseff_dre;
            dseff_ds_l = dseff_ds;
            dseff_dre_l = dseff_dre;
            dr_dre_l = dr_dre;

            mat dETM_ds = dlTF_ds;
            vec dETM_dF = zeros(6), dETM_dre = zeros(6);
            if (!ETMean_is_flow) {
                const double xi_before = xi_n + Dl_cur(0);
                dETM_ds = (dev_op*(Dl_cur(0)*dlTF_ds + (Dl_cur(2)*xi_start)*dr_ds))/xi_before;
                dETM_dF = (dev(lambdaTF) - ETMean)/xi_before;
                dETM_dre = (dev_op*(lambdaReo + (Dl_cur(2)*xi_start)*dr_dre))/xi_before;
            }

            // slopes of the clamped hardening functions and of the penalties
            const SMAHardening hs = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
            const double dHfF_dxi = hs.dF, dHfR_dxi = hs.dR;
            const double dl1_dxi = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
            const double dl0_dxi = -1.*dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);

            const vec DMsig = DM*sig;
            dPhiF_ds = Dr*dHcur + Hcur*g + DMsig + Dalpha*thp - D*(sig_eq*dHcur + Hcur*eta);
            dPhiR_ds = -1.*(1. + D)*(ETMean + dETM_ds.t()*sig) - DMsig - Dalpha*thp;
            dPhiReo_ds = dseff_ds.t()*r;
            dLamF_ds = dlTF_ds + DM;
            dLamR_ds = -1.*dETM_ds - DM;
            dLamReo_ds = xi_start*dr_ds;
            K_cpp.zeros();
            K_cpp(0, 0) = -1.*(dHfF_dxi + dl1_dxi);
            K_cpp(0, 1) = dHfF_dxi + dl1_dxi;
            K_cpp(1, 0) = -1.*(1. + D)*sum(sig%dETM_dF) + dHfR_dxi + dl0_dxi;
            K_cpp(1, 1) = -1.*(dHfR_dxi + dl0_dxi);
            K_cpp(1, 2) = -1.*(1. + D)*sum(sig%dETM_dre);
            K_cpp(2, 2) = sum(r%dseff_dre);
            c_F = -1.*Dl_cur(1)*dETM_dF;
            c_Reo = -1.*Dl_cur(1)*dETM_dre + (Dl_cur(2)*xi_start)*dr_dre;
        };
        std::vector<ReturnMechanism> mechs(3);
        mechs[0].Phi = PhiF_val;   mechs[1].Phi = PhiR_val;   mechs[2].Phi = PhiReo_val;
        mechs[0].Lambda = [&](const vec &sig, vec &Lam) { Lam = LamF_val(sig); };
        mechs[1].Lambda = [&](const vec &sig, vec &Lam) { Lam = LamR_val(sig); };
        mechs[2].Lambda = [&](const vec &sig, vec &Lam) { Lam = LamReo_val(sig); };
        mechs[0].dPhi_dsigma = [&](const vec &sig, vec &n) { linearise(sig); n = dPhiF_ds; };
        mechs[1].dPhi_dsigma = [&](const vec &sig, vec &n) { linearise(sig); n = dPhiR_ds; };
        mechs[2].dPhi_dsigma = [&](const vec &sig, vec &n) { linearise(sig); n = dPhiReo_ds; };
        mechs[0].dLambda_dsigma = [&](const vec &sig, mat &D) { linearise(sig); D = dLamF_ds; };
        mechs[1].dLambda_dsigma = [&](const vec &sig, mat &D) { linearise(sig); D = dLamR_ds; };
        mechs[2].dLambda_dsigma = [&](const vec &sig, mat &D) { linearise(sig); D = dLamReo_ds; };
        ReturnStateHooks hooks;
        hooks.update_state = refresh_state;
        hooks.K = [&](const vec &sig, const vec &, mat &K) { linearise(sig); K = K_cpp; };
        hooks.flow_state_coupling = [&](const vec &sig, const vec &, mat &dLdDs) { linearise(sig); dLdDs.col(0) = c_F; dLdDs.col(1).zeros(); dLdDs.col(2) = c_Reo; };
        // FB normalisation scales, floored at an absolute material scale (see unified_T).
        vec Ycrit_cpp(3);
        const double Yfloor = 0.1*std::max(Y0t + D*Hmax*sigmacaliber, 1.e-2);
        Ycrit_cpp(0) = std::max(fabs(Y0t) + fabs(D*Hcur*Mises_stress(sigma)), Yfloor);
        Ycrit_cpp(1) = std::max(fabs(Y0t) + fabs(D*sum(sigma%ETMean)), Yfloor);
        Ycrit_cpp(2) = std::max(YReo, Yfloor);

        ReturnMappingControl control;
        control.ncp = ncp_m;
        rm = closest_point_return_mapping(sigma, L_cpp, mechs, hooks, Ycrit_cpp, control);
        if (rm.converged) {
            error = rm.error;   // below precision_umat: the post-loop step-cut guard stays quiet
            sigma = rm.sigma;
            Ds_j = rm.Dlambda;
            s_j += Ds_j;
            xiF = s_j(0); xiR = s_j(1); pTR = s_j(2);
            DETF  = Ds_j(0)*lambdaTF;
            DETR  = Ds_j(1)*ETMean;
            DEReo = Ds_j(2)*lambdaReo;
            M_eff = xi*M_M + (1. - xi)*M_A;   // output L at the converged xi (the tangent keeps L_cpp)
            L = inv_cubic(M_eff);
            DM_sig = DM*sigma;
            A_xiF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sigma%DM_sig) + sum(sigma%Dalpha)*thp - HfF;
            A_xiR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sigma%DM_sig) - sum(sigma%Dalpha)*thp + HfR;
            dA_xiFdsigma = DM_sig + Dalpha*thp;
            dA_xiRdsigma = -1.*dA_xiFdsigma;
            lambdaTR = -1.*ETMean;
            // flows of the transformation strain: the mixture terms of the effective flows taken out
            dflow_dsigma = join_rows(dLamF_ds - DM, dLamR_ds + DM, dLamReo_ds);
            flow_chain = join_rows(c_F, zeros(6), c_Reo);
            // S = (1 + lambda1Reo) X and DEReo = Ds_2 lambdaReo as functions of (sigma, Ds) through areo
            const vec S_end = (1. + lambda1Reo)*X;
            const vec dS_dre = -1.*dseff_dre_l;
            const vec dDEReo_dre = lambdaReo + (Ds_j(2)*xi_start)*dr_dre_l;
            dWReo_dsigma = -0.5*((eye(6, 6) - dseff_ds_l).t()*DEReo + (Ds_j(2)*dLamReo_ds).t()*(S_start + S_end));
            dWReo_ds = {0., 0., -0.5*(sum(dS_dre%DEReo) + sum((S_start + S_end)%dDEReo_dre))};
        }
        else {
            // step cut (error stays 1, the guard below requests it), never a cutting-plane fallback
            sigma = sigma_start;
            refresh_state(sigma, zeros(3));
        }
    }
    else {
    const vec Ir05_v = Ir05();   // engineering -> tensorial shear
    // backtracking of the cutting plane: the state before the last step
    double bt_error = -1., bt_xi = xi;
    int n_bt = 0;
    vec bt_s = s_j, bt_ET = ET, bt_DETF = DETF, bt_DETR = DETR, bt_DEReo = DEReo, bt_EReo = EReo, bt_areo = areo, bt_Ds = Ds_j;
    vec bt_lF = zeros(6), bt_lR = zeros(6), bt_lReo = zeros(6), bt_la = zeros(6), bt_ds = zeros(3);
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        // Reuss mixing
        M_eff = xi*M_M + (1. - xi)*M_A;
        L = inv_cubic(M_eff);

        DM_sig   = DM*sigma;

        // equivalent stresses and gradients at the current sigma, shared below (Dalpha_T is loop-invariant)
        const double sig_eq = Mises_stress(sigma);
        const vec eta_sig = eta_stress(sigma);
        double Dr_sig;
        vec dDr_sig;
        if (aniso_criteria) {
            Dr_sig = Drucker_ani_stress(sigma, DFA_params, prager_b, prager_n);
            dDr_sig = dDrucker_ani_stress(sigma, DFA_params, prager_b, prager_n);
        } else {
            Dr_sig = Drucker_stress(sigma, prager_b, prager_n);
            dDr_sig = dDrucker_stress(sigma, prager_b, prager_n);
        }

        // Hcur explicit (depends on \sigma )
        if (sig_eq > sigmacrit) sigmastar = sig_eq - sigmacrit;
        else                    sigmastar = 0.;
        Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));

        // Re-evaluate flows at current \sigma , areo, \lambda_1Reo (build_lambdaTF at sigma)
        lambdaTF = Hcur * dDr_sig;
        lambdaTR = -1. * ETMean;

        X = HReo * (areo % Ir05_v);
        a_eq = Mises_strain(areo);
        // intrinsic saturation argument ||v^re||/ETRmax, no xi: the penalty would blow up as xi -> 0
        lambda1Reo = lagrange_pow_1(a_eq/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
        sigma_eff = sigma - (1. + lambda1Reo) * X;

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
        // only the macroscopic flow is xi-scaled (thesis §4 eq. 503, see header)
        lambdaReo   = xi_start * lambdaReo_raw;   // macroscopic flow into \varepsilon^T
        etaReo      = lambdaReo_raw;              // \partial Phi^re/\partial \sigma — un-scaled
        lambda_areo = lambdaReo_raw;              // back-strain v^re evolution — un-scaled

        kappa_j[0] = L*(lambdaTF  + DM_sig + Dalpha_T);
        kappa_j[1] = L*(lambdaTR  - DM_sig - Dalpha_T);
        kappa_j[2] = L*(lambdaReo);   // \xi -scaled — vanishes at \xi =0 as it should

        // Hardening
        compute_Hf(xi, HfF, HfR);

        // Forward thermodynamic force (build_PhihatF at sigma)
        PhihatF = Hcur * Dr_sig;
        A_xiF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sigma%DM_sig) + sum(sigma%Dalpha)*(T + DT - T_init) - HfF;
        lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        const double DHM = D*Hcur*sig_eq;
        double YtF = Y0t + DHM;
        Phi(0) = PhihatF + A_xiF - lambda1 - YtF;

        // Reverse thermodynamic force
        PhihatR = sum(sigma%ETMean);
        A_xiR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sigma%DM_sig) - sum(sigma%Dalpha)*(T + DT - T_init) + HfR;
        lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        double YtR = Y0t + D*PhihatR;
        Phi(1) = -1.*PhihatR + A_xiR + lambda0 - YtR;

        // Reorientation criterion, per unit martensite (thesis eq. 664, Form B)
        double Prager_sigma_eff = aniso_criteria
            ? Drucker_ani_stress(sigma_eff, DFA_params, prager_b, prager_n)
            : Drucker_stress(sigma_eff, prager_b, prager_n);
        Phi(2) = Prager_sigma_eff - YReo;

        // Hardening derivatives w.r.t. xi
        {
            const SMAHardening hs = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
            dHfF = hs.dF;
            dHfR = hs.dR;
        }

        // d(Hcur)/d\sigma 
        dHcurdsigma = k1*(Hmax - Hmin)*exp(-1.*k1*sigmastar)*eta_sig;

        // ---- Forward derivatives ----
        dPhihatFdsigma = dHcurdsigma * Dr_sig + Hcur * dDr_sig;
        dA_xiFdsigma = DM_sig + Dalpha_T;
        dlambda1dsigma = zeros(6);
        dYtFdsigma = D*(dHcurdsigma * sig_eq + Hcur * eta_sig);

        dPhiFdsigma = dPhihatFdsigma + dA_xiFdsigma - dlambda1dsigma - dYtFdsigma;
        double dlambda1dxi = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dPhiFdxiF = -dHfF - dlambda1dxi;
        dPhiFdxiR =  dHfF + dlambda1dxi;
        dPhiFdpTR = 0.;

        // ---- Reverse derivatives ----
        const double xi_div = xi_safe(xi);
        dPhihatRdsigma = ETMean;
        dA_xiRdsigma   = -1.*DM_sig - Dalpha_T;
        double dlambda0dxi = dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda0dsigma = zeros(6);
        dlambda0dxiF = -1.*dlambda0dxi;
        dlambda0dxiR =     dlambda0dxi;
        dYtRdsigma = D*ETMean;
        dPhiRdsigma = -1.*dPhihatRdsigma + dA_xiRdsigma + dlambda0dsigma - dYtRdsigma;
        dA_xiRdxiF =  dHfR;
        dA_xiRdxiR = -dHfR;
        dPhiRdpTR = 0.;

        // K(1,.) through ETMean: avoids a ~12-digit cancellation at small xi (see unified_T.hpp)
        const vec dPhiRdETMean = -(1. + D) * sigma;
        const double dPhiRdxi_finite_F =  dA_xiRdxiF + dlambda0dxiF;   //  +dHfR  - d\lambda_0/d\xi
        const double dPhiRdxi_finite_R =  dA_xiRdxiR + dlambda0dxiR;   //  -dHfR  + d\lambda_0/d\xi

        // ---- Reorientation derivatives (Form B) ----
        dlambda1Reo_arg = dlagrange_pow_1(a_eq/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo);
        const double eta_dot_X = sum(etaReo % X);                  // \eta_Reo : X (tensor double contraction)
        // eta_strain has engineering shears; the gradient of a_eq carries half shears, hence % Ir05
        const vec   eta_strain_areo = (a_eq > simcoon::iota) ? vec(eta_strain(areo) % Ir05_v) : vec(zeros(6));

        const vec dlambda1Reo_dareo = dlambda1Reo_arg * (1./ETRmax) * eta_strain_areo;

        dPhiReodsigma = etaReo;
        dPhiReodareo = -(1. + lambda1Reo) * HReo * (etaReo % Ir05_v)
                       - eta_dot_X * dlambda1Reo_dareo;
        dPhiReodpTR = 0.;

        // ---- K: state-coupling block, without the sigma-flow term ----
        K(0,0) = dPhiFdxiF;
        K(0,1) = dPhiFdxiR;
        K(0,2) = dPhiFdpTR;

        // K(1,.) via Lambda_ETMean^j (header); no division when the fallback ETMean = lambdaTF is active
        vec Lambda_ETMean_F = zeros(6);
        if (Mises_strain(lambdaTF - ETMean) > simcoon::iota) {
            Lambda_ETMean_F = (lambdaTF - ETMean) / xi_div;
        }
        const vec Lambda_ETMean_R = zeros(6);
        const vec Lambda_ETMean_re = lambdaReo_raw;
        K(1,0) = dPhiRdxi_finite_F + sum(dPhiRdETMean % Lambda_ETMean_F);
        K(1,1) = dPhiRdxi_finite_R + sum(dPhiRdETMean % Lambda_ETMean_R);
        K(1,2) = dPhiRdpTR        + sum(dPhiRdETMean % Lambda_ETMean_re);

        // Phi^re depends on the state only through areo: K(2,F) = K(2,R) = 0
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

        // Magnitude scales for the convergence measure |FB|/Y_crit, not the signed YtF/YtR
        // (they cross zero when D < 0): see the convergence-measure note in unified_T.hpp.
        Y_crit(0) = std::max(fabs(Y0t) + fabs(DHM), simcoon::iota);
        Y_crit(1) = std::max(fabs(Y0t) + fabs(D*PhihatR), simcoon::iota);
        Y_crit(2) = std::max(YReo, simcoon::iota);   // guard the FB /Y_crit normalization against a (degenerate) YReo = 0

        // active set: drop a mechanism with Phi < 0 and B(j,j) = 0 (singular FB row, see header)
        const arma::uvec active = arma::find((Phi >= 0.) + (arma::abs(B.diag()) >= simcoon::iota));
        vec Ds_active = Ds_j(active);   // the helper accumulates Dp += dp in place
        vec ds_active(active.n_elem);
        Fischer_Burmeister_m(Phi(active), Y_crit(active), B(active, active), Ds_active, ds_active, error);
        Ds_j(active) = Ds_active;
        ds_j.zeros();
        ds_j(active) = ds_active;
        if ((bt_error >= 0.) && (error > 2.*bt_error) && (n_bt < 5)) {
            // the merit grew: back to the state before the last step, half of it along its flows
            s_j = bt_s; ET = bt_ET; xi = bt_xi; DETF = bt_DETF; DETR = bt_DETR; DEReo = bt_DEReo;
            EReo = bt_EReo; areo = bt_areo; Ds_j = bt_Ds;
            lambdaTF = bt_lF; lambdaTR = bt_lR; lambdaReo = bt_lReo; lambda_areo = bt_la;
            ds_j = 0.5*bt_ds;
            bt_ds = ds_j;
            Ds_j += ds_j;
            L = inv_cubic(xi*M_M + (1. - xi)*M_A);
            n_bt++;
        }
        else {
            bt_error = error; n_bt = 0;
            bt_s = s_j; bt_ET = ET; bt_xi = xi; bt_DETF = DETF; bt_DETR = DETR; bt_DEReo = DEReo;
            bt_EReo = EReo; bt_areo = areo; bt_Ds = Ds_j - ds_j;
            bt_lF = lambdaTF; bt_lR = lambdaTR; bt_lReo = lambdaReo; bt_la = lambda_areo; bt_ds = ds_j;
        }

        s_j(0) += ds_j(0);
        s_j(1) += ds_j(1);
        s_j(2) += ds_j(2);

        // State updates
        ET   = ET   + ds_j(0)*lambdaTF + ds_j(1)*lambdaTR + ds_j(2)*lambdaReo;
        EReo = EReo + ds_j(2)*lambdaReo;
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
        sigma = el_pred(L, Eel);
    }
    }

    // closest point not solved: request a sub-step (the cutting plane commits its last iterate, as unified_T)
    if (use_cpp && error > simcoon::precision_umat) {
        tnew_dt = 0.5;
    }

    // forces of the mechanisms, start and end of the increment, with their slopes: the transformation
    // rows carry the bound penalties of xi (stored like the hardening), the reorientation row has none;
    // the back-strain channel has its own force A_Reo = -(1 + lambda1Reo) X, its work recovered (Wm_r)
    const uword N = 3;
    vec A_start = {A_xiF_start, A_xiR_start, 0.};
    vec A_end = {A_xiF, A_xiR, 0.};
    const mat dA_dsigma = join_rows(dA_xiFdsigma, dA_xiRdsigma, zeros(6));
    mat dA_ds;
    {
        const SMAHardening hs_end = sma_hardening(xi, a1, a2, a3, n1, n2, n3, n4);
        const double sF = hs_end.dF, sR = hs_end.dR;
        dA_ds = {{-sF, sF, 0.}, {sR, -sR, 0.}, {0., 0., 0.}};
        A_start(0) -= lagrange_pow_1(xi_start, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        A_start(1) -= lagrange_pow_0(xi_start, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        A_end(0) -= lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        A_end(1) -= lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        const double dl1 = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        const double dl0 = -1.*dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dA_ds(0, 0) -= dl1;
        dA_ds(0, 1) += dl1;
        dA_ds(1, 0) += dl0;
        dA_ds(1, 1) -= dl0;
    }
    const double dA_xiFdtheta = rhoDs0 + sum(sigma%Dalpha);
    const vec dA_dtheta = {dA_xiFdtheta, -1.*dA_xiFdtheta, 0.};
    const vec A_Reo_start = -1.*S_start;
    const vec A_Reo = -1.*(1. + lagrange_pow_1(Mises_strain(areo)/ETRmax, c_lambdaReo, p0_lambdaReo, n_lambdaReo, alpha_lambdaReo))
                      * (HReo * (areo % Ir05()));
    if (!use_cpp && Ds_j(2) > 0.) {
        dWReo_dsigma = zeros(6);
        dWReo_ds = {0., 0., -0.5*(sum((backstress_slope()*lambda_areo)%DEReo) - sum((A_Reo_start + A_Reo)%lambdaReo))};
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
        kappa = join_rows(kappa_j[0], kappa_j[1], kappa_j[2]);
        dPhi = join_rows(dPhiFdsigma, dPhiRdsigma, dPhiReodsigma);
        mat Bhat = zeros(3, 3);
        for (int l = 0; l < 3; ++l) {
            for (int j = 0; j < 3; ++j) {
                Bhat(l, j) = sum(dPhi.col(l) % kappa_j[j]) - K(l, j);
            }
        }
        ct = compute_tangent_operator(tangent_continuum, Bhat, kappa, dPhi, Ds_j, L);
    }
    dSdE = ct.Lt;
    // closest point without a solve (elastic trial, step cut): the elastic operator, zero sensitivities
    const bool solved = !use_cpp || (rm.converged && !rm.kappa_j.is_empty());
    const mat P_epsilon = solved ? ct.P_epsilon : mat(6, N, fill::zeros);
    const mat invBhat = solved ? ct.invBhat : mat(N, N, fill::zeros);

    // trial stress rate of a temperature move: the closest point freezes alpha at xi_n and folds
    // Dalpha into the flows, so its thermal strain is alpha(xi_n+1); the algorithmic operator acts through M^-1
    const bool algo = solved && !ct.kappa_tilde.is_empty();
    const vec alpha_q = use_cpp ? vec(alpha + (xi - xi_start)*Dalpha) : alpha;
    const vec L_alpha = algo ? vec(ct.L_tilde*alpha_q) : vec(L_tan*alpha_q);
    const mat kappa_q = !solved ? mat(6, N, fill::zeros) : (algo ? ct.kappa_tilde : kappa);
    const mat dPhi_q = solved ? dPhi : mat(6, N, fill::zeros);

    // P_theta = invBhat (dPhi/dtheta - dPhi/dsigma : L alpha), Bhat(l, j) = criterion l, mechanism j
    vec rhs(N), P_theta(N);
    for (uword l = 0; l < N; l++) rhs(l) = dA_dtheta(l) - sum(dPhi_q.col(l)%L_alpha);
    for (uword l = 0; l < N; l++) {
        P_theta(l) = invBhat(l, 0)*rhs(0) + invBhat(l, 1)*rhs(1);
        P_theta(l) += invBhat(l, 2)*rhs(2);
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
    dA_theta_ds(0) = Dc_0*log(T1/T_init) + sum(Dalpha%sigma) + rhoDs0;
    dA_theta_ds(1) = -1.*Dc_0*log(T1/T_init) - sum(Dalpha%sigma) - rhoDs0;

    const vec DET = ET - ET_start;
    const mat flows = join_rows(lambdaTF, lambdaTR, lambdaReo);
    const double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DET) + 0.5*(A_start(0) + A_end(0))*Ds_j(0) + 0.5*(A_start(1) + A_end(1))*Ds_j(1)
                            + 0.5*sum((A_Reo_start + A_Reo)%DEReo);

    // dissipation derivatives, flow directions held fixed over the increment
    vec dD = dSdE.t()*DET;
    vec g_T(N);
    for (uword j = 0; j < N; j++) {
        vec g = dSdE.t()*dA_dsigma.col(j);
        for (uword k = 0; k < N; k++) g += dA_ds(j, k)*P_epsilon.col(k);
        dD += g*Ds_j(j);
        g_T(j) = dA_dtheta(j) + sum(dA_dsigma.col(j)%dSdT) + dA_ds(j, 0)*P_theta(0) + dA_ds(j, 1)*P_theta(1);
        g_T(j) += dA_ds(j, 2)*P_theta(2);
    }
    double dD_T_turn = 0.;
    vec f(N);
    for (uword j = 0; j < N; j++) f(j) = 0.5*(sum((sigma_start+sigma)%flows.col(j)) + A_start(j) + A_end(j));
    // closest point: the turning of the flows and their multiplier chain, d DET = sum_k (Lambda_k dDs_k + Ds_k dLambda_k)
    if (!dflow_dsigma.is_empty()) {
        vec turn = zeros(6);
        for (uword k = 0; k < N; k++) turn += Ds_j(k)*(dflow_dsigma.cols(6*k, 6*k + 5).t()*(sigma_start+sigma));
        dD += dSdE.t()*turn;
        for (uword j = 0; j < N; j++) f(j) += 0.5*sum((sigma_start+sigma)%flow_chain.col(j));
        dD_T_turn = sum(dSdT%turn);
    }
    vec dDgamma_dE = 0.5*dD;
    for (uword j = 0; j < N; j++) dDgamma_dE += f(j)*P_epsilon.col(j);
    double dD_T = sum(dSdT%DET) + g_T(0)*Ds_j(0) + g_T(1)*Ds_j(1);
    dD_T += g_T(2)*Ds_j(2);
    dD_T += dD_T_turn;
    double dDgamma_dT = 0.5*dD_T + f(0)*P_theta(0) + f(1)*P_theta(1);
    dDgamma_dT += f(2)*P_theta(2);
    // the work of the back-strain force, d[1/2 (A_Reo_n + A_Reo) : DEReo]
    if (!dWReo_ds.is_empty()) {
        dDgamma_dE += dSdE.t()*dWReo_dsigma + P_epsilon*dWReo_ds;
        dDgamma_dT += sum(dWReo_dsigma%dSdT) + sum(dWReo_ds%P_theta);
    }

    // Entropy derivatives at the end of the increment
    const vec deta_dE = dSdE.t()*alpha + dA_theta_ds(0)*P_epsilon.col(0) + dA_theta_ds(1)*P_epsilon.col(1);
    const double deta_dT = c_0/T1 + sum(dSdT%alpha) + dA_theta_ds(0)*P_theta(0) + dA_theta_ds(1)*P_theta(1);

    // heat of the increment and its derivatives (heat_source_increment)
    heat_source_increment(Dgamma_loc, dDgamma_dE, dDgamma_dT, Deta, deta_dE, deta_dT, T1, r, drdE, drdT);

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DET)) - 0.5*(A_start(0) + A_end(0))*Ds_j(0) - 0.5*(A_start(1) + A_end(1))*Ds_j(1)
          - 0.5*sum((A_Reo_start + A_Reo)%DEReo);
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;

    Wt += T1*Deta;
    Wt_r += T1*Deta;
    Wt_ir += 0.;

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
