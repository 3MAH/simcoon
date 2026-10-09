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

///@file plastic_johnson_cook_ccp.cpp
///@brief Thermomechanical Johnson-Cook elastic-viscoplastic UMAT (EPJCK), CCP integration with
///       the exact heat source of the Gibbs framework
///@version 1.0

#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void umat_plasticity_johnson_cook_CCP_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);

    //From the props to the material properties
    double rho = props(0);
    double c_p = props(1);
    double E = props(2);
    double nu = props(3);
    double alpha_iso = props(4);
    double A_jc = props(5);
    double B_jc = props(6);
    double n_jc = props(7);
    double C_jc = props(8);
    double edot0 = props(9);
    double m_jc = props(10);
    double T_ref = props(11);
    double T_melt = props(12);

    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();

    //Elastic stiffness tensor
    mat L = L_iso(E, nu, "Enu");

    //Temperature initialization
    double T_init = statev(0);
    //From the statev to the internal variables
    double p = statev(1);
    vec EP(6);
    EP(0) = statev(2);
    EP(1) = statev(3);
    EP(2) = statev(4);
    EP(3) = statev(5);
    EP(4) = statev(6);
    EP(5) = statev(7);

    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);

    //Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        sigma = vide;
        EP = vide;
        p = 0.;

        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;

        Wt = 0.;
        Wt_r = 0.;
        Wt_ir = 0.;
    }

    //Additional parameters
    double c_0 = rho*c_p;

    // Thermal softening factor f_T = 1 - T*^m and its temperature derivative; T* clamped to
    // [0, 1) (defined below T_ref, never at the singular melting point), the derivative is zero
    // where the clamped factor is flat
    auto homologous = [&](const double &theta) {
        double Tstar = (theta - T_ref) / (T_melt - T_ref);
        if (Tstar < 0.) Tstar = 0.;
        if (Tstar >= 1.) Tstar = 1. - simcoon::iota;
        return Tstar;
    };
    auto dthermal_dT = [&](const double &Tstar) {
        return (Tstar > simcoon::iota && Tstar < 1. - simcoon::iota)
            ? -m_jc * pow(Tstar, m_jc - 1.) / (T_melt - T_ref) : 0.;
    };
    const double Tstar_start = homologous(T);
    const double Tstar = homologous(T + DT);
    const double thermal_factor_start = 1. - pow(Tstar_start, m_jc);
    const double thermal_factor = 1. - pow(Tstar, m_jc);
    const double dthermal_dT_start = dthermal_dT(Tstar_start);
    const double dthermal_dT_end = dthermal_dT(Tstar);

    // Strain hardening Hp = B p^n, C1-regularized at the onset for n < 1 (the modular
    // power-law block carries the blend; same law, no second implementation)
    PowerLawHardening hardening;
    int hardening_offset = 0;
    hardening.configure(vec{B_jc, n_jc}, hardening_offset);
    double Hp = hardening.R(p);
    double dHpdp = hardening.dR_dp(p);
    // Stored energy of the hardening, G^ir = f_T B p^(n+1)/(n+1) (the blend below p_reg
    // contributes O(p_reg^(n+1)), neglected in the integral)
    auto stored_hardening = [&](const double &pp) {
        return (pp > 0.) ? B_jc * pow(pp, n_jc + 1.) / (n_jc + 1.) : 0.;
    };

    // Rate factor and yield stress: evaluated in the CCP loop from the implicit Dp/DTime
    double rate_factor = 1.;
    double sigmaY_jc = 0.;

    //Variables values at the start of the increment; hardening force A_p = -f_T B p^n
    vec sigma_start = sigma;
    vec EP_start = EP;
    const double p_start = p;
    double A_p_start = -thermal_factor_start*Hp;

    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);

    //Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    sigma = el_pred(L, Eel, ndi);

    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B_mat = zeros(1,1);
    vec Y_crit = zeros(1);

    double dPhidp = 0.;
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;

    // Flow direction Lambda = dPhi/dsigma (associated J2) and kappa = L:Lambda, set in the loop
    vec Lambdap = zeros(6);
    std::vector<vec> kappa_j(1);
    mat K = zeros(1,1);

    // Branch selection at the rate kink x0 = edot0*DTime, where Phi(Dp) changes from the
    // quasi-static branch (rate factor 1) to the convex logarithmic one. A Newton started at
    // Dp = 0 with the quasi-static slope overshoots the kink and the log branch throws it back
    // below zero: a 2-cycle to maxiter_umat at small DTime. So the branch is decided first from
    // Phi at the kink (radial estimate), the loop starts AT the kink on the log branch with the
    // right derivative of the rate factor, and Newton is then monotone (convex branch).
    bool log_branch = false;
    if (DTime > simcoon::iota) {
        const double x0 = edot0*DTime;
        const vec Lambda_trial = eta_stress(sigma);
        const double Phi_x0 = Mises_stress(sigma) - x0*sum(Lambda_trial%(L*Lambda_trial))
                              - (A_jc + hardening.R(p + x0))*thermal_factor;
        if (Phi_x0 > 0.) {
            log_branch = true;
            Ds_j(0) = x0;
            s_j(0) += x0;
            EP = EP + x0*Lambda_trial;
            Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
            sigma = el_pred(L, Eel, ndi);
        }
    }

    //Loop parameters
    int compteur = 0;
    double error = 1.;

    //CCP Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        p = s_j(0);

        Hp = hardening.R(p);
        dHpdp = hardening.dR_dp(p);

        // Strain rate, fully implicit: pdot = Dp/DTime, clamped at edot0 (rate_factor >= 1)
        const double Dp_j = Ds_j(0);
        double drate_dDp = 0.;
        if (DTime > simcoon::iota) {
            const double edot_eff = std::max(Dp_j / DTime, edot0);
            rate_factor = 1. + C_jc * log(edot_eff / edot0);
            if (edot_eff > edot0 || log_branch) {
                drate_dDp = C_jc / (edot_eff * DTime);
            }
        }
        else {
            rate_factor = 1.;
        }

        sigmaY_jc = (A_jc + Hp) * rate_factor * thermal_factor;

        dPhidsigma = eta_stress(sigma);
        // K = dPhi/dp + dPhi/dDp: hardening and rate sensitivity
        dPhidp = -dHpdp * rate_factor * thermal_factor
                 - (A_jc + Hp) * drate_dDp * thermal_factor;

        //compute Phi and the derivatives
        Phi(0) = Mises_stress(sigma) - sigmaY_jc;

        Lambdap = dPhidsigma;
        kappa_j[0] = L*Lambdap;

        K(0,0) = dPhidp;
        B_mat(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = std::max(sigmaY_jc, simcoon::precision_umat);

        Fischer_Burmeister_m(Phi, Y_crit, B_mat, Ds_j, ds_j, error);

        s_j(0) += ds_j(0);
        EP = EP + ds_j(0)*Lambdap;

        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        sigma = el_pred(L, Eel, ndi);
    }

    // The loop reads p at the top and corrects s_j at the bottom: bring p and the hardening to
    // the converged iterate (EP and sigma already are)
    p = s_j(0);
    Hp = hardening.R(p);
    dHpdp = hardening.dR_dp(p);

    // A loop that leaves at maxiter_umat has not met the yield surface: the state would be
    // committed off-surface, so ask the solver for a step cut instead
    if (error > simcoon::precision_umat) {
        tnew_dt = 0.5;
    }

    //Computation of the increments of variables
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];

    // Effective plastic strain rate over the increment (diagnostic output)
    double edot_p_out = 0.;
    if (DTime > simcoon::iota && Dp > simcoon::iota) {
        edot_p_out = Dp / DTime;
    }

    // dPhi/dT = -(A + Hp) f_rate df_T/dT > 0 (softening); zero where the clamped factor is flat
    dPhidtheta = -(A_jc + Hp) * rate_factor * dthermal_dT_end;

    //Computation of the tangent modulus via the shared leading-mechanism helper (doc 7.4).
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    const std::vector<vec> dPhidsigma_l = { dPhidsigma };
    // tangent_none must NOT zero P_epsilon/invBhat here: these sensitivities feed the
    // PHYSICAL heat source r and its linearization, not just the Newton operator
    // (same rule as the other thermomechanical kernels).
    const int tangent_mode_eff = (tangent_mode == tangent_none)
        ? tangent_continuum : tangent_mode;
    const ContinuumTangent ct = compute_tangent_operator(
        tangent_mode_eff, Bhat, kappa_j, dPhidsigma_l, Ds_j, L,
        [&]() -> std::vector<mat> {  // lazy: evaluated only in algorithmic mode
            // J2 associated flow: dLambda/dsigma = deta_stress(sigma). Only dSdE is
            // algorithmically corrected; the thermal cross-tangents keep the continuum form.
            const std::vector<mat> dLambda_dsigma_l = { deta_stress(sigma) };
            return dLambda_dsigma_l;
        });
    dSdE = ct.Lt;
    const std::vector<vec>& P_epsilon = ct.P_epsilon;

    std::vector<double> P_theta(1);
    // P_theta = invBhat (dPhi/dtheta - dPhi/dsigma:L:alpha); invBhat is active-set masked,
    // so elastic steps give P_theta = 0
    P_theta[0] = ct.invBhat(0, 0) * (dPhidtheta - sum(dPhidsigma%(L*alpha)));

    dSdT = -1.*L*alpha - (kappa_j[0]*P_theta[0]);

    //computation of the internal energy production
    double eta_r = c_0*log((T+DT)/T_init) + sum(alpha%sigma);
    double eta_r_start = c_0*log(T/T_init) + sum(alpha%sigma_start);

    // Irreversible entropy of the temperature-dependent stored hardening, eta_ir = -dG^ir/dT
    double eta_ir = -dthermal_dT_end*stored_hardening(p);
    double eta_ir_start = -dthermal_dT_start*stored_hardening(p_start);
    // d eta_ir / dp at the end of the increment: the reversible coupling -theta deta_ir/dp pdot
    // enters the heat source like the thermoelastic -theta alpha:dsigma/dt
    const double deta_ir_dp = -dthermal_dT_end*Hp;

    double eta = eta_r + eta_ir;
    double eta_start = eta_r_start + eta_ir_start;

    double Deta = eta - eta_start;
    double Deta_r = eta_r - eta_r_start;
    double Deta_ir = eta_ir - eta_ir_start;

    vec Gamma_epsilon = zeros(6);
    double Gamma_theta = 0.;

    vec N_epsilon = zeros(6);
    double N_theta = 0.;

    // Hardening force A_p = -f_T B p^n (the dissipation f_T [A f_rate + B p^n (f_rate - 1)] Dp
    // is then non-negative at every temperature and rate) and its derivatives
    double A_p = -thermal_factor*Hp;
    double dA_pdp = -thermal_factor*dHpdp;
    double dA_pdtheta = -dthermal_dT_end*Hp;

    // Dissipation increment over the step (the quantity Wm_d accumulates)
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP) + 0.5*(A_p_start + A_p)*Dp;

    if(DTime < 1.E-12) {
        r = 0.;
        drdE = zeros(6);
        drdT = 0.;
    }
    else {
        Gamma_epsilon = dA_pdp*P_epsilon[0]*(Dp/DTime) + A_p/DTime*P_epsilon[0] + (dSdE*DEP)*(1./DTime) + sum(sigma%Lambdap)*P_epsilon[0]/DTime;
        Gamma_theta = (dA_pdp*P_theta[0] + dA_pdtheta)*(Dp/DTime) + A_p/DTime*P_theta[0] + sum(dSdT%DEP)*(1./DTime) + sum(sigma%Lambdap)*P_theta[0]/DTime;

        // Thermoelastic coupling and the coupling of the stored hardening, -theta deta_ir/dp pdot
        N_epsilon = -1./DTime*(T + DT)*(dSdE*alpha) - (T + DT)/DTime*deta_ir_dp*P_epsilon[0];
        N_theta = -1./DTime*(T + DT)*sum(dSdT%alpha) -1.*Deta/DTime - rho*c_p*(1./DTime)
                  - (T + DT)/DTime*deta_ir_dp*P_theta[0] - deta_ir_dp*Dp/DTime;

        drdE = N_epsilon + Gamma_epsilon;
        drdT = N_theta + Gamma_theta;

        // The dissipation enters r through its converged increment (the form of the viscous
        // kernels, viscous_heat_source), not its linearization Gamma_epsilon:DEtot +
        // Gamma_theta*DT (EPICP_T): at fixed DTime the rate-dependent Dp is a logarithmic
        // function of DEtot, so the first-order expansion built on P_epsilon = invBhat
        // L:dPhi/dsigma -- with the rate term C (A + Hp)/Dp in Bhat -- underestimates it and
        // vanishes under step refinement. Gamma_epsilon / Gamma_theta are its first-order
        // derivatives, so drdE / drdT stay the Newton linearizations of this r.
        // Thermoelastic and hardening-entropy couplings on the converged increment
        const double N_coupling = -(T + DT)*(sum(alpha%(sigma - sigma_start)) + deta_ir_dp*Dp)/DTime
                                  - Deta*DT/DTime - rho*c_p*DT/DTime;
        r = N_coupling + Dgamma_loc/DTime;
    }

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;

    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir += (T+0.5*DT)*Deta_ir;

    //statev evolving variables
    statev(0) = T_init;
    statev(1) = p;

    statev(2) = EP(0);
    statev(3) = EP(1);
    statev(4) = EP(2);
    statev(5) = EP(3);
    statev(6) = EP(4);
    statev(7) = EP(5);

    statev(8) = edot_p_out;
}

} //namespace simcoon
