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
///@brief Johnson-Cook elastic-viscoplastic UMAT (EPJCK), convex cutting plane integration
///@version 1.0

#include <string>
#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

void umat_plasticity_johnson_cook_CCP(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);

    //From the props to the material properties
    double E = props(0);
    double nu = props(1);
    double alpha_iso = props(2);
    double A_jc = props(3);
    double B_jc = props(4);
    double n_jc = props(5);
    double C_jc = props(6);
    double edot0 = props(7);
    double m_jc = props(8);
    double T_ref = props(9);
    double T_melt = props(10);

    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();

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

    //Elastic stiffness tensor
    L = L_iso(E, nu, "Enu");

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
    }

    // Thermal softening factor f_T = 1 - T*^m, T* clamped to [0, 1) (defined below T_ref, never
    // at the singular melting point), at the start and end-of-increment temperatures
    auto homologous = [&](const double &theta) {
        double Tstar = (theta - T_ref) / (T_melt - T_ref);
        if (Tstar < 0.) Tstar = 0.;
        if (Tstar >= 1.) Tstar = 1. - simcoon::iota;
        return Tstar;
    };
    const double thermal_factor_start = 1. - pow(homologous(T), m_jc);
    const double thermal_factor = 1. - pow(homologous(T + DT), m_jc);

    // Strain hardening Hp = B p^n, C1-regularized at the onset for n < 1 (the modular
    // power-law block carries the blend; same law, no second implementation)
    PowerLawHardening hardening;
    int hardening_offset = 0;
    hardening.configure(vec{B_jc, n_jc}, hardening_offset);
    double Hp = hardening.R(p);
    double dHpdp = hardening.dR_dp(p);

    // Rate factor and yield stress: evaluated in the CCP loop from the implicit Dp/DTime
    double rate_factor = 1.;
    double sigmaY_jc = 0.;

    //Variables values at the start of the increment; hardening force A_p = -f_T B p^n
    vec sigma_start = sigma;
    vec EP_start = EP;
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

    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        p = s_j(0);

        Hp = hardening.R(p);
        dHpdp = hardening.dR_dp(p);

        // Strain rate, fully implicit: pdot = Dp/DTime, clamped at edot0 so that
        // rate_factor >= 1 (no softening below the reference rate, no ln(0) at yield onset)
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
            // DTime = 0: rate-independent form
            rate_factor = 1.;
        }

        sigmaY_jc = (A_jc + Hp) * rate_factor * thermal_factor;

        dPhidsigma = eta_stress(sigma);
        // K = dPhi/dp + dPhi/dDp: hardening and rate sensitivity (dp = dDp in the CCP correction)
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

    //Computation of the tangent modulus via the shared leading-mechanism helper (doc 7.4).
    //The rate term sits in K, hence in Bhat: the algorithmic operator is the exact Jacobian of
    //the discrete map at fixed DTime. tangent_none returns L (explicit integration).
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    const std::vector<vec> dPhidsigma_l = { dPhidsigma };
    const ContinuumTangent ct = compute_tangent_operator(
        tangent_mode, Bhat, kappa_j, dPhidsigma_l, Ds_j, L,
        [&]() -> std::vector<mat> {  // lazy: evaluated only in algorithmic mode
            // J2 associated flow: dLambda/dsigma = deta_stress(sigma), the complete
            // Simo-Hughes correction for a single isotropic mechanism
            const std::vector<mat> dLambda_dsigma_l = { deta_stress(sigma) };
            return dLambda_dsigma_l;
        });
    Lt = ct.Lt;

    // Energy split: stored G^ir = f_T B p^(n+1)/(n+1), so A_p = -f_T B p^n and the dissipation
    // f_T [A f_rate + B p^n (f_rate - 1)] Dp stays non-negative at every temperature and rate
    double A_p = -thermal_factor*Hp;
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP) + 0.5*(A_p_start + A_p)*Dp;

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;

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
