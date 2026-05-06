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
///@brief User subroutine for elastic-plastic materials with Johnson-Cook hardening in 1D-2D-3D case
///@brief This subroutine uses a convex cutting plane algorithm
///@brief Johnson-Cook rate-dependent and temperature-dependent hardening is considered
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_johnson_cook_ccp.hpp>


using namespace std;
using namespace arma;

namespace simcoon {

///@brief The Johnson-Cook UMAT requires 11 constants for a mechanical coupling:

///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE
///@brief props[3] : JC initial yield stress A
///@brief props[4] : JC hardening coefficient B
///@brief props[5] : JC hardening exponent n
///@brief props[6] : JC strain rate sensitivity C
///@brief props[7] : JC reference strain rate edot0
///@brief props[8] : JC thermal softening exponent m_JC
///@brief props[9] : JC reference temperature T_ref
///@brief props[10] : JC melting temperature T_melt

///@brief The Johnson-Cook UMAT requires 9 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : Accumulative plastic parameter: p
///@brief statev[2] : Plastic strain 11: EP(0,0)
///@brief statev[3] : Plastic strain 22: EP(1,1)
///@brief statev[4] : Plastic strain 33: EP(2,2)
///@brief statev[5] : Plastic strain 12: EP(0,1) (*2)
///@brief statev[6] : Plastic strain 13: EP(0,2) (*2)
///@brief statev[7] : Plastic strain 23: EP(1,2) (*2)
///@brief statev[8] : Plastic strain rate edot_p (diagnostic output)

void umat_plasticity_johnson_cook_CCP(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);

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

    ///@brief Temperature initialization
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

    ///@brief Initialization
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

    // Current temperature for thermal softening
    double T_cur = T + DT;

    // Compute the homologous temperature T*
    // Clamp to [0, 1) to avoid negative or melting
    double Tstar = (T_cur - T_ref) / (T_melt - T_ref);
    if (Tstar < 0.) Tstar = 0.;
    if (Tstar >= 1.) Tstar = 1. - simcoon::iota;

    // Thermal softening factor: (1 - T*^m)
    double thermal_factor = 1. - pow(Tstar, m_jc);

    // Strain hardening: Hp = B * p^n
    double Hp = 0.;
    double dHpdp = 0.;

    if (p > simcoon::iota) {
        dHpdp = n_jc * B_jc * pow(p, n_jc - 1.);
        Hp = B_jc * pow(p, n_jc);
    }
    else {
        dHpdp = 0.;
        Hp = 0.;
    }

    // Rate factor: will be computed in the CCP loop
    // For the initial evaluation, use reference rate (rate_factor = 1)
    double rate_factor = 1.;

    // Total JC yield stress
    double sigmaY_jc = (A_jc + Hp) * rate_factor * thermal_factor;

    //Variables values at the start of the increment
    vec sigma_start = sigma;
    vec EP_start = EP;
    double A_p_start = -Hp;

    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);

    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    sigma = el_pred(L, Eel, ndi);

    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B_mat = zeros(1,1);
    vec Y_crit = zeros(1);

    double dPhidp = 0.;
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;

    //Compute the explicit flow direction
    vec Lambdap = eta_stress(sigma);
    std::vector<vec> kappa_j(1);
    kappa_j[0] = L*Lambdap;
    mat K = zeros(1,1);

    //Loop parameters
    int compteur = 0;
    double error = 1.;

    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        p = s_j(0);

        // Strain hardening
        if (p > simcoon::iota) {
            dHpdp = n_jc * B_jc * pow(p, n_jc - 1.);
            Hp = B_jc * pow(p, n_jc);
        }
        else {
            dHpdp = 0.;
            Hp = 0.;
        }

        // Strain rate: fully implicit, ε̇ = Δp / ΔTime
        // Clamp to edot0 to avoid ln(0) and ensure rate_factor >= 1
        double Dp = Ds_j(0);
        double edot_eff = edot0; // default: reference rate
        double drate_dDp = 0.;   // derivative of rate_factor w.r.t. Dp

        if (DTime > simcoon::iota) {
            edot_eff = std::max(Dp / DTime, edot0);

            // rate_factor = 1 + C * max(ln(edot_eff/edot0), 0)
            double ln_ratio = log(edot_eff / edot0);
            if (ln_ratio < 0.) ln_ratio = 0.;
            rate_factor = 1. + C_jc * ln_ratio;

            // Derivative: d(rate_factor)/d(Dp) = C / (edot_eff * DTime) when Dp/DTime > edot0
            if (Dp / DTime > edot0) {
                drate_dDp = C_jc / (edot_eff * DTime);
            }
            else {
                drate_dDp = 0.;
            }
        }
        else {
            // DTime ~ 0: quasi-static, no rate effect
            rate_factor = 1.;
            drate_dDp = 0.;
        }

        // Total JC yield stress
        sigmaY_jc = (A_jc + Hp) * rate_factor * thermal_factor;

        dPhidsigma = eta_stress(sigma);

        // dΦ/dp = -dσ_Y/dp = -(dHpdp * rate_factor * thermal_factor)
        // dΦ/dDp (rate contribution) = -(A_jc + Hp) * drate_dDp * thermal_factor
        // Total K(0,0) = dΦ/dp + dΦ/dDp (since dp = dDp in the CCP correction)
        dPhidp = -dHpdp * rate_factor * thermal_factor
                 - (A_jc + Hp) * drate_dDp * thermal_factor;

        //compute Phi and the derivatives
        Phi(0) = Mises_stress(sigma) - sigmaY_jc;

        Lambdap = eta_stress(sigma);
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

    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];

    // Store the effective plastic strain rate for output
    double edot_p_out = 0.;
    if (DTime > simcoon::iota && Dp > simcoon::iota) {
        edot_p_out = Dp / DTime;
    }

    //Computation of the tangent modulus
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    vec op = zeros(1);
    mat delta = eye(1,1);

    for (int i=0; i<1; i++) {
        if(Ds_j[i] > simcoon::iota)
            op(i) = 1.;
    }

    mat Bbar = zeros(1,1);
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
        }
    }

    mat invBbar = zeros(1, 1);
    mat invBhat = zeros(1, 1);
    invBbar = inv(Bbar);
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }

    std::vector<vec> P_epsilon(1);
    P_epsilon[0] = invBhat(0, 0)*(L*dPhidsigma);
    std::vector<double> P_theta(1);
    P_theta[0] = dPhidtheta - sum(dPhidsigma%(L*alpha));

    Lt = L - (kappa_j[0]*P_epsilon[0].t());

    double A_p = -Hp;
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP) + 0.5*(A_p_start + A_p)*Dp;

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;

    ///@brief statev evolving variables
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
