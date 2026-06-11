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

///@file fibre_Fp.cpp
///@brief Single-fibre plasticity with explicit Fp (multiplicative, exponential map)
///@author Chemisky

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/fibre_Fp.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

///@brief props[0]  : Young modulus (Hencky potential)
///@brief props[1]  : Poisson ratio
///@brief props[2]  : CTE
///@brief props[3]  : J2 equivalent yield stress limit : sigmaY
///@brief props[4]  : hardening parameter k
///@brief props[5]  : exponent m
///@brief props[6]  : c1 (fibre-axial deviatoric weight)
///@brief props[7]  : c2 (fibre-transverse-shear weight)
///@brief props[8-10): a0 (reference / relaxed-configuration fibre direction)
///@brief props[11] : anchor (0 = embedded line, rides F; 1 = lattice, rides Fe)

///@brief statev[0]    : T_init
///@brief statev[1]    : p (accumulated plastic strain)
///@brief statev[2-10] : Fp, row-major

void umat_plasticity_fibre_Fp(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const natural_basis &nb, const int &tangent_mode)
{
    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    UNUSED(DR);
    UNUSED(tangent_mode);
    UNUSED(Etot);

    double E = props(0);
    double nu = props(1);
    double alpha_iso = props(2);
    double sigmaY = props(3);
    double k = props(4);
    double m = props(5);
    double c1 = props(6);
    double c2 = props(7);
    vec a0 = {props(8), props(9), props(10)};
    a0 = a0 / norm(a0, 2);
    const bool lattice = (props(11) > 0.5);

    vec alpha = alpha_iso * Ith();
    L = L_iso(E, nu, "Enu");

    double T_init = statev(0);
    double p = statev(1);
    // Fp is anchored in the reference configuration: it is NOT rotated by DR.
    mat Fp = {{statev(2), statev(3), statev(4)},
              {statev(5), statev(6), statev(7)},
              {statev(8), statev(9), statev(10)}};

    if (start) {
        T_init = T;
        sigma = zeros(6);
        p = 0.;
        Fp = eye(3, 3);
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }

    // The solver's state (etot, sigma) is carried in LABORATORY components
    // (verified: ELISO under corate 3/5 returns sigma = L:lnV in lab axes), so
    // the UMAT works directly with the lab deformation gradient carried by the
    // natural basis. nb holds F only under corate_type = 5 (log_F); under any
    // other corate F = I and the multiplicative kinematics below degenerate.
    mat F_tot = nb.F();
    double J = det(F_tot);

    // Trial elastic state: Fe = F Fp^{-1} (lab components), the stress an exact
    // function of (F, Fp) through the Hencky potential on ln(Ve).
    mat Fp_inv = inv(Fp);
    mat Fe_tr = F_tot * Fp_inv;
    mat Be_tr = Fe_tr * Fe_tr.t();
    mat ee_tr = 0.5 * logmat_sympd(Be_tr);
    vec ee_tr_v = t2v_strain(ee_tr);

    vec sigma_start = sigma;

    vec Eth = alpha * (T + DT - T_init);
    sigma = el_pred(L, ee_tr_v - Eth, ndi) / J;

    // Fibre in lab components, explicit over the increment:
    //   embedded: a = F a0 / |F a0|   (material line, rides Fe Fp)
    //   lattice : a = Fe a0 / |Fe a0| (invariant under plastic flow)
    vec av = lattice ? vec(Fe_tr * a0) : vec(F_tot * a0);
    vec a_fibre = av / norm(av, 2);

    double Hp = 0.;
    double dHpdp = 0.;
    if (p > simcoon::iota) {
        dHpdp = m * k * pow(p, m - 1);
        Hp = k * pow(p, m);
    }
    double A_p_start = -Hp;

    // Return mapping in ln(Ve) space: the small-strain algorithm applies to the
    // logarithmic elastic strain unchanged (consistency O(Dp*||ee||) for
    // non-coaxial flow).
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);

    vec Phi = zeros(1);
    mat B = zeros(1, 1);
    vec Y_crit = zeros(1);

    vec dPhidsigma = zeros(6);
    double dPhidp = 0.;

    vec DEP = zeros(6);
    vec Lambdap = dEq_stress_TI(sigma, a_fibre, c1, c2);
    std::vector<vec> kappa_j(1);
    kappa_j[0] = (L * Lambdap) / J;
    mat K = zeros(1, 1);

    int compteur = 0;
    double error = 1.;

    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {

        p = s_j(0);
        if (p > simcoon::iota) {
            dHpdp = m * k * pow(p, m - 1);
            Hp = k * pow(p, m);
        }
        else {
            dHpdp = 0.;
            Hp = 0.;
        }
        dPhidsigma = dEq_stress_TI(sigma, a_fibre, c1, c2);
        dPhidp = -1. * dHpdp;

        Phi(0) = Eq_stress_TI(sigma, a_fibre, c1, c2) - Hp - sigmaY;

        Lambdap = dEq_stress_TI(sigma, a_fibre, c1, c2);
        kappa_j[0] = (L * Lambdap) / J;

        K(0, 0) = dPhidp;
        B(0, 0) = -1. * sum(dPhidsigma % kappa_j[0]) + K(0, 0);
        Y_crit(0) = sigmaY;

        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);

        s_j(0) += ds_j(0);
        DEP = DEP + ds_j(0) * Lambdap;

        sigma = el_pred(L, ee_tr_v - Eth - DEP, ndi) / J;
    }

    double Dp = Ds_j(0);

    if (Dp > 0.) {
        // Exponential-map update of Fp in the relaxed configuration: the flow
        // increment is rotated by the elastic rotation Re (zero plastic spin),
        // and tr(DEP) = 0 keeps det Fp = 1 exactly.
        mat Ve_tr = sqrtmat_sympd(Be_tr);
        mat Re_tr = inv(Ve_tr) * Fe_tr;
        mat DEp_t = v2t_strain(DEP);
        mat DEp_bar = Re_tr.t() * DEp_t * Re_tr;
        Fp = expmat_sym(0.5 * (DEp_bar + DEp_bar.t())) * Fp;

        // Final state recomputed exactly from (U, Fp): the stress remains a state
        // function of (F, Fp) after plastic flow.
        Fp_inv = inv(Fp);
        mat Fe = F_tot * Fp_inv;
        mat Be = Fe * Fe.t();
        vec ee_v = t2v_strain(0.5 * logmat_sympd(Be));
        sigma = el_pred(L, ee_v - Eth, ndi) / J;
    }

    // Tangent: shared leading-mechanism assembly on the Hencky stiffness.
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma % kappa_j[0]) - K(0, 0);
    const std::vector<vec> dPhidsigma_l = {dPhidsigma};
    const ContinuumTangent ct = assemble_continuum_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L / J);
    Lt = ct.Lt;

    double A_p = -Hp;
    double Dgamma_loc = 0.5 * sum((sigma_start + sigma) % DEP) + 0.5 * (A_p_start + A_p) * Dp;

    Wm += 0.5 * sum((sigma_start + sigma) % DEtot);
    Wm_r += 0.5 * sum((sigma_start + sigma) % (DEtot - DEP));
    Wm_ir += -0.5 * (A_p_start + A_p) * Dp;
    Wm_d += Dgamma_loc;

    statev(0) = T_init;
    statev(1) = p;
    statev(2) = Fp(0, 0);
    statev(3) = Fp(0, 1);
    statev(4) = Fp(0, 2);
    statev(5) = Fp(1, 0);
    statev(6) = Fp(1, 1);
    statev(7) = Fp(1, 2);
    statev(8) = Fp(2, 0);
    statev(9) = Fp(2, 1);
    statev(10) = Fp(2, 2);
}

} //namespace simcoon
