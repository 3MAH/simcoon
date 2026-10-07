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

/**
 * @file umat_tutorial_J2.cpp
 * @brief TUTORIAL UMAT — J2 plasticity, linear isotropic hardening, typed
 * tensors. Read alongside docs/simulation/umat_tutorial.rst.
 */

#include "umat_tutorial_J2.hpp"

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace arma;

namespace simcoon {

void umat_tutorial_J2(const std::string & /*umat_name*/, const vec &Etot,
                      const vec &DEtot, vec &sigma, mat &Lt, mat &L,
                      const mat &DR, const int & /*nprops*/, const vec &props,
                      const int & /*nstatev*/, vec &statev, const double &T,
                      const double &DT, const double & /*Time*/,
                      const double & /*DTime*/, double &Wm, double &Wm_r,
                      double &Wm_ir, double &Wm_d, const int & /*ndi*/,
                      const int & /*nshr*/, const bool &start,
                      double & /*tnew_dt*/, const int &tangent_mode) {

    // ------------------------------------------------------------------ 1.
    // Material properties. Units: MPa; Voigt order [11, 22, 33, 12, 13, 23].
    const double E = props(0);
    const double nu = props(1);
    const double alpha = props(2);   // CTE (isotropic)
    const double sigma_Y = props(3);
    const double H = props(4);       // linear hardening modulus: R(p) = H p

    const double mu = E / (2. * (1. + nu));

    // ------------------------------------------------------------------ 2.
    // Elasticity as a TYPED fourth-order tensor: the Tensor4Type carries the
    // stiffness convention, so contractions below produce stress-typed
    // results automatically.
    const tensor4 L_el(L_iso(E, nu, "Enu"), Tensor4Type::stiffness);
    L = L_iso(E, nu, "Enu");

    // ------------------------------------------------------------------ 3.
    // State variables: statev = [T_init, p, EP(6, strain Voigt)] — the EPICP
    // layout. On the start increment the state is virgin.
    double T_init = statev(0);
    double p = statev(1);
    tensor2 EP = strain(vec(statev.subvec(2, 7)));

    if (start) {
        T_init = T;
        p = 0.;
        EP = tensor2::zeros(Tensor2Type::strain);
        sigma = zeros(6);
        Wm = 0.; Wm_r = 0.; Wm_ir = 0.; Wm_d = 0.;
    }

    // Objectivity: the stored plastic strain co-rotates with the material
    // frame. The typed rotate() picks the strain rotation kernel (factor-2
    // shear convention) from the tensor's own Tensor2Type.
    EP = EP.rotate(Rotation::from_matrix(DR));

    const vec sigma_start = sigma;

    // ------------------------------------------------------------------ 4.
    // Elastic prediction (total form): the trial state freezes plastic flow.
    //   eps_el^trial = Etot + DEtot - eps_th - EP
    //   sigma^trial  = L : eps_el^trial
    const tensor2 eps_th =
        strain(vec(alpha * (T + DT - T_init) * Ith()));
    const tensor2 eps_el_trial =
        strain(vec(Etot + DEtot)) - eps_th - EP;
    const tensor2 sig_trial = L_el.contract(eps_el_trial);  // stress-typed

    // ------------------------------------------------------------------ 5.
    // Yield check:  f^trial = sigma_eq(sigma^trial) - (sigma_Y + H p).
    // Mises() dispatches on the type tag (sqrt(3/2 s:s) for a stress).
    const double f_trial = Mises(sig_trial) - (sigma_Y + H * p);

    double Dp = 0.;
    tensor2 Lam = tensor2::zeros(Tensor2Type::strain);

    if (f_trial > 0.) {
        // -------------------------------------------------------------- 6.
        // RADIAL RETURN, closed form. The flow direction
        //   Lambda = (3/2) dev(sigma)/sigma_eq
        // is preserved by the return (the correction is radial in deviatoric
        // space), and with LINEAR hardening the consistency condition
        //   f(Dp) = f^trial - (3 mu + H) Dp = 0
        // is linear in Dp — no local Newton iteration is needed:
        Lam = flow_normal(sig_trial);            // strain-typed
        Dp = f_trial / (3. * mu + H);

        p += Dp;
        EP += Dp * Lam;
    }

    // Converged stress from the total form (elastic when Dp = 0).
    const tensor2 sig = L_el.contract(eps_el_trial - Dp * Lam);
    sigma = sig.to_arma_voigt();

    // ------------------------------------------------------------------ 7.
    // TANGENT MODULUS. For the continuum operator everything is analytic:
    //
    //   L_t = L - (L:Lambda) o (L:Lambda) / (Lambda:L:Lambda + H)
    //
    // and for isotropic L with a deviatoric Lambda:
    //   L:Lambda        = 2 mu Lambda           (kappa)
    //   Lambda:L:Lambda = 3 mu
    //   =>  L_t = L - 4 mu^2 / (3 mu + H) * Lambda o Lambda.
    //
    // The code keeps the GENERIC kappa-form so it survives a change of yield
    // criterion; the specialization above is what the docs derive by hand.
    // NB: auto_dyadic is the FULL outer product kappa (x) kappa — sym_dyadic
    // is the (ik)(jl)-symmetrized product, a different operator.
    if (Dp <= simcoon::iota || tangent_mode == tangent_none) {
        Lt = L;  // elastic step, or explicit integration
    } else if (tangent_mode == tangent_continuum) {
        const tensor2 kappa = L_el.contract(Lam);            // stress-typed
        const double Lam_L_Lam = accu(Lam.mat() % kappa.mat());  // = 3 mu
        const tensor4 Lt_t =
            L_el - auto_dyadic(kappa) * (1. / (Lam_L_Lam + H));
        Lt = mat(Lt_t.mat());
    } else {
        // Algorithmic (Simo-Hughes) and beyond: hand the physical inputs to
        // the shared dispatch — Bhat = Lambda:kappa + H, the flow Hessian is
        // dLambda/dsigma = deta_stress. This is what the production kernels do.
        const tensor2 kappa = L_el.contract(Lam);
        const double Bhat = accu(Lam.mat() % kappa.mat()) + H;
        const vec sig_v = sigma;
        const ContinuumTangent ct = compute_tangent_operator(
            tangent_mode, Bhat, kappa.to_arma_voigt(), Lam.to_arma_voigt(),
            Dp, L, [&sig_v]() { return deta_stress(sig_v); });
        Lt = ct.Lt;
    }

    // ------------------------------------------------------------------ 8.
    // Work quantities — CUMULATIVE (the solver reports path integrals):
    // trapezoidal total work, plastic dissipation sigma:Dp*Lambda, the
    // recoverable part closes the balance.
    const double Wm_inc = 0.5 * dot(sigma_start + sigma, DEtot);
    const double Wm_d_inc =
        0.5 * accu((stress(vec(sigma_start)).mat() + sig.mat()) % (Dp * Lam).mat());
    Wm += Wm_inc;
    Wm_d += Wm_d_inc;
    Wm_r += Wm_inc - Wm_d_inc;
    Wm_ir += 0.;

    // ------------------------------------------------------------------ 9.
    // Pack the state back.
    statev(0) = T_init;
    statev(1) = p;
    statev.subvec(2, 7) = EP.to_arma_voigt();
}

} // namespace simcoon
