/**
 * EPICP external UMAT plugin using Tensor2/Tensor4
 * ==================================================
 *
 * Elastic-plastic with isotropic hardening (von Mises, CCP return mapping),
 * written almost entirely with the typed Tensor2/Tensor4 API.
 *
 * Design rule — typed tensors carry the *physical* quantities, primitives carry
 * the *scalar* return map:
 *   - stress, plastic strain, elastic strain, flow normal -> ``tensor2``
 *   - stiffness and consistent tangent                    -> ``tensor4``
 *   - p, Phi, the Fischer-Burmeister 1x1 system           -> ``double`` / small vec
 *
 * The result reads like the continuum equations:
 *   ``sigma = L.contract(Eel)``            stress = L : eps_e
 *   ``Mises(sigma)``                       von Mises equivalent stress
 *   ``flow_normal(sigma)``                 yield-surface normal  N = (3/2) dev(s)/seq
 *   ``sum(N % kappa)``                     double contraction  N : kappa
 *   ``L - dyadic(kappa, kappa/Bhat)``      L - (kappa (x) kappa) / Bhat
 *
 * Overhead: the typed objects are fixed-size (stack, no heap). Construction uses
 * the ``VoigtType`` *enum* (no per-call string parsing), and the only Voigt->tensor
 * conversion sits *before* the Newton loop — the hot loop never touches a string
 * tag or rebuilds a tensor from a raw vector.
 *
 * Compare with the built-in:
 *   ``src/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.cpp``
 *
 * Material properties (6): E, nu, alpha, sigmaY, k, m
 * State variables (8): T_init, p, EP[6]
 *
 * Build (macOS example):
 *   g++ -shared -fPIC -O3 -std=c++20 \
 *       -I<simcoon_root>/include \
 *       -isystem <simcoon_build>/_deps/fastor-src \
 *       -I/opt/homebrew/include \
 *       -L/opt/homebrew/lib -L<simcoon_build> -lsimcoon -larmadillo \
 *       -o umat_plugin_ext.dylib umat_plugin_ext_EPICP.cpp
 */

#include <iostream>
#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>

using namespace arma;
using namespace simcoon;

// LIB_EXPORT comes from umat_plugin_api.hpp

class LIB_EXPORT umat_epicp_tensor : public umat_plugin_ext_api {
public:

    std::string name() const override { return "epicp_tensor"; }

    void umat_external_M(
        const std::string &umat_name,
        const vec &Etot, const vec &DEtot,
        vec &sigma, mat &Lt_out, mat &L_out,
        const mat &DR,
        const int &nprops, const vec &props,
        const int &nstatev, vec &statev,
        const double &T, const double &DT,
        const double &Time, const double &DTime,
        double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
        const int &ndi, const int &nshr,
        const bool &start, double &tnew_dt) override
    {
        UNUSED(umat_name); UNUSED(Time); UNUSED(DTime); UNUSED(DR);
        UNUSED(ndi); UNUSED(nshr); UNUSED(tnew_dt); UNUSED(nprops); UNUSED(nstatev);

        // ---- Material properties ----
        double E      = props(0);
        double nu     = props(1);
        double alpha  = props(2);
        double sigmaY = props(3);
        double k      = props(4);
        double m      = props(5);

        // ---- Initialization ----
        if (start) {
            statev(0) = T;                 // T_init
            statev(1) = 0.;                // p
            statev.subvec(2, 7).zeros();   // EP
            sigma.zeros();
            Wm = Wm_r = Wm_ir = Wm_d = 0.;
        }

        // ---- State variables ----
        double T_init = statev(0);
        double p_old  = statev(1);
        double p      = p_old;

        // ---- Typed stiffness (built once, enum tag -> no string parse) ----
        tensor4 L(L_iso(E, nu, "Enu"), Tensor4Type::stiffness);
        L_out = L.mat();

        // ---- Wrap incoming state as tensors (the only Voigt->tensor conversions) ----
        tensor2 sig       = tensor2::from_voigt(sigma, VoigtType::stress);
        tensor2 EP        = tensor2::from_voigt(vec(statev.subvec(2, 7)), VoigtType::strain);
        tensor2 sig_start = sig;
        tensor2 EP_start  = EP;

        // ---- Trial mechanical strain (built once, reused every iteration) ----
        vec alpha_v = {alpha, alpha, alpha, 0., 0., 0.};
        tensor2 eps_tot = tensor2::from_voigt(
            vec(Etot + DEtot - alpha_v * (T + DT - T_init)), VoigtType::strain);

        // ---- Hardening at the start of the step ----
        double Hp = 0., dHpdp = 0.;
        if (p > simcoon::iota) {
            dHpdp = m * k * std::pow(p, m - 1);
            Hp    = k * std::pow(p, m);
        }
        double A_p_start = -Hp;

        // ---- Elastic prediction ----
        sig = L.contract(eps_tot - EP);

        // ---- Return mapping (CCP) — scalars stay primitive, tensors stay typed ----
        tensor2 N(VoigtType::strain), kappa(VoigtType::stress);
        double error = 1.;
        for (int iter = 0;
             iter < simcoon::maxiter_umat && error > simcoon::precision_umat;
             iter++)
        {
            if (p > simcoon::iota) {
                dHpdp = m * k * std::pow(p, m - 1);
                Hp    = k * std::pow(p, m);
            }

            double Phi = Mises(sig) - Hp - sigmaY;        // von Mises yield function

            // Associative von Mises: flow direction = yield-surface normal (3/2) dev(sigma)/seq
            N = flow_normal(sig);
            kappa = L.contract(N);                        // kappa = L : N

            double B_val = -sum(N % kappa) - dHpdp;       // - N:kappa + dPhi/dp

            vec Phi_v = {Phi};
            vec Y_v   = {sigmaY};
            mat B_m(1, 1); B_m(0, 0) = B_val;
            vec Ds_j  = zeros(1), ds_j = zeros(1);
            Fischer_Burmeister_m(Phi_v, Y_v, B_m, Ds_j, ds_j, error);

            p  += ds_j(0);
            EP += ds_j(0) * N;

            sig = L.contract(eps_tot - EP);               // stress update
        }

        // ---- Consistent tangent modulus:  Lt = L - (kappa (x) kappa) / Bhat ----
        if (p > simcoon::iota && (p - p_old) > simcoon::iota) {
            double Bhat = sum(N % kappa) + dHpdp;
            Lt_out = (L - dyadic(kappa, kappa * (1.0 / Bhat))).mat();
        } else {
            Lt_out = L_out;
        }

        // ---- Work (typed double contractions  sigma : d_eps) ----
        tensor2 sig_sum = sig_start + sig;
        tensor2 DEP     = EP - EP_start;
        tensor2 DEtot_t = tensor2::from_voigt(DEtot, VoigtType::strain);
        double  Dp      = p - p_old;
        double  A_p     = -Hp;

        Wm    += 0.5 * sum(sig_sum % DEtot_t);
        Wm_r  += 0.5 * sum(sig_sum % (DEtot_t - DEP));
        Wm_ir += -0.5 * (A_p_start + A_p) * Dp;
        Wm_d  += 0.5 * sum(sig_sum % DEP) + 0.5 * (A_p_start + A_p) * Dp;

        // ---- Write back ----
        sigma = sig.voigt();
        statev(0) = T_init;
        statev(1) = p;
        statev.subvec(2, 7) = EP.voigt();
    }
};

extern "C" LIB_EXPORT umat_plugin_ext_api* create_api() {
    return new umat_epicp_tensor();
}

extern "C" LIB_EXPORT void destroy_api(umat_plugin_ext_api* p) {
    delete p;
}
