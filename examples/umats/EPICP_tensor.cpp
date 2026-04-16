/**
 * EPICP External UMAT Plugin using Tensor2/Tensor4
 * ==================================================
 *
 * Same elastic-plastic (isotropic hardening, von Mises, CCP return mapping)
 * as the built-in EPICP, but using the Tensor2/Tensor4 API at the interface
 * level for type safety.
 *
 * The return mapping inner loop stays in raw Voigt vectors — that's where
 * the algorithm naturally lives and where raw arma::vec is fastest.
 * Tensor2/Tensor4 add value at the boundaries:
 *   - Typed stiffness: tensor4 L("stiffness") ensures correct contraction
 *   - Typed stress/strain: from_voigt(v, "stress") makes convention explicit
 *   - Mises/dev dispatch correctly based on type tag
 *   - Tangent contraction: L @ eps -> sigma with automatic type inference
 *
 * Compare with the built-in version in
 *   src/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.cpp
 *
 * Material properties (6): E, nu, alpha, sigmaY, k, m
 * State variables (8): T_init, p, EP[6]
 *
 * Build as plugin:
 *   g++ -shared -fPIC -O3 -std=c++20 \
 *       -I<simcoon_include> -I<armadillo_include> -I<fastor_include> \
 *       -L<simcoon_lib> -lsimcoon \
 *       -o external/umat_plugin_ext.so EPICP_tensor.cpp
 */

#include <iostream>
#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>

using namespace arma;
using namespace simcoon;

#if defined(_WIN32) || defined(_WIN64)
    #define LIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
    #if __GNUC__ >= 4
        #define LIB_EXPORT __attribute__((visibility("default")))
    #else
        #define LIB_EXPORT
    #endif
#else
    #define LIB_EXPORT
#endif

class LIB_EXPORT umat_epicp_tensor : public umat_plugin_ext_api {
public:

    std::string name() const { return "epicp_tensor"; }

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
        UNUSED(umat_name); UNUSED(Time); UNUSED(DTime);
        UNUSED(nshr); UNUSED(tnew_dt); UNUSED(nprops); UNUSED(nstatev);

        // ---- Material properties ----
        double E      = props(0);
        double nu     = props(1);
        double alpha  = props(2);
        double sigmaY = props(3);
        double k      = props(4);
        double m      = props(5);

        // ---- State variables ----
        double T_init = statev(0);
        double p      = statev(1);
        vec EP        = statev.subvec(2, 7);

        // ---- Typed stiffness tensor ----
        tensor4 L(L_iso(E, nu, "Enu"), "stiffness");
        L_out = L.mat();

        // ---- Initialization ----
        if (start) {
            T_init = T;
            sigma = zeros(6);
            EP = zeros(6);
            p = 0.;
            Wm = Wm_r = Wm_ir = Wm_d = 0.;
        }

        // ---- Hardening ----
        double Hp = 0., dHpdp = 0.;
        if (p > simcoon::iota) {
            dHpdp = m * k * std::pow(p, m - 1);
            Hp = k * std::pow(p, m);
        }

        // Save start values
        vec sigma_start = sigma;
        vec EP_start = EP;
        double A_p_start = -Hp;

        // ---- Elastic prediction ----
        // Build elastic strain as typed tensor2, then contract with stiffness
        vec alpha_v = zeros(6);
        alpha_v(0) = alpha_v(1) = alpha_v(2) = alpha;

        tensor2 Eel = tensor2::from_voigt(Etot + DEtot - alpha_v * (T + DT - T_init) - EP, "strain");
        tensor2 sigma_t = L.contract(Eel);  // stiffness * strain -> stress (automatic)
        sigma = sigma_t.voigt();

        // ---- Yield check ----
        double sigma_eq = Mises(sigma_t);     // type-aware: uses stress convention
        double Phi = sigma_eq - Hp - sigmaY;

        // ---- Return mapping (CCP) — raw Voigt arithmetic ----
        // The inner loop stays in vec space: this is pure numerical optimization,
        // not continuum mechanics. Raw arma::vec is the right tool here.
        vec Lambdap = eta_stress(sigma);
        vec kappa = L_out * Lambdap;

        double error = 1.;
        for (int iter = 0; iter < simcoon::maxiter_umat && error > simcoon::precision_umat; iter++) {

            if (p > simcoon::iota) {
                dHpdp = m * k * std::pow(p, m - 1);
                Hp = k * std::pow(p, m);
            }

            vec dPhidsigma = eta_stress(sigma);
            Phi = Mises_stress(sigma) - Hp - sigmaY;

            Lambdap = eta_stress(sigma);
            kappa = L_out * Lambdap;

            double dPhidp = -dHpdp;
            double B_val = -sum(dPhidsigma % kappa) + dPhidp;
            double Y_crit = sigmaY;

            // Fischer-Burmeister scalar update
            vec Phi_v = {Phi};
            vec Y_v = {Y_crit};
            mat B_m = {{B_val}};
            vec Ds_j = zeros(1), ds_j = zeros(1);
            Fischer_Burmeister_m(Phi_v, Y_v, B_m, Ds_j, ds_j, error);

            p += ds_j(0);
            EP = EP + ds_j(0) * Lambdap;

            // Recompute stress via typed contraction
            Eel = tensor2::from_voigt(Etot + DEtot - alpha_v * (T + DT - T_init) - EP, "strain");
            sigma_t = L.contract(Eel);
            sigma = sigma_t.voigt();
        }

        // ---- Consistent tangent modulus ----
        if (p > simcoon::iota && (p - statev(1)) > simcoon::iota) {
            vec dPhidsigma = eta_stress(sigma);
            kappa = L_out * dPhidsigma;
            double Bhat = sum(dPhidsigma % kappa) - (-dHpdp);
            vec P_eps = (1.0 / Bhat) * kappa;
            Lt_out = L_out - kappa * P_eps.t();
        } else {
            Lt_out = L_out;
        }

        // ---- Work computation ----
        vec DEP = EP - EP_start;
        double Dp = p - statev(1);
        double A_p = -Hp;

        Wm   += 0.5 * sum((sigma_start + sigma) % DEtot);
        Wm_r += 0.5 * sum((sigma_start + sigma) % (DEtot - DEP));
        Wm_ir += -0.5 * (A_p_start + A_p) * Dp;
        Wm_d += 0.5 * sum((sigma_start + sigma) % DEP) + 0.5 * (A_p_start + A_p) * Dp;

        // ---- Update state variables ----
        statev(0) = T_init;
        statev(1) = p;
        statev.subvec(2, 7) = EP;
    }
};

extern "C" LIB_EXPORT umat_plugin_ext_api* create_api() {
    return new umat_epicp_tensor();
}

extern "C" LIB_EXPORT void destroy_api(umat_plugin_ext_api* p) {
    delete p;
}
