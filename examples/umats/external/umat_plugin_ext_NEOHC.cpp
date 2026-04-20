/**
 * Neo-Hookean external UMAT plugin using Tensor2/Tensor4
 * ========================================================
 *
 * This plugin implements the compressible Neo-Hookean hyperelastic law
 * using the Tensor2/Tensor4 API.  It can be loaded by simcoon's solver
 * as an external UMAT (umat name ``"UMEXT"``).
 *
 * Strain energy (compressible Neo-Hookean, isochoric + volumetric split):
 *
 * .. math::
 *    W = \\frac{\\mu}{2}(I_1 - 3) - \\mu\\ln J + \\frac{\\lambda}{2}(\\ln J)^2
 *
 * PKII stress:  S = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1}
 * Material tangent: L_IJKL = lambda*C^{-1}_IJ C^{-1}_KL
 *                          + 2*(mu - lambda*ln J) * d(C^{-1})/dC_sym
 * Cauchy stress and spatial tangent obtained via Tensor2/Tensor4 push_forward.
 *
 * Material properties (3):
 *   props[0] = E       (Young's modulus)
 *   props[1] = nu      (Poisson's ratio)
 *   props[2] = alpha   (CTE)
 *
 * State variables (1):  statev[0] = T_init
 *
 * Build (macOS example):
 *   g++ -shared -fPIC -O3 -std=c++20 \
 *       -I<simcoon_root>/include \
 *       -isystem <simcoon_build>/_deps/fastor-src \
 *       -I/opt/homebrew/include \
 *       -L/opt/homebrew/lib -L<simcoon_build> -lsimcoon -larmadillo \
 *       -o umat_plugin_ext.dylib umat_plugin_ext.cpp
 *
 * Then run from the examples/umats/ directory (the solver looks for
 * ``external/umat_plugin_ext.<ext>``).
 */

#include <iostream>
#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
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

class LIB_EXPORT umat_neohookean_tensor : public umat_plugin_ext_api {
public:

    std::string name() const { return "neohookean_tensor"; }

    void umat_external_M(
        const std::string &umat_name,
        const vec &Etot, const vec &DEtot,
        vec &sigma_v, mat &Lt_v, mat &L_v,
        const mat &DR,
        const int &nprops, const vec &props,
        const int &nstatev, vec &statev,
        const double &T, const double &DT,
        const double &Time, const double &DTime,
        double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
        const int &ndi, const int &nshr,
        const bool &start, double &tnew_dt) override
    {
        UNUSED(umat_name); UNUSED(Etot); UNUSED(DR);
        UNUSED(nprops); UNUSED(nstatev);
        UNUSED(Time); UNUSED(DTime); UNUSED(nshr); UNUSED(tnew_dt);

        double T_init = statev(0);
        double E     = props(0);
        double nu    = props(1);
        double alpha = props(2);
        double mu     = E / (2.0 * (1.0 + nu));
        double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

        if (start) {
            T_init = T;
            sigma_v = zeros(6);
            Wm = Wm_r = Wm_ir = Wm_d = 0.;
        }

        vec sigma_start = sigma_v;

        // ---- Kinematics ----
        // NOTE: The generic external API does not provide F0/F1 directly.
        // For a real finite-strain plugin, store F in statev or use the
        // Abaqus-compatible API (umat_plugin_aba_api) which provides F
        // through phase_characteristics.
        //
        // Here we reconstruct F from the Green-Lagrange strain E = 0.5*(C-I):
        // C = I + 2*(Etot+DEtot),  F = sqrt(C) via eigendecomposition.
        vec E_GL = Etot + DEtot;
        mat C_mat = eye(3, 3);
        C_mat(0, 0) += 2.0 * E_GL(0);
        C_mat(1, 1) += 2.0 * E_GL(1);
        C_mat(2, 2) += 2.0 * E_GL(2);
        C_mat(0, 1) += E_GL(3); C_mat(1, 0) += E_GL(3);
        C_mat(0, 2) += E_GL(4); C_mat(2, 0) += E_GL(4);
        C_mat(1, 2) += E_GL(5); C_mat(2, 1) += E_GL(5);

        double J = std::sqrt(det(C_mat));
        double lnJ = std::log(J);
        mat invC_mat = inv(C_mat);

        // ---- PKII stress using Tensor2 ----
        tensor2 invC_t = tensor2::from_voigt(t2v_stress(invC_mat), "stress");
        tensor2 I_t    = tensor2::identity("stress");
        tensor2 S = mu * (I_t - invC_t) + lambda * lnJ * invC_t;

        // ---- Cauchy stress via push-forward ----
        // Reconstruct F from C via eigendecomposition: F = R * sqrt(D) * V^T
        vec eigval;
        mat eigvec;
        eig_sym(eigval, eigvec, C_mat);
        mat F1 = eigvec * diagmat(sqrt(eigval)) * eigvec.t();

        tensor2 sigma_cauchy = S.push_forward(mat::fixed<3,3>(F1));
        sigma_v = sigma_cauchy.voigt();

        // ---- Material tangent using Tensor4 ----
        tensor4 L_ref(
            lambda * auto_dyadic(invC_mat)
          + 2.0 * (mu - lambda * lnJ) * dinvSdSsym(C_mat),
            "stiffness");
        L_v = L_ref.mat();

        // ---- Spatial tangent via push_forward ----
        tensor4 Lt = L_ref.push_forward(mat::fixed<3,3>(F1));
        Lt_v = Lt.mat();

        // ---- Work ----
        Wm   += 0.5 * sum((sigma_start + sigma_v) % DEtot);
        Wm_r += 0.5 * sum((sigma_start + sigma_v) % DEtot);
        Wm_ir = Wm_d = 0.;

        statev(0) = T_init;
    }
};

extern "C" LIB_EXPORT umat_plugin_ext_api* create_api() {
    return new umat_neohookean_tensor();
}

extern "C" LIB_EXPORT void destroy_api(umat_plugin_ext_api* p) {
    delete p;
}
