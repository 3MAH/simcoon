/*
 * PROTOTYPE: Elastic-plastic UMAT using tensor2 / tensor4
 * ========================================================
 *
 * Side-by-side comparison with the current raw arma::vec / arma::mat style.
 * NOT meant to compile or replace anything — just to show how the code reads.
 */

#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace arma;
using namespace simcoon;


// ============================================================================
//  Isotropic elasticity — 3 lines of mechanics
// ============================================================================
//
//  BEFORE (raw arma):
//      L = L_iso(E, nu, "Enu");                         // what type? who knows
//      vec Eel = Etot + DEtot - alpha*(T+DT-T_init);    // stress voigt? strain voigt?
//      sigma = el_pred(L, Eel, ndi);                     // hope you got the convention right
//      Lt = L;
//
//  AFTER (typed):

void umat_elasticity_iso_typed(
    const vec &Etot, const vec &DEtot,
    vec &sigma_out, mat &Lt_out, mat &L_out,
    const vec &props, double T_init)
{
    // --- The actual mechanics: 3 lines ---
    tensor4 L(L_iso(props(0), props(1), "Enu"), Tensor4Type::stiffness);
    tensor2 sigma = L * strain(Etot + DEtot - props(2)*(T_init)));
    Lt_out = L_out = L.mat();

    sigma_out = sigma.voigt();
}


// ============================================================================
//  NEW STYLE — elastic-plastic with isotropic hardening (CCP)
// ============================================================================

void umat_plasticity_iso_CCP_typed(
    const vec &Etot, const vec &DEtot,
    vec &sigma_out, mat &Lt_out, mat &L_out,
    vec &sigma_in,
    const mat &DR,
    const vec &props, vec &statev,
    const double &T, const double &DT,
    const bool &start, const int &solver_type)
{
    double E        = props(0);
    double nu       = props(1);
    double alpha_v  = props(2);
    double sigmaY   = props(3);
    double k        = props(4);
    double m        = props(5);

    double T_init = statev(0);
    double p      = statev(1);

    // --- Typed stiffness ---
    tensor4 L(L_iso(E, nu, "Enu"), Tensor4Type::stiffness);

    // --- Plastic strain: tagged as VoigtType::strain ---
    // Can't accidentally rotate it with rotate_stress anymore.
    tensor2 EP = tensor2::from_voigt(
        vec::fixed<6>(statev.memptr() + 2), VoigtType::strain
    );

    // --- Rotation of internal variables ---
    // EP.rotate() dispatches to apply_strain automatically.
    // With raw code you'd write:  EP = rotate_strain(EP, DR);
    //   ... and if you wrote rotate_stress(EP, DR) by mistake: silent bug.
    Rotation dR = Rotation::from_matrix(DR);
    EP = EP.rotate(dR);

    if (start) {
        T_init = T;
        EP = tensor2::zeros(VoigtType::strain);
        p = 0.0;
    }

    // --- CTE as a strain tensor ---
    tensor2 alpha_tensor = tensor2::identity(VoigtType::strain) * alpha_v;

    // --- Elastic predictor ---
    tensor2 Etot_t = strain(Etot + DEtot);

    tensor2 thermal = alpha_tensor * (T + DT - T_init);
    tensor2 Eel = Etot_t - thermal - EP;

    // --- sigma = L : Eel ---
    // Output is automatically VoigtType::stress. No ambiguity.
    tensor2 sigma = L * Eel;

    // --- Yield function ---
    double Hp = (p > 1e-12) ? k * std::pow(p, m) : 0.0;
    double dHpdp = (p > 1e-12) ? m * k * std::pow(p, m - 1) : 0.0;

    double phi = Mises(sigma) - Hp - sigmaY;

    // --- Return mapping (CCP loop) ---
    tensor2 sigma_start = sigma;
    tensor2 EP_start = EP;

    if (phi > 0.0) {
        // Flow direction: dev(sigma) / ||dev(sigma)||
        // This is a stress-like quantity — dev() returns stress Voigt
        vec::fixed<6> n_v = dev(sigma);
        double norm_n = arma::norm(n_v);
        if (norm_n > 1e-14) n_v /= norm_n;

        // Flow direction as stress tensor
        tensor2 N = tensor2::from_voigt(n_v, VoigtType::stress);

        // kappa = L : N  → gives stress (stiffness contracts stress-normal → stress)
        tensor2 kappa = L.contract(
            tensor2::from_voigt(n_v, VoigtType::strain)  // n as strain for L:n
        );

        // Newton iterations on the plastic multiplier dp
        double dp = 0.0;
        for (int iter = 0; iter < 25; ++iter) {
            Hp = (p + dp > 1e-12) ? k * std::pow(p + dp, m) : 0.0;
            dHpdp = (p + dp > 1e-12) ? m * k * std::pow(p + dp, m - 1) : 0.0;

            phi = Mises(sigma) - Hp - sigmaY;
            if (std::abs(phi) < 1e-10) break;

            double denom = arma::dot(vec(dev(sigma)) / Mises(sigma), vec(kappa.voigt())) + dHpdp;
            double ddp = phi / denom;
            dp += ddp;

            // Update plastic strain and stress
            // EP += ddp * N  — both are strain-tagged, arithmetic preserves tag
            EP = EP + tensor2::from_voigt(vec::fixed<6>(n_v * ddp), VoigtType::strain);

            Eel = Etot_t - thermal - EP;
            sigma = L * Eel;
        }
        p += dp;
    }

    // --- Tangent modulus ---
    if (solver_type == 0 || solver_type == 2) {
        if (phi > -1e-10 && Mises(sigma) > 1e-10) {
            // Elastoplastic tangent
            vec::fixed<6> n_v = dev(sigma);
            double norm_n = arma::norm(n_v);
            if (norm_n > 1e-14) n_v /= norm_n;

            // kappa = L : n (as strain)
            vec::fixed<6> kappa_v = L.mat() * n_v;

            // Denominator: n : L : n + dH/dp
            double denom = arma::dot(n_v, kappa_v) + dHpdp;

            // Lt = L - (kappa * kappa^T) / denom
            // Type is preserved: stiffness - stiffness = stiffness
            tensor4 correction(
                mat::fixed<6,6>(kappa_v * kappa_v.t() / denom),
                Tensor4Type::stiffness
            );
            tensor4 Lt = L - correction;
            Lt_out = Lt.mat();
        } else {
            Lt_out = L.mat();
        }
    } else if (solver_type == 1) {
        sigma_in = -L.mat() * EP.voigt();
    }

    // --- Export ---
    sigma_out = sigma.voigt();
    L_out = L.mat();

    statev(0) = T_init;
    statev(1) = p;
    vec::fixed<6> ep_v = EP.voigt();
    for (int i = 0; i < 6; ++i) statev(2 + i) = ep_v(i);
}


// ============================================================================
//  Neo-Hookean (finite strain) — 5 lines of mechanics
// ============================================================================
//
//  BEFORE: 15 lines of FTensor boilerplate + manual /J everywhere
//  AFTER:  push_forward(F) includes the 1/J metric by default

void umat_neohookean_typed(
    const mat::fixed<3,3> &F,
    vec &sigma_out, mat &c_out,
    const vec &props)
{
    double mu  = props(0) / (2.0 * (1.0 + props(1)));
    double lam = props(0) * props(1) / ((1.0 + props(1)) * (1.0 - 2.0 * props(1)));
    double J   = arma::det(F);

    // --- 5 lines: stress + tangent ---
    tensor2 I = tensor2::identity(VoigtType::stress);
    tensor2 b = stress(F * F.t());
    tensor2 sigma = (mu * (b - I) + (lam * std::log(J)) * I) * (1.0 / J);

    tensor4 C_ref = lam * tensor4::volumetric() + 2.0 * mu * tensor4::identity();
    c_out = C_ref.push_forward(F).mat();       // metric=true: includes 1/J

    sigma_out = sigma.voigt();
}


// ============================================================================
//  Rotated anisotropic elasticity — rotation dispatch is automatic
// ============================================================================
//
//  BEFORE: you must remember QS for stiffness, QE for compliance
//  AFTER:  .rotate(R) picks the right rule from the type tag

void umat_elastic_aniso_typed(
    const vec &Etot, const vec &DEtot,
    vec &sigma_out, mat &Lt_out, mat &L_out,
    const vec &props, const double &theta_deg)
{
    tensor4 L_mat(L_isotrans(props(0), props(1), props(2), props(3), props(4), 1),
                  Tensor4Type::stiffness);

    Rotation R = Rotation::from_axis_angle(theta_deg * datum::pi / 180.0, 3);
    tensor4 L_glob = L_mat.rotate(R);   // dispatches to QS * L * QS^T

    // compliance rotation uses QE automatically — no manual dispatch
    tensor4 M_glob = L_mat.inverse().rotate(R);  // QE * M * QE^T

    tensor2 sigma = L_glob * strain(Etot + DEtot));
    sigma_out = sigma.voigt();
    Lt_out = L_out = L_glob.mat();
}
