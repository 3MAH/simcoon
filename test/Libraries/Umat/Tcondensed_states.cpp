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

///@file Tcondensed_states.cpp
///@brief Condensed stress states (ndi, nshr) driven through a real kernel (EPICP): truss (1, 0)
///       and beam fibre (1, 2). The kernel only sees the kept strain components; the condensed
///       ones are stress-free and their plastic strain is free, so the converged state is the
///       exact uniaxial / tension-torsion J2 state.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <cmath>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// EPICP: props = [E, nu, alpha, sigmaY, k, m]; statev = [T_init, p, EP(6)].
const double E_ = 70000., nu_ = 0.3, sigmaY_ = 200., k_ = 1000., m_ = 0.3;
const double G_ = E_ / (2. * (1. + nu_));

struct EpicpState {
    vec Etot = zeros(6), sigma = zeros(6), statev = zeros(8);
    mat Lt = zeros(6, 6), L = zeros(6, 6);
    bool start = true;
};

// One increment of EPICP under the condensed state (ndi, nshr). Only the kept components of
// DEtot are meaningful; the condensed ones are passed as zero, as an element would.
void step(EpicpState &s, const vec &DEtot, int ndi, int nshr, int tangent_mode = tangent_default) {
    const vec props = {E_, nu_, 0., sigmaY_, k_, m_};
    mat DR = eye(3, 3);
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0., tnew_dt = 1.;
    umat_plasticity_iso("EPICP", s.Etot, DEtot, s.sigma, s.Lt, s.L, DR, props.n_elem, props,
                        s.statev.n_elem, s.statev, 293.15, 0., 0., 1., Wm, Wm_r, Wm_ir, Wm_d,
                        ndi, nshr, s.start, tnew_dt, tangent_mode);
    ASSERT_DOUBLE_EQ(tnew_dt, 1.);
    s.Etot += DEtot;
    s.start = false;
}

double yield_stress(double p) { return sigmaY_ + k_ * pow(p, m_); }

}  // namespace

// Truss (1, 0): elastic sigma = E eps exactly, then the uniaxial J2 state — sigma on the
// hardening curve, eps_p11 = p, incompressible transverse plastic strain, every other stress zero.
TEST(Tcondensed_states, truss_EPICP_uniaxial)
{
    EpicpState s;
    vec DE = zeros(6);
    DE(0) = 1.e-3;   // elastic: sigma_Y / E = 2.86e-3
    step(s, DE, 1, 0);
    EXPECT_NEAR(s.sigma(0), E_ * 1.e-3, 1.e-9);
    EXPECT_LT(norm(s.sigma.subvec(1, 5), 2), 1.e-12);

    DE(0) = 1.e-3;
    for (int i = 0; i < 20; i++) step(s, DE, 1, 0);   // eps11 = 2.1 %

    const double p = s.statev(1);
    const vec EP = s.statev.subvec(2, 7);
    EXPECT_GT(p, 1.e-3);
    EXPECT_NEAR(s.sigma(0), yield_stress(p), 1.e-3);
    EXPECT_NEAR(s.sigma(0), E_ * (s.Etot(0) - EP(0)), 1.e-6);
    EXPECT_NEAR(EP(0), p, 1.e-9);
    EXPECT_NEAR(EP(1), -0.5 * p, 1.e-9);
    EXPECT_NEAR(EP(2), -0.5 * p, 1.e-9);
    EXPECT_LT(norm(s.sigma.subvec(1, 5), 2), 1.e-12);
}

// Beam fibre (1, 2): proportional tension-torsion. The kept stresses are sigma11, sigma12,
// sigma13; the J2 state has sqrt(sigma^2 + 3 tau^2) on the hardening curve and each kept
// stress given by its own elastic modulus (E, G) on the elastic strain.
TEST(Tcondensed_states, beam_EPICP_tension_torsion)
{
    for (int tangent_mode : {tangent_default, tangent_closest_point}) {
        EpicpState s;
        vec DE = zeros(6);
        DE(0) = 1.e-3;
        DE(3) = 1.e-3;   // engineering shear gamma_12
        step(s, DE, 1, 2, tangent_mode);   // elastic
        EXPECT_NEAR(s.sigma(0), E_ * 1.e-3, 1.e-9);
        EXPECT_NEAR(s.sigma(3), G_ * 1.e-3, 1.e-9);

        for (int i = 0; i < 20; i++) step(s, DE, 1, 2, tangent_mode);

        const double p = s.statev(1);
        const vec EP = s.statev.subvec(2, 7);
        const double sig = s.sigma(0), tau = s.sigma(3);
        EXPECT_GT(p, 1.e-3);
        EXPECT_NEAR(sqrt(sig * sig + 3. * tau * tau), yield_stress(p), 1.e-3);
        EXPECT_NEAR(Mises_stress(s.sigma), yield_stress(p), 1.e-3);
        EXPECT_NEAR(sig, E_ * (s.Etot(0) - EP(0)), 1.e-6);
        EXPECT_NEAR(tau, G_ * (s.Etot(3) - EP(3)), 1.e-6);
        EXPECT_NEAR(EP(1), EP(2), 1.e-12);   // transverse isotropy of the fibre
        EXPECT_NEAR(s.sigma(4), 0., 1.e-12);   // no gamma_13 imposed
        EXPECT_NEAR(s.sigma(1), 0., 1.e-12);
        EXPECT_NEAR(s.sigma(2), 0., 1.e-12);
        EXPECT_NEAR(s.sigma(5), 0., 1.e-12);
    }
}
