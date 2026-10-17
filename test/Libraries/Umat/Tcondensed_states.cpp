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
///@brief Condensed stress states (ndi, nshr) through the Abaqus route: abaqus2smart, select_umat_M
///       / select_umat_T (the condensed-state driver around the 3D law), smart2abaqus. The element
///       passes only its ntens components and stores no strain outside them; the converged point
///       is the 3D law at the strain that makes the other components stress-free, and DDSDDE
///       (DDSDDT, DRPLDE, DRPLDT) its exact derivatives.
///@version 2.0

#include <gtest/gtest.h>
#include <armadillo>
#include <cmath>
#include <string>
#include <vector>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>
#include <simcoon/Continuum_mechanics/Umat/condensed_state.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// EPICP: props = [E, nu, alpha, sigmaY, k, m]; statev = [T_init, p, EP(6)].
const double E_ = 70000., nu_ = 0.3, sigmaY_ = 200., k_ = 1000., m_ = 0.3;
const double G_ = E_ / (2. * (1. + nu_));
const vec EPICP_props = {E_, nu_, 0., sigmaY_, k_, m_};
// ZENER (Zener_fast): props = [E0, nu0, alpha, E1, nu1, etaB1, etaS1]; 8 statev
const vec ZENER_props = {3000., 0.35, 0., 1500., 0.35, 3000., 1200.};

// What an Abaqus element keeps at a material point: its ntens components and STATEV (Wm first)
struct element_point {
    int ndi, nshr, ntens;
    vector<double> stress, stran, statev;
    element_point(int ndi_, int nshr_, int nstatev_law)
        : ndi(ndi_), nshr(nshr_), ntens(ndi_ + nshr_), stress(ndi_ + nshr_, 0.), stran(ndi_ + nshr_, 0.),
          statev(nstatev_law + 4, 0.) {}
};

// One mechanical increment as umat_singleM runs it. commit = false leaves the point untouched
// (finite-difference probes); ddsdde (ntens x ntens, column-major) and the stress are returned.
void increment_M(const string &name, const vec &props, element_point &pt, const vector<double> &dstran,
                 const double &time, const double &dtime, vector<double> &ddsdde, vector<double> &stress,
                 const bool &commit = true) {
    const int nstatev = pt.statev.size();
    vector<double> statev = pt.statev;
    stress = pt.stress;
    ddsdde.assign(pt.ntens*pt.ntens, 0.);
    const double time_arr[2] = {time, time};
    const double drot[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    phase_characteristics rve;
    rve.construct(0, 1);
    auto sv = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    sv->resize(nstatev - 4);
    rve.sptr_matprops->resize(props.n_elem);
    vec props_smart = zeros(props.n_elem);
    double Time = 0., DTime = 0.;
    mat DR = eye(3, 3);
    bool start = false;
    abaqus2smart_M(stress.data(), ddsdde.data(), pt.stran.data(), dstran.data(), time_arr, dtime, 293.15, 0.,
                   props.n_elem, props.memptr(), nstatev, statev.data(), pt.ndi, pt.nshr, drot, sv->sigma, sv->Lt,
                   sv->Etot, sv->DEtot, sv->T, sv->DT, Time, DTime, props_smart, sv->Wm, sv->statev, DR, start);
    rve.sptr_matprops->update(0, name, 1., 0., 0., 0., props.n_elem, props_smart);
    double pnewdt = 1.;
    select_umat_M(rve, DR, Time, DTime, pt.ndi, pt.nshr, start, 0, pnewdt);
    EXPECT_GE(pnewdt, 1.);
    smart2abaqus_M(stress.data(), ddsdde.data(), statev.data(), pt.ndi, pt.nshr, sv->sigma, sv->statev, sv->Wm, sv->Lt, false);
    if (commit) {
        pt.stress = stress;
        pt.statev = statev;
        for (int i = 0; i < pt.ntens; i++) pt.stran[i] += dstran[i];
    }
}

void increment_M(const string &name, const vec &props, element_point &pt, const vector<double> &dstran,
                 const double &time, const double &dtime, vector<double> &ddsdde) {
    vector<double> stress;
    increment_M(name, props, pt, dstran, time, dtime, ddsdde, stress, true);
}

// The 3D law driven at the stress-free strain: the reference the condensed route must reproduce.
// The test keeps the full strain (ndi = 3) and solves sigma_C = 0 itself with the 3D DDSDDE.
void increment_reference(const string &name, const vec &props, element_point &ref, const condensed_state &cs,
                         const vec &DE_K, const double &time, const double &dtime, vec &DE_C) {
    vector<double> dstran(6, 0.), ddsdde, stress;
    for (uword i = 0; i < cs.K.n_elem; i++) dstran[cs.K(i)] = DE_K(i);
    for (int it = 0; it < 50; it++) {
        for (uword i = 0; i < cs.C.n_elem; i++) dstran[cs.C(i)] = DE_C(i);
        increment_M(name, props, ref, dstran, time, dtime, ddsdde, stress, false);
        const mat L(ddsdde.data(), 6, 6);
        const vec s(stress.data(), 6);
        if (abs(s(cs.C)).max() < 1.e-11 * abs(s).max()) break;
        DE_C -= solve(L(cs.C, cs.C), s(cs.C));
    }
    for (uword i = 0; i < cs.C.n_elem; i++) dstran[cs.C(i)] = DE_C(i);
    increment_M(name, props, ref, dstran, time, dtime, ddsdde, stress, true);
}

double yield_stress(double p) { return sigmaY_ + k_ * pow(p, m_); }

}  // namespace

TEST(Tcondensed_states, kept_and_free_sets)
{
    const condensed_state truss = make_condensed_state(1, 0);
    EXPECT_TRUE(all(truss.K == uvec({0})));
    EXPECT_TRUE(all(truss.C == uvec({1, 2, 3, 4, 5})));
    const condensed_state beam = make_condensed_state(1, 2);
    EXPECT_TRUE(all(beam.K == uvec({0, 3, 4})));
    EXPECT_TRUE(all(beam.C == uvec({1, 2, 5})));
    const condensed_state ps = make_condensed_state(2);   // Abaqus default: (2, 1)
    EXPECT_TRUE(all(ps.K == uvec({0, 1, 3})));
    EXPECT_TRUE(all(ps.C == uvec({2, 4, 5})));
    const condensed_state pe = make_condensed_state(3, 1);   // plane strain: nothing condensed
    EXPECT_TRUE(all(pe.K == uvec({0, 1, 2, 3})));
    EXPECT_FALSE(pe.active());
    EXPECT_THROW(make_condensed_state(0, 0), std::invalid_argument);
    EXPECT_THROW(make_condensed_state(1, 4), std::invalid_argument);
}

// Truss (1, 0): elastic sigma = E eps exactly, then the uniaxial J2 state — sigma on the hardening
// curve, eps_p11 = p, incompressible transverse plastic strain, DDSDDE the consistent modulus.
TEST(Tcondensed_states, truss_EPICP_uniaxial)
{
    element_point pt(1, 0, 8);
    vector<double> ddsdde;
    increment_M("EPICP", EPICP_props, pt, {1.e-3}, 0., 1., ddsdde);   // elastic: sigma_Y / E = 2.86e-3
    EXPECT_NEAR(pt.stress[0], E_ * 1.e-3, 1.e-9);
    EXPECT_NEAR(ddsdde[0], E_, 1.e-6);
    for (int i = 0; i < 20; i++) increment_M("EPICP", EPICP_props, pt, {1.e-3}, 1. + i, 1., ddsdde);

    const double p = pt.statev[4 + 1];
    EXPECT_GT(p, 1.e-3);
    EXPECT_NEAR(pt.stress[0], yield_stress(p), 1.e-6);
    EXPECT_NEAR(pt.statev[4 + 2], p, 1.e-9);           // EP11
    EXPECT_NEAR(pt.statev[4 + 3], -0.5 * p, 1.e-9);    // EP22
    EXPECT_NEAR(pt.statev[4 + 4], -0.5 * p, 1.e-9);    // EP33
    // consistent uniaxial modulus of the backward-Euler J2 map: E H / (E + H), H = dR/dp
    const double H = m_ * k_ * pow(p, m_ - 1.);
    EXPECT_NEAR(ddsdde[0], E_ * H / (E_ + H), 1.e-6 * E_);
}

// Beam fibre (1, 2): proportional tension-torsion, sqrt(sigma^2 + 3 tau^2) on the hardening
// curve, every condensed stress zero.
TEST(Tcondensed_states, beam_EPICP_tension_torsion)
{
    element_point pt(1, 2, 8);
    vector<double> ddsdde;
    increment_M("EPICP", EPICP_props, pt, {1.e-3, 1.e-3, 0.}, 0., 1., ddsdde);   // elastic
    EXPECT_NEAR(pt.stress[0], E_ * 1.e-3, 1.e-9);
    EXPECT_NEAR(pt.stress[1], G_ * 1.e-3, 1.e-9);
    for (int i = 0; i < 20; i++) increment_M("EPICP", EPICP_props, pt, {1.e-3, 1.e-3, 0.}, 1. + i, 1., ddsdde);
    const double p = pt.statev[4 + 1];
    const double sig = pt.stress[0], tau = pt.stress[1];
    EXPECT_GT(p, 1.e-3);
    EXPECT_NEAR(sqrt(sig * sig + 3. * tau * tau), yield_stress(p), 1.e-6);
    EXPECT_NEAR(pt.stress[2], 0., 1.e-9);   // no gamma_13 imposed
}

// Plane stress, the beam fibre and the truss, on a non-proportional path: the condensed route is
// the 3D law at the stress-free strain (EPICP: plasticity, ZENER: closed-form viscoelasticity),
// Wm and the internal variables included, and DDSDDE the derivative of its stress.
TEST(Tcondensed_states, matches_the_3D_law_and_its_tangent)
{
    const vector<pair<int, int>> states = {{1, 0}, {1, 2}, {2, 1}, {2, 3}};
    const vector<pair<string, vec>> laws = {{"EPICP", EPICP_props}, {"ZENER", ZENER_props}};
    const vec path = {3e-4, -1e-4, 0., 2e-4, -1.5e-4, 1e-4};
    for (const auto &law : laws) {
        const double scale = (law.first == "EPICP") ? 1. : 10.;
        for (const auto &st : states) {
            const condensed_state cs = make_condensed_state(st.first, st.second);
            element_point pt(st.first, st.second, 8);
            element_point ref(3, 3, 8);
            vec DE_C = zeros(cs.C.n_elem);
            const vec DE_K = scale * path(cs.K);
            vector<double> dstran(pt.ntens), ddsdde, stress;
            for (int i = 0; i < pt.ntens; i++) dstran[i] = DE_K(i);
            for (int inc = 0; inc < 12; inc++) {
                const double time = 0.05 * inc;
                if (inc == 11) {   // finite differences of the condensed route at the last increment
                    increment_M(law.first, law.second, pt, dstran, time, 0.05, ddsdde, stress, false);
                    const double h = 1.e-8;
                    for (int j = 0; j < pt.ntens; j++) {
                        vector<double> dp = dstran, dm = dstran, sp, sm, dummy;
                        dp[j] += h;
                        dm[j] -= h;
                        increment_M(law.first, law.second, pt, dp, time, 0.05, dummy, sp, false);
                        increment_M(law.first, law.second, pt, dm, time, 0.05, dummy, sm, false);
                        for (int i = 0; i < pt.ntens; i++) {
                            EXPECT_NEAR(ddsdde[j*pt.ntens + i], (sp[i] - sm[i]) / (2.*h), 1.e-6 * E_)
                                << law.first << " (" << st.first << ", " << st.second << ") DDSDDE(" << i << ", " << j << ")";
                        }
                    }
                }
                increment_M(law.first, law.second, pt, dstran, time, 0.05, ddsdde);
                increment_reference(law.first, law.second, ref, cs, DE_K, time, 0.05, DE_C);
            }
            for (int i = 0; i < pt.ntens; i++) {
                EXPECT_NEAR(pt.stress[i], ref.stress[cs.K(i)], 1.e-8 * (1. + fabs(ref.stress[cs.K(i)])))
                    << law.first << " (" << st.first << ", " << st.second << ") stress " << i;
            }
            for (size_t i = 0; i < pt.statev.size(); i++) {   // Wm and the internal variables
                EXPECT_NEAR(pt.statev[i], ref.statev[i], 1.e-8 * (1. + fabs(ref.statev[i])))
                    << law.first << " (" << st.first << ", " << st.second << ") statev " << i;
            }
        }
    }
}

// Thermomechanical plane stress (ELISO, thermoelastic): with sigma_33 = 0 the stress and its
// temperature derivative follow the reduced (plane-stress) moduli, and DDSDDE, DDSDDT, DRPLDE,
// DRPLDT are the derivatives of the condensed route.
TEST(Tcondensed_states, thermoelastic_plane_stress)
{
    const double rho = 7.8e-9, c_p = 4.6e8, alpha = 1.2e-5, T0 = 293.15;
    const vec props = {rho, c_p, E_, nu_, alpha};
    struct out_T { vector<double> stress, ddsdde, ddsddt, drplde; double drpldt = 0., rpl = 0.; };
    auto run = [&](const vector<double> &dstran, const double &DT) {
        const int ndi = 2, nshr = 1, ntens = 3, nstatev = 1 + 7;
        vector<double> statev(nstatev, 0.), stran(ntens, 0.);
        out_T o;
        o.stress.assign(ntens, 0.);
        o.ddsdde.assign(ntens*ntens, 0.);
        o.ddsddt.assign(ntens, 0.);
        o.drplde.assign(ntens, 0.);
        const double time_arr[2] = {0., 0.};
        const double drot[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
        phase_characteristics rve;
        rve.construct(0, 2);
        auto sv = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
        sv->resize(nstatev - 7);
        rve.sptr_matprops->resize(props.n_elem);
        vec props_smart = zeros(props.n_elem);
        double Time = 0., DTime = 0.;
        mat DR = eye(3, 3);
        bool start = false;
        abaqus2smart_T(o.stress.data(), o.ddsdde.data(), o.ddsddt.data(), o.drplde.data(), o.drpldt, stran.data(),
                       dstran.data(), time_arr, 1., T0, DT, props.n_elem, props.memptr(), nstatev, statev.data(), ndi,
                       nshr, drot, sv->sigma, sv->dSdE, sv->dSdT, sv->drdE, sv->drdT, sv->Etot, sv->DEtot, sv->T,
                       sv->DT, Time, DTime, props_smart, sv->Wm, sv->Wt, sv->statev, DR, start);
        start = true;   // the first increment of the point: T_init = T
        rve.sptr_matprops->update(0, "ELISO", 1., 0., 0., 0., props.n_elem, props_smart);
        double pnewdt = 1.;
        select_umat_T(rve, DR, Time, DTime, ndi, nshr, start, 0, pnewdt);
        EXPECT_GE(pnewdt, 1.);
        smart2abaqus_T(o.stress.data(), o.ddsdde.data(), o.ddsddt.data(), o.drplde.data(), o.drpldt, o.rpl,
                       statev.data(), ndi, nshr, sv->sigma, sv->statev, sv->r, sv->Wm, sv->Wt, sv->dSdE, sv->dSdT,
                       sv->drdE, sv->drdT, false);
        return o;
    };
    const vector<double> de = {1.e-3, -2.e-4, 5.e-4};
    const double DT = 2.;
    const out_T o = run(de, DT);

    // plane-stress thermoelasticity: sigma = Q (eps - alpha DT), Q = E/(1-nu^2) [[1, nu], [nu, 1]], G
    const double Q11 = E_ / (1. - nu_*nu_), Q12 = nu_ * Q11;
    EXPECT_NEAR(o.stress[0], Q11*(de[0] - alpha*DT) + Q12*(de[1] - alpha*DT), 1.e-8);
    EXPECT_NEAR(o.stress[1], Q12*(de[0] - alpha*DT) + Q11*(de[1] - alpha*DT), 1.e-8);
    EXPECT_NEAR(o.stress[2], G_*de[2], 1.e-8);
    EXPECT_NEAR(o.ddsdde[0], Q11, 1.e-6);
    EXPECT_NEAR(o.ddsdde[1], Q12, 1.e-6);
    EXPECT_NEAR(o.ddsdde[8], G_, 1.e-6);
    EXPECT_NEAR(o.ddsddt[0], -(Q11 + Q12)*alpha, 1.e-9);
    EXPECT_NEAR(o.ddsddt[2], 0., 1.e-12);   // the shear entry (12), not sigma_33

    // DDSDDE, DRPLDE and DDSDDT, DRPLDT (rates at dtime = 1) against finite differences
    const double h = 1.e-8, hT = 1.e-5;
    for (int j = 0; j < 3; j++) {
        vector<double> dp = de, dm = de;
        dp[j] += h;
        dm[j] -= h;
        const out_T op = run(dp, DT), om = run(dm, DT);
        for (int i = 0; i < 3; i++) EXPECT_NEAR(o.ddsdde[j*3 + i], (op.stress[i] - om.stress[i]) / (2.*h), 1.e-5 * E_);
        EXPECT_NEAR(o.drplde[j], (op.rpl - om.rpl) / (2.*h), 1.e-6 * fabs(o.drplde[0]) + 1.e-9);
    }
    const out_T op = run(de, DT + hT), om = run(de, DT - hT);
    for (int i = 0; i < 3; i++) EXPECT_NEAR(o.ddsddt[i], (op.stress[i] - om.stress[i]) / (2.*hT), 1.e-6);
    EXPECT_NEAR(o.drpldt, (op.rpl - om.rpl) / (2.*hT), 1.e-6 * fabs(o.drpldt));
}

// The bridge maps the element's ntens components in Abaqus order whatever (ndi, nshr): the beam
// fibre (1, 2) carries 11, 12, 13 and DDSDDE is the K x K block of the (condensed) tangent.
TEST(Tcondensed_states, abaqus_component_mapping)
{
    mat Lt = zeros(6, 6);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) Lt(i, j) = 10.*i + j;   // non-symmetric: catches a transpose
    const vec sigma = {1., 2., 3., 4., 5., 6.};
    const vec statev_smart = zeros(1), Wm = zeros(4);
    double stress[6], ddsdde[36], statev[5];
    smart2abaqus_M(stress, ddsdde, statev, 1, 2, sigma, statev_smart, Wm, Lt, false);
    const uvec K = {0, 3, 4};
    for (int i = 0; i < 3; i++) {
        EXPECT_DOUBLE_EQ(stress[i], sigma(K(i)));
        for (int j = 0; j < 3; j++) EXPECT_DOUBLE_EQ(ddsdde[j*3 + i], Lt(K(i), K(j)));
    }
    // plane strain (3, 1): 11, 22, 33, 12, column-major
    smart2abaqus_M(stress, ddsdde, statev, 3, 1, sigma, statev_smart, Wm, Lt, false);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) EXPECT_DOUBLE_EQ(ddsdde[j*4 + i], Lt(i, j));
}
