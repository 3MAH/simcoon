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

///@file Tabaqus_finite.cpp
///@brief The Abaqus bridges under NLGEOM: the law on DFGRD0/DFGRD1, Cauchy stress, Jaumann
///       Jacobian and the heat per unit current volume (fea_transfer.hpp)
///@version 1.0

#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// One Abaqus material point, 3D, carried across increments as Abaqus carries it: STRESS rotated
// by DROT before each call, STATEV as returned
struct point {
    string cmname;
    vec props;
    vec stress = zeros(6);
    vector<double> statev;
    mat ddsdde = zeros(6, 6);
    vec ddsddt = zeros(6), drplde = zeros(6);
    double rpl = 0., drpldt = 0., pnewdt = 1.;
    double time = 0.;
    double temperature;
    bool thermo;

    point(const string &name, const vec &p, const int &nstatev_law, const bool &thermo_, const double &T0 = 293.15)
        : cmname(name), props(p), temperature(T0), thermo(thermo_) {
        statev.assign((thermo ? 7 : 4) + nstatev_law + abaqus_kinematic_statev(abaqus_kinematic_form(cmname.c_str())), 0.);
    }

    // the body of umat_singleM / umat_singleT on (DFGRD0, DFGRD1, DROT, DTEMP)
    void call(const mat &F0, const mat &F1, const mat &DROT, const double &dtime, const double &dtemp = 0.) {
        const int ndi = 3, nshr = 3, nstatev = static_cast<int>(statev.size()), nprops = static_cast<int>(props.n_elem);
        const vec stress_rot = rotate_stress(stress, DROT);
        double stress_a[6], ddsdde_a[36], ddsddt_a[6], drplde_a[6];
        std::copy(stress_rot.begin(), stress_rot.end(), stress_a);
        std::fill(ddsdde_a, ddsdde_a + 36, 0.);
        std::fill(ddsddt_a, ddsddt_a + 6, 0.);
        std::fill(drplde_a, drplde_a + 6, 0.);
        const double stran[6] = {0., 0., 0., 0., 0., 0.}, dstran[6] = {0., 0., 0., 0., 0., 0.};
        const double time_a[2] = {time, time};
        const int corate_type = abaqus_kinematic_form(cmname.c_str());
        const int nstatev_kin = abaqus_kinematic_statev(corate_type);
        double *kin_statev = statev.data() + nstatev - nstatev_kin;
        const string umat_name = cmname.substr(0, 5);
        const vec props_v = props;
        bool start = false;
        double Time = 0., DTime = 0.;
        mat DR = zeros(3, 3);
        vec props_smart = zeros(nprops);
        phase_characteristics rve;
        if (!thermo) {
            rve.construct(0, 1);
            auto sv = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
            sv->resize(nstatev - 4 - nstatev_kin);
            rve.sptr_matprops->resize(nprops);
            abaqus2smart_M(stress_a, ddsdde_a, stran, dstran, time_a, dtime, temperature, dtemp, nprops, props_v.memptr(), nstatev - nstatev_kin, statev.data(), ndi, nshr, DROT.memptr(), sv->sigma, sv->Lt, sv->Etot, sv->DEtot, sv->T, sv->DT, Time, DTime, props_smart, sv->Wm, sv->statev, DR, start);
            rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
            ASSERT_TRUE(abaqus_finite_route(umat_name, ndi, corate_type, false, F1.memptr()));
            abaqus2smart_finite(*sv, DROT.memptr(), F0.memptr(), F1.memptr(), kin_statev, corate_type, DTime);
            run_umat_M(rve, sv->DR, Time, DTime, ndi, nshr, start, 0, 5, corate_type, pnewdt);
            vec s;
            mat C;
            smart2abaqus_finite_M(*sv, sv->Lt, corate_type, s, C, kin_statev);
            smart2abaqus_M(stress_a, ddsdde_a, statev.data(), ndi, nshr, s, sv->statev, sv->Wm, C, false);
        }
        else {
            rve.construct(0, 2);
            auto sv = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
            sv->resize(nstatev - 7 - nstatev_kin);
            rve.sptr_matprops->resize(nprops);
            const double drpldt_in = 0.;
            abaqus2smart_T(stress_a, ddsdde_a, ddsddt_a, drplde_a, drpldt_in, stran, dstran, time_a, dtime, temperature, dtemp, nprops, props_v.memptr(), nstatev - nstatev_kin, statev.data(), ndi, nshr, DROT.memptr(), sv->sigma, sv->dSdE, sv->dSdT, sv->drdE, sv->drdT, sv->Etot, sv->DEtot, sv->T, sv->DT, Time, DTime, props_smart, sv->Wm, sv->Wt, sv->statev, DR, start);
            rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
            ASSERT_TRUE(abaqus_finite_route(umat_name, ndi, corate_type, true, F1.memptr()));
            abaqus2smart_finite(*sv, DROT.memptr(), F0.memptr(), F1.memptr(), kin_statev, corate_type, DTime);
            run_umat_T(rve, sv->DR, Time, DTime, ndi, nshr, start, 0, 5, corate_type, pnewdt);
            vec s, dsdt, drde;
            mat C;
            double q = 0., dqdt = 0.;
            smart2abaqus_finite_T(*sv, corate_type, DTime, s, C, dsdt, q, drde, dqdt, kin_statev);
            smart2abaqus_T(stress_a, ddsdde_a, ddsddt_a, drplde_a, drpldt, rpl, statev.data(), ndi, nshr, s, sv->statev, q, sv->Wm, sv->Wt, C, dsdt, drde, mat({dqdt}), false);
            ddsddt = vec(ddsddt_a, 6);
            drplde = vec(drplde_a, 6);
        }
        stress = vec(stress_a, 6);
        ddsdde = mat(ddsdde_a, 6, 6);
        time += dtime;
        temperature += dtemp;
    }
};

const double E_ = 210000., nu_ = 0.3;
const vec props_EPICP_M = {E_, nu_, 1.2e-5, 300., 500., 0.4};
const vec props_EPICP_T = {7.8e-9, 4.6e8, E_, nu_, 1.2e-5, 300., 500., 0.4};

// Stretch with shear and a finite rotation
mat F_test(const double &s) {
    const mat U = {{1. + 0.12*s, 0.05*s, 0.}, {0.05*s, 1. - 0.03*s, 0.02*s}, {0., 0.02*s, 1. + 0.04*s}};
    return Rotation::from_axis_angle(0.4*s, 3).as_matrix()*U;
}

// Abaqus's increment of the strain: F1' = exp(h e) F1, e the Voigt unit (engineering shear)
mat perturbed(const mat &F1, const int &j, const double &h) {
    vec e = zeros(6);
    e(j) = h;
    return expmat_sym(v2t_strain(e))*F1;
}

}  // namespace

TEST(Tabaqus_finite, kinematic_form_and_route)
{
    EXPECT_EQ(abaqus_kinematic_form("EPICP"), corate_type_total);
    EXPECT_EQ(abaqus_kinematic_form("EPICP                         "), corate_type_total);
    EXPECT_EQ(abaqus_kinematic_form("EPICP_R                       "), 3);
    EXPECT_EQ(abaqus_kinematic_statev(corate_type_total), 0);
    EXPECT_EQ(abaqus_kinematic_statev(3), 15);

    const mat F = F_test(1.), I = eye(3, 3);
    EXPECT_TRUE(abaqus_finite_route("EPICP", 3, corate_type_total, false, F.memptr()));
    EXPECT_FALSE(abaqus_finite_route("EPICP", 3, corate_type_total, false, I.memptr()));   // NLGEOM = NO
    EXPECT_FALSE(abaqus_finite_route("EPICP", 2, corate_type_total, false, F.memptr()));   // plane stress
    EXPECT_FALSE(abaqus_finite_route("SMADI", 3, corate_type_total, false, F.memptr()));   // small-strain-only mechanical law
    EXPECT_FALSE(abaqus_finite_route("HYPOO", 3, corate_type_total, false, F.memptr()));   // not a total-form law...
    EXPECT_TRUE(abaqus_finite_route("HYPOO", 3, 3, false, F.memptr()));                    // ...but a rate one
    EXPECT_FALSE(abaqus_finite_route("UMABA", 3, 3, false, F.memptr()));
    EXPECT_TRUE(abaqus_finite_route("SMADI", 3, corate_type_total, true, F.memptr()));      // every thermomechanical law
}

TEST(Tabaqus_finite, total_ELISO_is_hencky)
{
    // Isotropic Hencky: T = L:ln U, tau = R T R^T, STRESS = tau/J, in one increment from rest
    point pt("ELISO", {E_, nu_, 1.2e-5}, 1, false);
    const mat F1 = F_test(1.);
    pt.call(eye(3, 3), F1, eye(3, 3), 1.);
    mat R, U;
    RU_decomposition(R, U, F1);
    const vec T = L_iso(E_, nu_, "Enu")*t2v_strain(logmat_sympd(U));
    const vec sigma = t2v_stress(R*v2t_stress(T)*R.t())/det(F1);
    EXPECT_LT(norm(pt.stress - sigma), 1.e-9*norm(sigma));
}

TEST(Tabaqus_finite, objectivity_is_exact_in_total_form)
{
    // Superposed rigid rotation Q (constant: DROT = I for both): STRESS' = Q STRESS Q^T, at every
    // increment of a plastic path, the Jacobian rotated alike
    const mat Q = Rotation::from_axis_angle(1.1, 1).as_matrix()*Rotation::from_axis_angle(0.7, 3).as_matrix();
    point a("EPICP", props_EPICP_M, 8, false), b("EPICP", props_EPICP_M, 8, false);
    const int n = 4;
    for (int k = 0; k < n; k++) {
        const mat F0 = F_test(double(k)/n), F1 = F_test(double(k + 1)/n);
        a.call(F0, F1, eye(3, 3), 0.25);
        b.call(Q*F0, Q*F1, eye(3, 3), 0.25);
        EXPECT_LT(norm(b.stress - rotate_stress(a.stress, Q)), 1.e-10*norm(a.stress));
        EXPECT_LT(norm(b.ddsdde - rotate_stiffness(a.ddsdde, Q)), 1.e-10*norm(a.ddsdde));
    }
    EXPECT_GT(a.statev[5], 1.e-3);   // p: the path is plastic
}

TEST(Tabaqus_finite, rate_form_keeps_its_strain_and_frame)
{
    // _R under a coaxial stretch: the stored strain is ln V (to the log_R accuracy, A^R:D
    // integrated, ~1e-4) and the frame I; ELISO matches the total form to the same accuracy
    point r("ELISO_R", {E_, nu_, 1.2e-5}, 1, false), t("ELISO", {E_, nu_, 1.2e-5}, 1, false);
    ASSERT_EQ(r.statev.size(), 4u + 1u + 15u);
    auto F = [](const double &s) { return mat(diagmat(vec({1. + 0.3*s, 1. - 0.1*s, 1. + 0.05*s}))); };
    for (int k = 0; k < 5; k++) {
        r.call(F(0.2*k), F(0.2*(k + 1)), eye(3, 3), 0.2);
        t.call(F(0.2*k), F(0.2*(k + 1)), eye(3, 3), 0.2);
    }
    const vec etot(r.statev.data() + 5, 6);
    const mat R(r.statev.data() + 11, 3, 3);
    const vec lnV = t2v_strain(logmat_sympd(F(1.)));
    EXPECT_LT(norm(etot - lnV), 5.e-4*norm(lnV));
    EXPECT_LT(norm(R - eye(3, 3)), 1.e-12);
    EXPECT_LT(norm(r.stress - t.stress), 5.e-4*norm(t.stress));
}

TEST(Tabaqus_finite, ddsdde_is_the_jaumann_jacobian)
{
    // DDSDDE = (1/J) d(J sigma)/d eps for an increment of Abaqus's strain without spin, by central
    // differences on DFGRD1, on a plastic increment: exact in the total form; the rate form's box
    // tangent ignores the variation of the log_R frame increment (~3e-3)
    for (const string name : {"EPICP", "EPICP_R"}) {
        point base(name, props_EPICP_M, 8, false);
        base.call(eye(3, 3), F_test(0.5), eye(3, 3), 0.5);
        const mat F0 = F_test(0.5), F1 = F_test(0.6);
        point ref = base;
        ref.call(F0, F1, eye(3, 3), 0.1);
        const double h = 1.e-6;
        mat C_fd(6, 6);
        for (int j = 0; j < 6; j++) {
            point p = base, m = base;
            const mat Fp = perturbed(F1, j, h), Fm = perturbed(F1, j, -h);
            p.call(F0, Fp, eye(3, 3), 0.1);
            m.call(F0, Fm, eye(3, 3), 0.1);
            C_fd.col(j) = (det(Fp)*p.stress - det(Fm)*m.stress)/(2.*h*det(F1));
        }
        EXPECT_LT(norm(ref.ddsdde - C_fd), (name == "EPICP" ? 1.e-6 : 1.e-2)*norm(C_fd)) << name;
    }
}

TEST(Tabaqus_finite, thermal_outputs_are_the_derivatives_of_rpl_and_stress)
{
    // RPL = r/(J dt): DRPLDE, DRPLDT and DDSDDT by central differences on DFGRD1 and DTEMP, on a
    // plastic increment (dissipation and thermoelastic coupling), both forms. DRPLDE carries the
    // accuracy of the kernel's own dr/deps (6e-4 here for EPICP; the conversion alone is exact,
    // drplde_is_exact_for_thermoelasticity)
    for (const string name : {"EPICP", "EPICP_R"}) {
        point base(name, props_EPICP_T, 8, true);
        base.call(eye(3, 3), F_test(0.5), eye(3, 3), 0.5, 2.);
        const mat F0 = F_test(0.5), F1 = F_test(0.6);
        const double dt = 0.1, dT = 1.5;
        point ref = base;
        ref.call(F0, F1, eye(3, 3), dt, dT);
        EXPECT_GT(fabs(ref.rpl), 0.) << name;

        const double h = 1.e-6;
        vec drplde_fd(6);
        for (int j = 0; j < 6; j++) {
            point p = base, m = base;
            p.call(F0, perturbed(F1, j, h), eye(3, 3), dt, dT);
            m.call(F0, perturbed(F1, j, -h), eye(3, 3), dt, dT);
            drplde_fd(j) = (p.rpl - m.rpl)/(2.*h);
        }
        EXPECT_LT(norm(ref.drplde - drplde_fd), 2.e-3*norm(drplde_fd)) << name;

        const double hT = 1.e-4;
        point p = base, m = base;
        p.call(F0, F1, eye(3, 3), dt, dT + hT);
        m.call(F0, F1, eye(3, 3), dt, dT - hT);
        EXPECT_NEAR(ref.drpldt, (p.rpl - m.rpl)/(2.*hT), 1.e-4*fabs(ref.drpldt)) << name;
        const vec ddsddt_fd = (p.stress - m.stress)/(2.*hT);
        EXPECT_LT(norm(ref.ddsddt - ddsddt_fd), 1.e-4*norm(ddsddt_fd)) << name;
    }
}


TEST(Tabaqus_finite, drplde_is_exact_for_thermoelasticity)
{
    // ELISO's dr/deps is exact: DRPLDE matches the central differences to rounding in the total
    // form, to the log_R frame-increment accuracy in the rate form
    const vec props_ELISO_T = {7.8e-9, 4.6e8, E_, nu_, 1.2e-5};
    for (const string name : {"ELISO", "ELISO_R"}) {
        point base(name, props_ELISO_T, 1, true);
        base.call(eye(3, 3), F_test(0.5), eye(3, 3), 0.5, 2.);
        const mat F0 = F_test(0.5), F1 = F_test(0.6);
        point ref = base;
        ref.call(F0, F1, eye(3, 3), 0.1, 1.5);
        const double h = 1.e-6;
        vec drplde_fd(6);
        for (int j = 0; j < 6; j++) {
            point p = base, m = base;
            p.call(F0, perturbed(F1, j, h), eye(3, 3), 0.1, 1.5);
            m.call(F0, perturbed(F1, j, -h), eye(3, 3), 0.1, 1.5);
            drplde_fd(j) = (p.rpl - m.rpl)/(2.*h);
        }
        EXPECT_LT(norm(ref.drplde - drplde_fd), (name == "ELISO" ? 1.e-8 : 1.e-3)*norm(drplde_fd)) << name;
    }
}
