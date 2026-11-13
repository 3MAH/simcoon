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

///@file Tumat_conventions.cpp
///@brief The registry of the law names: one row per name, consistent across the faces of the
///       dispatch, and the finite conventions every finite kernel declares through it.
///@version 1.0

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <stdexcept>

#include <simcoon/Continuum_mechanics/Umat/umat_registry.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>

using namespace simcoon;
using K = UmatKernel;

namespace {

bool is_plugin(const K &k) {
    return k == K::plugin_ext || k == K::plugin_abaqus || k == K::plugin_ansys;
}

}  // namespace

// A name the registry holds is served by at least one face; the mechanical faces agree on the
// kernel (sim.umat feeds the one kernel the small strain or the box strain); a kernel served
// through the modular adapter has its props translator, and every translator has its name.
TEST(UmatRegistry, EveryRowServesAFaceAndTheMechanicalFacesAgree) {
    ASSERT_FALSE(umat_registry().empty());
    for (const auto &[name, entry] : umat_registry()) {
        EXPECT_TRUE(entry.mechanical != K::none || entry.finite != K::none || entry.thermomechanical != K::none)
            << name << " is served by no face";
        if (entry.mechanical != K::none && entry.finite != K::none)
            EXPECT_EQ(entry.mechanical, entry.finite) << name << ": the mechanical faces call different kernels";
        const bool adapter = entry.mechanical == K::legacy_modular || entry.finite == K::legacy_modular
                             || entry.thermomechanical == K::legacy_modular;
        EXPECT_EQ(adapter, has_legacy_adapter(name)) << name << ": the registry and legacy_adapters.cpp disagree";
        EXPECT_EQ(umat_entry_of(name), &entry);
    }
    EXPECT_EQ(umat_entry_of("NOPE1"), nullptr);
}

// A _closed name keeps the dedicated kernel of a name the modular engine serves: the plain name
// exists and runs the adapter, the _closed one never does.
TEST(UmatRegistry, ClosedNamesPairWithTheirAdapterName) {
    int closed = 0;
    for (const auto &[name, entry] : umat_registry()) {
        if (name.size() <= 7 || name.compare(name.size() - 7, 7, "_closed") != 0) continue;
        ++closed;
        const umat_entry *plain = umat_entry_of(name.substr(0, name.size() - 7));
        ASSERT_NE(plain, nullptr) << name;
        for (const K k : {entry.mechanical, entry.finite, entry.thermomechanical})
            EXPECT_NE(k, K::legacy_modular) << name;
        EXPECT_TRUE(plain->mechanical == K::legacy_modular || plain->finite == K::legacy_modular
                    || plain->thermomechanical == K::legacy_modular) << name;
        if (entry.finite != K::none && entry.convention.material_frame)   // two boxes: the same statev layout (EPICP_closed = EPICP)
            EXPECT_EQ(entry.convention.material_frame, plain->convention.material_frame) << name;
    }
    EXPECT_EQ(closed, 7);   // EPICP, EPCHA, PRONK, ZENER, EPKCP (thermomechanical), NEOHC and NEOHI (built on F)
}

// The hyperelastic names run the HYPER_INVARIANTS block on the log-strain box: adapter on the finite
// face only, a declared (empty) statev layout, served under corate 3 and the total form (umat_smart).
TEST(UmatRegistry, HyperelasticNamesRunTheModularBlock) {
    for (const char *name : {"NEOHC", "NEOHI", "MOORI", "YEOHH", "ISHAH", "GETHH", "SWANH", "HOLZA"}) {
        EXPECT_TRUE(adapter_is_hyperelastic(name)) << name;
        const umat_entry *e = umat_entry_of(name);
        ASSERT_NE(e, nullptr) << name;
        EXPECT_EQ(e->finite, K::legacy_modular) << name;
        EXPECT_EQ(e->mechanical, K::none) << name;
        EXPECT_TRUE(e->convention.material_frame && e->convention.layout_declared) << name;
    }
    for (const char *name : {"NEOHC_closed", "NEOHI_closed", "OGDEN", "SNTVE"}) {
        EXPECT_FALSE(adapter_is_hyperelastic(name)) << name;
        EXPECT_TRUE(output_convention_of(name).material_tangent) << name;   // built on F
    }
    EXPECT_FALSE(adapter_is_hyperelastic("EPICP"));
}

// The one thing the convention table exists to prevent: a kernel on the finite route whose stress
// and tangent meanings are not stated. Every finite kernel is Kirchhoff-native (HYPOO included: it
// integrates the Kirchhoff rate), a box kernel runs in the material frame, a kernel built on F
// has the material tangent, and the external plugins stay off the finite route.
TEST(UmatRegistry, EveryFiniteKernelDeclaresKirchhoffConventions) {
    int finite = 0;
    for (const auto &[name, entry] : umat_registry()) {
        if (is_plugin(entry.mechanical)) {
            EXPECT_EQ(entry.finite, K::none) << name << " must stay off the finite route";
            EXPECT_FALSE(stress_output_is_kirchhoff(name)) << name;
            continue;
        }
        if (entry.finite == K::none) {
            EXPECT_THROW(output_convention_of(name), std::invalid_argument) << name;
            EXPECT_FALSE(stress_output_is_kirchhoff(name)) << name;
            continue;
        }
        ++finite;
        const umat_convention conv = output_convention_of(name);
        EXPECT_EQ(conv.stress, StressMeasure::kirchhoff) << name << " does not declare a Kirchhoff stress";
        EXPECT_TRUE(stress_output_is_kirchhoff(name)) << name;
        EXPECT_NE(conv.material_frame, conv.material_tangent) << name << ": a box kernel or a kernel built on F";
        if (!conv.material_frame) EXPECT_TRUE(conv.statev_tensors.empty()) << name;
        if (!conv.layout_declared) EXPECT_TRUE(conv.statev_tensors.empty()) << name;
        EXPECT_EQ(material_tangent_serves(name), conv.material_tangent) << name;
    }
    EXPECT_GE(finite, 20);
}

TEST(UmatRegistry, AnUnknownNameThrowsRatherThanDefaulting) {
    EXPECT_THROW(output_convention_of("NOPE1"), std::invalid_argument);
    EXPECT_THROW(output_convention_of("SMADI"), std::invalid_argument);   // small-strain only
    // ...and the message must point at the fix, not just say "not found".
    try {
        output_convention_of("NOPE1");
        FAIL() << "expected a throw";
    } catch (const std::invalid_argument &e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("umat_registry.cpp"), std::string::npos) << msg;
    }
    EXPECT_FALSE(stress_output_is_kirchhoff("NOPE1"));   // unknown: leave the stress alone
    EXPECT_FALSE(total_formulation_serves("NOPE1"));
    EXPECT_FALSE(total_formulation_serves_T("NOPE1"));
    EXPECT_FALSE(is_external_plugin("NOPE1"));
}

// The tangent RATE is deliberately not a declared property: every kernel is handed the
// solver's corate_type and must answer in it. If a `tangent` field ever comes back, it means
// some kernel has again been allowed to pick its own rate.
TEST(UmatRegistry, TheNativeKernelsDeclareOnlyTheirStressMeasure) {
    for (const char *name : {"MODUL", "HOLZA", "NEOHC", "SNTVE", "EPICP", "HYPOO"})
        EXPECT_EQ(output_convention_of(name).stress, StressMeasure::kirchhoff) << name;
    EXPECT_TRUE(is_external_plugin("UMEXT"));
    EXPECT_TRUE(is_external_plugin("UMABA"));
    EXPECT_TRUE(is_external_plugin("UMANS"));
    EXPECT_FALSE(is_external_plugin("EPICP"));
}

// The capabilities read off the kernels: the complementarity option and the closest point on
// each face, the mean-field scheme codes and the sub-phase layout they imply.
TEST(UmatRegistry, CapabilitiesFollowTheKernels) {
    EXPECT_NO_THROW(check_ncp_serves("MODUL", 1, false));
    EXPECT_NO_THROW(check_ncp_serves("EPICP", 1, false));
    EXPECT_NO_THROW(check_ncp_serves("EPICP_closed", 1, false));
    EXPECT_NO_THROW(check_ncp_serves("EPKCP", 1, true));
    EXPECT_THROW(check_ncp_serves("EPICP_closed", 1, true), std::invalid_argument);
    EXPECT_THROW(check_ncp_serves("SMADI", 1, false), std::invalid_argument);
    EXPECT_THROW(check_ncp_serves("NOPE1", 1, false), std::invalid_argument);
    EXPECT_NO_THROW(check_integrator_serves("EPCHA_closed", tangent_closest_point, false));
    EXPECT_NO_THROW(check_integrator_serves("EPJCK", tangent_closest_point, false));
    EXPECT_NO_THROW(check_integrator_serves("SMRDI", tangent_closest_point, false));
    EXPECT_NO_THROW(check_integrator_serves("EPKCP_closed", tangent_closest_point, true));
    EXPECT_THROW(check_integrator_serves("EPJCK", tangent_closest_point, true), std::invalid_argument);
    EXPECT_THROW(check_integrator_serves("ZENER_closed", tangent_closest_point, false), std::invalid_argument);   // the Kelvin kernel; ZENER runs the modular engine
    EXPECT_NO_THROW(check_integrator_serves("ZENER", tangent_closest_point, false));
    EXPECT_NO_THROW(check_integrator_serves("HOLZA", tangent_closest_point, false));   // the modular block (no mechanism)
    EXPECT_THROW(check_integrator_serves("OGDEN", tangent_closest_point, false), std::invalid_argument);

    const std::map<std::string, int> schemes = {{"MIHEN", 100}, {"MIMTN", 101}, {"MISCN", 103}, {"MIPLN", 104}};
    for (const auto &[name, entry] : umat_registry()) {
        const auto it = schemes.find(name);
        const int scheme = (it != schemes.end()) ? it->second : 0;
        EXPECT_EQ(mean_field_scheme(entry.mechanical), scheme) << name;
        EXPECT_EQ(sub_phase_shape(name), scheme == 104 ? 1 : (scheme > 0 ? 2 : 0)) << name;
        if (scheme > 0) EXPECT_EQ(entry.elastic_stiffness, entry.mechanical) << name;
    }
    EXPECT_EQ(sub_phase_shape("NOPE1"), 0);
}

// Every thermomechanical kernel is a log-strain box: the thermomechanical face is the admission
// to the total formulation.
TEST(UmatRegistry, TheThermomechanicalFaceServesTheTotalForm) {
    for (const auto &[name, entry] : umat_registry())
        EXPECT_EQ(total_formulation_serves_T(name), entry.thermomechanical != K::none) << name;
}
