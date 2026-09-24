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
///@brief Every finite kernel must DECLARE what its raw outputs are expressed in.
///@version 1.0

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <stdexcept>

#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>

using namespace simcoon;

// The one thing this file exists to prevent: a kernel added to the finite dispatch without
// saying what its stress and tangent mean. Before the conventions table those two facts lived
// apart -- a public predicate and a function-local set 190 lines away -- so a new kernel had to
// be remembered twice, and being forgotten in either was SILENT: a missing stress measure is an
// error of exactly J, a missing tangent rate is a wrong rate. Both are invisible at J ~ 1.
TEST(UmatConventions, EveryFiniteKernelDeclaresItsOutputConventions) {
    ASSERT_FALSE(finite_umat_names().empty()) << "the finite dispatch map must not be empty";
    for (const auto &entry : finite_umat_names()) {
        EXPECT_NO_THROW({ output_convention_of(entry.first); })
            << "the umat '" << entry.first << "' is served by the finite dispatch but has not "
            << "declared its output conventions (see output_convention_of in umat_smart.cpp)";
    }
}

TEST(UmatConventions, AnUndeclaredKernelThrowsRatherThanDefaulting) {
    EXPECT_THROW(output_convention_of("NOPE1"), std::invalid_argument);
    // ...and the message must point at the fix, not just say "not found".
    try {
        output_convention_of("NOPE1");
        FAIL() << "expected a throw";
    } catch (const std::invalid_argument &e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("umat_smart.cpp"), std::string::npos) << msg;
    }
}

// The tangent RATE is deliberately not a declared property: every kernel is handed the
// solver's corate_type and must answer in it. It was one briefly, while the finite kernels
// baked the log box and the dispatcher re-expressed it. If a `tangent` field ever comes back,
// it means some kernel has again been allowed to pick its own rate -- which is what made the
// old two-set arrangement able to disagree with itself.
TEST(UmatConventions, EveryNativeKernelDeclaresOnlyItsStressMeasure) {
    EXPECT_EQ(output_convention_of("MODUL").stress, StressMeasure::kirchhoff);
    EXPECT_EQ(output_convention_of("HOLZA").stress, StressMeasure::kirchhoff);
    EXPECT_EQ(output_convention_of("NEOHC").stress, StressMeasure::kirchhoff);
    EXPECT_EQ(output_convention_of("SNTVE").stress, StressMeasure::kirchhoff);
    EXPECT_EQ(output_convention_of("EPICP").stress, StressMeasure::kirchhoff);
}

// Every NATIVE kernel is Kirchhoff since the Kirchhoff-native refactor. The only Cauchy
// declarations left are conventions simcoon does not own, and they are declared EXPLICITLY --
// which is the whole point: HYPOO used to sit in neither set and inherit "Cauchy + in-rate" by
// omission, with a comment reading "kept in the Cauchy group for now".
TEST(UmatConventions, OnlyForeignConventionsAreCauchyAndTheyAreDeclared) {
    for (const auto &entry : finite_umat_names()) {
        const umat_convention conv = output_convention_of(entry.first);
        if (conv.stress == StressMeasure::cauchy) {
            const bool foreign = (entry.first == "HYPOO"     // Cauchy-rate hypoelastic law
                               || entry.first == "UMEXT"     // external dylib plugin
                               || entry.first == "UMABA");   // Abaqus wrapper
            EXPECT_TRUE(foreign)
                << "'" << entry.first << "' declares Cauchy, but every NATIVE kernel should be "
                << "Kirchhoff-native. If this is deliberate, add it to the foreign list here "
                << "with the reason; if not, the kernel needs converting.";
        }
    }
    EXPECT_EQ(output_convention_of("HYPOO").stress, StressMeasure::cauchy);
}

// stress_output_is_kirchhoff is consumed by the python wrapper for EVERY name it serves,
// including small-strain-only ones the finite dispatch never sees. It must stay total where
// output_convention_of throws, and it must not pay a thrown exception to do so.
TEST(UmatConventions, TheKirchhoffPredicateIsTotal) {
    EXPECT_TRUE(stress_output_is_kirchhoff("MODUL"));
    EXPECT_TRUE(stress_output_is_kirchhoff("HOLZA"));
    EXPECT_FALSE(stress_output_is_kirchhoff("HYPOO"));
    EXPECT_FALSE(stress_output_is_kirchhoff("SMADI"));   // small-strain only: no finite convention
    EXPECT_FALSE(stress_output_is_kirchhoff("NOPE1"));   // unknown: leave the stress alone
}
