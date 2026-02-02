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

///@file Tdamage.cpp
///@brief Test for damage evolution functions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/damage.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tdamage, damage_weibull_vonmises)
{
    double alpha = 2.0;
    double beta = 200.;
    double damage = 0.;
    double DTime = 1.0;

    // Zero stress -> zero damage rate
    vec sigma_zero = zeros(6);
    double dD = damage_weibull(sigma_zero, damage, alpha, beta, DTime, "Mises");
    EXPECT_LT(fabs(dD), 1.E-9);

    // Non-zero uniaxial stress
    vec sigma = {400., 0., 0., 0., 0., 0.};
    dD = damage_weibull(sigma, damage, alpha, beta, DTime, "Mises");
    // d = (1-D)/DTime * (1 - exp(-(sigma_eq/beta)^alpha))
    // sigma_eq = 400, so d = 1*(1 - exp(-(400/200)^2)) = 1 - exp(-4)
    double expected = 1.0 * (1.0 - exp(-pow(400. / 200., 2.)));
    EXPECT_LT(fabs(dD - expected), 1.E-9);

    // With existing damage, rate should be reduced
    double dD_damaged = damage_weibull(sigma, 0.5, alpha, beta, DTime, "Mises");
    EXPECT_LT(dD_damaged, dD);
    double expected_damaged = 0.5 * (1.0 - exp(-pow(400. / 200., 2.)));
    EXPECT_LT(fabs(dD_damaged - expected_damaged), 1.E-9);
}

TEST(Tdamage, damage_weibull_hydro)
{
    double alpha = 2.0;
    double beta = 200.;
    double damage = 0.;
    double DTime = 1.0;

    // Hydrostatic stress state
    vec sigma = {300., 300., 300., 0., 0., 0.};
    double dD = damage_weibull(sigma, damage, alpha, beta, DTime, "hydro");
    // tr(sigma) = 900, hydro = 900/(3*200) = 1.5
    double expected = 1.0 * (1.0 - exp(-pow(900. / (3. * 200.), 2.)));
    EXPECT_LT(fabs(dD - expected), 1.E-9);
}

TEST(Tdamage, damage_weibull_J3)
{
    double alpha = 2.0;
    double beta = 200.;
    double damage = 0.;
    double DTime = 1.0;

    // Simple stress state for J3
    vec sigma = {400., 0., 0., 0., 0., 0.};
    double dD = damage_weibull(sigma, damage, alpha, beta, DTime, "J3");
    // Just check it returns a valid positive value
    EXPECT_GE(dD, 0.);
}

TEST(Tdamage, damage_kachanov_vonmises)
{
    double A0 = 500.;
    double r = 2.;
    double damage = 0.;

    // Kachanov with vonmises criterion
    vec stress = {100., 0., 0., 0., 0., 0.};
    vec strain = {0.001, -0.0003, -0.0003, 0., 0., 0.};

    double dD = damage_kachanov(stress, strain, damage, A0, r, "vonmises");
    // Result should be positive for non-zero stress
    EXPECT_GT(dD, 0.);

    // Higher stress should give higher damage rate
    vec stress2 = {200., 0., 0., 0., 0., 0.};
    double dD2 = damage_kachanov(stress2, strain, damage, A0, r, "vonmises");
    EXPECT_GT(dD2, dD);
}

TEST(Tdamage, damage_kachanov_hydro)
{
    double A0 = 500.;
    double r = 2.;
    double damage = 0.;

    vec stress = {100., 100., 100., 0., 0., 0.};
    vec strain = {0.001, 0.001, 0.001, 0., 0., 0.};

    double dD = damage_kachanov(stress, strain, damage, A0, r, "hydro");
    EXPECT_GT(dD, 0.);
}

TEST(Tdamage, damage_miner)
{
    double S_max = 300.;
    double S_mean = 100.;
    double S_ult = 500.;
    double b = 0.;
    double B0 = 400.;
    double beta = 1.;

    // Coffin-Manson-like: damage per cycle
    double dN = damage_miner(S_max, S_mean, S_ult, b, B0, beta);
    EXPECT_GT(dN, 0.);

    // Higher S_max should increase damage
    double dN2 = damage_miner(400., S_mean, S_ult, b, B0, beta);
    EXPECT_GT(dN2, dN);
}

TEST(Tdamage, damage_manson)
{
    double S_amp = 200.;
    double C2 = 400.;
    double gamma2 = 2.;

    // Coffin-Manson: d = (S_amp/C2)^gamma2
    double dN = damage_manson(S_amp, C2, gamma2);
    double expected = pow(200. / 400., 2.);
    EXPECT_LT(fabs(dN - expected), 1.E-9);

    // Zero amplitude should give zero damage
    double dN_zero = damage_manson(0., C2, gamma2);
    EXPECT_LT(fabs(dN_zero), 1.E-9);

    // Larger amplitude should give more damage
    double dN_large = damage_manson(300., C2, gamma2);
    EXPECT_GT(dN_large, dN);
}
