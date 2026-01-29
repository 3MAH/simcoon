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

///@file Tstats.cpp
///@brief Test for statistical functions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/stats.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tstats, normal_distrib)
{
    // CDF of standard normal at 0 should be 0.5
    double cdf_at_zero = normal_distrib(0., 0., 1.);
    EXPECT_LT(fabs(cdf_at_zero - 0.5), 1.E-6);

    // CDF at large positive should be close to 1
    double cdf_at_large = normal_distrib(5., 0., 1.);
    EXPECT_LT(fabs(cdf_at_large - 1.0), 1.E-4);

    // CDF at large negative should be close to 0
    double cdf_at_neg = normal_distrib(-5., 0., 1.);
    EXPECT_LT(cdf_at_neg, 1.E-4);

    // Symmetry: CDF(x) + CDF(-x) = 1
    double x = 1.5;
    double sum_sym = normal_distrib(x, 0., 1.) + normal_distrib(-x, 0., 1.);
    EXPECT_LT(fabs(sum_sym - 1.0), 1.E-6);

    // Non-zero mean: CDF(mu, mu, sigma) = 0.5
    double cdf_at_mean = normal_distrib(3.0, 3.0, 2.0);
    EXPECT_LT(fabs(cdf_at_mean - 0.5), 1.E-6);
}

TEST(Tstats, weibull_distributions)
{
    double alpha = 2.0;
    double beta = 1.0;

    // PDF at x=0 should be 0 (Macaulay bracket makes it 0 for x<=0)
    double pdf_at_zero = proba_distrib_weibull(0., alpha, beta);
    EXPECT_LT(fabs(pdf_at_zero), 1.E-9);

    // PDF at negative should be 0
    double pdf_neg = proba_distrib_weibull(-1., alpha, beta);
    EXPECT_LT(fabs(pdf_neg), 1.E-9);

    // PDF at positive x should be positive
    double pdf_pos = proba_distrib_weibull(0.5, alpha, beta);
    EXPECT_GT(pdf_pos, 0.);

    // CDF at x=0 should be 0
    double cdf_at_zero = cumul_distrib_weibull(0., alpha, beta);
    EXPECT_LT(fabs(cdf_at_zero), 1.E-9);

    // CDF at negative should be 0
    double cdf_neg = cumul_distrib_weibull(-1., alpha, beta);
    EXPECT_LT(fabs(cdf_neg), 1.E-9);

    // CDF should be monotonically increasing
    double cdf1 = cumul_distrib_weibull(0.5, alpha, beta);
    double cdf2 = cumul_distrib_weibull(1.0, alpha, beta);
    double cdf3 = cumul_distrib_weibull(2.0, alpha, beta);
    EXPECT_LT(cdf1, cdf2);
    EXPECT_LT(cdf2, cdf3);

    // CDF at large x should approach 1
    double cdf_large = cumul_distrib_weibull(10., alpha, beta);
    EXPECT_LT(fabs(cdf_large - 1.0), 1.E-6);

    // Exponential distribution case: alpha=1
    // CDF = 1 - exp(-x/beta) for x >= 0
    double x_test = 1.0;
    double beta_exp = 2.0;
    double cdf_exp = cumul_distrib_weibull(x_test, 1.0, beta_exp);
    double cdf_ref = 1.0 - exp(-x_test / beta_exp);
    EXPECT_LT(fabs(cdf_exp - cdf_ref), 1.E-9);
}

TEST(Tstats, tri_sum)
{
    // tri_sum(a, b) = b*(2*a + (b-1))/2
    // For a=1, b=1: 1*(2+0)/2 = 1
    EXPECT_EQ(tri_sum(1, 1), 1);

    // For a=1, b=3: 3*(2+2)/2 = 6
    EXPECT_EQ(tri_sum(1, 3), 6);

    // For a=3, b=2: 2*(6+1)/2 = 7
    EXPECT_EQ(tri_sum(3, 2), 7);

    // For a=5, b=4: 4*(10+3)/2 = 26
    EXPECT_EQ(tri_sum(5, 4), 26);
}

TEST(Tstats, ODF_functions)
{
    // ODF_sd at Theta=mean should return alpha1
    vec params = {1.0, 0.5, 2.0, 1.0};
    double odf_at_mean = ODF_sd(0., 0., params);
    EXPECT_LT(fabs(odf_at_mean - 1.0), 1.E-6);

    // ODF_sd at Theta = mean + pi/2 should return 0
    double odf_at_perp = ODF_sd(0.5 * sim_pi, 0., params);
    EXPECT_LT(fabs(odf_at_perp), 1.E-6);

    // ODF_hard should be positive at the mean
    double odf_hard = ODF_hard(0., 0., 1.0, 1.0);
    EXPECT_GT(odf_hard, 0.);

    // ODF_hard peak is at Theta = mean, value = ampl
    EXPECT_LT(fabs(odf_hard - 1.0), 1.E-9);

    // ODF_hard should decay away from mean
    double odf_hard_away = ODF_hard(2.0, 0., 1.0, 1.0);
    EXPECT_LT(odf_hard_away, odf_hard);
}

TEST(Tstats, Gaussian)
{
    double mu = 0.;
    double sigma = 1.;
    double ampl = 1.;

    // Peak at x = mu
    double g_peak = Gaussian(mu, mu, sigma, ampl);
    double g_away = Gaussian(mu + 1., mu, sigma, ampl);
    EXPECT_GT(g_peak, g_away);

    // Known value: Gaussian(0,0,1,1) = 1/(sqrt(2*pi)) ~ 0.3989
    double expected = 1.0 / (sigma * sqrt(2.0 * sim_pi));
    EXPECT_LT(fabs(g_peak - expected), 1.E-9);

    // Symmetry: G(x) = G(-x) for mean=0
    double g_pos = Gaussian(1.0, 0., 1., 1.);
    double g_neg = Gaussian(-1.0, 0., 1., 1.);
    EXPECT_LT(fabs(g_pos - g_neg), 1.E-9);
}

TEST(Tstats, Lorentzian)
{
    double x0 = 0.;
    double gamma = 1.;
    double ampl = 1.;

    // Peak at x = x0
    double l_peak = Lorentzian(x0, x0, gamma, ampl);
    double l_away = Lorentzian(x0 + 1., x0, gamma, ampl);
    EXPECT_GT(l_peak, l_away);

    // Known value: L(0,0,1,1) = 1*1/(2*pi*(0 + 0.25)) = 1/(0.5*pi) = 2/pi
    double expected = ampl * gamma / (2.0 * sim_pi * pow(gamma / 2.0, 2.0));
    EXPECT_LT(fabs(l_peak - expected), 1.E-9);

    // Symmetry
    double l_pos = Lorentzian(1.0, 0., 1., 1.);
    double l_neg = Lorentzian(-1.0, 0., 1., 1.);
    EXPECT_LT(fabs(l_pos - l_neg), 1.E-9);
}

TEST(Tstats, PseudoVoigt)
{
    double x0 = 0.;
    double sigma = 1.;
    double gamma = 1.;
    double ampl = 1.;

    // eta=0 should give pure Gaussian
    vec params_gauss = {0.0};
    double pv_gauss = PseudoVoigt(0.5, x0, sigma, gamma, ampl, params_gauss);
    double g_ref = Gaussian(0.5, x0, sigma, ampl);
    EXPECT_LT(fabs(pv_gauss - g_ref), 1.E-9);

    // eta=1 should give pure Lorentzian
    vec params_lor = {1.0};
    double pv_lor = PseudoVoigt(0.5, x0, sigma, gamma, ampl, params_lor);
    double l_ref = Lorentzian(0.5, x0, gamma, ampl);
    EXPECT_LT(fabs(pv_lor - l_ref), 1.E-9);

    // eta=0.5 should give average
    vec params_mix = {0.5};
    double pv_mix = PseudoVoigt(0.5, x0, sigma, gamma, ampl, params_mix);
    double expected = 0.5 * l_ref + 0.5 * g_ref;
    EXPECT_LT(fabs(pv_mix - expected), 1.E-9);
}

TEST(Tstats, Pearson7)
{
    double x0 = 0.;
    double inv_width = 1.0;
    vec params = {2.0, 3.0}; // max=2.0, shape=3.0

    // At center x=x0: result = max * (1 + 0)^(-shape) = max = 2.0
    double p7 = Pearson7(x0, x0, inv_width, params);
    EXPECT_LT(fabs(p7 - 2.0), 1.E-9);

    // Away from center, value should be less
    double p7_away = Pearson7(1.0, x0, inv_width, params);
    EXPECT_LT(p7_away, p7);
}
