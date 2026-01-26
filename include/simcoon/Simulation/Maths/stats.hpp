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

///@file stats.hpp
///@brief Usefull statistical functions
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file stats.hpp
 * @brief Mathematical utility functions.
 */

/** @addtogroup maths
 *  @{
 */

/**
 * @brief Computes the value of a normal (Gaussian) distribution at a given point.
 * @param x The point at which to evaluate the distribution (double)
 * @param mu The mean of the distribution (double)
 * @param sigma The standard deviation of the distribution (double)
 * @return The probability density value at x (double)
 * @details The normal distribution is given by:
 * \f[
 *     f(x) = \frac{1}{\sigma \sqrt{2\pi}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
 * \f]
 * @code
 *      double p = normal_distrib(0.0, 0.0, 1.0);  // Standard normal at x=0
 * @endcode
 */
double normal_distrib(const double &x, const double &mu, const double &sigma);

/**
 * @brief Computes the probability density of a Weibull distribution.
 * @param x The point at which to evaluate the distribution (double)
 * @param alpha The shape parameter \f$ \alpha > 0 \f$ (double)
 * @param beta The scale parameter \f$ \beta > 0 \f$ (double)
 * @return The probability density value at x (double)
 * @details The Weibull probability density function is:
 * \f[
 *     f(x) = \frac{\alpha}{\beta} \left(\frac{x}{\beta}\right)^{\alpha-1} \exp\left(-\left(\frac{x}{\beta}\right)^\alpha\right)
 * \f]
 * @code
 *      double p = proba_distrib_weibull(1.5, 2.0, 1.0);
 * @endcode
 */
double proba_distrib_weibull(const double &x, const double &alpha, const double &beta);

/**
 * @brief Computes the cumulative distribution function of a Weibull distribution.
 * @param x The point at which to evaluate the CDF (double)
 * @param alpha The shape parameter \f$ \alpha > 0 \f$ (double)
 * @param beta The scale parameter \f$ \beta > 0 \f$ (double)
 * @return The cumulative probability at x (double)
 * @details The Weibull cumulative distribution function is:
 * \f[
 *     F(x) = 1 - \exp\left(-\left(\frac{x}{\beta}\right)^\alpha\right)
 * \f]
 * @code
 *      double F = cumul_distrib_weibull(1.5, 2.0, 1.0);
 * @endcode
 */
double cumul_distrib_weibull(const double &x, const double &alpha, const double &beta);
    
/**
 * @brief Computes the triangular sum of two integers.
 * @param a First integer (int)
 * @param b Second integer (int)
 * @return The triangular sum (int)
 * @details Returns \f$ \frac{(a+b)(a+b+1)}{2} + b \f$, which provides a unique mapping from pairs of integers to a single integer.
 */
int tri_sum(const int &a, const int &b);

/**
 * @brief Computes a classic Orientation Distribution Function (ODF).
 * @param Theta The angle \f$ \Theta \f$ in radians (double)
 * @param Phi The angle \f$ \Phi \f$ in radians (double)
 * @param params Vector of parameters: \f$ [a_1, p_1, a_2, p_2] \f$ (arma::vec)
 * @return The ODF value at the given angles (double)
 * @details The ODF is given by:
 * \f[
 *     \textrm{ODF}(\Theta) = a_1 \cos^{2p_1}(\Theta) + a_2 \cos^{2p_2+1}(\Theta) \sin^{2p_2}(\Theta)
 * \f]
 */
double ODF_sd(const double &Theta, const double &Phi, const arma::vec &params);

/**
 * @brief Computes a hardening-type Orientation Distribution Function.
 * @param Theta The angle \f$ \Theta \f$ in radians (double)
 * @param Phi The angle \f$ \Phi \f$ in radians (double)
 * @param alpha Parameter \f$ \alpha \f$ (double)
 * @param beta Parameter \f$ \beta \f$ (double)
 * @return The ODF value at the given angles (double)
 */
double ODF_hard(const double &Theta, const double &Phi, const double &alpha, const double &beta);

/**
 * @brief Computes a Gaussian peak function.
 * @param x The position at which to evaluate (double)
 * @param mu The center of the peak (double)
 * @param sigma The width parameter (standard deviation) (double)
 * @param ampl The amplitude (default = 1.0) (double)
 * @return The Gaussian value at x (double)
 * @details The Gaussian is given by:
 * \f[
 *     G(x) = A \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
 * \f]
 * @code
 *      double g = Gaussian(0.5, 0.0, 1.0, 2.0);
 * @endcode
 */
double Gaussian(const double &x, const double &mu, const double &sigma, const double &ampl = 1.);

/**
 * @brief Computes a Lorentzian (Cauchy) peak function.
 * @param x The position at which to evaluate (double)
 * @param x0 The center of the peak (double)
 * @param gamma The half-width at half-maximum (double)
 * @param ampl The amplitude (default = 1.0) (double)
 * @return The Lorentzian value at x (double)
 * @details The Lorentzian is given by:
 * \f[
 *     L(x) = \frac{A}{\pi} \frac{\gamma}{(x-x_0)^2 + \gamma^2}
 * \f]
 * @code
 *      double l = Lorentzian(0.5, 0.0, 1.0, 2.0);
 * @endcode
 */
double Lorentzian(const double &x, const double &x0, const double &gamma, const double &ampl = 1.);

/**
 * @brief Computes a Pseudo-Voigt peak function (weighted sum of Gaussian and Lorentzian).
 * @param x The position at which to evaluate (double)
 * @param x0 The center of the peak (double)
 * @param sigma The Gaussian width parameter (double)
 * @param gamma The Lorentzian width parameter (double)
 * @param eta The mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian) (double)
 * @param params Optional additional parameters (arma::vec)
 * @return The Pseudo-Voigt value at x (double)
 * @details The Pseudo-Voigt is given by:
 * \f[
 *     V(x) = \eta \cdot L(x) + (1-\eta) \cdot G(x)
 * \f]
 */
double PseudoVoigt(const double &x, const double &x0, const double &sigma, const double &gamma, const double &eta, const arma::vec &params = arma::ones(1));

/**
 * @brief Computes a Pearson VII peak function.
 * @param x The position at which to evaluate (double)
 * @param x0 The center of the peak (double)
 * @param w The width parameter (double)
 * @param params Vector containing additional shape parameters (arma::vec)
 * @return The Pearson VII value at x (double)
 * @details The Pearson VII function provides flexible peak shapes that can approximate both Gaussian and Lorentzian profiles.
 */
double Pearson7(const double &x, const double &x0, const double &w, arma::vec &params);


/** @} */ // end of maths group

} //namespace simcoon
