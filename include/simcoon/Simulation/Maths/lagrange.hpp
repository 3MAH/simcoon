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

///@file lagrange.hpp
///@brief Function that are utilized in constrained problems
///@version 1.0

#pragma once

namespace simcoon{

/**
 * @file lagrange.hpp
 * @brief Mathematical utility functions.
 */

/** @addtogroup maths
 *  @{
 */


/**
 * @brief Exponential Lagrange multiplier (penalty), as the exponential contact pressure in Abaqus.
 * @details \f[ \lambda(h) = \frac{p_0}{e - 1}\left(\frac{h}{c} + 1\right)\left(e^{h/c + 1} - 1\right) \f]
 * for \f$ h > -c \f$, and \f$ \lambda = 0 \f$ for \f$ h \leq -c \f$.
 * @param h constraint function (gap)
 * @param c width of the transition zone
 * @param p0 multiplier value at \f$ h = 0 \f$
 * @return \f$ \lambda(h) \f$
 */
double lagrange_exp(const double &h, const double &c, const double &p0);

/**
 * @brief First derivative \f$ \mathrm{d}\lambda/\mathrm{d}h \f$ of lagrange_exp.
 * @param h constraint function (gap)
 * @param c width of the transition zone
 * @param p0 multiplier value at \f$ h = 0 \f$
 * @return \f$ \mathrm{d}\lambda/\mathrm{d}h \f$ (0 for \f$ h \leq -c \f$)
 */
double dlagrange_exp(const double &h, const double &c, const double &p0);

/**
 * @brief Power-law Lagrange multiplier enforcing \f$ x \geq 0 \f$.
 * @details \f$ \lambda(x) = p_0 (1 - x)\, x^{-n} \f$ for \f$ x \geq c \f$; below \f$ c \f$ the law is
 * continued by its tangent at \f$ c \f$ plus \f$ \alpha (x - c)^2 \f$, so it stays finite at \f$ x = 0 \f$.
 * @param x constrained variable
 * @param c regularisation threshold
 * @param p0 multiplier scale
 * @param n power-law exponent
 * @param alpha coefficient of the quadratic continuation below \f$ c \f$
 * @return \f$ \lambda(x) \f$
 */
double lagrange_pow_0(const double &x, const double &c, const double &p0, const double &n, const double &alpha);

/**
 * @brief First derivative \f$ \mathrm{d}\lambda/\mathrm{d}x \f$ of lagrange_pow_0.
 * @param x constrained variable
 * @param c regularisation threshold
 * @param p0 multiplier scale
 * @param n power-law exponent
 * @param alpha coefficient of the quadratic continuation below \f$ c \f$
 * @return \f$ \mathrm{d}\lambda/\mathrm{d}x \f$
 */
double dlagrange_pow_0(const double &x, const double &c, const double &p0, const double &n, const double &alpha);

/**
 * @brief Power-law Lagrange multiplier enforcing \f$ x \leq 1 \f$.
 * @details \f$ \lambda(x) = p_0\, x\, (1 - x)^{-n} \f$ for \f$ x \leq 1 - c \f$; above \f$ 1 - c \f$ the law is
 * continued by its tangent at \f$ 1 - c \f$ plus \f$ \alpha (x - 1 + c)^2 \f$, so it stays finite at \f$ x = 1 \f$.
 * @param x constrained variable
 * @param c regularisation threshold (distance to 1)
 * @param p0 multiplier scale
 * @param n power-law exponent
 * @param alpha coefficient of the quadratic continuation above \f$ 1 - c \f$
 * @return \f$ \lambda(x) \f$
 */
double lagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha);

/**
 * @brief First derivative \f$ \mathrm{d}\lambda/\mathrm{d}x \f$ of lagrange_pow_1.
 * @param x constrained variable
 * @param c regularisation threshold (distance to 1)
 * @param p0 multiplier scale
 * @param n power-law exponent
 * @param alpha coefficient of the quadratic continuation above \f$ 1 - c \f$
 * @return \f$ \mathrm{d}\lambda/\mathrm{d}x \f$
 */
double dlagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha);

/**
 * @brief Second derivative \f$ \mathrm{d}^2\lambda/\mathrm{d}x^2 \f$ of lagrange_pow_1 (\f$ 2\alpha \f$ above \f$ 1 - c \f$).
 * @param x constrained variable
 * @param c regularisation threshold (distance to 1)
 * @param p0 multiplier scale
 * @param n power-law exponent
 * @param alpha coefficient of the quadratic continuation above \f$ 1 - c \f$
 * @return \f$ \mathrm{d}^2\lambda/\mathrm{d}x^2 \f$
 */
double d2lagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha);
    

/** @} */ // end of maths group

} //namespace simcoon

