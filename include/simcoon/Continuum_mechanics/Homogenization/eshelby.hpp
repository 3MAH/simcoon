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

///@file eshelby.hpp
///@brief Definition of the Eshelby tensor for ellipsoidal inclusions with
// Parts of this methods are copyrighted by Gavazzi & Lagoudas 1992 - Fair use only
///@version 1.0

#pragma once

#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/**
 * @file eshelby.hpp
 * @brief Eshelby tensor computation functions for homogenization.
 */

/** @addtogroup eshelby
 *  @{
 */


/**
 * @brief Provides the Eshelby tensor of a spherical inclusion for isotropic linear elasticity.
 * @param nu the Poisson ratio
 * @return Provides the Eshelby tensor of a spherical inclusion for isotropic linear elasticity in the Simcoon
 * formalism. Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor, as a function
 * of the Poisson ratio \f$ \nu \f$
\f[
    \mathbf{S}= \begin{array}{cccccc}
    \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
    \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
    \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
    0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
    0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
    0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{array} \right)
\f]
 * @details Example: 
 * @code 
       mat S = Eshelby_sphere(nu);
 * @endcode
*/
arma::mat Eshelby_sphere(const double &);

/**
 * @brief Provides the Eshelby tensor of a cylindrical inclusion for isotropic linear elasticity
 * @param nu the Poisson ratio
 * @return  Provides the Eshelby tensor of a cylindrical inclusion for isotropic linear elasticity in the Simcoon formalism,
 * as a function of the Poisson ratio \f$ \nu \f$. The cylinder is oriented such as the longitudinal axis is the axis \f$ 1 \f$.
 * Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor.
\f[
    \mathbf{S}= \begin{array}{cccccc}
        \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{array} \right)
\f]
 * @details Example: 
 * @code 
       mat S = Eshelby_cylinder(nu);
 * @endcode
*/
arma::mat Eshelby_cylinder(const double &);

/**
 * @brief Provides the Eshelby tensor of a prolate inclusion for isotropic linear elasticity
 * @param nu the Poisson ratio
 * @param a_r Aspect ratio defined as \f$ a_r = \frac{a1}{a2} = \frac{a1}{a3} \f$
 * @return Provides the Eshelby tensor of a prolate inclusion for isotropic linear elasticity in the Simcoon formalism, 
 * as a function of the Poisson ratio \f$ \nu \f$ and the aspect ratio \f$ a_r = \frac{a1}{a2} = \frac{a1}{a3} \f$. 
 * The prolate inclusion is oriented such as the axis of rotation is the axis \f$ 1 \f$.
 * Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor.
\f[
    \mathbf{S}= \begin{array}{cccccc}
        S_{11} & S_{12} & S_{12} & 0 & 0 & 0 \\ 
        S_{21} & S_{22} & S_{23} & 0 & 0 & 0 \\
        S_{21} & S_{23} & S_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & S_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & S_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & S_{66} \end{array} \right)
\f]
 * with
\f[     \begin{array}{cc}
        S_{11} & = \frac{1}{2(1-\nu)}\left(1-2\nu+\frac{3a_r^2-1}{a_r^2-1}-g\left(1-2\nu+\frac{3a_r^2}{a_r^2-1}\right)\right) \\
        S_{12} & = \frac{-1}{2(1-\nu)}\left(1-2\nu+\frac{1}{a_r^2-1}+g\left(1-2\nu+\frac{3}{a_r^2-1}\right)\right) \\
        S_{21} & = \frac{-a_r^2}{2(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(\frac{3a_r^2}{a_r^2-1}-\left(1-2\nu\right)\right) \\
        S_{22} & = \frac{3a_r^2}{8(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(1-2\nu-\frac{9}{4\left(a_r^2-1\right)}\right) \\
        S_{23} & = \frac{1}{4(1-\nu)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}-g\left(1-2\nu+\frac{3}{4\left(a_r^2-1\right)}\right)\right) \\
        S_{44} & = \frac{2}{4\left(1-\nu\right)}\left(1-2\nu-\frac{a_r^2+1}{a_r^2-1}-\frac{g}{2}\left(1-2\nu-\frac{3a_r^2+1}{a_r^2-1}\right)\right) \\
        S_{66} & = \frac{2}{4\left(1-\nu\right)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}+g\left(1-2\nu-\frac{3}{4\left(a_r^2-1\right)}\right)\right)
        \end{array}        
\f]
 * considering \f$ g = a_r\frac{a_r\sqrt{a_r^2-1}}{\left(a_r^2-1\right)^{\frac{3}{2}}} - acos(a_r) \f$
 *
 * @details Example: 
 * @code 
       mat S = Eshelby_prolate(nu,a_r);
 * @endcode
*/
arma::mat Eshelby_prolate(const double &, const double &);

/**
 * @brief Provides the Eshelby tensor of a oblate inclusion for isotropic linear elasticity
 * @param nu the Poisson ratio
 * @param a_r Aspect ratio defined as \f$ a_r = \frac{a1}{a2} = \frac{a1}{a3} \f$
 * @return Provides the Eshelby tensor of a oblate inclusion for isotropic linear elasticity in the Simcoon formalism,
 * as a function of the Poisson ratio \f$ \nu \f$ and the aspect ratio \f$ a_r = \frac{a1}{a2} = \frac{a1}{a3} \f$.
 * The oblate inclusion is oriented such as the axis of rotation is the axis \f$ 1 \f$.
\f[
    \mathbf{S}= \begin{array}{cccccc}
        S_{11} & S_{12} & S_{12} & 0 & 0 & 0 \\ 
        S_{21} & S_{22} & S_{23} & 0 & 0 & 0 \\
        S_{21} & S_{23} & S_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & S_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & S_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & S_{66} \end{array} \right)
\f]
 * with
\f[     \begin{array}{cc}
        S_{11} &= \frac{1}{2(1-\nu)}\left(1-2\nu+\frac{3a_r^2-1}{a_r^2-1}-g\left(1-2\nu+\frac{3a_r^2}{a_r^2-1}\right)\right) \\
        S_{12} &= \frac{-1}{2(1-\nu)}\left(1-2\nu+\frac{1}{a_r^2-1}+g\left(1-2\nu+\frac{3}{a_r^2-1}\right)\right) \\
        S_{21} &= \frac{-a_r^2}{2(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(\frac{3a_r^2}{a_r^2-1}-\left(1-2\nu\right)\right) \\
        S_{22} &= \frac{3a_r^2}{8(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(1-2\nu-\frac{9}{4\left(a_r^2-1\right)}\right) \\
        S_{23} &= \frac{1}{4(1-\nu)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}-g\left(1-2\nu+\frac{3}{4\left(a_r^2-1\right)}\right)\right) \\
        S_{44} &= \frac{2}{4\left(1-\nu\right)}\left(1-2\nu-\frac{a_r^2+1}{a_r^2-1}-\frac{g}{2}\left(1-2\nu-\frac{3a_r^2+1}{a_r^2-1}\right)\right) \\
        S_{66} &= \frac{2}{4\left(1-\nu\right)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}+g\left(1-2\nu-\frac{3}{4\left(a_r^2-1\right)}\right)\right)
        \end{array}        
\f]
 * considering \f$ g = a_r\frac{-a_r\sqrt{1-a_r^2}}{\left(1-a_r^2\right)^{\frac{3}{2}}} - acos(a_r) \f$
 *
 * @details Example: 
 * @code 
       mat S = Eshelby_oblate(nu,a_r);
 * @endcode
*/
arma::mat Eshelby_oblate(const double &, const double &);

/**
 * @brief Provides the Eshelby tensor of a penny-shaped (crack) inclusion for isotropic linear elasticity
 * @param nu the Poisson ratio
 * @return Provides the Eshelby tensor of a penny-shaped (crack) inclusion for isotropic linear elasticity in the Simcoon formalism,
 * as a function of the Poisson ratio \f$ \nu \f$. This corresponds to the limit case of an oblate spheroid when the aspect ratio
 * \f$ a_r = \frac{a1}{a2} = \frac{a1}{a3} \rightarrow 0 \f$.
 * The penny-shaped inclusion is oriented such as the normal to the crack plane is the axis \f$ 1 \f$.
 * Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor
 * (Voigt order [11, 22, 33, 12, 13, 23]): \f$ S_{1212} = S_{1313} = 1/2 \f$ (1 with the engineering
 * shear), the in-plane shear \f$ S_{2323} = 0 \f$.
\f[
    \mathbf{S}= \left( \begin{array}{cccccc}
        1 & \frac{\nu}{1-\nu} & \frac{\nu}{1-\nu} & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 1 & 0 & 0 \\
        0 & 0 & 0 & 0 & 1 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \end{array} \right)
\f]
 * @details Example: 
 * @code 
       mat S = Eshelby_penny(nu);
 * @endcode
*/
arma::mat Eshelby_penny(const double &);

/// Uses the Voigt notation for the tensors.
void calG(const double &, const double &, const double &, const double &, const double &, const arma::Mat<int> &, const arma::mat &, arma::mat &);

//Weighted Gauss integration over a sphere to represent the integration over the ellipsoid
void Gauss(arma::Mat<int> &, const arma::mat &, arma::mat &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const int &, const int &);

/**
 * @brief Provides a numerical estimation of the Eshelby tensor of an ellipsoisal, possibly anisotropic inclusion, 
 * from the full parametrization of integration points
 * @param L the stiffness tensor of the considered material
 * @param a1 the aspect ratio in the direction \f$ 1 \f$
 * @param a2 the aspect ratio in the direction \f$ 2 \f$
 * @param a3 the aspect ratio in the direction \f$ 3 \f$
 * @param x the position vector of integration points in the direction \f$ 1 \f$
 * @param wx the weight coefficient of integration points in the direction \f$ 1 \f$
 * @param y the position vector of integration points in the direction \f$ 2 \f$
 * @param wy the weight coefficient of integration points in the direction \f$ 2 \f$
 * @param mp the number of integration points in the direction \f$ 1 \f$
 * @param np the number of integration points in the direction \f$ 2 \f$
 * @return Provides the numerical estimation of the Eshelby tensor of an ellispoid in the general case of anisotropic media,
 * as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
 * \f$ 1 \f$, \f$ 2 \f$ and \f$ 3 \f$, respectively. It also requires the list of integration points and their respective weight
 * for the numerical integration, as well as the number of integration points in the \f$ 1 \f$ and \f$ 2 \f$ directions.
 * The points and weights are calculated using the *point* function that require to be called previously. 
 * Returns the Eshelby tensor as a mat, according to the convention of a localisation tensor.
 * @details Example: 
 * @code 
        mat L = L_iso(210000., 0.3, "Enu");
        double a1 = 2.0;
        double a2 = 1.0;
        double a3 = 1.0;
        int mp = 4;
        int np = 4;
        vec x, wx, y, wy;
        points(x, wx, y, wy, mp, np);
        mat S = Eshelby(L, a1, a2, a3, x, wx, y, wy, mp, np);

 * @endcode
*/
arma::mat Eshelby(const arma::mat &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const int &mp, const int &np);

/**
 * @brief Provides a numerical estimation of the Eshelby tensor of an ellipsoisal, possibly anisotropic inclusion.
 * @param L the stiffness tensor of the considered material
 * @param a1 the aspect ratio in the direction \f$ 1 \f1
 * @param a2 the aspect ratio in the direction \f$ 2 \f1
 * @param a3 the aspect ratio in the direction \f$ 3 \f1
 * @param mp the number of integration points in the direction \f$ 1 \f$
 * @param np the number of integration points in the direction \f$ 2 \f$
 * @return Provides the numerical estimation of the Eshelby tensor of an ellispoid in the general case of anisotropic media,
 * as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
 * \f$ 1 \f$, \f$ 2 \f$ and \f$ 3 \f$, respectively. The list of integration and their respective weight
 * for the numerical integration is automatically determined from the number of integration points 
 * in the \f$ 1 \f$ and \f$ 2 \f$ directions. (The points and weights are calculated once per (mp, np), see eshelby_quadrature_points).
 * Returns the Eshelby tensor as a mat, according to the convention of a localisation tensor.
 * @details Example: 
 * @code 
        mat L = L_iso(210000., 0.3, "Enu");
        double a1 = 2.0;
        double a2 = 1.0;
        double a3 = 1.0;
        int mp = 4;
        int np = 4;
        mat S = Eshelby(L, a1, a2, a3, mp, np);
 * @endcode
*/
arma::mat Eshelby(const arma::mat &, const double &, const double &, const double &, const int &, const int &);

/**
 * @brief Provides a numerical estimation of the Hill interaction tensor of an ellipsoisal, possibly anisotropic inclusion, 
 * from the full parametrization of integration points
 * @param L the stiffness tensor of the considered material
 * @param a1 the aspect ratio in the direction \f$ 1 \f1
 * @param a2 the aspect ratio in the direction \f$ 2 \f1
 * @param a3 the aspect ratio in the direction \f$ 3 \f1
 * @param x the position vector of integration points in the direction \f$ 1 \f$
 * @param wx the weight coefficient of integration points in the direction \f$ 1 \f$
 * @param y the position vector of integration points in the direction \f$ 2 \f$
 * @param wy the weight coefficient of integration points in the direction \f$ 2 \f$
 * @param mp the number of integration points in the direction \f$ 1 \f$
 * @param np the number of integration points in the direction \f$ 2 \f$
 * @return Provides the numerical estimation of the Hill interaction tensor of an ellispoid in the general case of anisotropic media,
 * as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
 * \f$ 1 \f$, \f$ 2 \f$ and \f$ 3 \f$, respectively. It also requires the list of integration points and their respective weight
 * for the numerical integration, as well as the number of integration points in the \f$ 1 \f$ and \f$ 2 \f$ directions.
 * The points and weights are calculated using the *point* function that require to be called previously. 
 * Returns the Hill interaction tensor as a mat, according to the convention of a localisation tensor.
 * @details Example: 
 * @code 
        mat L = L_iso(210000., 0.3, "Enu");
        double a1 = 2.0;
        double a2 = 1.0;
        double a3 = 1.0;
        int mp = 4;
        int np = 4;
        vec x, wx, y, wy;
        points(x, wx, y, wy, mp, np);
        mat S = T_II(L, a1, a2, a3, x, wx, y, wy, mp, np);

 * @endcode
*/
arma::mat T_II(const arma::mat &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const int &, const int &);

/**
 * @brief Provides a numerical estimation of the Hill interaction tensor of an ellipsoisal, possibly anisotropic inclusion.
 * @param L the stiffness tensor of the considered material
 * @param a1 the aspect ratio in the direction \f$ 1 \f$
 * @param a2 the aspect ratio in the direction \f$ 2 \f$
 * @param a3 the aspect ratio in the direction \f$ 3 \f$
 * @param mp the number of integration points in the direction \f$ 1 \f$
 * @param np the number of integration points in the direction \f$ 2 \f$
 * @return Provides the numerical estimation of the Hill interaction tensor of an ellispoid in the general case of anisotropic media,
 * as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
 * \f$ 1 \f$, \f$ 2 \f$ and \f$ 3 \f$, respectively. The list of integration and their respective weight
 * for the numerical integration is automatically determined from the number of integration points 
 * in the \f$ 1 \f$ and \f$ 2 \f$ directions. (The points and weights are calculated once per (mp, np), see eshelby_quadrature_points).
 * Returns the Hill interaction tensor as a mat, according to the convention of a localisation tensor.
 * @details Example: 
 * @code 
        mat L = L_iso(210000., 0.3, "Enu");
        double a1 = 2.0;
        double a2 = 1.0;
        double a3 = 1.0;
        int mp = 4;
        int np = 4;
        mat S = T_II(L, a1, a2, a3, mp, np);
 * @endcode
*/
arma::mat T_II(const arma::mat &, const double &, const double &, const double &, const int &, const int &);
    
/**
 * @brief Provides a numerical estimation of the Hill interaction tensor of an ellipsoisal, possibly anisotropic inclusion.
 * @param x the position vector of integration points in the direction \f$ 1 \f$
 * @param wx the weight coefficient of integration points in the direction \f$ 1 \f$
 * @param y the position vector of integration points in the direction \f$ 2 \f$
 * @param wy the weight coefficient of integration points in the direction \f$ 2 \f$
 * @param mp the number of integration points in the direction \f$ 1 \f$
 * @param np the number of integration points in the direction \f$ 2 \f$
 * @details 
 * This methods computes the list of integration and their respective weight for the numerical integration,
 * as a function of the list of integration points and their respective weight
 * for the numerical integration, as well as the number of integration points in the \f$ 1 \f$ and \f$ 2 \f$ directions.
 * The methods actually update *x*, *wx*, *y* and *wy* according to *mp* and *np*. 
 * Note that *x*, *wx*, *y*, *wy* have to be initialized first with the size of *mp* and *np*, respectively.
 * 
 * Example: 
 * @code 
        vec x(mp);
        vec wx(mp);
        vec y(np);
        vec wy(np);
        points(x, wx, y, wy, mp, np);
 * @endcode
*/
void points(arma::vec &, arma::vec &, arma::vec &, arma::vec &, const int &, const int &);

/**
 * @brief How the Eshelby tensors of a mean-field scheme are evaluated: the quadrature tables of
 * the integral for one (mp, np, rule) and whether closed forms may replace it. Shared and immutable.
 *
 * @details The only way to obtain one is eshelby_quadrature_points(), which builds each
 * (mp, np, rule, closed_form) once and hands out a const reference to it: never modified nor
 * freed afterwards, so any thread can use it while other threads integrate with other settings.
 * A setting identifies its tables exactly (the object address can serve as a cache key).
 *
 * Rules (the integrand is centrally symmetric, G(-xi) = G(xi)):
 * - 0: Gauss-Legendre on the whole sphere, mp points in cos(theta) over [-1, 1] and np points in
 *   phi over [0, 2 pi] (the historical rule, mp x np evaluations);
 * - 1: half sphere, (mp+1)/2 Gauss-Legendre points in cos(theta) over [0, 1] (weights doubled)
 *   and np equally spaced points in phi with the trapezoidal rule, which converges spectrally on
 *   this periodic integrand: half the evaluations for an accuracy as good or better.
 */
struct eshelby_quadrature {
    int mp;             ///< number of polar points in the tables (mp, or (mp+1)/2 for rule 1)
    int np;             ///< number of azimuthal points in the tables
    int rule;           ///< 0 Gauss-Legendre on the sphere, 1 half sphere with a trapezoidal rule in phi
    bool closed_form;   ///< use the closed forms (isotropic medium, sphere or spheroid) when they apply
    arma::vec x;        ///< polar points (cos(theta))
    arma::vec wx;       ///< polar weights
    arma::vec y;        ///< azimuthal points
    arma::vec wy;       ///< azimuthal weights

    eshelby_quadrature(const eshelby_quadrature &) = delete;
    eshelby_quadrature &operator=(const eshelby_quadrature &) = delete;

private:
    int mp_asked;       ///< the mp of the request (the key)
    eshelby_quadrature(const int &, const int &, const int &, const bool &);
    friend const eshelby_quadrature &eshelby_quadrature_points(const int &, const int &, const int &, const bool &);
};

/**
 * @brief The shared Eshelby evaluation settings and quadrature tables for (mp, np, rule, closed_form).
 * @param mp the number of integration points in the polar direction (at least 1)
 * @param np the number of integration points in the azimuthal direction (at least 1)
 * @param rule 0 (default) Gauss-Legendre on the sphere, 1 half sphere with a trapezoidal rule in phi
 * @param closed_form false (default): always integrate; true: the schemes use the closed forms
 *        (Eshelby_sphere, Eshelby_prolate, Eshelby_oblate) when the medium is isotropic and the
 *        inclusion a sphere or a spheroid of axis 1, and integrate otherwise
 * @return A const reference, valid until the process exits.
 * @details Thread-safe: each new setting is built once under a lock, and each thread remembers
 * the last one it asked for. Rule 0 tables are those of *points*.
 * @code
        const eshelby_quadrature &q = eshelby_quadrature_points(50, 50);
        mat S = Eshelby(L, a1, a2, a3, q.x, q.wx, q.y, q.wy, q.mp, q.np);
 * @endcode
*/
const eshelby_quadrature &eshelby_quadrature_points(const int &mp, const int &np, const int &rule = 0, const bool &closed_form = false);

/**
 * @brief Closed-form Eshelby tensor of an inclusion in an isotropic medium, when one applies.
 * @param L stiffness of the medium (6x6), in the inclusion frame
 * @param a1, a2, a3 semi-axes of the inclusion
 * @param S the Eshelby tensor, set when the function returns true
 * @return true when L is isotropic (every entry within 1e-10 of its largest one from the
 * isotropic tensor of lambda = L(0,1), mu = L(3,3)) and the inclusion is a sphere or a spheroid
 * of axis 1 (a2 = a3) whose aspect ratio a1/a2 is 1 or at least 1% away from 1 (closer to 1 the
 * spheroid formulas lose accuracy by cancellation); false otherwise, S untouched.
 */
bool Eshelby_closed_form(const arma::mat &L, const double &a1, const double &a2, const double &a3, arma::mat &S);


/** @} */ // end of eshelby group

} //namespace simcoon
