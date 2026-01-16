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
 * Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor.
\f[
    \mathbf{S}= \begin{array}{cccccc}
        1 & \frac{\nu}{1-\nu} & \frac{\nu}{1-\nu} & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 1 \end{array} \right)
\f]
 * @details Example: 
 * @code 
       mat S = Eshelby_penny(nu);
 * @endcode
*/
arma::mat Eshelby_penny(const double &);

//This methods is using the Voigt notations for the tensors.
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
 * in the \f$ 1 \f$ and \f$ 2 \f$ directions. (The points and weights are calculated using the internal *point* function).
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
 * @param a1 the aspect ratio in the direction \f$ 1 \f1
 * @param a2 the aspect ratio in the direction \f$ 2 \f1
 * @param a3 the aspect ratio in the direction \f$ 3 \f1
 * @param mp the number of integration points in the direction \f$ 1 \f$
 * @param np the number of integration points in the direction \f$ 2 \f$
 * @return Provides the numerical estimation of the Hill interaction tensor of an ellispoid in the general case of anisotropic media,
 * as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction
 * \f$ 1 \f$, \f$ 2 \f$ and \f$ 3 \f$, respectively. The list of integration and their respective weight
 * for the numerical integration is automatically determined from the number of integration points 
 * in the \f$ 1 \f$ and \f$ 2 \f$ directions. (The points and weights are calculated using the internal *point* function).
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

} //namespace simcoon
