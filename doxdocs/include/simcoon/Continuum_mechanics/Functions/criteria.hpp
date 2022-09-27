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

///@file criteria.hpp
///@brief Provide function for yield surfaces
///@version 1.0

#pragma once
#include <string.h>
#include <armadillo>

namespace simcoon{
    
/**
 * @brief Provides the Prager equivalent stress, given its vector representation
 * @param v
 * @return The Prager equivalent stress (double)
 * @details Returns the Prager equivalent stress \f$ \boldsymbol{\sigma}^{P} \f$, considering
\f[
    \sigma^{P} = \sigma^{VM} \left(\frac{1 + b \cdot J_3 \left(\boldsymbol{\sigma} \right)}{\left(J_2 \left(\boldsymbol{\sigma} \right) \right)^{3/2} } \right)^{m}
\f]
    considering the input stress \f$ \boldsymbol{\sigma} \f$, \f$ \boldsymbol{\sigma}^{VM} \f$ is the Von Mises computed equivalent stress, and \f$ b \f$ and \f$ m \f$ are parameter that define the equivalent stress.
    Note that if n > 10, the Prager criteria is sufficiently close to the Mises that the Mises norm is used to avoid numerical instabilities from high-power computations
 * @code 
        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        double sigma_Prager = Prager_stress(sigma, b, n);
 * @endcode
*/
double Prager_stress(const arma::vec &v, const double &b, const double &n);

/**
 * @brief Provides derivative of the Prager equivalent stress, given its vector representation
 * @param v
 * @return The derivative of the Prager equivalent stress (arma::vec)
 * @details Returns the derivative of the Prager equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a convex yield surface
\f[
    \frac{\partial \sigma^{P}}{\partial \mathbf{\sigma}} = \sqrt{3} \left( \frac{1. + b * J_3}{J_2^{3/2}} \right)^{1/n-1} \, \frac{1}{2} \frac{\sqrt{J_2}}{\mathbf{\sigma}'} 
    + b \, m \left( 6 J_2^2 \right) \left( 6 J_2 \mathbf{\sigma}' \cdot \mathbf{\sigma}' - 4 J_2^2 \mathbf{I} + \frac{3}{m-9} J_3 \mathbf{\sigma}' \right)
\f]
    considering the input stress \f$ \boldsymbol{\sigma} \f$, \f$ \boldsymbol{\sigma}^{VM} \f$ is the Von Mises computed equivalent stress, and \f$ b \f$ and \f$ m \f$ are parameter that define the equivalent stress.
    Note that if n > 10, the Prager criteria is sufficiently close to the Mises that the derivative of the Mises norm is used to avoid numerical instabilities from high-power computations
 * @code 
        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        vec dsigma_Pragerdsigma = dPrager_stress(sigma, b, n);
 * @endcode
*/
arma::vec dPrager_stress(const arma::vec &, const double &, const double &);

/**
 * @brief Provides the Tresca equivalent stress, given its vector representation
 * @param v
 * @return The Prager equivalent stress (double)
 * @details Returns the Tresca equivalent stress \f$ \boldsymbol{\sigma}^{T} \f$, considering
\f[
    \sigma^{T} = \sigma_{I} - \sigma_{III},
\f]
    considering the input stress \f$ \boldsymbol{\sigma} \f$.
    Note that the principal stress are classified such that \f$ \sigma_{I} \geq \sigma_{II} \geq \sigma_{III}.
 * @code 
        vec sigma = randu(6);
        double sigma_Tresca = Tresca_stress(sigma);
 * @endcode
*/
double Tresca_stress(const arma::vec &);

//This function returns the the derivative of the Tresca equivalent stress.
arma::vec dTresca_stress(const arma::vec &);
    
//This function returns the Full anisotropic tensor from components
arma::mat P_ani(const arma::vec &);

//This function returns the Hill anisotropic tensor, providing F,G,H,L,M,N
arma::mat P_hill(const arma::vec &);

//This function returns the anisotropic equivalent stress, given a matrix H
double Ani_stress(const arma::vec &, const arma::mat &);

//This function returns the derivative of the anisotropic equivalent stress, given a matrix H
arma::vec dAni_stress(const arma::vec &, const arma::mat &H);
    
//This function returns the anisotropic equivalent stress, providing the anisotropic tensor H components
double Ani_stress(const arma::vec &, const arma::vec &);

//This function returns the derivative of the anisotropic equivalent stress, providing the anisotropic tensor H components
arma::vec dAni_stress(const arma::vec &, const arma::vec &);
    
//This function returns the anisotropic equivalent stress, providing F,G,H,L,M,N
double Hill_stress(const arma::vec &v, const arma::vec &);

//This function returns the derivative of the anisotropic equivalent stress, providing F,G,H,L,M,N
arma::vec dHill_stress(const arma::vec &, const arma::vec &);
    
//This function computes the selected equivalent stress function
double Eq_stress(const arma::vec &, const std::string &, const arma::vec &);

//This function computes the derivative of the selected equivalent stress function
arma::vec dEq_stress(const arma::vec &, const std::string &, const arma::vec &);
    
} //namespace simcoon
