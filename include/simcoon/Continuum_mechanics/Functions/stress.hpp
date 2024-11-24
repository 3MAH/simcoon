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

#pragma once
#include <armadillo>


namespace simcoon{

/**
* @file stress.hpp
* @author Yves Chemisky
* @brief A set of functions that transforms stress measures into anothers (considering finite transformations)
*/

/**
 * @brief Provides the first Piola Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 *
 * Returns a matrix, that is the first Piola-Kirchoff stress tensor from
 * - the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param sigma (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the first Piola Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat sigma = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat PKI = Cauchy2PKI(sigma, F);
 * @endcode
*/
arma::mat Cauchy2PKI(const arma::mat &sigma, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Biot stress tensor \f$ \mathbf{T} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 *
 * Returns a matrix, that is the Biot stress tensor from
 * - the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 * Note that the first Piola Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$ is computed, and the Biot stress is obtained through
 * 
 * \f[
 *      \mathbf{T} = \frac{1}{2} \left( \mathbf{R}^T \cdot mathbf{\Sigma} + mathbf{\Sigma}^T \cdot \mathbf{R} \right)
 * \f] 
 * 
 *
 * @param sigma (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param R (3x3 arma::mat, optional) rotation part of the transformation gradient \f$ \mathbf{R} \f$ 
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the Biot stress tensor \f$ \mathbf{T} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat sigma = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat Biot = Cauchy2Biot(sigma, F);
 * @endcode
*/
arma::mat Cauchy2Biot(const arma::mat &sigma, const arma::mat &F, const arma::mat &R = arma::zeros(3,3), const double &J = 0.);

/**
 * @brief Provides the second Piola Kirchoff stress tensor \f$ \mathbf{S} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 *
 * Returns a matrix, that is the first Piola-Kirchoff stress tensor from
 * - the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param sigma (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the second Piola Kirchoff stress tensor \f$ \mathbf{S} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat sigma = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat PKII = Cauchy2PKII(sigma, F);
 * @endcode
*/
arma::mat Cauchy2PKII(const arma::mat &sigma, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma}\f$
 *
 * Returns a matrix, that is the Kirchoff stress tensor from
 * - the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param sigma (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat sigma = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat tau = Cauchy2Kirchoff(sigma, F);
 * @endcode
*/
arma::mat Cauchy2Kirchoff(const arma::mat &sigma, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma}\f$
 *
 * Returns a vector (considering Voigt Notation), that is the Kirchoff stress tensor from
 * - the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param sigma (arma::vec size=6) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (arma::vec size=6) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      vec sigma = m2t_stress(0.5*(matrand + matrand.t()));
 *      mat F = eye(3,3);
 *      vec tau = Cauchy2Kirchoff(sigma, F);
 * @endcode
*/
arma::vec Cauchy2Kirchoff(const arma::vec &sigma, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Cauchy stress tensor \f$ \mathbf{\sigma}\f$ from the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 *
 * Returns a matrix, that is the Cauchy stress tensor from
 * - the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param tau (3x3 arma::mat) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat tau = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat sigma = Kirchoff2Cauchy(tau, F);
 * @endcode
*/
arma::mat Kirchoff2Cauchy(const arma::mat &tau, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Cauchy stress tensor \f$ \mathbf{\sigma}\f$ from the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 *
 * Returns a vecto, that is the Cauchy stress tensor from
 * - the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param tau (arma::vec size=6) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (arma::vec size=6) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      vec tau = m2t_stress(0.5*(matrand + matrand.t()));
 *      mat F = eye(3,3);
 *      mat sigma = Kirchoff2Cauchy(tau, F);
 * @endcode
*/
arma::vec Kirchoff2Cauchy(const arma::vec &tau, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma}\f$
 *
 * Returns a matrix, that is the first Piola-Kirchoff stress tensor from
 * - the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param tau (3x3 arma::mat) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat tau = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat PKI = Kirchoff2PKI(tau, F);
 * @endcode
*/
arma::mat Kirchoff2PKI(const arma::mat &tau, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma}\f$
 *
 * Returns a matrix, that is the second Piola-Kirchoff stress tensor from
 * - the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ \mathbf{J} \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 *
 * @param tau (3x3 arma::mat) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \mathbf{F} \f$
 * @return (3x3 arma::mat) the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 *
 * @details Example:
 * @code
 *      mat matrand = randu(3,3);
 *      mat tau = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat PKII = Kirchoff2PKII(tau, F);
 * @endcode
*/
arma::mat Kirchoff2PKII(const arma::mat &tau, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$ from the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * 
 * Returns a vector, that is the second Piola-Kirchoff stress tensor from
 * - the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ J \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 * 
 * @param tau (arma::vec size=6)) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ 
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \textrm{det}\,\mathbf{F} \f$ 
 * @return (arma::vec size=6) the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * 
 * @details Example: 
 * @code
 *      mat matrand = randu(3,3);
 *      vec tau = m2t_stress(0.5*(matrand + matrand.t()));
 *      mat F = eye(3,3);
 *      vec PKII = Kirchoff2PKII(tau, F);
 * @endcode
*/
arma::vec Kirchoff2PKII(const arma::vec &tau, const arma::mat &F, const double &J= 0.);

/**
 * @brief Provides the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ from the first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 * 
 * Returns a matrix, that is the Kirchoff stress tensor from
 * - the first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ J \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 * 
 * @param Sigma (3x3 arma::mat) first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \textrm{det}\,\mathbf{F} \f$ 
 * @return (3x3 arma::mat) the Kirchoff stress tensor \f$ \mathbf{\tau} \f$
 * 
 * @details Example: 
 * @code
 *      mat PKI = randu(3,3);
 *      mat F = eye(3,3);
 *      mat tau = PKI2Kirchoff(PKI, F);
 * @endcode
*/
arma::mat PKI2Kirchoff(const arma::mat &Sigma, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ from the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * 
 * Returns a matrix, that is the Kirchoff stress tensor from
 * - the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ J \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 * 
 * @param S (3x3 arma::mat) the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \textrm{det}\,\mathbf{F} \f$ 
 * @return (3x3 arma::mat) the Kirchoff stress tensor \f$ \mathbf{\tau}
 * 
 * @details Example: 
 * @code
 *      mat matrand = randu(3,3);
 *      mat PKII = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat tau = PKI2Kirchoff(PKII, F);
 * @endcode
*/
arma::mat PKII2Kirchoff(const arma::mat &S, const arma::mat &F, const double &J = 0.);


/**
 * @brief Provides the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ from the first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 * 
 * Returns a matrix, that is the Cauchy stress tensor from
 * - the first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ J \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 * 
 * @param Sigma (3x3 arma::mat) first Piola-Kirchoff stress tensor \f$ \mathbf{\Sigma} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \textrm{det}\,\mathbf{F} \f$ 
 * @return (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * 
 * @details Example: 
 * @code
 *      mat PKI = randu(3,3);
 *      mat F = eye(3,3);
 *      mat sigma = PKI2Cauchy(PKI, F);
 * @endcode
*/
arma::mat PKI2Cauchy(const arma::mat &Sigma, const arma::mat &F, const double &J = 0.);

/**
 * @brief Provides the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ from the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * 
 * Returns a matrix, that is the Cauchy stress tensor from
 * - the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * - the transformation gradient \f$ \mathbf{F} \f$ and its determinant \f$ J \f$ (this last one is optional, if not indicated, it will be computed from \f$ \mathbf{F} \f$)
 * 
 * @param S (3x3 arma::mat) the second Piola-Kirchoff stress tensor \f$ \mathbf{S} \f$
 * @param F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param J (double, optional) determinant of the transformation gradient \f$ \textrm{det}\,\mathbf{F} \f$ 
 * @return (3x3 arma::mat) the Cauchy stress tensor \f$ \mathbf{\sigma} \f$
 * 
 * @details Example: 
 * @code
 *      mat matrand = randu(3,3);
 *      mat PKII = 0.5*(matrand + matrand.t());
 *      mat F = eye(3,3);
 *      mat sigma = PKII2Cauchy(PKII, F);
 * @endcode
*/
arma::mat PKII2Cauchy(const arma::mat &S, const arma::mat &F, const double &J = 0.);
    
} //namespace simcoon
