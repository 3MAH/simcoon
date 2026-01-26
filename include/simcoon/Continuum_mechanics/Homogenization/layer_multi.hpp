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

///@file layer_multi.hpp
///@brief Characteristics of a layer, similar to a phase in this version for homogenization purposes
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

/**
 * @file layer_multi.hpp
 * @brief Layered composite homogenization
 * @author Yves Chemisky
 * @version 1.0
 *
 * This class extends phase_multi for layered (laminated) composites.
 * It implements the mechanics for homogenization of stratified media
 * with parallel layers of different materials.
 */

/** @addtogroup homogenization
 *  @{
 */

/**
 * @brief Layer class for laminated composite homogenization
 *
 * @details This class implements the mechanics of layered composites for
 * mean-field homogenization. It is used for materials with a stratified
 * microstructure where layers are stacked in a specific direction.
 *
 * **Applications:**
 * - Laminated composites (e.g., carbon fiber laminates)
 * - Layered geological media
 * - Multilayer coatings
 * - Sandwich structures
 *
 * **Coordinate System:**
 * - The stacking direction is aligned with the normal direction (n)
 * - The in-plane directions are tangential (t)
 *
 * **Stress/Strain Decomposition:**
 *
 * For layered media, it is convenient to decompose stress and strain into
 * normal (n) and tangential (t) components:
 * - Normal stress components: \f$ \sigma_{nn}, \sigma_{nt} \f$
 * - Tangential strain components: \f$ \varepsilon_{tt} \f$
 *
 * The continuity conditions at layer interfaces impose:
 * - Continuous normal stress: \f$ \sigma_n^{(1)} = \sigma_n^{(2)} \f$
 * - Continuous tangential strain: \f$ \varepsilon_t^{(1)} = \varepsilon_t^{(2)} \f$
 *
 * @see phase_multi for base class documentation
 */
class layer_multi : public phase_multi
{
	private:

	protected:

	public :

    /**
     * @brief Normal-normal block of the tangent modulus (3×3)
     *
     * Submatrix of the stiffness tensor relating normal stress to normal strain:
     * \f$ \mathbf{D}_{nn} = \frac{\partial \boldsymbol{\sigma}_n}{\partial \boldsymbol{\varepsilon}_n} \f$
     */
    arma::mat Dnn;

    /**
     * @brief Normal-tangential block of the tangent modulus (3×3)
     *
     * Submatrix coupling normal stress to tangential strain:
     * \f$ \mathbf{D}_{nt} = \frac{\partial \boldsymbol{\sigma}_n}{\partial \boldsymbol{\varepsilon}_t} \f$
     */
    arma::mat Dnt;

    /**
     * @brief Derivative of deformation gradient with respect to normal position (3×3)
     *
     * Used in incremental formulations for geometric nonlinearity
     */
    arma::mat dXn;

    /**
     * @brief Derivative of deformation gradient with respect to tangential position (3×3)
     *
     * Used in incremental formulations for geometric nonlinearity
     */
    arma::mat dXt;

    /**
     * @brief Stress-like quantity for interface conditions (6×1)
     *
     * Intermediate variable used in the layer homogenization algorithm
     */
    arma::vec sigma_hat;

    /**
     * @brief Derivative of local field with respect to layer position (6×1)
     *
     * Used for computing gradients through the layer thickness
     */
    arma::vec dzdx1;

    /** @brief Default constructor */
    layer_multi();

    /**
     * @brief Full constructor with all parameters
     *
     * @param A Strain concentration tensor
     * @param A_start Initial strain concentration tensor
     * @param B Stress concentration tensor
     * @param B_start Initial stress concentration tensor
     * @param A_in Inelastic concentration vector
     * @param Dnn Normal-normal stiffness block
     * @param Dnt Normal-tangential stiffness block
     * @param dXn Deformation gradient derivative (normal)
     * @param dXt Deformation gradient derivative (tangential)
     * @param sigma_hat Interface stress variable
     * @param dzdx1 Field gradient through thickness
     */
    layer_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::vec&);

    /**
     * @brief Copy constructor
     * @param lm layer_multi object to copy
     */
    layer_multi(const layer_multi&);

    /** @brief Destructor */
    ~layer_multi();

    /**
     * @brief Assignment operator
     * @param lm layer_multi object to assign from
     * @return Reference to this object
     */
    virtual layer_multi& operator = (const layer_multi&);

    /**
     * @brief Output stream operator for debugging
     * @param os Output stream
     * @param lm layer_multi object to output
     * @return Reference to output stream
     */
    friend std::ostream& operator << (std::ostream&, const layer_multi&);

};


/** @} */ // end of homogenization group

} //namespace simcoon
