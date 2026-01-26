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

///@file cylinder_multi.hpp
///@brief characteristics of an cylinder for homogenization purposes
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

/**
 * @file cylinder_multi.hpp
 * @brief Cylindrical inclusion for homogenization
 * @author Yves Chemisky
 * @version 1.0
 *
 * This class extends phase_multi for cylindrical (fiber) inclusions.
 * It is typically used for fiber-reinforced composites where fibers
 * can be modeled as infinite cylinders.
 */

/** @addtogroup homogenization
 *  @{
 */

/**
 * @brief Cylindrical inclusion class for fiber-reinforced composite homogenization
 *
 * @details This class implements the mechanics of cylindrical (fiber) inclusions
 * for mean-field homogenization of fiber-reinforced composites. Cylindrical
 * inclusions are a limiting case of ellipsoids with one infinite axis.
 *
 * **Applications:**
 * - Unidirectional fiber composites
 * - Long fiber reinforced polymers
 * - Carbon fiber composites
 *
 * **Coordinate System:**
 * - The cylinder axis is aligned with the local x₃ direction
 * - The cross-section is circular in the x₁-x₂ plane
 *
 * **Concentration Tensors:**
 *
 * The strain concentration relates macroscopic to local strain:
 * \f[
 * \boldsymbol{\varepsilon}^{fiber} = \mathbf{A} : \bar{\boldsymbol{\varepsilon}}
 * \f]
 *
 * Local and global concentration tensors are related by rotation:
 * \f[
 * \mathbf{A} = \mathbf{R}^T : \mathbf{A}_{loc} : \mathbf{R}
 * \f]
 *
 * @see phase_multi for base class documentation
 * @see ellipsoid_multi for general ellipsoidal inclusions
 */
class cylinder_multi : public phase_multi
{
private:

protected:

	public :

    /**
     * @brief Interaction tensor in local coordinates (6×6)
     *
     * Computed for cylinder geometry in the local frame where
     * the cylinder axis is aligned with x₃.
     */
    arma::mat T_loc;

    /**
     * @brief Interaction tensor in global coordinates (6×6)
     *
     * Rotated from local to global using cylinder orientation
     */
    arma::mat T;

    /**
     * @brief Strain concentration tensor in local coordinates (6×6)
     *
     * Local frame concentration tensor before rotation to global frame
     */
    arma::mat A_loc;

    /**
     * @brief Stress concentration tensor in local coordinates (6×6)
     *
     * Local frame concentration tensor before rotation to global frame
     */
    arma::mat B_loc;

    /** @brief Default constructor */
    cylinder_multi();

    /**
     * @brief Full constructor with all parameters
     *
     * @param A Strain concentration tensor (global)
     * @param A_start Initial strain concentration tensor
     * @param B Stress concentration tensor (global)
     * @param B_start Initial stress concentration tensor
     * @param A_in Inelastic concentration vector
     * @param T_loc Local interaction tensor
     * @param T Global interaction tensor
     * @param A_loc Local strain concentration tensor
     * @param B_loc Local stress concentration tensor
     */
    cylinder_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&);

    /**
     * @brief Copy constructor
     * @param cm cylinder_multi object to copy
     */
    cylinder_multi(const cylinder_multi&);

    /** @brief Destructor */
    ~cylinder_multi();

    /**
     * @brief Assignment operator
     * @param cm cylinder_multi object to assign from
     * @return Reference to this object
     */
	virtual cylinder_multi& operator = (const cylinder_multi&);

    /**
     * @brief Output stream operator for debugging
     * @param os Output stream
     * @param cm cylinder_multi object to output
     * @return Reference to output stream
     */
    friend std::ostream& operator << (std::ostream&, const cylinder_multi&);

};


/** @} */ // end of homogenization group

} //namespace simcoon
