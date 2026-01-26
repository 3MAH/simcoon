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

///@file phase_multi.hpp
///@brief Characteristics of a phase for multiphase simulations, the parent class of:
// - ellipsoid_multi
// - layer_multi
// - cylinder_multi
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>

namespace simcoon{

/**
 * @file phase_multi.hpp
 * @brief Base class for multiphase homogenization
 * @author Yves Chemisky
 * @version 1.0
 *
 * This file defines the base class for phases in mean-field homogenization schemes.
 * The phase_multi class stores concentration tensors that relate macroscopic fields
 * to local (phase-level) fields.
 */

/** @addtogroup homogenization
 *  @{
 */

/**
 * @brief Base class for phases in multiphase homogenization schemes
 *
 * @details The phase_multi class serves as the base class for all phase types
 * used in mean-field homogenization methods. It stores the concentration tensors
 * that relate macroscopic strain/stress to local phase-level quantities.
 *
 * **Concentration Tensor Relationships:**
 *
 * The strain concentration tensor \f$ \mathbf{A} \f$ relates macroscopic strain to local strain:
 * \f[
 * \boldsymbol{\varepsilon}^{(r)} = \mathbf{A}^{(r)} : \bar{\boldsymbol{\varepsilon}}
 * \f]
 *
 * The stress concentration tensor \f$ \mathbf{B} \f$ relates macroscopic stress to local stress:
 * \f[
 * \boldsymbol{\sigma}^{(r)} = \mathbf{B}^{(r)} : \bar{\boldsymbol{\sigma}}
 * \f]
 *
 * For inelastic behavior, an additional concentration vector handles eigenstrain contributions:
 * \f[
 * \boldsymbol{\varepsilon}^{(r)} = \mathbf{A}^{(r)} : \bar{\boldsymbol{\varepsilon}} + \mathbf{A}_{in}^{(r)}
 * \f]
 *
 * **Derived Classes:**
 * - ellipsoid_multi: Ellipsoidal inclusions (Eshelby-based methods)
 * - cylinder_multi: Cylindrical inclusions
 * - layer_multi: Layered composites
 *
 * @see ellipsoid_multi, cylinder_multi, layer_multi
 */
class phase_multi
{
	private:

	protected:

	public :

    /**
     * @brief Strain concentration tensor (6×6 matrix in Voigt notation)
     *
     * Relates macroscopic strain to local phase strain:
     * \f$ \boldsymbol{\varepsilon}^{local} = \mathbf{A} : \bar{\boldsymbol{\varepsilon}} \f$
     */
    arma::mat A;

    /**
     * @brief Strain concentration tensor at the start of increment
     *
     * Stored for incremental schemes and convergence checks
     */
    arma::mat A_start;

    /**
     * @brief Stress concentration tensor (6×6 matrix in Voigt notation)
     *
     * Relates macroscopic stress to local phase stress:
     * \f$ \boldsymbol{\sigma}^{local} = \mathbf{B} : \bar{\boldsymbol{\sigma}} \f$
     */
    arma::mat B;

    /**
     * @brief Stress concentration tensor at the start of increment
     *
     * Stored for incremental schemes and convergence checks
     */
    arma::mat B_start;

    /**
     * @brief Inelastic strain concentration vector (6×1 in Voigt notation)
     *
     * Accounts for eigenstrain/transformation strain contributions:
     * \f$ \boldsymbol{\varepsilon}^{local} = \mathbf{A} : \bar{\boldsymbol{\varepsilon}} + \mathbf{A}_{in} \f$
     */
    arma::vec A_in;

    /**
     * @brief Default constructor
     *
     * Initializes concentration tensors to zero matrices/vectors of appropriate size
     */
    phase_multi();

    /**
     * @brief Constructor with parameters
     *
     * @param A_init Initial strain concentration tensor (6×6)
     * @param A_start_init Initial strain concentration tensor for start state (6×6)
     * @param B_init Initial stress concentration tensor (6×6)
     * @param B_start_init Initial stress concentration tensor for start state (6×6)
     * @param A_in_init Initial inelastic concentration vector (6×1)
     */
    phase_multi(const arma::mat& A_init, const arma::mat& A_start_init,
                const arma::mat& B_init, const arma::mat& B_start_init,
                const arma::vec& A_in_init);

    /**
     * @brief Copy constructor
     * @param pm phase_multi object to copy
     */
    phase_multi(const phase_multi& pm);

    /**
     * @brief Virtual destructor
     */
    virtual ~phase_multi();

    /**
     * @brief Reset current tensors to start-of-increment values
     *
     * Sets A = A_start and B = B_start for iterative schemes
     */
    void to_start();

    /**
     * @brief Store current tensors as start-of-increment values
     *
     * Sets A_start = A and B_start = B after convergence
     */
    void set_start();

    /**
     * @brief Assignment operator
     * @param pm phase_multi object to assign from
     * @return Reference to this object
     */
	virtual phase_multi& operator = (const phase_multi& pm);

    /**
     * @brief Output stream operator for debugging
     * @param os Output stream
     * @param pm phase_multi object to output
     * @return Reference to output stream
     */
    friend std::ostream& operator << (std::ostream& os, const phase_multi& pm);
};


/** @} */ // end of homogenization group

} //namespace simcoon
