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

///@file ellipsoid_multi.hpp
///@brief characteristics of an ellipsoid for homogenization purposes
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

/**
 * @file ellipsoid_multi.hpp
 * @brief Ellipsoidal inclusion for Eshelby-based homogenization
 * @author Yves Chemisky
 * @version 1.0
 *
 * This class extends phase_multi for ellipsoidal inclusions using Eshelby's
 * equivalent inclusion method. It computes the Eshelby tensor, Hill's polarization
 * tensor, and interaction tensors for mean-field homogenization schemes.
 */

/** @addtogroup homogenization
 *  @{
 */

/**
 * @brief Ellipsoidal inclusion class for Eshelby-based homogenization
 *
 * @details This class implements the mechanics of ellipsoidal inclusions embedded
 * in an infinite matrix, based on Eshelby's classical solution (1957). It is used
 * in mean-field homogenization methods such as:
 * - Mori-Tanaka scheme
 * - Self-consistent scheme
 * - Differential scheme
 *
 * **Key Tensors:**
 *
 * The Eshelby tensor \f$ \mathbf{S} \f$ relates the constrained strain in an inclusion
 * to its eigenstrain:
 * \f[
 * \boldsymbol{\varepsilon}^c = \mathbf{S} : \boldsymbol{\varepsilon}^*
 * \f]
 *
 * Hill's polarization tensor \f$ \mathbf{P} \f$ is related to the Eshelby tensor by:
 * \f[
 * \mathbf{P} = \mathbf{S} : \mathbf{L}_0^{-1}
 * \f]
 *
 * The interaction tensor \f$ \mathbf{T} \f$ (Wu's tensor) relates inclusion stress
 * to the strain difference:
 * \f[
 * \mathbf{T} = \mathbf{L}_0 : (\mathbf{I} - \mathbf{S})
 * \f]
 *
 * @note All tensors are computed in the local coordinate system of the ellipsoid
 *       and can be rotated to global coordinates using the ellipsoid orientation.
 *
 * @see phase_multi for base class documentation
 * @see Eshelby, J.D. (1957). The determination of the elastic field of an ellipsoidal
 *      inclusion, and related problems. Proc. R. Soc. Lond. A 241, 376-396.
 */
class ellipsoid_multi : public phase_multi
{
private:

protected:

	public :

    /**
     * @brief Eshelby tensor in local coordinates (6×6 in Voigt notation)
     *
     * The Eshelby tensor depends only on the ellipsoid aspect ratios and
     * the matrix Poisson ratio (for isotropic matrix).
     */
    arma::mat S_loc;

    /**
     * @brief Hill's polarization tensor in local coordinates (6×6 in Voigt notation)
     *
     * Related to Eshelby tensor: \f$ \mathbf{P} = \mathbf{S} : \mathbf{L}_0^{-1} \f$
     */
    arma::mat P_loc;

    /**
     * @brief Interaction tensor (Wu's tensor) in local coordinates (6×6)
     *
     * \f$ \mathbf{T}_{loc} = \mathbf{L}_0 : (\mathbf{I} - \mathbf{S}) \f$
     */
    arma::mat T_loc;

    /**
     * @brief Interaction tensor in global coordinates (6×6)
     *
     * Rotated from local to global using ellipsoid orientation
     */
    arma::mat T;

    /**
     * @brief Inelastic interaction tensor in local coordinates (6×6)
     *
     * Used for eigenstrain/transformation strain contributions
     */
    arma::mat T_in_loc;

    /**
     * @brief Inelastic interaction tensor in global coordinates (6×6)
     */
    arma::mat T_in;

    /** @brief Number of integration points in polar direction for numerical integration */
    static int mp;

    /** @brief Number of integration points in azimuthal direction */
    static int np;

    /** @brief Gauss points for polar integration */
    static arma::vec x;

    /** @brief Gauss weights for polar integration */
    static arma::vec wx;

    /** @brief Gauss points for azimuthal integration */
    static arma::vec y;

    /** @brief Gauss weights for azimuthal integration */
    static arma::vec wy;

    /** @brief Default constructor */
    ellipsoid_multi();

    /**
     * @brief Full constructor with all parameters
     *
     * @param A Strain concentration tensor
     * @param A_start Initial strain concentration tensor
     * @param B Stress concentration tensor
     * @param B_start Initial stress concentration tensor
     * @param A_in Inelastic concentration vector
     * @param S_loc Local Eshelby tensor
     * @param P_loc Local polarization tensor
     * @param T_loc Local interaction tensor
     * @param T Global interaction tensor
     * @param T_in_loc Local inelastic interaction tensor
     * @param T_in Global inelastic interaction tensor
     */
    ellipsoid_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&);

    /**
     * @brief Copy constructor
     * @param em ellipsoid_multi object to copy
     */
    ellipsoid_multi(const ellipsoid_multi&);

    /** @brief Destructor */
    ~ellipsoid_multi();

    /**
     * @brief Compute the Eshelby tensor in local coordinates
     *
     * @details Computes \f$ \mathbf{S} \f$ based on the matrix stiffness and
     * ellipsoid geometry using numerical integration over the unit sphere.
     *
     * @param L0 Matrix (reference medium) stiffness tensor (6×6)
     * @param ell Ellipsoid geometry (aspect ratios and orientation)
     */
    virtual void fillS_loc(const arma::mat&, const ellipsoid &);

    /**
     * @brief Compute Hill's polarization tensor in local coordinates
     *
     * @details Computes \f$ \mathbf{P} = \mathbf{S} : \mathbf{L}_0^{-1} \f$
     *
     * @param L0 Matrix stiffness tensor (6×6)
     * @param ell Ellipsoid geometry
     */
    virtual void fillP_loc(const arma::mat&, const ellipsoid &);

    /**
     * @brief Compute the interaction tensor
     *
     * @details Computes the interaction tensor \f$ \mathbf{T} \f$ that relates
     * the stress in the inclusion to the strain mismatch with the matrix.
     *
     * @param L0 Matrix stiffness tensor (6×6)
     * @param Lt Inclusion tangent stiffness tensor (6×6)
     * @param ell Ellipsoid geometry
     */
    virtual void fillT(const arma::mat&, const arma::mat&, const ellipsoid &);

    /**
     * @brief Compute the interaction tensor for isotropic matrix
     *
     * @details Optimized computation when the matrix is isotropic, using
     * analytical expressions for the Eshelby tensor.
     *
     * @param L0 Isotropic matrix stiffness tensor (6×6)
     * @param Lt Inclusion tangent stiffness tensor (6×6)
     * @param ell Ellipsoid geometry
     */
    virtual void fillT_iso(const arma::mat&, const arma::mat&, const ellipsoid &);

    /**
     * @brief Compute the inelastic interaction tensor
     *
     * @details Computes the tensor relating eigenstrain/transformation strain
     * to the resulting stress perturbation.
     *
     * @param L0 Matrix stiffness tensor (6×6)
     * @param Lt Inclusion tangent stiffness tensor (6×6)
     * @param ell Ellipsoid geometry
     */
    virtual void fillT_mec_in(const arma::mat&, const arma::mat&, const ellipsoid &);

    /**
     * @brief Assignment operator
     * @param em ellipsoid_multi object to assign from
     * @return Reference to this object
     */
	virtual ellipsoid_multi& operator = (const ellipsoid_multi&);

    /**
     * @brief Output stream operator for debugging
     * @param os Output stream
     * @param em ellipsoid_multi object to output
     * @return Reference to output stream
     */
    friend std::ostream& operator << (std::ostream&, const ellipsoid_multi&);

};


/** @} */ // end of homogenization group

} //namespace simcoon
