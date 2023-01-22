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

/**
 * @file recovery_props.hpp
 * @author Yves Chemisky 
 * @section A set of function that allow to get the Elastic properties from stiffness/compliance tensors
*/

#pragma once
#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @brief Check the symmetries of a stiffness matrix \f$ \mathbf{L} \f$, and fill the vector of material properties. 
 * 
 * Depending on the symmetry found, the string umat_type, the axis of symmetry (if applicable) the vector of material properties, and the  major symmetry maj_sym \f$ L_{ij} = L_{ji} \f$ ?.
 * If the major symmetry condition is not fulfilled, the check of symmetries if performed on the symmetric part of \f$ \mathbf{L} \f$. For fully anisotropic and monoclinic symmetries, the vector of parameters is not returned (the full stiffness tensor is generally directly utilized)
 *
 * @param[in] L : The stiffness matrix
 * @param[out] umat_type : The material symmetry type
 * @param[out] axis : The axis of symmetry (if applicable)
 * @param[out] props : The material properties vector
 * @param[out] maj_sym : The major symmetry condition (L_ij = L_ji ?).
 *
 * Material Symmetries considered:
 *
 *  | Symmetry        | umat_type   | axis       | size of props |
 *  |------------------|------------|------------|---------------|
 *  | Fully anisotropic | ELANI  | 0          | N/A           |
 *  | Monoclinic        | ELMON  | 1,2 or 3   | N/A           |
 *  | Orthotropic       | ELORT  | 0          | 9             |
 *  | Cubic             | ELCUB  | 0          | 3             |
 *  | Transversely isotropic | ELITR | 1,2 or 3 | 5             |
 *  | Isotropic         | ELISO  | 0          | 2             |
 *
 * @details Example: 
 *  @code
 *     mat L;
 *     std::string umat_name;
 *     int axis;
 *     vec props;
 *     int maj_sym;
 *     check_symetries(L, umat_name, axis, props, maj_sym);
 *  @endcode
 */
void check_symetries(const arma::mat &L, std::string &umat_type, int &axis, arma::vec &props, int &maj_sym);

/**
 * @brief Provides material parameters of a linear isotropic material from a stiffness matrix
 *
 * Returns a vector containing the Young's modulus and the Poisson ratio :\f$ \left(E, \nu \right) \f$  of a linear elastic isotropic material, 
 * providing the stiffness matrix :\f$ \mathbf{L} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param L (6x6 arma::mat) A stiffness tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E, \nu \right) \f$
 * 
 * @details Example: 
 * @code
 *     mat L = L_iso(70000., 0.3, 'Enu');
 *     vec eliso_props = L_iso_props(L);
 * @endcode
*/
arma::vec L_iso_props(const arma::mat &L);

/**
 * @brief Provides material parameters of a linear isotropic material from a compliance matrix
 *
 * Returns a vector containing the Young's modulus and the Poisson ratio :\f$ \left(E, \nu \right) \f$  of a linear elastic isotropic material, 
 * providing the compliance matrix :\f$ \mathbf{M} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param M (6x6 arma::mat) A compliance tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E, \nu \right) \f$
 * 
 * @details Example: 
 * @code
 *     mat M = M_iso(70000., 0.3, 'Enu');
 *     vec eliso_props = L_iso_props(L);
 * @endcode
*/
arma::vec M_iso_props(const arma::mat &);
    
/**
 * @brief Provides material parameters of a linear transversely isotropic material from a stiffness matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right) \f$ of a linear elastic transversely isotropic material, 
 * providing the stiffness matrix :\f$ \mathbf{L} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param L (6x6 arma::mat) A stiffness tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right) \f$
 * 
 * @details Example: 
 * @code
    int axis = 1;
    double E_L = 4500;
    double E_T = 2300;
    double nu_TL = 0.05;
    double nu_TT = 0.3;
    double G_LT = 2700;
    mat L = L_isotrans(E_L, E_T, nu_TL, nu_TT, G_LT., axis);
    vec isotrans_props = L_isotrans_props(L, axis);
 * @endcode
*/
arma::vec L_isotrans_props(const arma::mat &, const int &);

/**
 * @brief Provides material parameters of a linear transversely isotropic material from a compliance matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right) \f$ of a linear elastic transversely isotropic material, 
 * providing the compliance matrix :\f$ \mathbf{M} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param M (6x6 arma::mat) A compliance tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right) \f$
 * 
 * @details Example: 
 * @code
    int axis = 1;
    double E_L = 4500;
    double E_T = 2300;
    double nu_TL = 0.05;
    double nu_TT = 0.3;
    double G_LT = 2700;
    mat M = M_isotrans(E_L, E_T, nu_TL, nu_TT, G_LT., axis);
    vec isotrans_props = M_isotrans_props(M, axis);
 * @endcode
*/
arma::vec M_isotrans_props(const arma::mat &, const int &);

//This function recovers the properties of a cubic isotropic stiffness tensor
arma::vec L_cubic_props(const arma::mat &);

//This function recovers the properties of a cubic isotropic compliance tensor
arma::vec M_cubic_props(const arma::mat &);
    
//This function recovers the properties of an orthotropic compliance tensor
arma::vec L_ortho_props(const arma::mat &);
    
//This function recovers the properties of an orthotropic compliance tensor
arma::vec M_ortho_props(const arma::mat &);
    
//This function recovers the properties of a fully anisotropic compliance tensor
arma::vec M_aniso_props(const arma::mat &);

} //namespace simcoon
