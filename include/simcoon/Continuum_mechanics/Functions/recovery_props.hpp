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
 * @param[in] tol : The tolerance utilized to check the symetries. If less than the global sim_limit (1.E-8), sim_limit is utilized. Default is 0. (so sim_limit is utilized by default)
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
void check_symetries(const arma::mat &L, std::string &umat_type, int &axis, arma::vec &props, int &maj_sym, const double &tol = 0.);

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
arma::vec M_iso_props(const arma::mat &M);
    
/**
 * @brief Provides material parameters of a linear transversely isotropic material from a stiffness matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right) \f$ of a linear elastic transversely isotropic material, 
 * providing the stiffness matrix :\f$ \mathbf{L} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param L (6x6 arma::mat) A stiffness tensor
 * @param axis axis of symetry for the linear transversely isotropic material
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
arma::vec L_isotrans_props(const arma::mat &L, const int &axis);

/**
 * @brief Provides material parameters of a linear transversely isotropic material from a compliance matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \right) \f$ of a linear elastic transversely isotropic material, 
 * providing the compliance matrix :\f$ \mathbf{M} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param M (6x6 arma::mat) A compliance tensor
 * @param axis axis of symetry for the linear transversely isotropic material
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
arma::vec M_isotrans_props(const arma::mat &M, const int &axis);

/**
 * @brief Provides material parameters of a cubic material from a stiffness matrix
 *
 * Returns a vector containing the parameters \f$ \left(E, \nu, G \right) \f$ of a linear cubic material, 
 * providing the stiffness matrix :\f$ \mathbf{L} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param L (6x6 arma::mat) A stiffness tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E, \nu, G \right) \f$
 * 
 * @details Example: 
 * @code
    mat L = L_cubic(185000., 158000., 39700., 'Cii') //C11, C12, C44
    vec cubic_props = L_cubic_props(L);
 * @endcode
*/
arma::vec L_cubic_props(const arma::mat &L);

/**
 * @brief Provides material parameters of a cubic material from a compliance matrix
 *
 * Returns a vector containing the parameters \f$ \left(E, \nu, G \right) \f$ of a linear cubic material, 
 * providing the compliance matrix :\f$ \mathbf{M} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param M (6x6 arma::mat) A compliance tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E, \nu, G \right) \f$
 * 
 * @details Example: 
 * @code
    mat M = M_cubic(185000., 158000., 39700., 'Cii') //C11, C12, C44
    vec cubic_props = M_cubic_props(M);
 * @endcode
*/
arma::vec M_cubic_props(const arma::mat &M);
    
/**
 * @brief Provides material parameters of an orthotropic material from a stiffness matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23} \right) \f$ of a linear orthotropic material, 
 * providing the stiffness matrix :\f$ \mathbf{L} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param L (6x6 arma::mat) A stiffness tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23} \right) \f$
 * 
 * @details Example: 
 * @code
    double E_1 = 4500;
    double E_2 = 2300;
    double E_3 = 2700;
    double nu_12 = 0.06;
    double nu_13 = 0.08;
    double nu_23 = 0.3;
    double G_12 = 2200;
    double G_13 = 2100;
    double G_23 = 2400;
    mat L = L_ortho(E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23);
    vec ortho_props = L_ortho_props(L);
 * @endcode
*/
arma::vec L_ortho_props(const arma::mat &L);
    
/**
 * @brief Provides material parameters of an orthotropic material from a compliance matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23} \right) \f$ of a linear orthotropic material, 
 * providing the compliance matrix :\f$ \mathbf{M} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param M (6x6 arma::mat) A compliance tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23} \right) \f$
 * 
 * @details Example: 
 * @code
    double E_1 = 4500;
    double E_2 = 2300;
    double E_3 = 2700;
    double nu_12 = 0.06;
    double nu_13 = 0.08;
    double nu_23 = 0.3;
    double G_12 = 2200;
    double G_13 = 2100;
    double G_23 = 2400;
    mat M = M_ortho(E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23);
    vec ortho_props = M_ortho_props(M);
 * @endcode
*/
arma::vec M_ortho_props(const arma::mat &M);
    
/**
 * @brief Provides material parameters of an anisotropic material from a compliance matrix
 *
 * Returns a vector containing the parameters \f$ \left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23}, \eta_{14}, \eta_{15}, \eta_{16}, \eta_{24}, \eta_{25}, \eta_{26}, \eta_{34}, \eta_{35}, \eta_{36}, \eta_{45}, \eta_{46}, \eta_{56} \right) \f$ of a linear anisotropic material, 
 * providing the compliance matrix :\f$ \mathbf{M} \f$. Note that an averaging over the component is operated (usefull when the provided matrix 
 * do not exactly correspond to an isotropic material)
 * 
 * @param M (6x6 arma::mat) A compliance tensor
 * @return (arma::vec) the vector of parameters \f$ \left(E_1, E_1, E_3, \nu_{12} \nu_{13}, \nu_{23}, G_{12}, G_{13}, G_{23}, \eta_{14}, \eta_{15}, \eta_{16}, \eta_{24}, \eta_{25}, \eta_{26}, \eta_{34}, \eta_{35}, \eta_{36}, \eta_{45}, \eta_{46}, \eta_{56} \right) \f$
 * 
 * \f{eqnarray*}{
	   E_1 & = & \frac{1}{M_{11}} \\
	   E_2 & = & \frac{1}{M_{22}} \\
	   E_3 & = & \frac{1}{M_{33}} \\
	
      \nu_{12} & = & - \frac{1}{2} E_1 \left( M_{12} + M_{21} \right) \\
      \nu_{13} & = & - \frac{1}{2} E_1 \left( M_{13} + M_{31} \right) \\
      \nu_{23} & = & - \frac{1}{2} E_2 \left( M_{23} + M_{32} \right) \\
    
      G_{12} & = & \frac{1}{M_{44}} \\ 
      G_{13} & = & \frac{1}{M_{55}} \\ 
      G_{23} & = & \frac{1}{M_{66}} \\ 

      \eta_{14} & = & \frac{1}{2} E_1 \left( M_{14} + M_{41} \right) \\
      \eta_{25} & = & \frac{1}{2} E_1 \left( M_{15} + M_{51} \right) \\
      \eta_{26} & = & \frac{1}{2} E_1 \left( M_{16} + M_{61} \right) \\

      \eta_{24} & = & \frac{1}{2} E_2 \left( M_{24} + M_{42} \right) \\
      \eta_{25} & = & \frac{1}{2} E_2 \left( M_{25} + M_{52} \right) \\
      \eta_{26} & = & \frac{1}{2} E_2 \left( M_{26} + M_{62} \right) \\

      \eta_{34} & = & \frac{1}{2} E_3 \left( M_{34} + M_{43} \right) \\
      \eta_{35} & = & \frac{1}{2} E_3 \left( M_{35} + M_{53} \right) \\
      \eta_{36} & = & \frac{1}{2} E_3 \left( M_{36} + M_{63} \right) \\

      \eta_{45} & = & \frac{1}{2} G_{12} \left( M_{45} + M_{54} \right) \\
      \eta_{46} & = & \frac{1}{2} G_{12} \left( M_{46} + M_{64} \right) \\
      \eta_{56} & = & \frac{1}{2} G_{12} \left( M_{56} + M_{65} \right)

 * \f}
 * 
 * @details Example: 
 * @code
    mat A(6, 6, fill::randu);
    mat M = A.t() * A;
    vec aniso_props = M_aniso_props(M);
 * @endcode
*/
arma::vec M_aniso_props(const arma::mat &M);

} //namespace simcoon
