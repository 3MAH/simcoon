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
#include <string.h>
#include <armadillo>

namespace simcoon{

/**
* @file constitutive.hpp
* @author Yves Chemisky 
* @section Constitutive library contains all the function helpful to write a constitutive relation of constitutive model
*/

/**
 * @brief Returns the fourth order identity tensor Ireal written in Voigt notation
 * @param None
 * @return The following 6x6 mat (arma::mat)
\f[
    I_{real} = \left( \begin{array}{cccccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)
\f]
 * @details Example: 
 * @code 
        mat Ir = Ireal();
 * @endcode
*/
arma::mat Ireal();

/**
 * @brief Returns the volumic of the identity tensor Ivol written in Voigt notation
 * @param None
 * @return The following 6x6 mat (arma::mat) 
\f[
    I_{vol} = \left( \begin{array}{cccccc}
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \end{array} \right)
\f]
 * @details Example: 
 * @code 
        mat Iv = Ivol();
 * @endcode
*/
arma::mat Ivol();

/**
 * @brief Returns the deviatoric of the identity tensor Idev written in Voigt notation
 * @param None
 * @return The following 6x6 mat (arma::mat)
\f[ 
    I_{dev} = I_{real} - I_{vol} = \left( \begin{array}{cccccc}
        2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)
\f]
 * @details Example: 
 * @code 
        mat Id = Idev();
 * @endcode
 */
arma::mat Idev();

/**
 * @brief Returns the fourth order identity tensor \f$ \widehat{I} \f$ written in Voigt notation
 * @param None
 * @return The following 6x6 mat (arma::mat)
\f[ 
    \widehat{I} = \left( \begin{array}{cccccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)
\f]
 * @details For example, this tensor allows to obtain : \f$ L*\widehat{M}=I \f$ or \f$ \widehat{L}*M=I \f$,
 * where a matrix \f$ \widehat{A} \f$ is set by \f$ \widehat{A}=\widehat{I}\,A\,\widehat{I} \f$
 * @code 
        mat Ir2 = Ireal2();
 * @endcode
 */  
 arma::mat Ireal2();

/**
 * @brief Returns the deviatoric part of the identity tensor, in the form of \f$ \widehat{I} \f$ considering the Voigt notation
 * @param None
 * @return The following 6x6 mat (arma::mat) 
\f[ 
    I_{dev2} = \left( \begin{array}{cccccc}
        2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
        -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)
\f]
 * @details Example: 
 * @code 
 *  mat Id2 = Idev2();
 * @endcode
 */
arma::mat Idev2();

/**
 * @brief Returns the expansion vector
 * @param None
 * @return The following 6 vec (arma::vec) 
\f[ 
    I_{th} = \left( \begin{array}{c}
        1 \\
        1 \\
        1 \\
        0 \\
        0 \\
        0 \end{array} \right)
\f]
 * @details Example: 
 * @code 
        vec It = Ith();
 * @endcode
*/
arma::vec Ith();

/**
 * @brief Returns the operator from a stress to strain Voigt convention
 * @param None
 * @return The following 6 vec (arma::vec) 
\f[ 
    I_{r2} = \left( \begin{array}{ccc}
        1 \\
        1 \\
        1 \\
        2 \\
        2 \\
        2 \end{array} \right)
\f]
 * @details Example: 
 * @code 
        vec I2 = Ir2();
 * @endcode
*/
arma::vec Ir2();

/**
 * @brief Returns the operator from a strain to stress Voigt convention
 * @param None
 * @return The following 6 vec (arma::vec)
\f[ 
    I_{r05} = \left( \begin{array}{ccc}
        1 \\
        1 \\
        1 \\
        0.5 \\
        0.5 \\
        0.5 \end{array} \right)
\f] 
 * @details Example: 
 * @code 
        vec I05 = Ir05();
 * @endcode
*/
arma::vec Ir05();

/**
 * @brief Provides the elastic stiffness tensor for an isotropic material.
 * @details The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
 * Exhaustive list of possible third argument :
 * ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
 * @param C1, C2, conv
 * @return The 6x6 stiffness matrix considering a Voigt notation (arma::mat)
 * @details Example: 
 * @code 
        double E = 210000;
        double nu = 0.3;
        mat Liso = L_iso(E, nu, "Enu");
 * @endcode
 */
arma::mat L_iso(const double &C1, const double &C2, const std::string &conv="Enu");

/**
 * @brief Provides the elastic compliance tensor for an isotropic material.
 * @param C1, C2, conv
 * @return The 6x6 compliance matrix considering a Voigt notation (arma::mat) 
 * @details The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
 * Exhaustive list of possible third argument :
 * ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
 * @code 
        double E = 210000;
        double nu = 0.3;
        mat Miso = M_iso(E, nu, "Enu");
 * @endcode
*/
arma::mat M_iso(const double &C1, const double &C2, const std::string& conv="Enu");

/**
 * @brief Provides the elastic stiffness tensor for a cubic material.
 * @param C1, C2, C3, conv
 * @return The 6x6 stiffness matrix considering a Voigt notation (arma::mat) 
 * @details Arguments are the stiffness coefficients C11, C12 and C44 or E, nu and G. 
 * The fourth argument specify which couple has been provided and the order of coefficients: ‘EnuG’ or ’Cii’.
 * @code 
        double C11 = alead(10000., 100000.);
        double C12 = alead(10000., 100000.);
        double C44 = alead(10000., 100000.);
        mat Lcubic = L_cubic(C11, C12, C44, "Cii");
 * @endcode
*/
arma::mat L_cubic(const double &C1, const double &C2, const double &C3, const std::string& conv="EnuG");

/**
 * @brief Provides the elastic compliance tensor for a cubic material.
 * @param C1, C2, C3, conv
 * @return The 6x6 compliance matrix considering a Voigt notation (arma::mat)
 * @details Arguments are the stiffness coefficients C11, C12 and C44 or E, nu and G. 
 * The fourth argument specify which couple has been provided and the order of coefficients: ‘EnuG’ or ’Cii’.
 * @code 
        double C11 = alead(10000., 100000.);
        double C12 = alead(10000., 100000.);
        double C44 = alead(10000., 100000.);
        mat Mcubic = M_cubic(C11,C12,C44);
 * @endcode 
*/
arma::mat M_cubic(const double &C1, const double &C2, const double &C3, const std::string& = "EnuG");

/**
 * @brief Provides the elastic stiffness tensor for an orthotropic material.
 * @param C11, C12, C22, C23, C33, C44, C55, C66, conv
 * @return The 6x6 stiffness matrix considering a Voigt notation (arma::mat)
 * @details Arguments could be all the stiffness coefficients (\f$ C_{ij} \f$) or the material parameter. In this last case for an orthotropic material the material parameters should be : \f$ E_x, E_y, E_z, \nu_{xy}, \nu_{yz}, \nu_{xz}, G_{xy}, G_{yz}, G_{xz} \f$.
 * The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.
 * @code 
        double C11 = alead(10000., 100000.);
        double C12 = alead(10000., 100000.);
        double C13 = alead(10000., 100000.);
        double C22 = alead(10000., 100000.);
        double C23 = alead(10000., 100000.);
        double C33 = alead(10000., 100000.);
        double C44 = alead(10000., 100000.);
        double C55 = alead(10000., 100000.);
        double C66 = alead(10000., 100000.);
        mat Lortho = L_ortho(C11, C12, C13, C22, C23, C33, C44, C55, C66,"Cii");
 * @endcode 
 */
arma::mat L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const std::string& conv="EnuG");

/**
 * @brief Provides the elastic compliance tensor for an orthotropic material.
 * @param C11, C12, C22, C23, C33, C44, C55, C66, conv
 * @return The 6x6 compliance matrix considering a Voigt notation (arma::mat)
 * @details Arguments could be all the stiffness coefficients (\f$ C_{ij} \f$) or the material parameter. In this last case for an orthotropic material the material parameters should be : \f$ E_x, E_y, E_z, \nu_{xy}, \nu_{yz}, \nu_{xz}, G_{xy}, G_{yz}, G_{xz} \f$.
 * The last argument must be set to “Cii” if the inputs are the stiffness coefficients or to “EnuG” if the inputs are the material parameters.
 * @code 
       double C11 = alead(10000., 100000.);
       double C12 = alead(10000., 100000.);
       double C13 = alead(10000., 100000.);
       double C22 = alead(10000., 100000.);
       double C23 = alead(10000., 100000.);
       double C33 = alead(10000., 100000.);
       double C44 = alead(10000., 100000.);
       double C55 = alead(10000., 100000.);
       double C66 = alead(10000., 100000.);
       mat Mortho = M_ortho(C11, C12, C13, C22, C23, C33, C44, C55, C66,"Cii");
 * @endcode 
 */
arma::mat M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const std::string& conv="EnuG");

/**
 * @brief Provides the elastic stiffness tensor for an isotropic transverse material.
 * @param EL, ET, nuTL, nuTT, GLT, axis
 * @return The 6x6 stiffness matrix considering a Voigt notation (arma::mat)
 * @details Arguments are longitudinal Young modulus \f$ E_L \f$ (EL), transverse young modulus \f$ E_T \f$ (ET), Poisson’s ratio for loading along the longitudinal axis \f$ nu_{TL} \f$ (nuTL), Poisson’s ratio for loading along the transverse axis \f$ \nu_{TT} \f$ (nuTT), shear modulus \f$ G_{LT} \f$ (GLT) and the axis of symmetry.
 * @code 
        double EL = alead(10000., 100000.);
        double ET = alead(10000., 100000.);
        double nuTL = alead(0., 0.5);
        double nuTT = alead(0.5, 0.5);
        double GLT = alead(10000., 100000.);
        double axis = 1;
        mat Lisotrans = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);
 * @endcode 
 */
arma::mat L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis);

/**
 * @brief Provides the elastic compliance tensor for an isotropic transverse material.
 * @param EL, ET, nuTL, nuTT, GLT, axis
 * @return The 6x6 compliance matrix considering a Voigt notation (arma::mat)
 * @details Arguments are longitudinal Young modulus \f$ E_L \f$ (EL), transverse young modulus \f$ E_T \f$ (ET), Poisson’s ratio for loading along the longitudinal axis \f$ nu_{TL} \f$ (nuTL), Poisson’s ratio for loading along the transverse axis \f$ \nu_{TT} \f$ (nuTT), shear modulus \f$ G_{LT} \f$ (GLT) and the axis of symmetry.
 * @code 
        double EL = alead(10000., 100000.);
        double ET = alead(10000., 100000.);
        double nuTL = alead(0., 0.5);
        double nuTT = alead(0.5, 0.5);
        double GLT = alead(10000., 100000.);
        double axis = 1;
        mat Misotrans = M_isotrans(EL, ET, nuTL, nuTT, GLT, axis);
 * @endcode 
 */
arma::mat M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis);

/**
 * @brief Provides the viscous tensor H an isotropic material.
 * @param etaB, etaS. 
 * @return The 6x6 viscous matrix considering a Voigt notation (arma::mat)
\f[ 
    H_{iso} = \left( \begin{array}{cccccc}
    \eta_B + 4/3 etaS & \eta_B - 2/3 etaS & \eta_B - 2/3 etaS & 0 & 0 & 0 \\
    \eta_B - 2/3 etaS & \eta_B + 4/3 etaS & \eta_B - 2/3 etaS & 0 & 0 & 0 \\
    \eta_B - 2/3 etaS & \eta_B - 2/3 etaS & \eta_B + 4/3 etaS & 0 & 0 & 0 \\
    0 & 0 & 0 & \eta_S & 0 & 0 \\
    0 & 0 & 0 & 0 & \eta_S & 0 \\
    0 & 0 & 0 & 0 & 0 & \eta_S \end{array} \right)
\f]
 * @details Arguments are the Bulk viscosity \f$ \eta_B \f$ (etaB) and shear viscosity \f$ \eta_S \f$ (etaS). 
 * @code 
        double etaB = alead(0., 0.1);
        double etaS = alead(0., 0.1);
        mat Hiso = H_iso(etaB, etaS);
 * @endcode 
 */
arma::mat H_iso(const double &etaB, const double &etaS);

//Update the elastic prediction, providing the stiffness tensor and the trial elastic strain

/**
 * @brief Provides the stress tensor (elastic prediction), from the previous stress increment, providing the elastic stiffness tensor and the trial elastic strain increment:
 * @param sigma_start, L, Eel, ndi
 * @return The predicted 6 vec stress vector (arma::vec)
 * @details The input parameters are : sigma_start: The previous stress, L : Stiffness matrix; Eel : elastic strain vector, ndi (optional, default = 3): number of dimensions
 * @code 
        vec sigma_start = zeros(6);
        sigma_start.randu(6);
        mat L = L_iso(70000, 0.3,"Enu");
        vec Eel;
        Eel.randu(6);
        int ndi = 3;
        vec sigma =  el_pred(sigma_start,L, Eel, ndi);
 * @endcode 
 */
arma::vec el_pred(const arma::vec &sigma_start, const arma::mat &L, const arma::vec &Eel, const int &ndi=3);
    
/**
 * @brief Provides the stress tensor (elastic prediction), from the elastic stiffness tensor and the trial elastic strain:
 * @param L, Eel, ndi
 * @return The predicted 6 vec stress vector (arma::vec)
 * @details The input parameters are : L : Stiffness matrix; Eel ; elastic strain vector, ndi (optional, default = 3): number of dimensions
 * @code 
        mat L = L_iso(70000, 0.3,"Enu");
        vec Eel;
        Eel.randu(6);
        int ndi = 3;
        vec sigma =  el_pred(L, Eel, ndi);
 * @endcode 
*/
arma::vec el_pred(const arma::mat &L, const arma::vec &Eel, const int &ndi=3);

//Return 
/**
 * @brief Provides the isotropized tangent modulus from the spectral decomposition (see \cite Bornet.etal.2001)
 * @param Lt
 * @return The 6x6 isotropic matrix considering a Voigt notation (arma::mat)
 * @details The returned isotropic tensor is called consistent since for any given strain it return the same stress as the anisotropic version.
 * @code 
        double EL = (double)rand();
        double ET = (double)rand();
        double nuTL = (double)rand();
        double nuTT = (double)rand();
        double GLT = (double)rand();
        double axis = 1;
        mat L_isotrans = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);
        mat L_iso = Isotropize(Lisotrans);
 * @endcode 
*/
arma::mat Isotropize(const arma::mat &Lt);

} //namespace simcoon
