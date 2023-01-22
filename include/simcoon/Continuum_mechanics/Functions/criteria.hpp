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
 * @param v, b, n
 * @return The Prager equivalent stress (double)
 * @details Returns the Prager equivalent stress \f$ \mathbf{\sigma}^{P} \f$, considering
\f[
    \sigma^{P} = \sigma^{VM} \left(\frac{1 + b \cdot J_3 \left(\mathbf{\sigma} \right)}{\left(J_2 \left(\mathbf{\sigma} \right) \right)^{3/2} } \right)^{m}
\f]
    considering the input stress \f$ \mathbf{\sigma} \f$, \f$ \mathbf{\sigma}^{VM} \f$ is the Von Mises computed equivalent stress, and \f$ b \f$ and \f$ m \f$ are parameter that define the equivalent stress.
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
 * @param v, b, n
 * @return The derivative of the Prager equivalent stress (arma::vec)
 * @details Returns the derivative of the Prager equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a convex yield surface
\f[
    \frac{\partial \sigma^{P}}{\partial \mathbf{\sigma}} = \sqrt{3} \left( \frac{1. + b * J_3}{J_2^{3/2}} \right)^{1/n-1} \, \frac{1}{2} \frac{\sqrt{J_2}}{\mathbf{\sigma}'} 
    + b \, m \left( 6 J_2^2 \right) \left( 6 J_2 \mathbf{\sigma}' \cdot \mathbf{\sigma}' - 4 J_2^2 \mathbf{I} + \frac{3}{m-9} J_3 \mathbf{\sigma}' \right)
\f]
    considering the input stress \f$ \mathbf{\sigma} \f$, \f$ \mathbf{\sigma}^{VM} \f$ is the Von Mises computed equivalent stress, and \f$ b \f$ and \f$ m \f$ are parameter that define the equivalent stress.
    Note that if n > 10, the Prager criteria is sufficiently close to the Mises that the derivative of the Mises norm is used to avoid numerical instabilities from high-power computations
 * @code 
        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        vec dsigma_Pragerdsigma = dPrager_stress(sigma, b, n);
 * @endcode
*/
arma::vec dPrager_stress(const arma::vec &v, const double &b, const double &n);

/**
 * @brief Provides the Tresca equivalent stress, given its vector representation
 * @param v
 * @return The Prager equivalent stress (double)
 * @details Returns the Tresca equivalent stress \f$ \mathbf{\sigma}^{T} \f$, considering
\f[
    \sigma^{T} = \sigma_{I} - \sigma_{III},
\f]
    considering the input stress \f$ \mathbf{\sigma} \f$.
    Note that the principal stress are classified such that \f$ \sigma_{I} \geq \sigma_{II} \geq \sigma_{III}.
 * @code 
        vec sigma = randu(6);
        double sigma_Tresca = Tresca_stress(sigma);
 * @endcode
*/
double Tresca_stress(const arma::vec &v);

/**
 * @brief Provides derivative of the Tresca equivalent stress, given its vector representation
 * @param v
 * @return The derivative of the Tresca equivalent stress (arma::vec)
 * @details Returns the derivative of the Tresca equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a Tresca convex but not smooth yield surface
 * @warning Note that so far that the correct derivative it is not implemented! Only stress flow \f$ \eta_{stress}=\frac{3/2\sigma_{dev}}{\sigma_{Mises}} \f$ is returned
 * @code 
        vec sigma = randu(6);
        double b = 1.2;
        double m = 0.5;
        vec dsigma_Trescadsigma = dTresca_stress(sigma, b, n);
 * @endcode
*/
arma::vec dTresca_stress(const arma::vec &v);
    
/**
 * @brief an anisotropic configurational tensor \f$ P \f$ in the Voigt format (6x6 matrix), given its vector representation
 * @param P_params
 * @return The anisotropic configurational tensor (arma::mat)
 * @details The vector of parameters must be constituted of 9 values, respectively:
    \f$ P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44}=P_{1212},P_{55}=P_{1313},P_{66}=P_{2323} \f$
\f[
    P_{ani} = \left( \begin{array}{cccccc}
        P_{11} & P_{12} & P_{13} & 0 & 0 & 0 \\
        P_{12} & P_{22} & P_{23} & 0 & 0 & 0 \\
        P_{13} & P_{23} & P_{33} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 \, P_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 \, P_{55} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \, P_{66} \end{array} \right)
\f]
    considering the input stress \f$ \mathbf{\sigma} \f$.
    Note that the equivalent anisotropic tensor is written as : \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$

    It reduces to : 
\f[ 
    \begin{align}
    \sigma^{ani} & = \left( P_{11}\,\sigma_{11}^2 + P_{22}\,\sigma_{22}^2 + P_{33}\,\sigma_{33}^2 \right. \\
     & + 2\,P_{12}\,\sigma_{11}\,\sigma_{22} + 2\,P_{13}\,\sigma_{11}\,\sigma_{33} + 2\,P_{23}\,\sigma_{22} \sigma_{33} \\
     & + \left. 2\,P_{44}\,\sigma_{12}^2 + 2\,P_{55}\,\sigma_{13}^2 + 2\,P_{66}\,\sigma_{23}^2 \right)^{1/2}
    \end{align}
\f]
    Considering the Mises equivalent strain 
\f[ 
    \begin{align}    
    \sigma^{M} & = \left( \frac{1}{2} \left[ \left( \sigma_{11} - \sigma_{22} \right)^2 + \left( \sigma_{11} - \sigma_{33} \right)^2 + \left( \sigma_{22} - \sigma_{33} \right)^2 \right. \\
     & + \left. 6\,\sigma_{12}^2 + 6\,\sigma_{13}^2 + 6\,\sigma_{23}^2 \right)^{1/2}
    \end{align}
\f]
    In that case, \f$ P \f$ reduces to:
\f[
    P_{ani} = \left( \begin{array}{cccccc}
        1 & -1/2 & -1/2 & 0 & 0 & 0 \\
        -1/2 & 1 & -1/2 & 0 & 0 & 0 \\
        -1/2 & -1/2 & -1/2 & 0 & 0 & 0 \\
        0 & 0 & 0 & 3 & 0 & 0 \\
        0 & 0 & 0 & 0 & 3 & 0 \\
        0 & 0 & 0 & 0 & 0 & 3 \end{array} \right)
\f]
 * @code 
        vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
        mat P = P_Ani(P_params);
 * @endcode
*/
arma::mat P_ani(const arma::vec &P_params);

/**
 * @brief Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion \cite Hill.48 in the Voigt format (6x6 matrix), given its vector representation
 * @param P_params
 * @return The anisotropic Hill48 configurational tensor (arma::mat)
 * @details The vector of parameters must be constituted of 5 values, respectively:
    \f$ F^*,G^*,H^*,L,M,N \f$, such that
\f[
    P_{Hill48} = \left( \begin{array}{cccccc}
        G + H & -H & -G & 0 & 0 & 0 \\
        -H & F + H & -F & 0 & 0 & 0 \\
        -G & -F & F + G & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 \, L & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 \, M & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \, N \end{array} \right)
\f]
    considering the input stress \f$ \mathbf{\sigma} \f$.
    Note that the equivalent anisotropic Hill 1948) tensor is written as : \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$
    which reduces to : 
\f[ 
    \begin{align}    
    \sigma^{H48} & = \left( F\, \left( \sigma_{11} - \sigma_{22} \right)^2 + G\, \left( \sigma_{11} - \sigma_{33} \right)^2 + H\, \left( \sigma_{22} - \sigma_{33} \right)^2 \right. \\
     & + \left. 2\,L\,\sigma_{12}^2 + 2\,M\,\sigma_{13}^2 + 2\,N\,\sigma_{23}^2 \right)^{1/2}
    \end{align}
\f]
    Considering the Mises equivalent strain 
\f[ 
    \sigma^{M} & = \left( \frac{1}{2} \left[ \left( \sigma_{11} - \sigma_{22} \right)^2 + \left( \sigma_{11} - \sigma_{33} \right)^2 + \left( \sigma_{22} - \sigma_{33} \right)^2 \right. \\
     & + \left. 6\,\sigma_{12}^2 + 6\,\sigma_{13}^2 + 6\,\sigma_{23}^2 \right)^{1/2}
    \end{align}
\f]
    In that case, \f$ H \f$ reduces to:
\f[
    P_{Hill48} = \left( \begin{array}{cccccc}
        1 & -1/2 & -1/2 & 0 & 0 & 0 \\
        -1/2 & 1 & -1/2 & 0 & 0 & 0 \\
        -1/2 & -1/2 & -1/2 & 0 & 0 & 0 \\
        0 & 0 & 0 & 3 & 0 & 0 \\
        0 & 0 & 0 & 0 & 3 & 0 \\
        0 & 0 & 0 & 0 & 0 & 3 \end{array} \right)
\f]
So that \f$ F = H = G = 1/2 \f$, \f$ = L = M = N = 3/2 \f$ 
 * @code 
        vec P_params = {0.5,0.6,0.7,3.,3.,3.2};
        mat P = P_Hill(P_params);
 * @endcode
*/
arma::mat P_hill(const arma::vec &P_params);

/**
 * @brief Provides the anisotropic equivalent stress, given the stress in a vector format and given a configurational tensor P
 * @param sigma, H
 * @return The anisotropic equivalent stress (double)
 * @details Returns anisotropic equivalent stress \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$.
 * @code 
        vec P_params = {0.5,0.6,0.7,3.,3.,3.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        double sigma_ani = Eq_stress_P(sigma,P_Hill);
 * @endcode
*/
double Eq_stress_P(const arma::vec &sigma, const arma::mat &H);

/**
 * @brief Provides the derivative of the anisotropic equivalent stress, given the stress in a vector format and given a configurational tensor P
 * @param sigma, H
 * @return The derivative of the anisotropic equivalent stress (mat::vec)
 * @details Returns the derivative of the anisotropic equivalent stress \f$ \frac{\partial \mathbf{\sigma}^{ani}}{\partial \mathbf{\sigma}} \f$, considering
\f[
    \sigma^{eq}_{ani} = \frac{\mathbf{H} : \mathbf{\sigma}}{ \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } }
\f]
 * @code 
        vec P_params = {0.5,0.6,0.7,3.,3.,3.2};
        mat P = P_Hill(P_params);
        vec sigma = randu(6);
        vec dsigma_ani = dEq_stress_P(sigma,P_Hill);
 * @endcode
*/
arma::vec dEq_stress_P(const arma::vec &sigma, const arma::mat &H);
    
/**
 * @brief Provides the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N)
 * @param sigma, P_params
 * @return The Hill48 anisotropic equivalent stress (double)
 * @details Returns the Hill48 anisotropic equivalent stress \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$.
 * see the function P_hill() for more details to obtain the tensor H from the set of parameters (F,G,H,L,M,N).
 * @code 
        vec P_params = {0.5,0.6,0.7,3.,3.,3.2};
        vec sigma = randu(6);
        double sigma_Hill = Hill_stress(sigma, P_params);
 * @endcode
*/
double Hill_stress(const arma::vec &sigma, const arma::vec &P_params);

/**
 * @brief Provides the derivative of the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N)
 * @param sigma, P_params
 * @return The derivative of the Hill48 anisotropic equivalent stress (arma::vec)
 * @details Returns the derivative of the anisotropic equivalent stress \f$ \frac{\partial \mathbf{\sigma}^{ani}}{\partial \mathbf{\sigma}} \f$, considering
\f[
    \sigma^{eq}_{ani} = \frac{\mathbf{H} : \mathbf{\sigma}}{ \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } }
\f]
see the function P_hill() for more details to obtain the tensor H from the set of parameters (F,G,H,L,M,N).
 * @code 
        vec P_params = {0.5,0.6,0.7,3.,3.,3.2};
        vec sigma = randu(6);
        double sigma_Hill = Hill_stress(sigma, P_params);
 * @endcode
*/
arma::vec dHill_stress(const arma::vec &, const arma::vec &);

/**
 * @brief Provides the anisotropic equivalent stress, given the stress in a vector format and a vector of parameters \f$ ( P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44},P_{55},P_{66} ) \f$
 * @param sigma, P_params
 * @return The anisotropic equivalent stress (double)
 * @details Returns the anisotropic equivalent stress \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$.
 * see the function P_Ani() for more details to obtain the tensor P from the set of parameters.
 * @code 
        vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
        vec sigma = randu(6);
        double sigma_Ani = Ani_stress(sigma, P_params);
 * @endcode
*/
double Ani_stress(const arma::vec &sigma, const arma::vec &P_params);

/**
 * @brief Provides the derivative of the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N)
 * @param sigma, P_params
 * @return The derivative of the Hill48 anisotropic equivalent stress (arma::vec)
 * @details Returns the derivative of the anisotropic equivalent stress \f$ \frac{\partial \mathbf{\sigma}^{ani}}{\partial \mathbf{\sigma}} \f$, considering
\f[
    \sigma^{eq}_{ani} = \frac{\mathbf{H} : \mathbf{\sigma}}{ \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } }
\f]
see the function P_Ani() for more details to obtain the tensor P from the set of parameters.
 * @code 
        vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
        vec sigma = randu(6);
        vec dsigma_Ani = dAni_stress(sigma, P_params);
 * @endcode
*/
arma::vec dAni_stress(const arma::vec &, const arma::vec &);
    
//This function computes the selected equivalent stress function

/**
 * @brief Provides an equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress and a vector of parameters
 * @param sigma, type_eq, P_params
 * @return The equivalent stress (double)
 * @details Returns the anisotropic equivalent stress, depending on the equivalen type provided : "Mises", "Tresca", "Prager", "Hill", "Ani"
See the detailed function for each equivalent stress to determine the number and type of parameters to provide, in the form of an arma::vec (so for the Prager criteria, P_params = {b,n} for instance)
 * @code 
        vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
        vec sigma = randu(6);
        double sigma_Ani = Ani_stress(sigma, P_params);
 * @endcode
*/
double Eq_stress(const arma::vec &sigma, const std::string &type_eq, const arma::vec &P_params = arma::zeros(1));

/**
 * @brief Provides an derivative of the equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress and a vector of parameters
 * @param sigma, type_eq, P_params
 * @return The deriative of the equivalent stress (arma::vec)
 * @details Returns the anisotropic equivalent stress, depending on the equivalen type provided : "Mises", "Tresca", "Prager", "Hill", "Ani"
See the detailed function for each equivalent stress to determine the number and type of parameters to provide, in the form of an arma::vec (so for the Prager criteria, P_params = {b,n} for instance)
 * @code 
        vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
        vec sigma = randu(6);
        double sigma_Ani = Ani_stress(sigma, P_params);
 * @endcode
*/
arma::vec dEq_stress(const arma::vec &sigma, const std::string &type_eq, const arma::vec &P_params = arma::zeros(1));
    
} //namespace simcoon
