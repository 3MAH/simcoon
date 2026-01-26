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
 * @file criteria.hpp
 * @brief Functions for yield surfaces and equivalent stress computations.
 * @version 1.0
 */

#pragma once
#include <string.h>
#include <armadillo>

namespace simcoon
{

/** @addtogroup criteria
 *  @{
 */

    /**
     * @brief Provides the Drucker equivalent stress, given its vector representation
     * @param v The stress vector
     * @param b Parameter that defines the equivalent stress
     * @param n Parameter that defines the equivalent stress
     * @return The Drucker equivalent stress (double)
     * @details Returns the Drucker equivalent stress \f$ \mathbf{\sigma}^{P} \f$, considering
        \f[
            \sigma^{P} = \sigma^{VM} \left(\frac{1 + b \cdot J_3 \left(\mathbf{\sigma} \right)}{\left(J_2 \left(\mathbf{\sigma} \right) \right)^{3/2} } \right)^{m}
        \f]
        considering the input stress \f$ \mathbf{\sigma} \f$, \f$ \mathbf{\sigma}^{VM} \f$ is the Von Mises computed equivalent stress, and \f$ b \f$ and \f$ m \f$ are parameter that define the equivalent stress.
        Note that if n > 10, the Drucker criteria is sufficiently close to the Mises that the Mises norm is used to avoid numerical instabilities from high-power computations
     * @code
            vec sigma = randu(6);
            double b = 1.2;
            double m = 0.5;
            double sigma_Drucker = Drucker_stress(sigma, b, n);
     * @endcode
    */
    double Drucker_stress(const arma::vec &v, const double &b, const double &n);

    /**
     * @brief Provides derivative of the Drucker equivalent stress, given its vector representation
     * @param v The stress vector
     * @param b Parameter that defines the equivalent stress
     * @param n Parameter that defines the equivalent stress
     * @return The derivative of the Drucker equivalent stress (arma::vec)
     * @details Returns the derivative of the Drucker equivalent stress with respect to stress. It main use is to define evolution equations for strain based on an associated rule of a convex yield surface
        \f[
            \frac{\partial \sigma^{P}}{\partial \mathbf{\sigma}} = \sqrt{3} \left( \frac{1. + b * J_3}{J_2^{3/2}} \right)^{1/n-1} \, \frac{1}{2} \frac{\sqrt{J_2}}{\mathbf{\sigma}'}
            + b \, m \left( 6 J_2^2 \right) \left( 6 J_2 \mathbf{\sigma}' \cdot \mathbf{\sigma}' - 4 J_2^2 \mathbf{I} + \frac{3}{m-9} J_3 \mathbf{\sigma}' \right)
        \f]
     *   considering the input stress \f$ \mathbf{\sigma} \f$, \f$ \mathbf{\sigma}^{VM} \f$ is the Von Mises computed equivalent stress, and \f$ b \f$ and \f$ m \f$ are parameter that define the equivalent stress.
     *   Note that if n > 10, the Drucker criteria is sufficiently close to the Mises that the derivative of the Mises norm is used to avoid numerical instabilities from high-power computations
     * @code
            vec sigma = randu(6);
            double b = 1.2;
            double m = 0.5;
            vec dsigma_Druckerdsigma = dDrucker_stress(sigma, b, n);
     * @endcode
    */
    arma::vec dDrucker_stress(const arma::vec &v, const double &b, const double &n);

    /**
     * @brief Provides the Tresca equivalent stress, given its vector representation
     * @param v The stress vector
     * @return The Tresca equivalent stress (double)
     * @details Returns the Tresca equivalent stress \f$ \mathbf{\sigma}^{T} \f$, considering
        \f[
            \sigma^{T} = \sigma_{I} - \sigma_{III},
        \f]
        considering the input stress \f$ \mathbf{\sigma} \f$.
        Note that the principal stress are classified such that \f$ \sigma_{I} \geq \sigma_{II} \geq \sigma_{III} \f$.
     * @code
            vec sigma = randu(6);
            double sigma_Tresca = Tresca_stress(sigma);
     * @endcode
    */
    double Tresca_stress(const arma::vec &v);

    /**
     * @brief Provides derivative of the Tresca equivalent stress, given its vector representation
     * @param v The stress vector
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
     * @brief Provides the derivative of the second stress invariant \f$ J_2 \f$ with respect to stress
     * @param v The stress vector in Voigt notation (6 components)
     * @return The derivative of \f$ J_2 \f$ with respect to stress (arma::vec)
     * @details Returns the derivative of the second deviatoric stress invariant:
    \f[
        \frac{\partial J_2}{\partial \mathbf{\sigma}} = \mathbf{\sigma}_{dev}
    \f]
        where \f$ \mathbf{\sigma}_{dev} \f$ is the deviatoric part of the stress tensor.
        The second invariant \f$ J_2 \f$ is defined as \f$ J_2 = \frac{1}{2} \mathbf{\sigma}_{dev} : \mathbf{\sigma}_{dev} \f$.
     * @code
            vec sigma = randu(6);
            vec dJ2 = dJ2_stress(sigma);
     * @endcode
    */
    arma::vec dJ2_stress(const arma::vec &v);

    /**
     * @brief Provides the derivative of the third stress invariant \f$ J_3 \f$ with respect to stress
     * @param v The stress vector in Voigt notation (6 components)
     * @return The derivative of \f$ J_3 \f$ with respect to stress (arma::vec)
     * @details Returns the derivative of the third deviatoric stress invariant:
    \f[
        \frac{\partial J_3}{\partial \mathbf{\sigma}} = \mathbf{S} \cdot \mathbf{S} - \frac{1}{3} \text{tr}(\mathbf{S}^2) \mathbf{I}
    \f]
        where \f$ \mathbf{S} \f$ is the deviatoric stress tensor and \f$ \mathbf{I} \f$ is the identity tensor.
        The third invariant \f$ J_3 \f$ is defined as \f$ J_3 = \det(\mathbf{\sigma}_{dev}) \f$.
     * @code
            vec sigma = randu(6);
            vec dJ3 = dJ3_stress(sigma);
     * @endcode
    */
    arma::vec dJ3_stress(const arma::vec &v);

    /**
     * @brief Provides an anisotropic Drucker-type equivalent stress combining Drucker criteria with DFA anisotropy
     * @param v The stress vector in Voigt notation (6 components)
     * @param params The DFA parameters vector (F, G, H, L, M, N, K) - see P_DFA() for details
     * @param b Parameter controlling the influence of the third invariant \f$ J_3 \f$
     * @param n Exponent parameter that defines the equivalent stress shape
     * @return The anisotropic Drucker equivalent stress (double)
     * @details Returns an anisotropic equivalent stress that combines the Drucker yield criterion
        with the Deshpande-Fleck-Ashby (DFA) anisotropic formulation:
    \f[
        \sigma^{P}_{ani} = \sigma^{DFA} \left( \frac{1 + b \cdot J_3(\mathbf{\sigma})}{(J_2(\mathbf{\sigma}))^{3/2}} \right)^{1/n}
    \f]
        where \f$ \sigma^{DFA} \f$ is the DFA anisotropic equivalent stress (see DFA_stress()).
        This formulation is particularly useful for modeling porous shape memory alloys (SMA).
        Note that if \f$ n > 10 \f$, the criterion reduces to the DFA stress to avoid numerical instabilities.
     * @code
            vec sigma = randu(6);
            vec params = {0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2};
            double b = 1.2;
            double n = 2.0;
            double sigma_ani = Drucker_anisotrope_stress(sigma, params, b, n);
     * @endcode
    */
    double Drucker_ani_stress(const arma::vec &v, const arma::vec &params, const double &b, const double &n);

    /**
     * @brief Provides the derivative of the anisotropic Drucker-type equivalent stress
     * @param v The stress vector in Voigt notation (6 components)
     * @param params The DFA parameters vector (F, G, H, L, M, N, K) - see P_DFA() for details
     * @param b Parameter controlling the influence of the third invariant \f$ J_3 \f$
     * @param n Exponent parameter that defines the equivalent stress shape
     * @return The derivative of the anisotropic Drucker equivalent stress (arma::vec)
     * @details Returns the derivative of the anisotropic Drucker-DFA equivalent stress with respect to stress.
        This derivative is computed using the chain rule combining the derivatives of the DFA stress,
        \f$ J_2 \f$, and \f$ J_3 \f$ invariants. It is primarily used to define evolution equations
        for strain based on an associated flow rule.
        Note that if \f$ n > 10 \f$, the derivative reduces to the standard stress flow direction
        \f$ \eta_{stress} \f$ to avoid numerical instabilities.
     * @code
            vec sigma = randu(6);
            vec params = {0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2};
            double b = 1.2;
            double n = 2.0;
            vec dsigma_ani = dDrucker_anisotrope_stress(sigma, params, b, n);
     * @endcode
    */
    arma::vec dDrucker_ani_stress(const arma::vec &v, const arma::vec &params, const double &b, const double &n);

    /**
     * @brief Returns an anisotropic configurational tensor \f$ P \f$ in the Voigt format (6x6 matrix)
     * @param P_params The vector of parameters (9 components)
     * @return The anisotropic configurational tensor (arma::mat)
     * @details The vector of parameters must be constituted of 9 values, respectively:
        \f$ P_{11}, P_{22}, P_{33}, P_{12}, P_{13}, P_{23}, P_{44}=P_{1212}, P_{55}=P_{1313}, P_{66}=P_{2323} \f$
    \f[
        P_{ani} = \left( \begin{array}{cccccc}
            P_{11} & P_{12} & P_{13} & 0 & 0 & 0 \\
            P_{12} & P_{22} & P_{23} & 0 & 0 & 0 \\
            P_{13} & P_{23} & P_{33} & 0 & 0 & 0 \\
            0 & 0 & 0 & 2 \, P_{44} & 0 & 0 \\
            0 & 0 & 0 & 0 & 2 \, P_{55} & 0 \\
            0 & 0 & 0 & 0 & 0 & 2 \, P_{66} \end{array} \right)
    \f]
        The equivalent anisotropic stress is written as: \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$

        It reduces to:
    \f[
        \begin{align}
        \sigma^{ani} & = \left( P_{11}\,\sigma_{11}^2 + P_{22}\,\sigma_{22}^2 + P_{33}\,\sigma_{33}^2 \right. \\
         & + 2\,P_{12}\,\sigma_{11}\,\sigma_{22} + 2\,P_{13}\,\sigma_{11}\,\sigma_{33} + 2\,P_{23}\,\sigma_{22} \sigma_{33} \\
         & + \left. 2\,P_{44}\,\sigma_{12}^2 + 2\,P_{55}\,\sigma_{13}^2 + 2\,P_{66}\,\sigma_{23}^2 \right)^{1/2}
        \end{align}
    \f]
        For the isotropic Mises equivalent stress, \f$ P \f$ reduces to:
    \f[
        P_{Mises} = \left( \begin{array}{cccccc}
            1 & -1/2 & -1/2 & 0 & 0 & 0 \\
            -1/2 & 1 & -1/2 & 0 & 0 & 0 \\
            -1/2 & -1/2 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 3 & 0 & 0 \\
            0 & 0 & 0 & 0 & 3 & 0 \\
            0 & 0 & 0 & 0 & 0 & 3 \end{array} \right)
    \f]
     * @code
            vec P_params = {1., 1.2, 1.3, -0.2, -0.2, -0.33, 1., 1., 1.4};
            mat P = P_Ani(P_params);
     * @endcode
    */
    arma::mat P_Ani(const arma::vec &P_params);

    /**
     * @brief Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion (Hill, 1948) in the Voigt format (6x6 matrix), given its vector representation
     * @param P_params The vector of parameters
     * @return The anisotropic Hill48 configurational tensor (arma::mat)
     * @details The vector of parameters must be constituted of 6 values, respectively:
        \f$ F,G,H,L,M,N \f$, such that
    \f[
        P_{Hill48} = \left( \begin{array}{cccccc}
            G + H & -H & -G & 0 & 0 & 0 \\
            -H & F + H & -F & 0 & 0 & 0 \\
            -G & -F & F + G & 0 & 0 & 0 \\
            0 & 0 & 0 & 2 \, N & 0 & 0 \\
            0 & 0 & 0 & 0 & 2 \, M & 0 \\
            0 & 0 & 0 & 0 & 0 & 2 \, L \end{array} \right)
    \f]
        considering the input stress \f$ \mathbf{\sigma} \f$.
        Note that the equivalent anisotropic Hill 1948) tensor is written as : \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$
        which reduces to :
    \f[
        \begin{align}
        \sigma^{H48} & = \left( H\, \left( \sigma_{11} - \sigma_{22} \right)^2 + G\, \left( \sigma_{11} - \sigma_{33} \right)^2 + F\, \left( \sigma_{22} - \sigma_{33} \right)^2 \right. \\
         & + \left. 2\,N\,\sigma_{12}^2 + 2\,M\,\sigma_{13}^2 + 2\,L\,\sigma_{23}^2 \right)^{1/2}
        \end{align}
    \f]
        Considering the Mises equivalent strain
    \f[
        \begin{align}
        \sigma^{M} & = \frac{1}{2} \left( \left( \sigma_{11} - \sigma_{22} \right)^2 + \left( \sigma_{11} - \sigma_{33} \right)^2 + \left( \sigma_{22} - \sigma_{33} \right)^2 \right. \\
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
    So that \f$ F = H = G = 1/2 \f$ and \f$ L = M = N = 3/2 \f$
     * @code
            vec P_params = {0.5,0.6,0.7,3.,3.,3.2};
            mat P = P_Hill(P_params);
     * @endcode
    */
    arma::mat P_Hill(const arma::vec &P_params);

    /**
     * @brief Provides an anisotropic configurational tensor considering the quadratic Deshpande–Fleck–Ashby yield criterion (Deshpande, Fleck & Ashby, 2001) in the Voigt format (6x6 matrix), given its vector representation
     * @param P_params The vector of parameters
     * @return The anisotropic DFA configurational tensor (arma::mat)
     * @details The vector of parameters must be constituted of 7 values, respectively:
        \f$ F,G,H,L,M,N,K \f$, such that
    \f[
        P_{DFA} = \left( \begin{array}{cccccc}
            G + H + K/9 & -H + K/9 & -G + K/9 & 0 & 0 & 0 \\
            -H + K/9 & F + H + K/9 & -F + K/9 & 0 & 0 & 0 \\
            -G + K/9 & -F + K/9 & F + G + K/9 & 0 & 0 & 0 \\
            0 & 0 & 0 & 2 \, N & 0 & 0 \\
            0 & 0 & 0 & 0 & 2 \, M & 0 \\
            0 & 0 & 0 & 0 & 0 & 2 \, L \end{array} \right)
    \f]
        considering the input stress \f$ \mathbf{\sigma} \f$.
        Note that the equivalent anisotropic tensor is written as : \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$
        which reduces to :
    \f[
        \begin{align}
        \sigma^{DFA} & = \left( H\, \left( \sigma_{11} - \sigma_{22} \right)^2 + G\, \left( \sigma_{11} - \sigma_{33} \right)^2 + F\, \left( \sigma_{22} - \sigma_{33} \right)^2 \right. \\
         & + \left. 2\,L\,\sigma_{12}^2 + 2\,M\,\sigma_{13}^2 + 2\,N\,\sigma_{23}^2 \right)^{1/2} + K \left( \left( \sigma_{11} + \sigma_{22} + \sigma_{33} \right) /9 \right)^2
        \end{align}
    \f]
        Considering the full anisotric formulation:
    \f[
        \begin{align}
        \sigma^{ani} & = \left( P_{11}\,\sigma_{11}^2 + P_{22}\,\sigma_{22}^2 + P_{33}\,\sigma_{33}^2 \right. \\
         & + 2\,P_{12}\,\sigma_{11}\,\sigma_{22} + 2\,P_{13}\,\sigma_{11}\,\sigma_{33} + 2\,P_{23}\,\sigma_{22} \sigma_{33} \\
         & + \left. 2\,P_{44}\,\sigma_{12}^2 + 2\,P_{55}\,\sigma_{13}^2 + 2\,P_{66}\,\sigma_{23}^2 \right)^{1/2}
        \end{align}
    \f]
        the above matrix is identified.
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6,1.};
            mat P = P_Hill(P_params);
     * @endcode
    */
    arma::mat P_DFA(const arma::vec &P_params);

    /**
     * @brief Provides the anisotropic equivalent stress, given the stress in a vector format and given a configurational tensor P
     * @param sigma The stress vector
     * @param H The configurational tensor
     * @return The anisotropic equivalent stress (double)
     * @details Returns anisotropic equivalent stress \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$.
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6};
            mat P = P_Hill(P_params);
            vec sigma = randu(6);
            double sigma_ani = Eq_stress_P(sigma,P_Hill);
     * @endcode
    */
    double Eq_stress_P(const arma::vec &sigma, const arma::mat &H);

    /**
     * @brief Provides the derivative of the anisotropic equivalent stress, given the stress in a vector format and given a configurational tensor P
     * @param sigma The stress vector
     * @param H The configurational tensor
     * @return The derivative of the anisotropic equivalent stress (mat::vec)
     * @details Returns the derivative of the anisotropic equivalent stress \f$ \frac{\partial \mathbf{\sigma}^{ani}}{\partial \mathbf{\sigma}} \f$, considering
    \f[
        \sigma^{eq}_{ani} = \frac{\mathbf{H} : \mathbf{\sigma}}{ \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } }
    \f]
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6};
            mat P = P_Hill(P_params);
            vec sigma = randu(6);
            vec dsigma_ani = dEq_stress_P(sigma,P_Hill);
     * @endcode
    */
    arma::vec dEq_stress_P(const arma::vec &sigma, const arma::mat &H);

    /**
     * @brief Provides the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N)
     * @param sigma The stress vector
     * @param P_params The vector of parameters
     * @return The Hill48 anisotropic equivalent stress (double)
     * @details Returns the Hill48 anisotropic equivalent stress \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$.
     * see the function P_Hill() for more details to obtain the tensor H from the set of parameters (F,G,H,L,M,N).
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6};
            vec sigma = randu(6);
            double sigma_Hill = Hill_stress(sigma, P_params);
     * @endcode
    */
    double Hill_stress(const arma::vec &sigma, const arma::vec &P_params);

    /**
     * @brief Provides the derivative of the Hill48 anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N)
     * @param sigma The stress vector
     * @param P_params The vector of parameters
     * @return The derivative of the Hill48 anisotropic equivalent stress (arma::vec)
     * @details Returns the derivative of the anisotropic equivalent stress \f$ \frac{\partial \mathbf{\sigma}^{ani}}{\partial \mathbf{\sigma}} \f$, considering
    \f[
        \sigma^{eq}_{ani} = \frac{\mathbf{H} : \mathbf{\sigma}}{ \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } }
    \f]
    see the function P_Hill() for more details to obtain the tensor H from the set of parameters (F,G,H,L,M,N).
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6};
            vec sigma = randu(6);
            vec dsigma_Hill = dHill_stress(sigma, P_params);
     * @endcode
    */
    arma::vec dHill_stress(const arma::vec &sigma, const arma::vec &P_params);

    /**
     * @brief Provides the Deshpande–Fleck–Ashby (2001) anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N,K)
     * @param sigma The stress vector
     * @param P_params The vector of parameters
     * @return The DFA anisotropic equivalent stress (double)
     * @details Returns the DFA anisotropic equivalent stress \f$ \sigma^{eq}_{ani} = \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } \f$.
     * see the function P_DFA() for more details to obtain the tensor H from the set of parameters (F,G,H,L,M,N,K).
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6,1.2};
            vec sigma = randu(6);
            double sigma_DFA = DFA_stress(sigma, P_params);
     * @endcode
    */
    double DFA_stress(const arma::vec &sigma, const arma::vec &P_params);

    /**
     * @brief Provides the derivative of the Deshpande–Fleck–Ashby (2001) anisotropic equivalent stress, given the stress in a vector format and a vector of parameters (F,G,H,L,M,N,K)
     * @param sigma The stress vector
     * @param P_params The vector of parameters
     * @return The derivative of the DFA anisotropic equivalent stress (arma::vec)
     * @details Returns the derivative of the DFA anisotropic equivalent stress \f$ \frac{\partial \mathbf{\sigma}^{ani}}{\partial \mathbf{\sigma}} \f$, considering
    \f[
        \sigma^{eq}_{ani} = \frac{\mathbf{H} : \mathbf{\sigma}}{ \sqrt{ \mathbf{\sigma} : P : \mathbf{\sigma} } }
    \f]
    see the function P_DFA() for more details to obtain the tensor H from the set of parameters (F,G,H,L,M,N,K).
     * @code
            vec P_params = {0.5,0.6,0.7,1.5,1.5,1.6,1.2};
            vec sigma = randu(6);
            vec dsigma_DFA = dDFA_stress(sigma, P_params);
     * @endcode
    */
    arma::vec dDFA_stress(const arma::vec &sigma, const arma::vec &P_params);

    /**
     * @brief Provides the anisotropic equivalent stress, given the stress in a vector format and a vector of parameters \f$ ( P_{11},P_{22},P_{33},P_{12},P_{13},P_{23},P_{44},P_{55},P_{66} ) \f$
     * @param sigma The stress vector
     * @param P_params The vector of parameters
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
     * @param sigma The stress vector
     * @param P_params The vector of parameters
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
    arma::vec dAni_stress(const arma::vec &sigma, const arma::vec &P_params);

    // This function computes the selected equivalent stress function

    /**
     * @brief Provides an equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress and a vector of parameters
     * @param sigma The stress vector
     * @param type_eq The type of equivalent stress ("Mises", "Tresca", "Drucker", "Hill", "Ani")
     * @param P_params The vector of parameters
     * @return The equivalent stress (double)
     * @details Returns the anisotropic equivalent stress, depending on the equivalen type provided : "Mises", "Tresca", "Drucker", "Hill", "Ani"
    See the detailed function for each equivalent stress to determine the number and type of parameters to provide, in the form of an arma::vec (so for the Drucker criteria, P_params = {b,n} for instance)
     * @code
            vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
            vec sigma = randu(6);
            double sigma_Ani = Ani_stress(sigma, P_params);
     * @endcode
    */
    double Eq_stress(const arma::vec &sigma, const std::string &type_eq, const arma::vec &P_params = arma::zeros(1));

    /**
     * @brief Provides an derivative of the equivalent stress, given the stress in a vector format, a string indicating the type of equivalent stress and a vector of parameters
     * @param sigma The stress vector
     * @param type_eq The type of equivalent stress ("Mises", "Tresca", "Drucker", "Hill", "Ani")
     * @param P_params The vector of parameters
     * @return The deriative of the equivalent stress (arma::vec)
     * @details Returns the anisotropic equivalent stress, depending on the equivalen type provided : "Mises", "Tresca", "Drucker", "Hill", "Ani"
    See the detailed function for each equivalent stress to determine the number and type of parameters to provide, in the form of an arma::vec (so for the Drucker criteria, P_params = {b,n} for instance)
     * @code
            vec P_params = {1.,1.2,1.3,-0.2,-0.2,-0.33,1.,1.,1.4};
            vec sigma = randu(6);
            double sigma_Ani = Ani_stress(sigma, P_params);
     * @endcode
    */
    arma::vec dEq_stress(const arma::vec &sigma, const std::string &type_eq, const arma::vec &P_params = arma::zeros(1));

/** @} */ // end of criteria group

} // namespace simcoon
