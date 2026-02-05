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
 * @file contimech.hpp
 * @author Yves Chemisky 
 * @brief Contimech library contains functions that compute continuum mechanical quantities and operations on stress/strains, directions, etc.
 */

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup contimech
 *  @{
 */

/**
 * @brief Returns the deviatoric part of a 3x3 matrix
 * @param m The input 3x3 matrix (arma::mat)
 * @return The 3x3 deviatoric part of the matrix input (arma::mat)
 * @details The deviatoric part of a tensor \f$ \mathbf{m} \f$ is defined as:
\f[
    \mathbf{m}' = \mathbf{m} - \frac{1}{3} \textrm{tr}(\mathbf{m}) \mathbf{I}
\f]
 * where \f$ \textrm{tr}(\mathbf{m}) = m_{11} + m_{22} + m_{33} \f$ is the trace and \f$ \mathbf{I} \f$ is the identity tensor.
 * Example:
 * @code
        mat m = randu(3,3);
        mat m_dev = dev(m);
 * @endcode
*/
arma::mat dev(const arma::mat &m);

/**
 * @brief Returns the spherical part of a 3x3 matrix
 * @param m The input 3x3 matrix (arma::mat)
 * @return The 3x3 spherical part of the matrix input (arma::mat)
 * @details The spherical (or hydrostatic) part of a tensor \f$ \mathbf{m} \f$ is defined as:
\f[
    \mathbf{m}^{sph} = \frac{1}{3} \textrm{tr}(\mathbf{m}) \mathbf{I}
\f]
 * where \f$ \textrm{tr}(\mathbf{m}) = m_{11} + m_{22} + m_{33} \f$ is the trace and \f$ \mathbf{I} \f$ is the identity tensor.
 * Note that \f$ \mathbf{m} = \mathbf{m}' + \mathbf{m}^{sph} \f$.
 * Example:
 * @code
        mat m = randu(3,3);
        mat m_sph = sph(m);
 * @endcode
*/
arma::mat sph(const arma::mat &m);

/**
 * @brief Returns the trace of a tensor v expressed in Voigt notation
 * @param v The input tensor (arma::vec)
 * @return The trace of the input tensor (double)
 * @details The trace of a second-order tensor \f$ \mathbf{m} \f$ written in Voigt notation as \f$ \mathbf{v} = (v_1, v_2, v_3, v_4, v_5, v_6)^T \f$ is:
\f[
    \textrm{tr}(\mathbf{m}) = v_1 + v_2 + v_3 = m_{11} + m_{22} + m_{33}
\f]
 * Example:
 * @code
        vec v = randu(6);
        double trace = tr(v);
 * @endcode
*/
double tr(const arma::vec &v);

/**
 * @brief Returns the deviatoric part of a tensor v expressed in Voigt notation
 * @param v The input tensor (arma::vec)
 * @return The 6 vector, deviatoric part (arma::vec)
 * @details The deviatoric part of a tensor \f$ \mathbf{m} \f$ expressed in Voigt notation is:
\f[
    \mathbf{v}' = \mathbf{v} - \frac{1}{3} \textrm{tr}(\mathbf{v}) \mathbf{I}_{th}
\f]
 * where \f$ \mathbf{I}_{th} = (1, 1, 1, 0, 0, 0)^T \f$ is the expansion vector.
 * This gives \f$ v'_i = v_i - \frac{1}{3}(v_1 + v_2 + v_3) \f$ for \f$ i = 1,2,3 \f$ and \f$ v'_i = v_i \f$ for \f$ i = 4,5,6 \f$.
 * Example:
 * @code
        vec v = randu(6);
        vec v_dev = dev(v);
 * @endcode
*/
arma::vec dev(const arma::vec &v);

/**
 * @brief Provides the Von Mises stress \f$ \sigma^{Mises} \f$ of a second order stress tensor written as a vector v in Voigt notation
 * @param v The input stress tensor (arma::vec)
 * @return The Mises equivalent (double)
 * @details The Von Mises equivalent stress is defined as:
\f[
    \sigma^{Mises} = \sqrt{ \frac{3}{2} \mathbf{\sigma}' : \mathbf{\sigma}' } = \sqrt{3 \, J_2}
\f]
 * where \f$ \mathbf{\sigma}' = \textrm{dev}(\mathbf{\sigma}) \f$ is the deviatoric stress tensor and \f$ J_2 \f$ is the second deviatoric invariant.
 * Example:
 * @code
        vec v = randu(6);
        double Mises_sig = Mises_stress(v);
 * @endcode
*/
double Mises_stress(const arma::vec &v);

/**
 * @brief Provides the stress flow \f$ \boldsymbol{\eta}_{stress} \f$ of a second order stress tensor written as a vector v in Voigt notation
 * @param v The input stress tensor (arma::vec)
 * @return The 6 vector, flow (arma::vec)
 * @details The stress flow direction is defined as:
\f[
    \boldsymbol{\eta}_{stress} = \frac{\frac{3}{2} \mathbf{\sigma}'}{\sigma^{Mises}} = \frac{\frac{3}{2} \textrm{dev}(\mathbf{\sigma})}{\sigma^{Mises}}
\f]
 * This corresponds to the normal to the Von Mises yield surface and is used in associated flow rules.
 * Example:
 * @code
        vec v = randu(6);
        vec v_flow = eta_stress(v);
 * @endcode
*/
arma::vec eta_stress(const arma::vec &v);
    
/**
 * @brief Provides the strain flow (direction) from a stress tensor (Euclidian norm), according to the Voigt convention for strains
 * @param v The input stress tensor (arma::vec)
 * @return The 6 vector, flow (arma::vec)
 * @details The normalized flow direction based on the Euclidian norm is:
\f[
    \boldsymbol{\eta}_{norm} = \frac{\mathbf{\sigma}}{\|\mathbf{\sigma}\|} = \frac{\mathbf{\sigma}}{\sqrt{\mathbf{\sigma}:\mathbf{\sigma}}}
\f]
 * Note that the Euclidian norm is expressed as \f$ \|\mathbf{m}\| = \sqrt{\mathbf{m}:\mathbf{m}} \f$, considering that v is the vector representation of a stress matrix m.
 * @code
        vec v = randu(6);
        vec v_flow = eta_norm_stress(v);
 * @endcode
*/
arma::vec eta_norm_stress(const arma::vec &v);

/**
 * @brief Provides the strain flow (direction) from a strain tensor (Euclidian norm), according to the Voigt convention for strains
 * @param v The input strain tensor (arma::vec)
 * @return The 6 vector, flow (arma::vec)
 * @details The normalized flow direction based on the Euclidian norm is:
\f[
    \boldsymbol{\eta}_{norm} = \frac{\boldsymbol{\varepsilon}}{\|\boldsymbol{\varepsilon}\|} = \frac{\boldsymbol{\varepsilon}}{\sqrt{\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon}}}
\f]
 * Note that the Euclidian norm is expressed as \f$ \|\mathbf{m}\| = \sqrt{\mathbf{m}:\mathbf{m}} \f$, considering that v is the vector representation of a strain matrix m.
 * @code
        vec v = randu(6);
        vec v_flow = eta_norm_strain(v);
 * @endcode
*/
arma::vec eta_norm_strain(const arma::vec &v);
    
/**
 * @brief Provides the Euclidian norm of a stress tensor
 * @param v The input stress tensor (arma::vec)
 * @return The norm (double)
 * @details The Euclidian (Frobenius) norm of a stress tensor is:
\f[
    \|\mathbf{\sigma}\| = \sqrt{\mathbf{\sigma}:\mathbf{\sigma}} = \sqrt{\sigma_{ij}\,\sigma_{ij}}
\f]
 * In Voigt notation: \f$ \|\mathbf{\sigma}\| = \sqrt{v_1^2 + v_2^2 + v_3^2 + 2(v_4^2 + v_5^2 + v_6^2)} \f$
 * @code
        vec v = randu(6);
        double norm=norm_stress(v);
 * @endcode
*/
double norm_stress(const arma::vec &v);

/**
 * @brief Provides the Euclidian norm of a strain tensor
 * @param v The input strain tensor (arma::vec)
 * @return The norm (double)
 * @details The Euclidian (Frobenius) norm of a strain tensor is:
\f[
    \|\boldsymbol{\varepsilon}\| = \sqrt{\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon}} = \sqrt{\varepsilon_{ij}\,\varepsilon_{ij}}
\f]
 * In Voigt notation: \f$ \|\boldsymbol{\varepsilon}\| = \sqrt{v_1^2 + v_2^2 + v_3^2 + \frac{1}{2}(v_4^2 + v_5^2 + v_6^2)} \f$
 * @code
        vec v = randu(6);
        double norm=norm_strain(v);
 * @endcode
*/
double norm_strain(const arma::vec &v);
    
/**
 * @brief Provides the Von Mises strain \f$ \varepsilon^{Mises} \f$ of a second order strain tensor
 * @param v The input strain tensor (arma::vec)
 * @return The Mises norm (double)
 * @details The Von Mises equivalent strain is defined as:
\f[
    \varepsilon^{Mises} = \sqrt{ \frac{2}{3} \boldsymbol{\varepsilon}' : \boldsymbol{\varepsilon}' }
\f]
 * where \f$ \boldsymbol{\varepsilon}' = \textrm{dev}(\boldsymbol{\varepsilon}) \f$ is the deviatoric strain tensor.
 * This definition is conjugate to the Von Mises stress such that \f$ \sigma^{Mises} \cdot \varepsilon^{Mises} = \mathbf{\sigma}' : \boldsymbol{\varepsilon}' \f$.
 * @code
        vec v = randu(6);
        double Mises_eps = Mises_strain(v);
 * @endcode
*/
double Mises_strain(const arma::vec &v);

/**
 * @brief Provides the strain flow \f$ \boldsymbol{\eta}_{strain} \f$ of a second order strain tensor written as a vector v in Voigt notation
 * @param v The input strain tensor (arma::vec)
 * @return The 6 vector, flow (arma::vec)
 * @details The strain flow direction is defined as:
\f[
    \boldsymbol{\eta}_{strain} = \frac{\frac{2}{3} \boldsymbol{\varepsilon}'}{\varepsilon^{Mises}} = \frac{\frac{2}{3} \textrm{dev}(\boldsymbol{\varepsilon})}{\varepsilon^{Mises}}
\f]
 * This is conjugate to the stress flow \f$ \boldsymbol{\eta}_{stress} \f$ used in plasticity.
 * Example:
 * @code
        vec v = randu(6);
        vec eps_f = eta_strain(v);
 * @endcode
*/
arma::vec eta_strain(const arma::vec &v);

/**
 * @brief Provides the second invariant of the deviatoric part of a second order stress tensor, written as a vector v in Voigt notation
 * @param v The input stress tensor (arma::vec)
 * @return The second invariant (double)
 * @details The second deviatoric stress invariant \f$ J_2 \f$ is defined as:
\f[
    J_2 = \frac{1}{2} \mathbf{\sigma}' : \mathbf{\sigma}' = \frac{1}{2} \sigma'_{ij} \, \sigma'_{ij}
\f]
 * It is related to the Von Mises stress by \f$ \sigma^{Mises} = \sqrt{3 \, J_2} \f$.
 * @code
        vec v = randu(6);
        double J2 = J2_stress(v);
 * @endcode
*/
double J2_stress(const arma::vec &v);

/**
 * @brief Provides the second invariant of the deviatoric part of a second order strain tensor, written as a vector v in Voigt notation
 * @param v The input strain tensor (arma::vec)
 * @return The second invariant (double)
 * @details The second deviatoric strain invariant \f$ J_2 \f$ is defined as:
\f[
    J_2 = \frac{1}{2} \boldsymbol{\varepsilon}' : \boldsymbol{\varepsilon}' = \frac{1}{2} \varepsilon'_{ij} \, \varepsilon'_{ij}
\f]
 * where \f$ \boldsymbol{\varepsilon}' = \textrm{dev}(\boldsymbol{\varepsilon}) \f$ is the deviatoric strain tensor.
 * @code
        vec v = randu(6);
        double J2 = J2_strain(v);
 * @endcode
*/
double J2_strain(const arma::vec &v);

/**
 * @brief Provides the third invariant of the deviatoric part of a second order stress tensor, written as a vector v in Voigt notation
 * @param v The input stress tensor (arma::vec)
 * @return The third invariant (double)
 * @details The third deviatoric stress invariant \f$ J_3 \f$ is defined as:
\f[
    J_3 = \det(\mathbf{\sigma}') = \frac{1}{3} \sigma'_{ij} \, \sigma'_{jk} \, \sigma'_{ki}
\f]
 * This invariant characterizes the Lode angle and is used in pressure-dependent yield criteria.
 * @code
        vec v = randu(6);
        double J3 = J3_stress(v);
 * @endcode
*/
double J3_stress(const arma::vec &v);

/**
 * @brief Provides the third invariant of the deviatoric part of a second order strain tensor, written as a vector v in Voigt notation
 * @param v The input strain tensor (arma::vec)
 * @return The third invariant (double)
 * @details The third deviatoric strain invariant \f$ J_3 \f$ is defined as:
\f[
    J_3 = \det(\boldsymbol{\varepsilon}') = \frac{1}{3} \varepsilon'_{ij} \, \varepsilon'_{jk} \, \varepsilon'_{ki}
\f]
 * where \f$ \boldsymbol{\varepsilon}' = \textrm{dev}(\boldsymbol{\varepsilon}) \f$ is the deviatoric strain tensor.
 * @code
        vec v = randu(6);
        double J3 = J3_strain(v);
 * @endcode
*/
double J3_strain(const arma::vec &v);

/**
 * @brief Provides the results of the MacCaulay brackets operator \f$\langle d \rangle^+\f$
 * @param d The input value (double)
 * @return The value of \f$\langle d \rangle^+\f$ (double)
 * @details The positive Macaulay bracket (ramp function) is defined as:
\f[
    \langle d \rangle^+ = \begin{cases} d & \text{if } d > 0 \\ 0 & \text{if } d \leq 0 \end{cases} = \frac{d + |d|}{2}
\f]
 * @code
    constexpr int MIN = -10;
    constexpr int MAX = 100;
    double d = MIN + (double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)));
    double d_MC = Macaulay_p(d);
 * @endcode
*/
double Macaulay_p(const double &d);

/**
 * @brief Provides the results of the MacCaulay brackets operator \f$\langle d \rangle^-\f$
 * @param d The input value (double)
 * @return The value of \f$\langle d \rangle^-\f$ (double)
 * @details The negative Macaulay bracket is defined as:
\f[
    \langle d \rangle^- = \begin{cases} 0 & \text{if } d \geq 0 \\ d & \text{if } d < 0 \end{cases} = \frac{d - |d|}{2}
\f]
 * Note that \f$ d = \langle d \rangle^+ + \langle d \rangle^- \f$.
 * @code
    constexpr int MIN = -10;
    constexpr int MAX = 100;
    double d = MIN + (double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)));
    double d_MC = Macaulay_n(d);
 * @endcode
*/
double Macaulay_n(const double &d);

/**
 * @brief Provides the results of the sign operator
 * @param d The input value (double)
 * @return The sign of the input value (double)
 * @details The sign (signum) function is defined as:
\f[
    \textrm{sign}(d) = \begin{cases} +1 & \text{if } d > \iota \\ 0 & \text{if } |d| \leq \iota \\ -1 & \text{if } d < -\iota \end{cases}
\f]
 * where \f$ \iota \f$ is a small tolerance (set in the simcoon parameters, by default \f$ 10^{-12} \f$).
 * @code
    constexpr int MIN = -10;
    constexpr int MAX = 100;
    double d = MIN + (double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)));
    double d_sign = sign(d);
 * @endcode
*/
double sign(const double &d);

/**
 * @brief Provides the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u and v
 * @param u The angle in the plane defined by the center of the ellipsoid, a1 and a2 directions (double)
 * @param v The angle in the plane defined by the center of the ellipsoid, a1 and a3 directions (double)
 * @param a1 The length of the semi-principal axis a1 (double)
 * @param a2 The length of the semi-principal axis a2 (double)
 * @param a3 The length of the semi-principal axis a3 (double)
 * @return The normal vector (arma::vec)
 * @details Provides the normalized vector to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u and v. These 2 angles correspond to the rotations in the plan defined by the center of the ellipsoid, a1 and a2 directions for u, a1 and a3 ones for v. u = 0 corresponds to a1 direction and v = 0 correspond to a3 one. So the normal vector is set at the parametrized position :
\f[
    \begin{align}
    x & = a_{1} \cos(u) \sin(v) \\
    y & = a_{2} \sin(u) \sin(v) \\
    z & = a_{3} \cos(v)
    \end{align}
\f]
 * @code 
    const double Pi = 3.14159265358979323846;

    double u = (double)rand()/(double)(RAND_MAX) % (2*Pi) - (2*Pi);
    double v = (double)rand()/(double)(RAND_MAX) % Pi - Pi;
    double a1 = (double)rand();
    double a2 = (double)rand();
    double a3 = (double)rand();
    vec normal = normal_ellipsoid(u, v, a1, a2, a3);
 * @endcode
*/
arma::vec normal_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3);

/**
 * @brief Provides the curvature of an ellipsoid with semi-principal axes of length a1, a2, a3 at the angle u,v.
 * @param u The angle in the plane defined by the center of the ellipsoid, a1 and a2 directions (double)
 * @param v The angle in the plane defined by the center of the ellipsoid, a1 and a3 directions (double)
 * @param a1 The length of the semi-principal axis a1 (double)
 * @param a2 The length of the semi-principal axis a2 (double)
 * @param a3 The length of the semi-principal axis a3 (double)
 * @return The curvature (double)
 * @details Provides the normalized curvature of an ellipsoid with semi-principal axes of length a1, a2, a3. The position of the evaluated curvature is set by angles u and v. These 2 angles correspond to the rotations in the plan defined by the center of the ellipsoid, a1 and a2 directions for u, a1 and a3 ones for v. u = 0 corresponds to a1 direction and v = 0 correspond to a3 one. So the curvature is set at the parametrized position :
\f[
    \begin{align}
    x & = a_{1} \cos(u) \sin(v) \\
    y & = a_{2} \sin(u) \sin(v) \\
    z & = a_{3} \cos(v)
    \end{align}
    \Xi = \frac{a_1^2\,a_2^2\,a_3^2}{\left(a_1^2\,a_2^2\,\cos^2 v + a_3^2\,\sin^2 v +  a_2^2\,\cos^2 v + a_1^2\,\sin u \right)^2 }
\f]
 * @code 
    const double Pi = 3.14159265358979323846;

    double u = (double)rand()/(double)(RAND_MAX) % (2*Pi) - (2*Pi);
    double v = (double)rand()/(double)(RAND_MAX) % Pi - Pi;
    double a1 = (double)rand();
    double a2 = (double)rand();
    double a3 = (double)rand();
    double curvature = curvature_ellipsoid(u, v, a1, a2, a3);
 * @endcode
*/
double curvature_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3);

/**
 * @brief Provides the normal and tangent components of the traction vector in the normal direction n to an ellipsoid with axes a1, a2, a3 from an input inside stress \f$ \sigma \f$. The direction of the normalized vector is set by angles u and v
 * @param sigma_in The input stress tensor (arma::vec)
 * @param u The angle in the plane defined by the center of the ellipsoid, a1 and a2 directions (double)
 * @param v The angle in the plane defined by the center of the ellipsoid, a1 and a3 directions (double)
 * @param a1 The length of the semi-principal axis a1 (double)
 * @param a2 The length of the semi-principal axis a2 (double)
 * @param a3 The length of the semi-principal axis a3 (double)
 * @return The normal and tangent components of the traction vector (arma::vec)
 * @details Provides the normal and tangent components of a stress vector σin in accordance with the normal direction n to an ellipsoid with axes a1, a2, a3. The normal vector is set at the parametrized position :
\f[
    \begin{align}
    x & = a_{1} \cos(u) \sin(v) \\
    y & = a_{2} \sin(u) \sin(v) \\
    z & = a_{3} \cos(v)
    \end{align}
\f]
 * @code 
    const double Pi = 3.14159265358979323846;

    vec sigma_in = randu(6);
    double u = (double)rand()/(double)(RAND_MAX) % Pi - Pi/2;
    double v = (double)rand()/(double)(RAND_MAX) % (2*Pi) - Pi;
    double a1 = (double)rand();
    double a2 = (double)rand();
    double a3 = (double)rand();
    vec sigma_i = sigma_int(sigma_in, u, v, a1, a2, a3);
 * @endcode
*/
arma::vec sigma_int(const arma::vec &sigma_in, const double &u, const double &v, const double &a1, const double &a2, const double &a3);

/**
 * @brief Provides the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
 * @param a The input normal vector (arma::vec)
 * @return Hill interfacial operator 6x6 matrix (arma::mat)
 * @details The Hill interfacial operator \f$ \mathbf{P} \f$ is defined in index notation as:
\f[
    P_{ikjl} = \frac{1}{2} \left( a_i \, \delta_{kj} \, a_l + a_j \, \delta_{ik} \, a_l + a_i \, \delta_{lj} \, a_k + a_j \, \delta_{il} \, a_k \right)
\f]
 * where \f$ \mathbf{a} \f$ is the normal vector and \f$ \delta_{ij} \f$ is the Kronecker delta.
 * @code
        vec a = randu(3);
        mat H = p_ikjl(a);
 * @endcode
*/
arma::mat p_ikjl(const arma::vec &a);

/**
 * @brief Provides the dyadic product of a symmetric tensor with itself (auto dyadic product)
 * @param a The input symmetric tensor (arma::mat)
 * @details This function returns the dyadic (tensor) product of a symmetric tensor with itself:
\f[
    \mathbf{C} = \mathbf{a} \otimes \mathbf{a}, \quad C_{ijkl} = a_{ij} \, a_{kl}
\f]
 * The function returns a 6x6 matrix that corresponds to a 4th order tensor in Voigt notation (such as stiffness matrices).
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code
    mat a = randu(3,3);
    mat sym_a = 0.5*(a + a.t());
    mat c = auto_sym_dyadic(sym_a);
 * @endcode
*/
arma::mat auto_sym_dyadic(const arma::mat &a);

/**
 * @brief Provides the dyadic product of two symmetric tensors
 * @param a The first input symmetric tensor (arma::mat)
 * @param b The second input symmetric tensor (arma::mat)
 * @details This function returns the dyadic (tensor) product of two symmetric tensors:
\f[
    \mathbf{C} = \mathbf{a} \otimes \mathbf{b}, \quad C_{ijkl} = a_{ij} \, b_{kl}
\f]
 * The function returns a 6x6 matrix that corresponds to a 4th order tensor in Voigt notation (such as stiffness matrices).
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code
    mat a = randu(3,3);
    mat sym_a = 0.5*(a + a.t());
    mat b = randu(3,3);
    mat sym_b = 0.5*(b + b.t());
    mat c = sym_dyadic(sym_a,sym_b);
 * @endcode
*/
arma::mat sym_dyadic(const arma::mat &a, const arma::mat &b);

/**
 * @brief Provides the dyadic product of a tensor with itself (auto dyadic product)
 * @param a The input tensor (arma::mat)
 * @details This function returns the dyadic (tensor) product of a tensor with itself:
\f[
    \mathbf{C} = \mathbf{a} \otimes \mathbf{a}, \quad C_{ijkl} = a_{ij} \, a_{kl}
\f]
 * The function returns a 6x6 matrix that corresponds to a 4th order tensor in Voigt notation.
 * Unlike auto_sym_dyadic, this function does not assume the input tensor is symmetric.
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code
    mat a = randu(3,3);
    mat c = auto_dyadic(a);
 * @endcode
*/
arma::mat auto_dyadic(const arma::mat &a);

/**
 * @brief Provides the dyadic product of four vectors to provide a symmetric 4th order tensor (minor symmetry)
 * @param n_a The first input vector (arma::vec)
 * @param n_b The second input vector (arma::vec)
 * @param conv The convention for the dyadic product ("aabb" or "abab") (std::string)
 * @details This function returns a 4th order tensor from two vectors with minor symmetry.
 * For convention "aabb":
\f[
    C_{ijkl} = (n_a)_i \, (n_a)_j \, (n_b)_k \, (n_b)_l
\f]
 * For convention "abab":
\f[
    C_{ijkl} = \frac{1}{2} \left( (n_a)_i \, (n_b)_j \, (n_a)_k \, (n_b)_l + (n_a)_i \, (n_b)_j \, (n_b)_k \, (n_a)_l \right)
\f]
 * The function returns a 6x6 matrix corresponding to a 4th order tensor in Voigt notation.
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code
    vec n_a = randu(3);
    vec n_b = randu(3);
    mat c = dyadic_4vectors_sym(n_a, n_b, "aabb");
 * @endcode
*/
arma::mat dyadic_4vectors_sym(const arma::vec &n_a, const arma::vec &n_b, const std::string &conv);

/**
 * @brief Provides the dyadic product of two tensors
 * @param a The first input tensor (arma::mat)
 * @param b The second input tensor (arma::mat)
 * @details This function returns the dyadic (tensor) product of two tensors:
\f[
    \mathbf{C} = \mathbf{a} \otimes \mathbf{b}, \quad C_{ijkl} = a_{ij} \, b_{kl}
\f]
 * The function returns a 6x6 matrix corresponding to a 4th order tensor in Voigt notation (such as stiffness matrices).
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code
    mat a = randu(3,3);
    mat b = randu(3,3);
    mat c = dyadic(a,b);
 * @endcode
*/
arma::mat dyadic(const arma::mat &a, const arma::mat &b);

/**
 * @brief Provides the symmetric 4th-order dyadic product of a symmetric tensor with itself
 * @param a The input symmetric tensor (arma::mat)
 * @details This function returns the symmetric dyadic product \f$ \mathbf{C} = \mathbf{a} \odot \mathbf{a} \f$, defined as:
\f[
    \left(\mathbf{a} \odot \mathbf{a} \right)_{ijkl} = \frac{1}{2} \left( a_{ik} \, a_{jl} + a_{il} \, a_{jk} \right)
\f]
 * The function returns a 6x6 matrix corresponding to a 4th order tensor in Voigt notation (possessing both major and minor symmetries).
 * @returns The 6x6 matrix that represent the symmetric dyadic product (arma::mat)
 * @code
    mat a = randu(3,3);
    mat sym_a = 0.5*(a + a.t());
    mat c = auto_sym_dyadic_operator(sym_a);
 * @endcode
*/
arma::mat auto_sym_dyadic_operator(const arma::mat &a);

/**
 * @brief Provides the symmetric 4th-order dyadic product of two symmetric tensors
 * @param a The first input symmetric tensor (arma::mat)
 * @param b The second input symmetric tensor (arma::mat)
 * @details This function returns the symmetric dyadic product \f$ \mathbf{C} = \mathbf{a} \odot \mathbf{b} \f$, defined as:
\f[
    \left(\mathbf{a} \odot \mathbf{b} \right)_{ijkl} = \frac{1}{2} \left( a_{ik} \, b_{jl} + a_{il} \, b_{jk} \right)
\f]
 * The function returns a 6x6 matrix corresponding to a 4th order tensor in Voigt notation (possessing minor symmetry).
 * @returns The 6x6 matrix that represent the symmetric dyadic product (arma::mat)
 * @code
    mat a = randu(3,3);
    mat sym_a = 0.5*(a + a.t());
    mat b = randu(3,3);
    mat sym_b = 0.5*(b + b.t());
    mat c = sym_dyadic_operator(sym_a,sym_b);
 * @endcode
*/
arma::mat sym_dyadic_operator(const arma::mat &a, const arma::mat &b);

/**
 * @brief Provides a symmetric 4th-order tensor \f$ \mathbb{G}^{ij} \f$, defined such that \f$ \mathbf{B}_i : \mathbf{D} : \mathbf{B}_j = \mathbb{G}^{ij} : \mathbf{D} \f$, where \f$ \mathbf{B}_i \f$ and \f$ \mathbf{B}_j \f$ are eigenprojection tensors
 * @param b_i The first eigenvector (arma::vec)
 * @param b_j The second eigenvector (arma::vec)
 * @details The 4th-order tensor \f$ \mathbb{G}^{ij} \f$ is constructed from eigenvectors \f$ \mathbf{n}^i \f$ and \f$ \mathbf{n}^j \f$.
 * It is defined such that the triple product can be expressed as a single double contraction:
\f[
    \mathbf{B}_i : \mathbf{D} : \mathbf{B}_j = \mathbb{G}^{ij} : \mathbf{D}
\f]
 * where \f$ \mathbf{B}_i = \mathbf{n}^i \otimes \mathbf{n}^i \f$ and \f$ \mathbf{B}_j = \mathbf{n}^j \otimes \mathbf{n}^j \f$ are the eigenprojection tensors.
 * Defining \f$ \mathbf{N} = \mathbf{n}^i \otimes \mathbf{n}^j \f$, the tensor components are:
\f[
    \mathbb{G}^{ij}_{klmn} = \frac{1}{2} N_{kl} \left( N_{mn} + N_{nm} \right)
\f]
 * This tensor is useful in logarithmic strain computations, particularly for the derivative of the Hencky strain with respect to the Cauchy-Green tensor (see Xiao, Bruhns & Meyers, 1997).
 * The function returns a 6x6 matrix corresponding to a 4th order tensor in Voigt notation.
 * @returns The 6x6 matrix that represents the 4th-order tensor \f$ \mathbb{G}^{ij} \f$
 * @code
    mat B = L_Cauchy_Green(F1);

    vec bi = zeros(3);
    mat Bi;
    eig_sym(bi, Bi, B);
    mat BBBB = zeros(6,6);

    double f_z = 0.;
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                BBBB = BBBB + f_z*B_klmn(Bi.col(i),Bi.col(j));
            }
        }
    }
 * @endcode
*/
arma::mat linearop_eigsym(const arma::vec &b_i, const arma::vec &b_j);

/** @} */ // end of contimech group

} //namespace simcoon
