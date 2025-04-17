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
* @section Contimech library contains Functions that computes continuu mechanical quatities and operations on stress/strains, directions, etc
*/

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @brief Returns the deviatoric part of a 3x3 matrix
 * @param m The input 3x3 matrix (arma::mat)
 * @return The 3x3 deviatoric part of the matrix input (arma::mat)
 * @details Example: 
 * @code 
        mat m = randu(3,3);
        mat m_dev = dev(m);
 * @endcode
*/
arma::mat dev(const arma::mat &m);

/**
 * @brief Returns the the spherical part of a 3x3 matrix
 * @param m The input 3x3 matrix (arma::mat)
 * @return The 3x3 spherical part of the matrix input (arma::mat)
 * @details Example: 
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
 * @details Example: 
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
 * @details Example: 
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
 * @details \f$ \sigma^{Mises} = \sqrt{ \frac{3}/{2} \textrm{dev} (\sigma) : \textrm{dev} (\sigma) } \f$
 * Example: 
 * @code 
        vec v = randu(6);
        double Mises_sig = Mises_stress(v);
 * @endcode
*/
double Mises_stress(const arma::vec &v);

/**
 * @brief Provides the stress flow \f$ \eta_{stress}=\frac{3/2 \textrm{dev} (\sigma)}{\sigma_{Mises}} \f$ of a second order stress tensor written as a vector v in Voigt notation
 * @param v The input stress tensor (arma::vec)
 * @return The 6 vector, flow (arma::vec)
 * @details Example: 
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
 * @details Note that the Euclidian norm is expressed as \f$ \sqrt{m:m} \f$, considering that v is the vector representation of a stress matrix m
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
 * @details Note that the Euclidian norm is expressed as \f$ \sqrt{m:m} \f$, considering that v is the vector representation of a strain matrix m
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
 * @details Note that the Euclidian norm is expressed as \f$ \sqrt{m:m} \f$, considering that v is the vector representation of a stress matrix m
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
 * @details Note that the Euclidian norm is expressed as \f$ \sqrt{m:m} \f$, considering that v is the vector representation of a strain matrix m
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
 * @details \f$ \varepsilon^{Mises} = \sqrt{ \frac{2}/{3} \textrm{dev} (\varepsilon) : \textrm{dev} (\varepsilon) } \f$
 * Note that the Euclidian norm is expressed as \f$ \sqrt{m:m} \f$, considering that v is the vector representation of a strain matrix m
 * @code 
        vec v = randu(6);
        double Mises_eps = Mises_strain(v);
 * @endcode
*/
double Mises_strain(const arma::vec &v);

/**
 * @brief Provides the strain flow \f$ \eta_{strain}=\frac{3/2 \textrm{dev} (\varepsilon)}{\varepsilon_{Mises}} \f$ of a second order strain tensor written as a vector v in Voigt notation
 * @param v The input strain tensor (arma::vec)
 * @return The 6 vector, flow (arma::vec)
 * @details Example: 
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
 * @details Note that the second invariant norm is expressed as \f$ \frac{1}/{2} {m':m'} \f$, considering that v is the vector representation of a stress matrix m
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
 * @details Note that the second invariant norm is expressed as \f$ \frac{1}/{2} {m':m'} \f$, considering that v is the vector representation of a strain matrix m
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
 * @details Note that the second invariant norm is expressed as \f$ \frac{1}/{3} {m'_{ij}\,m'_{jk}\,m'_{ki}} \f$, considering that v is the vector representation of a stress matrix m
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
 * @details Note that the second invariant norm is expressed as \f$ \frac{1}/{3} {m'_{ij}\,m'_{jk}\,m'_{ki}} \f$, considering that v is the vector representation of a strain matrix m
 * @code 
        vec v = randu(6);
        double J3 = J3_strain(v);
 * @endcode
*/
double J3_strain(const arma::vec &v);

/**
 * @brief Provides the results of the MacCaulay brackets operator <>+
 * @param d The input value (double)
 * @return The value of <d>+ (double)
 * @details This function returns the value if d positive, zero if d is negative (Macaulay brackets <>+)
 * @code 
    constexpr int MIN = -10;
    constexpr int MAX = 100;
    double d = MIN + (double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)));
    double d_MC = Macaulay_p(d);
 * @endcode
*/
double Macaulay_p(const double &d);

/**
 * @brief Provides the results of the MacCaulay brackets operator <>-
 * @param d The input value (double)
 * @return The value of <d>- (double)
 * @details This function returns the value if d is negative, zero if d is positive (Macaulay brackets <>-)
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
 * @details This function returns the value 1 if d is positive, -1 if d is negative, and 0 if the absolute value of d is less than iota (set in the simcoon parameters, by default 1.E-12)
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
 * @code 
        vec a = randu(6);
        mat H = p_ikjl(a);
 * @endcode
*/
arma::mat p_ikjl(const arma::vec &a);

/**
 * @brief Provides the dyadic product of a symmetric tensor with itself (auto dyadic product)
 * @param a The input symmetric tensor (arma::mat)
 * @details This function returns the operation \f$ c = a \otimes a \f$. The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
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
 * @details This function returns the operation \f$ c = a \otimes b \f$. The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
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
 * @details This function returns the operation \f$ c = a \otimes a \f$. The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
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
 * @details This function returns the operation \f$ c = n_a \otimes n_a \otimes n_b \otimes  n_b \f$ if convention ('conv') is 'aabb' or \f$ c = 0.5 * \left( n_a \otimes n_b \otimes n_a \otimes  n_b + n_a \otimes n_b \otimes n_b \otimes n_a \right) if convention ('conv') is 'abab' \f$.
 *  The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code 
    vec n_a = randu(3);
    vec n_b = randu(3);    
    mat c = dyadic_4vectors_sym(n_a, n_b, "aabb");
 * @endcode
*/
arma::mat dyadic_4vectors_sym(const arma::vec &n_a, const arma::vec &n_b, const std::string &conv);

/**
 * @brief Provides the dyadic product of of two symmetric tensors
 * @param a The first input symmetric tensor (arma::mat)
 * @param b The second input symmetric tensor (arma::mat)
 * @details This function returns the operation \f$ c = a \otimes b \f$. The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
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
 * @details This function returns the operation \f$ c = a \odot a \f$. The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
\f[
    \begin{align}
    \left(\mathbf{a} \odot \mathbf{b} \right)_{ijkl} = \frac{1}/{2} \left( a_{ik} a_{jl} + a_{il} a_{jk} \right)
    \end{align}
\f]
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code 
    mat a = randu(3,3);
    mat b = randu(3,3);
    mat c = auto_sym_dyadic_operator(a);
 * @endcode
*/
arma::mat auto_sym_dyadic_operator(const arma::mat &a);

/**
 * @brief Provides the symmetric 4th-order dyadic product of two symmetric tensors
 * @param a The first input symmetric tensor (arma::mat)
 * @param b The second input symmetric tensor (arma::mat)
 * @details This function returns the operation \f$ c = a \odot b \f$. The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
\f[
    \begin{align}
    \left(\mathbf{a} \odot \mathbf{b} \right)_{ijkl} = \frac{1}/{2} \left( a_{ik} b_{jl} + a_{il} b_{jk} \right)
    \end{align}
\f]
 * @returns The 6x6 matrix that represent the dyadic product (arma::mat)
 * @code 
    mat a = randu(3,3);
    mat b = randu(3,3);
    mat c = sym_dyadic_operator(a,b);
 * @endcode
*/
arma::mat sym_dyadic_operator(const arma::mat &a, const arma::mat &b);

/**
 * @brief Provides the symmetric 4th-order tensor B_klmn, defined as B_i x D x B_j = B_klmn x D, considering B_i and B_j are projection tensors of the vectors b_i and b_j
 * @param b_i The first input vector (arma::vec)
 * @param b_j The second input vector (arma::vec)
 * @details  The function returns a 6x6 matrix that correspond to a 4th order tensor. Note that such conversion to 6x6 matrices product correspond to a conversion with the component of the 4th order tensor correspond to the component of the matrix (such as stiffness matrices)
 * @returns The 6x6 matrix that represent the tensor BBBB
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
arma::mat B_klmn(const arma::vec &b_i, const arma::vec &b_j);

} //namespace simcoon
