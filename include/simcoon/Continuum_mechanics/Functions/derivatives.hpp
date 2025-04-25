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

///@file derivatives.hpp
///@brief A set of functions that performs the derivative of specific tensors attributes (Invariants, determinant, inverse...)
///@version 1.0

#pragma once
#include <armadillo>
#include <simcoon/FTensor.hpp>

namespace simcoon{

// This function returns the derivative of the first invariant (trace) of a tensor
// @param S The input tensor
// @return The derivative of the first invariant
// @example arma::mat result = dI1DS(S);
arma::mat dI1DS(const arma::mat &);

// This function returns the derivative of the second invariant of a tensor : I_2 = 1/2 S_ij S_ij
// @param S The input tensor
// @return The derivative of the second invariant
// @example arma::mat result = dI2DS(S);
arma::mat dI2DS(const arma::mat &);

// This function returns the derivative of the third invariant of a tensor : I_3 = 1/3 S_ij S_jk S_ki
// @param S The input tensor
// @return The derivative of the third invariant
// @example arma::mat result = dI3DS(S);
arma::mat dI3DS(const arma::mat &);

// This function returns the derivative of the trace of a tensor
// @param S The input tensor
// @return The derivative of the trace
// @example arma::mat result = dtrSdS(S);
arma::mat dtrSdS(const arma::mat &);

// This function returns the derivative of the determinant of a tensor
// @param S The input tensor
// @return The derivative of the determinant
// @example arma::mat result = ddetSdS(S);
arma::mat ddetSdS(const arma::mat &);

// This function returns the derivative of the inverse of a symmetric tensor
// @param S The input tensor
// @return The derivative of the inverse
// @example arma::mat result = dinvSdSsym(S);
arma::mat dinvSdSsym(const arma::mat &);

} //namespace simcoon
