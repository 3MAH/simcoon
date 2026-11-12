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

///@file small_lu.hpp
///@brief LU factorization with partial pivoting of small dense systems, on the stack
///@version 1.0

#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <armadillo>

namespace simcoon{

/** @addtogroup maths
 *  @{
 */

/**
 * @brief LU factorization with partial pivoting of a small dense matrix (at most @p Nmax rows), on
 *        the stack.
 *
 * For the per-point linear systems of the integrators (the 6x6 / 7x7 global Newton Jacobian, the
 * 6x6 condensation of the closest-point return mapping): a LAPACK inverse (getrf + getri) and the
 * products with it cost 0.5 to 0.7 us on a 6x6, far more than the arithmetic, and allocate.
 * factor() once, then solve() for each right-hand side; both run their loops on the compile-time
 * size when n == Nmax (factor_n<Nmax>, solve_n<Nmax>) and on the runtime one otherwise. A pivot below
 * \f$ n\,\epsilon \max|A_{ij}| \f$, or a non-finite entry, reports a singular matrix.
 */
template <int Nmax>
class small_lu {
public:
    /// Factorizes the square matrix @p A (n <= Nmax); false when it is (numerically) singular.
    bool factor(const arma::mat &A) {
        n_ = static_cast<int>(A.n_rows);
        if (n_ > Nmax || static_cast<int>(A.n_cols) != n_) return false;
        return (n_ == Nmax) ? factor_n<Nmax>(A) : factor_n<0>(A);
    }

    /// Solves A x = b in place (b: n values), after a successful factor().
    void solve_in_place(double *b) const {
        if (n_ == Nmax) solve_n<Nmax>(b); else solve_n<0>(b);
    }

    /// X <- A^-1 X column by column, in place (X: n rows), after a successful factor().
    void solve_columns_in_place(arma::mat &X) const {
        for (arma::uword c = 0; c < X.n_cols; ++c) solve_in_place(X.colptr(c));
    }

    /// x = A^-1 b
    arma::vec solve(const arma::vec &b) const {
        arma::vec x(b);
        solve_in_place(x.memptr());
        return x;
    }

    /// X = A^-1 B, column by column
    arma::mat solve(const arma::mat &B) const {
        arma::mat X(B);
        for (arma::uword j = 0; j < X.n_cols; ++j) solve_in_place(X.colptr(j));
        return X;
    }

    int n_rows() const { return n_; }

private:
    /// NS > 0: loops of compile-time length (n == Nmax); NS = 0: the runtime n_
    template <int NS> bool factor_n(const arma::mat &A) {
        const int n = (NS > 0) ? NS : n_;
        double amax = 0.;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                LU_[i][j] = A(i, j);
                if (!std::isfinite(LU_[i][j])) return false;
                amax = std::max(amax, std::fabs(LU_[i][j]));
            }
        }
        const double tiny = amax*n*std::numeric_limits<double>::epsilon();
        for (int k = 0; k < n; ++k) {
            int p = k;
            for (int i = k + 1; i < n; ++i)
                if (std::fabs(LU_[i][k]) > std::fabs(LU_[p][k])) p = i;
            piv_[k] = p;
            if (!(std::fabs(LU_[p][k]) > tiny)) return false;
            if (p != k)
                for (int j = 0; j < n; ++j) std::swap(LU_[k][j], LU_[p][j]);
            const double inv_piv = 1./LU_[k][k];
            for (int i = k + 1; i < n; ++i) {
                const double f = LU_[i][k]*inv_piv;
                LU_[i][k] = f;
                if (f == 0.) continue;
                for (int j = k + 1; j < n; ++j) LU_[i][j] -= f*LU_[k][j];
            }
        }
        return true;
    }

    /// NS > 0: loops of compile-time length (n == Nmax); NS = 0: the runtime n_
    template <int NS> void solve_n(double *b) const {
        const int n = (NS > 0) ? NS : n_;
        for (int k = 0; k < n; ++k)
            if (piv_[k] != k) std::swap(b[k], b[piv_[k]]);
        for (int i = 1; i < n; ++i) {
            double v = b[i];
            for (int j = 0; j < i; ++j) v -= LU_[i][j]*b[j];
            b[i] = v;
        }
        for (int i = n - 1; i >= 0; --i) {
            double v = b[i];
            for (int j = i + 1; j < n; ++j) v -= LU_[i][j]*b[j];
            b[i] = v/LU_[i][i];
        }
    }

    double LU_[Nmax][Nmax];
    int piv_[Nmax];
    int n_ = 0;
};

/** @} */ // end of maths group

} //namespace simcoon
