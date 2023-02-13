/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2023 W. Liu                                                  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file matrix_arithmetic.h
 * @author W. Liu
 * @date 01 February 2023
 * @brief Implementation of sparse matrix arithmetic operations.
 */

#ifndef MUI_MATRIX_ARITHMETIC_H_
#define MUI_MATRIX_ARITHMETIC_H_

#include <cassert>
#include <math.h>

namespace mui {
namespace linalg {

// Overload addition operator to perform sparse matrix addition
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator+(sparse_matrix<ITYPE,VTYPE> &addend) const{

    if (rows_ != addend.rows_ || cols_ != addend.cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: matrix size mismatch during matrix addition" << std::endl;
        std::abort();
    }

    sparse_matrix<ITYPE,VTYPE> res(rows_, cols_);
    for (auto element : matrix_) {
        res.set_value(element.first.first, element.first.second, element.second);
    }
    for (auto element : addend.matrix_) {
        res.set_value(element.first.first, element.first.second, res.get_value(element.first.first, element.first.second) + element.second);
    }
    return res;
}

// Overload subtraction operator to perform sparse matrix subtraction
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator-(sparse_matrix<ITYPE,VTYPE> &subtrahend) const {
   if (rows_ != subtrahend.rows_ || cols_ != subtrahend.cols_) {
       std::cerr << "MUI Error [matrix_arithmetic.h]: matrix size mismatch during matrix subtraction" << std::endl;
       std::abort();
   }
   sparse_matrix<ITYPE,VTYPE> res(rows_, cols_);
   for (auto element : matrix_) {
       res.set_value(element.first.first, element.first.second, element.second);
   }
   for (auto element : subtrahend.matrix_) {
       res.set_value(element.first.first, element.first.second, res.get_value(element.first.first, element.first.second) - element.second);
   }
   return res;
}

// Overload multiplication operator to perform sparse matrix multiplication
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator*(sparse_matrix<ITYPE,VTYPE> &multiplicand) const {
    if (cols_ != multiplicand.rows_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: matrix size mismatch during matrix multiplication" << std::endl;
        std::abort();
    }
    sparse_matrix<ITYPE,VTYPE> res(rows_, multiplicand.cols_);
    for (auto element1 : matrix_) {
        for (auto element2 : multiplicand.matrix_) {
            if (element1.first.second == element2.first.first) {
                res.set_value(element1.first.first, element2.first.second, res.get_value(element1.first.first, element2.first.second) + element1.second * element2.second);
            }
        }
    }
    return res;
}

// Overload multiplication operator to perform scalar multiplication A*x
template <typename ITYPE, typename VTYPE>
template <typename STYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator*(const STYPE &scalar) const {
    static_assert(std::is_convertible<STYPE, VTYPE>::value,
            "MUI Error [matrix_arithmetic.h]: scalar type cannot be converted to matrix element type in scalar multiplication");
    sparse_matrix<ITYPE,VTYPE> res(rows_,cols_);
    for (const auto element : matrix_) {
        if (static_cast<VTYPE>(scalar) >= std::numeric_limits<VTYPE>::min())
            res.set_value(element.first.first, element.first.second, element.second * static_cast<VTYPE>(scalar));
   }
   return res;
}

// Overload multiplication operator to perform scalar multiplication x*A
template<typename ITYPE, typename VTYPE, typename STYPE>
sparse_matrix<ITYPE,VTYPE> operator*(const STYPE &scalar, const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
  return exist_mat * scalar;
}

// Member function of dot product
template <typename ITYPE, typename VTYPE>
VTYPE sparse_matrix<ITYPE,VTYPE>::dot_product(sparse_matrix<ITYPE,VTYPE> &exist_mat) const {
    assert(((cols_ == 1)&&(exist_mat.cols_ == 1)) &&
        "MUI Error [matrix_arithmetic.h]: dot_product function only works for column vectors");
    sparse_matrix<ITYPE,VTYPE> tempThis(*this);
    sparse_matrix<ITYPE,VTYPE> thisT(tempThis.transpose());
    sparse_matrix<ITYPE,VTYPE> tempMat(thisT * exist_mat);
    assert(((tempMat.get_rows() == 1)&&(tempMat.get_cols() == 1)) &&
                    "MUI Error [matrix_arithmetic.h]: result of dot_product function should be a scalar");
    return (tempMat.get_value(0,0));
}

// Member function of Hadamard product
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::hadamard_product(const sparse_matrix<ITYPE,VTYPE> &exist_mat) const {
    if (rows_ != exist_mat.rows_ || cols_ != exist_mat.cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: matrix size mismatch during matrix Hadamard product" << std::endl;
        std::abort();
    }
    sparse_matrix<ITYPE,VTYPE> res(rows_, cols_);
    for (auto element : matrix_) {
        res.matrix_[std::make_pair(element.first.first, element.first.second)] = element.second  * exist_mat.get_value(element.first.first, element.first.second);
    }
    return res;
}

// Member function to get transpose of matrix
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::transpose() const {
    sparse_matrix<ITYPE,VTYPE> res(cols_, rows_);
    for (auto element : matrix_)
        res.set_value(element.first.second, element.first.first, element.second);
    return res;
}

// Member function to perform LU decomposition
template <typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::lu_decomposition(sparse_matrix<ITYPE,VTYPE> &L, sparse_matrix<ITYPE,VTYPE> &U) const {
    if ((L.get_rows() != 0) || (U.get_rows() != 0) || (L.get_cols() != 0) || (U.get_cols() != 0)) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: L & U Matrices must be null in LU decomposition" << std::endl;
        std::abort();
    }
    if (rows_ != cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Only square matrix can perform LU decomposition" << std::endl;
        std::abort();
    }

    // Resize the lower triangular matrix
    L.resize_null(rows_, cols_);
    // Resize the upper triangular matrix
    U.resize_null(rows_, cols_);

    ITYPE n = rows_;
    for (ITYPE i = 0; i < rows_; ++i) {
        // Calculate the upper triangular matrix
        for (ITYPE k = i; k < cols_; ++k) {
            VTYPE sum = 0.0;
            for (ITYPE j = 0; j < i; ++j) {
                sum += L.get_value(i, j) * U.get_value(j, k);
            }
            U.set_value(i, k, (this->get_value(i, k) - sum));
        }

        // Calculate the lower triangular matrix
        for (ITYPE k = i; k < rows_; k++) {
            if (i == k) {
                L.set_value(i, i, static_cast<VTYPE>(1.0));
            } else {
                VTYPE sum = 0.0;
                for (ITYPE j = 0; j < i; ++j) {
                    sum += L.get_value(k, j) * U.get_value(j, i);
                }
                assert(std::abs(U.get_value(i, i)) >= std::numeric_limits<VTYPE>::min() &&
                                  "MUI Error [matrix_arithmetic.h]: Divide by zero assert for U.get_value(i, i)");
                L.set_value(k, i, (this->get_value(k, i) - sum) / U.get_value(i, i));
            }
        }
    }
}

// Member function to perform QR decomposition
template <typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::qr_decomposition(sparse_matrix<ITYPE,VTYPE> &Q, sparse_matrix<ITYPE,VTYPE> &R) const {
    if ((Q.get_rows() != 0) || (R.get_rows() != 0) || (Q.get_cols() != 0) || (R.get_cols() != 0)) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Q & R Matrices must be null in QR decomposition" << std::endl;
        std::abort();
    }
    assert((rows_ >= cols_) &&
          "MUI Error [matrix_arithmetic.h]: number of rows of matrix should larger or equals to number of columns in QR decomposition");
    // Resize the orthogonal matrix
    Q.resize_null(rows_, cols_);
    // Resize the upper triangular matrix
    R.resize_null(rows_, cols_);
    // Get a copy of the matrix
    sparse_matrix<ITYPE,VTYPE> mat_copy (*this);
    // Diagonal elements
    std::vector<VTYPE> r_diag (cols_);

    // Calculate the diagonal element values
    for (ITYPE c = 0; c <cols_; ++c)  {
        VTYPE  nrm (0.0);

       // Compute 2-norm of k-th column without under/overflow.
        for (ITYPE r = c; r < rows_; ++r)
            nrm = std::sqrt((nrm * nrm) + (mat_copy.matrix_[std::make_pair(r, c)] * mat_copy.matrix_[std::make_pair(r, c)]));

        if (nrm != static_cast<VTYPE>(0.0))  {

           // Form k-th Householder vector.
            if (mat_copy.matrix_[std::make_pair(c, c)] < static_cast<VTYPE>(0.0))
                nrm = -nrm;

            for (ITYPE r = c; r < rows_; ++r)
                mat_copy.matrix_[std::make_pair(r, c)] /= nrm;

            mat_copy.matrix_[std::make_pair(c, c)]  += static_cast<VTYPE>(1.0);

           // Apply transformation to remaining columns.
            for (ITYPE j = c + 1; j < cols_; ++j)  {
                VTYPE  s = 0.0;

                for (ITYPE r = c; r < rows_; ++r)
                    s += mat_copy.matrix_[std::make_pair(r, c)]  * mat_copy.matrix_[std::make_pair(r, j)];

                s /= -mat_copy.matrix_[std::make_pair(c, c)];
                for (ITYPE r = c; r < rows_; ++r)
                    mat_copy.matrix_[std::make_pair(r, j)]  += s * mat_copy.matrix_[std::make_pair(r, c)];
            }
        }
        r_diag[c] = -nrm;
    }

    // Calculate the orthogonal matrix
    for (ITYPE c = cols_ - 1; c >= 0; --c)  {
        Q.set_value(c, c, static_cast<VTYPE>(1.0));

        for (ITYPE cc = c; cc < cols_; ++cc)
            if (mat_copy.matrix_[std::make_pair(c, c)]  != static_cast<VTYPE>(0.0))  {
                VTYPE  s=0.0;

                for (ITYPE r = c; r < rows_; ++r)
                    s += mat_copy.matrix_[std::make_pair(r, c)]  * Q.get_value(r, cc);

                s /= -mat_copy.matrix_[std::make_pair(c, c)];
                for (ITYPE r = c; r < rows_; ++r)
                    Q.set_value(r, cc, (Q.get_value(r, cc) + s * mat_copy.matrix_[std::make_pair(r, c)]));
            }
    }

    // Calculate the upper triangular matrix
    for (ITYPE c = 0; c < cols_; ++c)
        for (ITYPE r = 0; r < rows_; ++r)
            if (c < r)
                R.set_value(c, r, mat_copy.matrix_[std::make_pair(c, r)]);
            else if (c == r)
                R.set_value(c, r, r_diag[c]);
}

// Member function to get the inverse of matrix by using Gaussian elimination
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::inverse() const {
    if (rows_ != cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Matrix must be square to find its inverse" << std::endl;
        std::abort();
    }

    sparse_matrix<ITYPE,VTYPE> mat_copy (*this);
    sparse_matrix<ITYPE,VTYPE> inverse_mat (rows_,"identity");

    for (ITYPE r = 0; r < rows_; ++r)  {

        ITYPE max_row = r;
        VTYPE max_value= static_cast<VTYPE>(-1.0);

        // Partial pivoting for Gaussian elimination
        ITYPE ppivot;
        for (ITYPE rb = r; rb < rows_; ++rb)  {
            const VTYPE tmp = std::abs(mat_copy.matrix_[std::make_pair(rb, r)]);

            if ((tmp > max_value) && (std::abs(tmp) >= std::numeric_limits<VTYPE>::min()))  {
                max_value = tmp;
                max_row = rb;
            }
        }

        assert(std::abs(mat_copy.matrix_[std::make_pair(max_row, r)]) >= std::numeric_limits<VTYPE>::min() &&
                          "MUI Error [matrix_arithmetic.h]: Divide by zero assert for mat_copy.matrix_[std::make_pair(max_row, r)]. Cannot perform matrix invert due to singular matrix.");

        if (max_row != r)  {
            for (ITYPE c = 0; c < cols_; ++c)
                std::swap (mat_copy.matrix_[std::make_pair(r, c)] , mat_copy.matrix_[std::make_pair(max_row, c)]);
            ppivot = max_row;
        } else {
            ppivot = 0;
        }

        const ITYPE indx = ppivot;

        if (indx != 0)
            for (ITYPE c = 0; c < cols_; ++c)
                std::swap (inverse_mat.matrix_[std::make_pair(r, c)] , inverse_mat.matrix_[std::make_pair(indx, c)]);

        const VTYPE diag = mat_copy.matrix_[std::make_pair(r, r)];

        for (ITYPE c = 0; c < cols_; ++c)  {
            mat_copy.matrix_[std::make_pair(r, c)]  /= diag;
            inverse_mat.matrix_[std::make_pair(r, c)] /= diag;
        }

        for (ITYPE rr = 0; rr < rows_; ++rr)
            if (rr != r)  {
                const VTYPE off_diag = mat_copy.matrix_[std::make_pair(rr, r)];

                for (ITYPE c = 0; c < cols_; ++c)  {
                    mat_copy.matrix_[std::make_pair(rr, c)]  -= off_diag * mat_copy.matrix_[std::make_pair(r, c)];
                    inverse_mat.matrix_[std::make_pair(rr, c)]  -= off_diag * inverse_mat.matrix_[std::make_pair(r, c)];
                }
            }
    }
    return inverse_mat;
}

} // linalg
} // mui

#endif /* MUI_MATRIX_ARITHMETIC_H_ */
