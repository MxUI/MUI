/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                    W. Liu                                                  *
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
 * @brief Implemantation of sparse matrix arithmetic operations.
 */

#ifndef MUI_MATRIX_ARITHMETIC_H_
#define MUI_MATRIX_ARITHMETIC_H_

#include <cassert>
#include <math.h>

namespace mui {
namespace linalg {

// Overload addition operator to perform sparse matrix addition
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator+(sparse_matrix<ITYPE,VTYPE> &addend) {

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
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator-(sparse_matrix<ITYPE,VTYPE> &subtrahend) {
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
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator*(sparse_matrix<ITYPE,VTYPE> &multiplicand) {
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

// Overloaded assignment operator
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>& sparse_matrix<ITYPE,VTYPE>::operator=(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    if (this != &exist_mat) { // check for self-assignment
        // copy the values from the other matrix to this matrix
        assert(matrix_.empty() &&
                  "MUI Error [matrix_arithmetic.h]: assignment operator '=' only works for empty (all zero elements) matrix");
        assert(((rows_ == exist_mat.rows_) && (cols_ == exist_mat.cols_)) &&
                  "MUI Error [matrix_arithmetic.h]: matrix size mismatch in assignment operator '='");
        (*this).copy(exist_mat);
    }
    return *this;
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
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::hadamard_product(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
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
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::transpose() {
    sparse_matrix<ITYPE,VTYPE> res(cols_, rows_);
    for (auto element : matrix_)
        res.set_value(element.first.second, element.first.first, element.second);
    return res;
}

// Member function to perform LU decomposition
template <typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::lu_decomposition(sparse_matrix<ITYPE,VTYPE> &L, sparse_matrix<ITYPE,VTYPE> &U) {
    if ((L.get_rows() != 0) || (U.get_rows() != 0) || (L.get_cols() != 0) || (U.get_cols() != 0)) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: L & U Matrices must be null in LU decomposition" << std::endl;
        std::abort();
    }
    if (rows_ != cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Only square matrix can perform LU decomposition" << std::endl;
        std::abort();
    }

    L.resize_null(rows_, cols_);
    U.resize_null(rows_, cols_);

    ITYPE n = rows_;
    for (ITYPE i = 0; i < rows_; ++i) {
        // calculate the upper triangular matrix
        for (ITYPE k = i; k < cols_; ++k) {
            VTYPE sum = 0.0;
            for (ITYPE j = 0; j < i; ++j) {
                sum += L.get_value(i, j) * U.get_value(j, k);
            }
            U.set_value(i, k, matrix_[std::make_pair(i, k)] - sum);
        }

        // calculate the lower triangular matrix
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
                L.set_value(k, i, (matrix_[std::make_pair(k, i)] - sum) / U.get_value(i, i));
            }
        }
    }
}

// Member function to perform QR decomposition
template <typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::qr_decomposition(sparse_matrix<ITYPE,VTYPE> &Q, sparse_matrix<ITYPE,VTYPE> &R) {
    if ((Q.get_rows() != 0) || (R.get_rows() != 0) || (Q.get_cols() != 0) || (R.get_cols() != 0)) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Q & R Matrices must be null in QR decomposition" << std::endl;
        std::abort();
    }
    assert((rows_ >= cols_) &&
          "MUI Error [matrix_arithmetic.h]: number of rows of matrix should larger or equals to number of columns in QR decomposition");

    Q.resize_null(rows_, cols_);
    R.resize_null(rows_, cols_);

    std::vector<VTYPE> r_diag (cols_);

    for (ITYPE c = 0; c <cols_; ++c)  {
        VTYPE  nrm (0.0);

       // Compute 2-norm of k-th column without under/overflow.
        for (ITYPE r = c; r < rows_; ++r)
            nrm = std::sqrt((nrm * nrm) + (matrix_[std::make_pair(r, c)] * matrix_[std::make_pair(r, c)]));

        if (nrm != static_cast<VTYPE>(0.0))  {

           // Form k-th Householder vector.
            if (matrix_[std::make_pair(c, c)] < static_cast<VTYPE>(0.0))
                nrm = -nrm;

            for (ITYPE r = c; r < rows_; ++r)
                matrix_[std::make_pair(r, c)] /= nrm;

            matrix_[std::make_pair(c, c)]  += static_cast<VTYPE>(1.0);

           // Apply transformation to remaining columns.
            for (ITYPE j = c + 1; j < cols_; ++j)  {
                VTYPE  s = 0.0;

                for (ITYPE r = c; r < rows_; ++r)
                    s += matrix_[std::make_pair(r, c)]  * matrix_[std::make_pair(r, j)];

                s /= -matrix_[std::make_pair(c, c)];
                for (ITYPE r = c; r < rows_; ++r)
                    matrix_[std::make_pair(r, j)]  += s * matrix_[std::make_pair(r, c)];
            }
        }
        r_diag[c] = -nrm;
    }

    for (ITYPE c = cols_ - 1; c >= 0; --c)  {
        Q.set_value(c, c, static_cast<VTYPE>(1.0));

        for (ITYPE cc = c; cc < cols_; ++cc)
            if (matrix_[std::make_pair(c, c)]  != static_cast<VTYPE>(0.0))  {
                VTYPE  s=0.0;

                for (ITYPE r = c; r < rows_; ++r)
                    s += matrix_[std::make_pair(r, c)]  * Q.get_value(r, cc);

                s /= -matrix_[std::make_pair(c, c)];
                for (ITYPE r = c; r < rows_; ++r)
                    Q.set_value(r, cc, (Q.get_value(r, cc) + s * matrix_[std::make_pair(r, c)]));
            }
    }

    for (ITYPE c = 0; c < cols_; ++c)
        for (ITYPE r = 0; r < rows_; ++r)
            if (c < r)
                R.set_value(c, r, matrix_[std::make_pair(c, r)]);
            else if (c == r)
                R.set_value(c, r, r_diag[c]);
}

// Member function to get the inverse of matrix by using Gaussian elimination
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::inverse() {
    if (rows_ != cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Matrix must be square to find its inverse" << std::endl;
        std::abort();
    }
    std::cout<< " ********* matrix: " << rows_ << " " << cols_ << std::endl;
    ITYPE n = rows_;
    sparse_matrix<ITYPE,VTYPE> inverse_mat(n, n);
    sparse_matrix<ITYPE,VTYPE> augmented_mat(n, n);
    std::cout<< " ********* matrix: " << n << std::endl;
    // Create a map to store the augmented matrix
    std::map<std::pair<ITYPE, ITYPE>, VTYPE> augmented_matrix;
    for (const auto& element : matrix_) {
      if (element.first.first == element.first.second) {
        augmented_mat.matrix_[std::make_pair(element.first.first,element.first.second)] = element.second;

      } else {
        augmented_mat.matrix_[std::make_pair(element.first.first,element.first.second)] = static_cast<VTYPE>(0.0);
      }
      augmented_mat.matrix_[std::make_pair(element.first.first,(n+element.first.second))] =
          (element.first.first == element.first.second) ? static_cast<VTYPE>(1.0) : static_cast<VTYPE>(0.0);
    }

    std::cout<< " ********* augmented_matrix: " << std::endl;
    for (const auto element : augmented_mat.matrix_) {
        std::cout<< "       " << element.first.first << " " << element.first.second << " " << element.second << std::endl;
    }
    std::cout<< " ********* augmented_matrix finished " << std::endl;

    // Gaussian elimination
    for (ITYPE i = 0; i < n; ++i) {
        // Find the maximum value in column i
        ITYPE max_row = i;
        VTYPE max_val = std::fabs(augmented_mat.matrix_[std::make_pair(i, i)]);
      for (ITYPE j = i+1; j < n; ++j) {
          VTYPE val = std::fabs(augmented_mat.matrix_[std::make_pair(j, i)]);
        if (val > max_val) {
          max_val = val;
          max_row = j;
        }
      }

      // Swap the current row with the row containing the maximum value
      if (max_row != i) {
        for (ITYPE j = i; j <= n; ++j) {
          std::swap(augmented_mat.matrix_[std::make_pair(i, j)], augmented_mat.matrix_[std::make_pair(max_row, j)]);
        }
      }

      // Zero out the entries below the pivot
      VTYPE pivot = augmented_mat.matrix_[std::make_pair(i, i)];
      for (ITYPE j = i+1; j < n; ++j) {
          assert(std::abs(pivot) >= std::numeric_limits<VTYPE>::min() &&
                  "MUI Error [matrix_arithmetic.h]: Divide by zero assert for pivot. Perhaps matrix is singular and inverse of the matrix cannot be computed");
          VTYPE ratio = augmented_mat.matrix_[std::make_pair(j, i)] / pivot;
        for (ITYPE k = i; k <= n; ++k) {
          augmented_mat.matrix_[std::make_pair(j, k)] -= ratio * augmented_mat.matrix_[std::make_pair(i, k)];
        }
      }
    }

    // Back substitution
    for (ITYPE i = n - 1; i >= 0; --i) {
        VTYPE pivot = augmented_mat.matrix_[std::make_pair(i, i)];
      for (ITYPE j = i - 1; j >= 0; --j) {
          assert(std::abs(pivot) >= std::numeric_limits<VTYPE>::min() &&
                  "MUI Error [matrix_arithmetic.h]: Divide by zero assert for pivot. Perhaps matrix is singular and inverse of the matrix cannot be computed");
          VTYPE ratio = augmented_mat.matrix_[std::make_pair(j, i)] / pivot;
        for (ITYPE k = i; k <= n; ++k) {
          augmented_mat.matrix_[std::make_pair(j, k)] -= ratio * augmented_mat.matrix_[std::make_pair(i, k)];
        }
      }
    }

    // Store the inverse in the sparse matrix
    for (ITYPE i = 0; i < n; ++i) {
        for (ITYPE j = 0; j < n; ++j) {
            inverse_mat.set_value(i, j, augmented_mat.matrix_[std::make_pair(i, j)]);
        }
    }

    return inverse_mat;
}

} // linalg
} // mui

#endif /* MUI_MATRIX_ARITHMETIC_H_ */
