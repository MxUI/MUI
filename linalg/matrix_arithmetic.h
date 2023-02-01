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
 * @date 29 January 2023
 * @brief Implemantation of sparse matrix arithmetic operations.
 */

#ifndef MUI_MATRIX_ARITHMETIC_H_
#define MUI_MATRIX_ARITHMETIC_H_

namespace mui {
namespace linalg {

// Overload addition operator to perform sparse matrix addition
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator+(sparse_matrix<ITYPE,VTYPE> &addend) {

    if (rows != addend.rows || cols != addend.cols) {
        std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix addition" << std::endl;
        std::abort();
    }

    sparse_matrix<ITYPE,VTYPE> res(rows, cols);
    for (auto elememt : matrix) {
        res.set_value(elememt.first.first, elememt.first.second, elememt.second);
    }
    for (auto elememt : addend.matrix) {
        res.set_value(elememt.first.first, elememt.first.second, res.get_value(elememt.first.first, elememt.first.second) + elememt.second);
    }
    return res;
}

// Overload subtraction operator to perform sparse matrix subtraction
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator-(sparse_matrix<ITYPE,VTYPE> &subtrahend) {
   if (rows != subtrahend.rows || cols != subtrahend.cols) {
       std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix subtraction" << std::endl;
       std::abort();
   }
   sparse_matrix<ITYPE,VTYPE> res(rows, cols);
   for (auto elememt : matrix) {
       res.set_value(elememt.first.first, elememt.first.second, elememt.second);
   }
   for (auto elememt : subtrahend.matrix) {
       res.set_value(elememt.first.first, elememt.first.second, res.get_value(elememt.first.first, elememt.first.second) - elememt.second);
   }
   return res;
}

// Overload multiplication operator to perform sparse matrix multiplication
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator*(sparse_matrix<ITYPE,VTYPE> &multiplicand) {
    if (cols != multiplicand.rows) {
        std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix multiplication" << std::endl;
        std::abort();
    }
    sparse_matrix<ITYPE,VTYPE> res(rows, multiplicand.cols);
    for (auto elememt1 : matrix) {
        for (auto elememt2 : multiplicand.matrix) {
            if (elememt1.first.second == elememt2.first.first) {
                res.set_value(elememt1.first.first, elememt2.first.second, res.get_value(elememt1.first.first, elememt2.first.second) + elememt1.second * elememt2.second);
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
            "MUI Error [matrix.h]: scalar type cannot be converted to matrix element type in scalar multiplication");
    sparse_matrix<ITYPE,VTYPE> res(rows,cols);
    for (const auto elememt : matrix) {
        if (static_cast<VTYPE>(scalar) >= std::numeric_limits<VTYPE>::min())
            res.set_value(elememt.first.first, elememt.first.second, elememt.second * static_cast<VTYPE>(scalar));
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
        assert(matrix.empty() &&
                  "MUI Error [matrix.h]: assignment operator '=' only works for empty (all zero elements) matrix");
        assert(((rows == exist_mat.rows) && (cols == exist_mat.cols)) &&
                  "MUI Error [matrix.h]: matrix size mismatch in assignment operator '='");
        (*this).copy(exist_mat);
    }
    return *this;
}

// Member function of dot product
template <typename ITYPE, typename VTYPE>
VTYPE sparse_matrix<ITYPE,VTYPE>::dot_product(sparse_matrix<ITYPE,VTYPE> &exist_mat) const {
    assert(((cols == 1)&&(exist_mat.cols == 1)) &&
        "MUI Error [matrix.h]: dot_product function only works for column vectors");
    sparse_matrix<ITYPE,VTYPE> tempThis(*this);
    sparse_matrix<ITYPE,VTYPE> thisT(tempThis.transpose());
    sparse_matrix<ITYPE,VTYPE> tempMat(thisT * exist_mat);
    assert(((tempMat.get_rows() == 1)&&(tempMat.get_cols() == 1)) &&
                    "MUI Error [matrix.h]: result of dot_product function should be a scalar");
    return (tempMat.get_value(0,0));
}

// Member function of Hadamard product
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::hadamard_product(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    if (rows != exist_mat.rows || cols != exist_mat.cols) {
        std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix Hadamard product" << std::endl;
        std::abort();
    }
    sparse_matrix<ITYPE,VTYPE> res(rows, cols);
    for (auto elememt : matrix) {
        res.matrix[std::make_pair(elememt.first.first, elememt.first.second)] = elememt.second  * exist_mat.get_value(elememt.first.first, elememt.first.second);
    }
    return res;
}

// Member function to get transpose of matrix
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::transpose() {
    sparse_matrix<ITYPE,VTYPE> res(cols, rows);
    for (auto elememt : matrix)
        res.set_value(elememt.first.second, elememt.first.first, elememt.second);
    return res;
}

} // linalg
} // mui

#endif /* MUI_MATRIX_ARITHMETIC_H_ */
