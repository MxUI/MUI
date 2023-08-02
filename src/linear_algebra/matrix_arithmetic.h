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

// **************************************************
// ************ Public member functions *************
// **************************************************

// Overload addition operator to perform sparse matrix addition
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator+(sparse_matrix<ITYPE,VTYPE> &addend) {

    if (rows_ != addend.rows_ || cols_ != addend.cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: matrix size mismatch during matrix addition" << std::endl;
        std::abort();
    }

    if (addend.matrix_format_ != matrix_format_) {
        addend.format_conversion(this->get_format(), true, true, "overwrite");
    } else {
        if (!addend.is_sorted_unique("matrix_arithmetic.h", "operator+()")){
            if (addend.matrix_format_ == format::COO) {
                addend.sort_coo(true, true, "overwrite");
            } else if (addend.matrix_format_ == format::CSR) {
                addend.sort_csr(true, "overwrite");
            } else if (addend.matrix_format_ == format::CSC) {
                addend.sort_csc(true, "overwrite");
            } else {
                std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised addend matrix format for matrix operator+()" << std::endl;
                std::cerr << "    Please set the addend matrix_format_ as:" << std::endl;
                std::cerr << "    format::COO: COOrdinate format" << std::endl;
                std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                std::abort();
            }
        }
    }

    if (!this->is_sorted_unique("matrix_arithmetic.h", "operator+()")){
        if (matrix_format_ == format::COO) {
            this->sort_coo(true, true, "overwrite");
        } else if (matrix_format_ == format::CSR) {
            this->sort_csr(true, "overwrite");
        } else if (matrix_format_ == format::CSC) {
            this->sort_csc(true, "overwrite");
        } else {
            std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix operator+()" << std::endl;
            std::cerr << "    Please set the matrix_format_ as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }
    }

    // Create a new sparse matrix object for the result
    sparse_matrix<ITYPE,VTYPE> res(rows_, cols_, this->get_format());

    if (matrix_format_ == format::COO) {

        // Perform element-wise addition of the COO vectors
        res.matrix_coo.values_.reserve(matrix_coo.values_.size() + addend.matrix_coo.values_.size());
        res.matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size() + addend.matrix_coo.row_indices_.size());
        res.matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size() + addend.matrix_coo.col_indices_.size());

        // Insert the COO vectors of the initial sparse matrix to the result sparse matrix
        res.matrix_coo.values_ = std::vector<VTYPE>(matrix_coo.values_.begin(), matrix_coo.values_.end());
        res.matrix_coo.row_indices_ = std::vector<ITYPE>(matrix_coo.row_indices_.begin(), matrix_coo.row_indices_.end());
        res.matrix_coo.col_indices_ = std::vector<ITYPE>(matrix_coo.col_indices_.begin(), matrix_coo.col_indices_.end());

        // Append the addend COO vectors to the result sparse matrix
        res.matrix_coo.values_.insert(res.matrix_coo.values_.end(), addend.matrix_coo.values_.begin(), addend.matrix_coo.values_.end());
        res.matrix_coo.row_indices_.insert(res.matrix_coo.row_indices_.end(), addend.matrix_coo.row_indices_.begin(), addend.matrix_coo.row_indices_.end());
        res.matrix_coo.col_indices_.insert(res.matrix_coo.col_indices_.end(), addend.matrix_coo.col_indices_.begin(), addend.matrix_coo.col_indices_.end());

        // Sort and deduplicate the result
        res.sort_coo(true, true, "plus");
        res.nnz_ = res.matrix_coo.values_.size();

    } else if (matrix_format_ == format::CSR) {

        // Perform element-wise addition of the CSR vectors
        res.matrix_csr.values_.reserve(matrix_csr.values_.size() + addend.matrix_csr.values_.size());
        res.matrix_csr.row_ptrs_.resize(rows_ + 1);
        res.matrix_csr.col_indices_.reserve(matrix_csr.col_indices_.size() + addend.matrix_csr.col_indices_.size());

        ITYPE row = 0;
        while (row < rows_) {
            ITYPE start = matrix_csr.row_ptrs_[row];
            ITYPE end = matrix_csr.row_ptrs_[row + 1];

            ITYPE addend_start = addend.matrix_csr.row_ptrs_[row];
            ITYPE addend_end = addend.matrix_csr.row_ptrs_[row + 1];

            res.matrix_csr.row_ptrs_[0] = 0;

            // Merge the values and column indices of the two rows
            ITYPE i = start;
            ITYPE j = addend_start;
            while (i < end && j < addend_end) {
                ITYPE col = matrix_csr.col_indices_[i];
                ITYPE addend_col = addend.matrix_csr.col_indices_[j];

                if (col == addend_col) {
                    // Add the corresponding values if the columns match
                    if (std::abs(matrix_csr.values_[i] + addend.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()){
                        res.matrix_csr.values_.emplace_back(matrix_csr.values_[i] + addend.matrix_csr.values_[j]);
                        res.matrix_csr.col_indices_.emplace_back(col);
                    }
                    i++;
                    j++;
                } else if (col < addend_col) {
                    // Add the current value from the initial matrix
                    if (std::abs(matrix_csr.values_[i]) >= std::numeric_limits<VTYPE>::min()){
                        res.matrix_csr.values_.emplace_back(matrix_csr.values_[i]);
                        res.matrix_csr.col_indices_.emplace_back(col);
                    }
                    i++;
                } else {
                    // Add the current value from the addend matrix
                    if (std::abs(addend.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()){
                        res.matrix_csr.values_.emplace_back(addend.matrix_csr.values_[j]);
                        res.matrix_csr.col_indices_.emplace_back(addend_col);
                    }
                    j++;
                }
            }

            // Add any remaining elements from the initial matrix
            for (; i < end; i++) {
                if (std::abs(matrix_csr.values_[i]) >= std::numeric_limits<VTYPE>::min()){
                    res.matrix_csr.values_.emplace_back(matrix_csr.values_[i]);
                    res.matrix_csr.col_indices_.emplace_back(matrix_csr.col_indices_[i]);
                }
            }

            // Add any remaining elements from the addend matrix
            for (; j < addend_end; j++) {
                if (std::abs(addend.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()){
                    res.matrix_csr.values_.emplace_back(addend.matrix_csr.values_[j]);
                    res.matrix_csr.col_indices_.emplace_back(addend.matrix_csr.col_indices_[j]);
                }
            }

            // Update the row pointer
            res.nnz_ = res.matrix_csr.col_indices_.size();
            res.matrix_csr.row_ptrs_[row + 1] = res.nnz_;

            row++;
        }

    } else if (matrix_format_ == format::CSC) {

        // Perform element-wise addition of the CSC vectors
        res.matrix_csc.values_.reserve(matrix_csc.values_.size() + addend.matrix_csc.values_.size());
        res.matrix_csc.row_indices_.reserve(matrix_csc.row_indices_.size() + addend.matrix_csc.row_indices_.size());
        res.matrix_csc.col_ptrs_.resize(cols_ + 1);

        ITYPE column = 0;
        while (column < cols_) {
            ITYPE start = matrix_csc.col_ptrs_[column];
            ITYPE end = matrix_csc.col_ptrs_[column + 1];

            ITYPE addend_start = addend.matrix_csc.col_ptrs_[column];
            ITYPE addend_end = addend.matrix_csc.col_ptrs_[column + 1];

            res.matrix_csc.col_ptrs_[0] = 0;

            // Merge the values and row indices of the two columns
            ITYPE i = start;
            ITYPE j = addend_start;
            while (i < end && j < addend_end) {
                ITYPE row = matrix_csc.row_indices_[i];
                ITYPE addend_row = addend.matrix_csc.row_indices_[j];

                if (row == addend_row) {
                    // Add the corresponding values if the columns match
                    if (std::abs(matrix_csc.values_[i] + addend.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csc.values_.emplace_back(matrix_csc.values_[i] + addend.matrix_csc.values_[j]);
                        res.matrix_csc.row_indices_.emplace_back(row);
                    }
                    i++;
                    j++;
                } else if (row < addend_row) {
                    // Add the current value from the initial matrix
                    if (std::abs(matrix_csc.values_[i]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csc.values_.emplace_back(matrix_csc.values_[i]);
                        res.matrix_csc.row_indices_.emplace_back(row);
                    }
                    i++;
                } else {
                    // Add the current value from the addend matrix
                    if (std::abs(addend.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csc.values_.emplace_back(addend.matrix_csc.values_[j]);
                        res.matrix_csc.row_indices_.emplace_back(addend_row);
                    }
                    j++;
                }
            }

            // Add any remaining elements from the initial matrix
            for (; i < end; i++) {
                if (std::abs(matrix_csc.values_[i]) >= std::numeric_limits<VTYPE>::min()){
                    res.matrix_csc.values_.emplace_back(matrix_csc.values_[i]);
                    res.matrix_csc.row_indices_.emplace_back(matrix_csc.row_indices_[i]);
                }
            }

            // Add any remaining elements from the addend matrix
            for (; j < addend_end; j++) {
                if (std::abs(addend.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()){
                    res.matrix_csc.values_.emplace_back(addend.matrix_csc.values_[j]);
                    res.matrix_csc.row_indices_.emplace_back(addend.matrix_csc.row_indices_[j]);
                }
            }

            // Update the column pointer
            res.nnz_ = res.matrix_csc.row_indices_.size();
            res.matrix_csc.col_ptrs_[column + 1] = res.nnz_;

            column++;
        }

    } else {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix operator+()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
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

   if (subtrahend.matrix_format_ != matrix_format_) {
       subtrahend.format_conversion(this->get_format(), true, true, "overwrite");
   } else {
       if (!subtrahend.is_sorted_unique("matrix_arithmetic.h", "operator-()")){
           if (subtrahend.matrix_format_ == format::COO) {
               subtrahend.sort_coo(true, true, "overwrite");
           } else if (subtrahend.matrix_format_ == format::CSR) {
               subtrahend.sort_csr(true, "overwrite");
           } else if (subtrahend.matrix_format_ == format::CSC) {
               subtrahend.sort_csc(true, "overwrite");
           } else {
               std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised subtrahend matrix format for matrix operator-()" << std::endl;
               std::cerr << "    Please set the subtrahend matrix_format_ as:" << std::endl;
               std::cerr << "    format::COO: COOrdinate format" << std::endl;
               std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
               std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
               std::abort();
           }
       }
   }

    if (!this->is_sorted_unique("matrix_arithmetic.h", "operator-()")){
        if (matrix_format_ == format::COO) {
            this->sort_coo(true, true, "overwrite");
        } else if (matrix_format_ == format::CSR) {
            this->sort_csr(true, "overwrite");
        } else if (matrix_format_ == format::CSC) {
            this->sort_csc(true, "overwrite");
        } else {
            std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix operator-()" << std::endl;
            std::cerr << "    Please set the matrix_format_ as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }
    }

   // Create a new sparse matrix object for the result
   sparse_matrix<ITYPE,VTYPE> res(rows_, cols_, this->get_format());

   if (matrix_format_ == format::COO) {

       std::vector<VTYPE> subtrahend_value;
       subtrahend_value.reserve(subtrahend.matrix_coo.values_.size());

        for (VTYPE &element : subtrahend.matrix_coo.values_) {
            subtrahend_value.emplace_back(element*(-1));
        }

        // Perform element-wise subtrahend of the COO vectors
        res.matrix_coo.values_.reserve(matrix_coo.values_.size() + subtrahend.matrix_coo.values_.size());
        res.matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size() + subtrahend.matrix_coo.row_indices_.size());
        res.matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size() + subtrahend.matrix_coo.col_indices_.size());

       // Insert the COO vectors of the initial sparse matrix to the result sparse matrix
       res.matrix_coo.values_ = std::vector<VTYPE>(matrix_coo.values_.begin(), matrix_coo.values_.end());
       res.matrix_coo.row_indices_ = std::vector<ITYPE>(matrix_coo.row_indices_.begin(), matrix_coo.row_indices_.end());
       res.matrix_coo.col_indices_ = std::vector<ITYPE>(matrix_coo.col_indices_.begin(), matrix_coo.col_indices_.end());

       // Append the subtrahend COO vectors to the result sparse matrix
       res.matrix_coo.values_.insert(res.matrix_coo.values_.end(), subtrahend_value.begin(), subtrahend_value.end());
       res.matrix_coo.row_indices_.insert(res.matrix_coo.row_indices_.end(), subtrahend.matrix_coo.row_indices_.begin(), subtrahend.matrix_coo.row_indices_.end());
       res.matrix_coo.col_indices_.insert(res.matrix_coo.col_indices_.end(), subtrahend.matrix_coo.col_indices_.begin(), subtrahend.matrix_coo.col_indices_.end());

       // Sort and deduplicate the result
       res.sort_coo(true, true, "plus");

    } else if (matrix_format_ == format::CSR) {

        // Perform element-wise subtraction of the CSR vectors
        res.matrix_csr.values_.reserve(matrix_csr.values_.size() + subtrahend.matrix_csr.values_.size());
        res.matrix_csr.row_ptrs_.resize(rows_ + 1);
        res.matrix_csr.col_indices_.reserve(matrix_csr.col_indices_.size() + subtrahend.matrix_csr.col_indices_.size());

        ITYPE row = 0;
        while (row < rows_) {
            ITYPE start = matrix_csr.row_ptrs_[row];
            ITYPE end = matrix_csr.row_ptrs_[row + 1];

            ITYPE subtrahend_start = subtrahend.matrix_csr.row_ptrs_[row];
            ITYPE subtrahend_end = subtrahend.matrix_csr.row_ptrs_[row + 1];

            res.matrix_csr.row_ptrs_[0] = 0;

            // Merge the values and column indices of the two rows
            ITYPE i = start;
            ITYPE j = subtrahend_start;
            while (i < end && j < subtrahend_end) {
                ITYPE col = matrix_csr.col_indices_[i];
                ITYPE subtrahend_col = subtrahend.matrix_csr.col_indices_[j];

                if (col == subtrahend_col) {
                    // Add the corresponding values if the columns match
                    if (std::abs(matrix_csr.values_[i] - subtrahend.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csr.values_.emplace_back(matrix_csr.values_[i] - subtrahend.matrix_csr.values_[j]);
                        res.matrix_csr.col_indices_.emplace_back(col);
                    }
                    i++;
                    j++;
                } else if (col < subtrahend_col) {
                    // Add the current value from the initial matrix
                    if (std::abs(matrix_csr.values_[i]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csr.values_.emplace_back(matrix_csr.values_[i]);
                        res.matrix_csr.col_indices_.emplace_back(col);
                    }
                    i++;
                } else {
                    // Add the current value from the subtrahend matrix
                    if (std::abs(-subtrahend.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csr.values_.emplace_back(-subtrahend.matrix_csr.values_[j]);
                        res.matrix_csr.col_indices_.emplace_back(subtrahend_col);
                    }
                    j++;
                }
            }

            // Add any remaining elements from the initial matrix
            for (; i < end; i++) {
                if (std::abs(matrix_csr.values_[i]) >= std::numeric_limits<VTYPE>::min()) {
                    res.matrix_csr.values_.emplace_back(matrix_csr.values_[i]);
                    res.matrix_csr.col_indices_.emplace_back(matrix_csr.col_indices_[i]);
                }
            }

            // Add any remaining elements from the subtrahend matrix
            for (; j < subtrahend_end; j++) {
                if (std::abs(-subtrahend.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                    res.matrix_csr.values_.emplace_back(-subtrahend.matrix_csr.values_[j]);
                    res.matrix_csr.col_indices_.emplace_back(subtrahend.matrix_csr.col_indices_[j]);
                }
            }

            // Update the row pointer
            res.nnz_ = res.matrix_csr.col_indices_.size();
            res.matrix_csr.row_ptrs_[row + 1] = res.nnz_;

            row++;
        }

    } else if (matrix_format_ == format::CSC) {

        // Perform element-wise subtraction of the CSC vectors
        res.matrix_csc.values_.reserve(matrix_csc.values_.size() + subtrahend.matrix_csc.values_.size());
        res.matrix_csc.row_indices_.reserve(matrix_csc.row_indices_.size() + subtrahend.matrix_csc.row_indices_.size());
        res.matrix_csc.col_ptrs_.resize(cols_ + 1);

        ITYPE column = 0;
        while (column < cols_) {
            ITYPE start = matrix_csc.col_ptrs_[column];
            ITYPE end = matrix_csc.col_ptrs_[column + 1];

            ITYPE subtrahend_start = subtrahend.matrix_csc.col_ptrs_[column];
            ITYPE subtrahend_end = subtrahend.matrix_csc.col_ptrs_[column + 1];

            res.matrix_csc.col_ptrs_[0] = 0;

            // Merge the values and row indices of the two columns
            ITYPE i = start;
            ITYPE j = subtrahend_start;
            while (i < end && j < subtrahend_end) {
                ITYPE row = matrix_csc.row_indices_[i];
                ITYPE subtrahend_row = subtrahend.matrix_csc.row_indices_[j];

                if (row == subtrahend_row) {
                    // Add the corresponding values if the columns match
                    if (std::abs(matrix_csc.values_[i] - subtrahend.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csc.values_.emplace_back(matrix_csc.values_[i] - subtrahend.matrix_csc.values_[j]);
                        res.matrix_csc.row_indices_.emplace_back(row);
                    }
                    i++;
                    j++;
                } else if (row < subtrahend_row) {
                    // Add the current value from the initial matrix
                    if (std::abs(matrix_csc.values_[i]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csc.values_.emplace_back(matrix_csc.values_[i]);
                        res.matrix_csc.row_indices_.emplace_back(row);
                    }
                    i++;
                } else {
                    // Add the current value from the subtrahend matrix
                    if (std::abs(-subtrahend.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csc.values_.emplace_back(-subtrahend.matrix_csc.values_[j]);
                        res.matrix_csc.row_indices_.emplace_back(subtrahend_row);
                    }
                    j++;
                }
            }

            // Add any remaining elements from the initial matrix
            for (; i < end; i++) {
                if (std::abs(matrix_csc.values_[i]) >= std::numeric_limits<VTYPE>::min()){
                    res.matrix_csc.values_.emplace_back(matrix_csc.values_[i]);
                    res.matrix_csc.row_indices_.emplace_back(matrix_csc.row_indices_[i]);
                }
            }

            // Add any remaining elements from the subtrahend matrix
            for (; j < subtrahend_end; j++) {
                if (std::abs(-subtrahend.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()){
                    res.matrix_csc.values_.emplace_back(-subtrahend.matrix_csc.values_[j]);
                    res.matrix_csc.row_indices_.emplace_back(subtrahend.matrix_csc.row_indices_[j]);
                }
            }

            // Update the column pointer
            res.nnz_ = res.matrix_csc.row_indices_.size();
            res.matrix_csc.col_ptrs_[column + 1] = res.nnz_;

            column++;
        }

    } else {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix operator-()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
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

    if (multiplicand.matrix_format_ != matrix_format_) {
        multiplicand.format_conversion(this->get_format(), true, true, "overwrite");
    } else {
        if (!multiplicand.is_sorted_unique("matrix_arithmetic.h", "operator*()")){
            if (multiplicand.matrix_format_ == format::COO) {
                multiplicand.sort_coo(true, true, "overwrite");
            } else if (multiplicand.matrix_format_ == format::CSR) {
                multiplicand.sort_csr(true, "overwrite");
            } else if (multiplicand.matrix_format_ == format::CSC) {
                multiplicand.sort_csc(true, "overwrite");
            } else {
                std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised multiplicand matrix format for matrix operator*()" << std::endl;
                std::cerr << "    Please set the multiplicand matrix_format_ as:" << std::endl;
                std::cerr << "    format::COO: COOrdinate format" << std::endl;
                std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                std::abort();
            }
        }
    }

     if (!this->is_sorted_unique("matrix_arithmetic.h", "operator*()")){
         if (matrix_format_ == format::COO) {
             this->sort_coo(true, true, "overwrite");
         } else if (matrix_format_ == format::CSR) {
             this->sort_csr(true, "overwrite");
         } else if (matrix_format_ == format::CSC) {
             this->sort_csc(true, "overwrite");
         } else {
             std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix operator*()" << std::endl;
             std::cerr << "    Please set the matrix_format_ as:" << std::endl;
             std::cerr << "    format::COO: COOrdinate format" << std::endl;
             std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
             std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
             std::abort();
         }
     }

    // Create a new sparse matrix object for the result
    sparse_matrix<ITYPE,VTYPE> res(rows_, multiplicand.cols_, this->get_format());

    if (matrix_format_ == format::COO) {

        // Perform element-wise multiplication of the COO vectors
        res.matrix_coo.values_.reserve((matrix_coo.values_.size() <= multiplicand.matrix_coo.values_.size()) ? multiplicand.matrix_coo.values_.size() : matrix_coo.values_.size());
        res.matrix_coo.row_indices_.reserve((matrix_coo.row_indices_.size() <= multiplicand.matrix_coo.row_indices_.size()) ? multiplicand.matrix_coo.row_indices_.size() : matrix_coo.row_indices_.size());
        res.matrix_coo.col_indices_.reserve((matrix_coo.col_indices_.size() <= multiplicand.matrix_coo.col_indices_.size()) ? multiplicand.matrix_coo.col_indices_.size() : matrix_coo.col_indices_.size());

        for (ITYPE i = 0; i < static_cast<ITYPE>(matrix_coo.row_indices_.size()); ++i) {
            for (ITYPE j = 0; j < static_cast<ITYPE>(multiplicand.matrix_coo.col_indices_.size()); ++j) {
                if (matrix_coo.col_indices_[i] == multiplicand.matrix_coo.row_indices_[j]) {
                    // Multiply the corresponding values if the columns match
                    VTYPE value = matrix_coo.values_[i] * multiplicand.matrix_coo.values_[j];
                    if (std::abs(value) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_coo.values_.emplace_back(value);
                        res.matrix_coo.row_indices_.emplace_back(matrix_coo.row_indices_[i]);
                        res.matrix_coo.col_indices_.emplace_back(multiplicand.matrix_coo.col_indices_[j]);
                    }
                }
            }
        }

        // Sort and deduplicate the result
        res.sort_coo(true, true, "plus");

        res.nnz_ = res.matrix_coo.values_.size();

    } else if (matrix_format_ == format::CSR) {

        // Perform element-wise multiplication of the CSR vectors
        res.matrix_csr.values_.reserve((matrix_csr.values_.size() <= multiplicand.matrix_csr.values_.size()) ? multiplicand.matrix_csr.values_.size() : matrix_csr.values_.size());
        res.matrix_csr.row_ptrs_.resize(rows_+1);
        res.matrix_csr.col_indices_.reserve((matrix_csr.col_indices_.size() <= multiplicand.matrix_csr.col_indices_.size()) ? multiplicand.matrix_csr.col_indices_.size() : matrix_csr.col_indices_.size());

        // Initialize a vector to store the intermediate results
        std::vector<VTYPE> intermediate(multiplicand.cols_, 0.0);

        res.matrix_csr.row_ptrs_[0] = 0;

        // Iterate over each row of the initial matrix
        for (ITYPE i = 0; i < rows_; ++i) {
            // Clear the intermediate results vector for each row
            std::fill(intermediate.begin(), intermediate.end(), 0.0);

            ITYPE start = matrix_csr.row_ptrs_[i];
            ITYPE end = matrix_csr.row_ptrs_[i + 1];

            // Iterate over the non-zero elements of the row
            for (ITYPE j = start; j < end; ++j) {
                // Get the column index and value of the element
                ITYPE col = matrix_csr.col_indices_[j];
                VTYPE value = matrix_csr.values_[j];

                ITYPE multiplicand_start = multiplicand.matrix_csr.row_ptrs_[col];
                ITYPE multiplicand_end = multiplicand.matrix_csr.row_ptrs_[col + 1];

                // Multiply the element with the corresponding column of the other matrix
                for (ITYPE k = multiplicand_start; k < multiplicand_end; ++k) {
                    ITYPE multiplicand_col = multiplicand.matrix_csr.col_indices_[k];
                    VTYPE multiplicand_value = multiplicand.matrix_csr.values_[k];
                    intermediate[multiplicand_col] += value * multiplicand_value;
                }
            }

            // Add the intermediate results to the result vectors
            for (ITYPE j = 0; j < multiplicand.cols_; ++j) {
                VTYPE result_value = intermediate[j];
                if (std::abs(result_value) >= std::numeric_limits<VTYPE>::min()) {
                    res.matrix_csr.values_.emplace_back(result_value);
                    res.matrix_csr.col_indices_.emplace_back(j);
                }
            }
            res.matrix_csr.row_ptrs_[i+1]=res.matrix_csr.values_.size();
        }

        res.nnz_ = res.matrix_csr.values_.size();

    } else if (matrix_format_ == format::CSC) {

        // Perform element-wise multiplication of the CSC vectors
        res.matrix_csc.values_.reserve((matrix_csc.values_.size() <= multiplicand.matrix_csc.values_.size()) ? multiplicand.matrix_csc.values_.size() : matrix_csc.values_.size());
        res.matrix_csc.row_indices_.reserve((matrix_csc.row_indices_.size() <= multiplicand.matrix_csc.row_indices_.size()) ? multiplicand.matrix_csc.row_indices_.size() : matrix_csc.row_indices_.size());
        res.matrix_csc.col_ptrs_.resize(cols_+1);

        // Initialize a vector to store the intermediate results
        std::vector<VTYPE> intermediate(rows_, 0.0);

        res.matrix_csc.col_ptrs_[0] = 0;

        // Iterate over each column of the initial matrix
        for (ITYPE j = 0; j < multiplicand.cols_; ++j) {
            // Clear the intermediate results vector for each row
            std::fill(intermediate.begin(), intermediate.end(), 0.0);

            ITYPE multiplicand_start = multiplicand.matrix_csc.col_ptrs_[j];
            ITYPE multiplicand_end = multiplicand.matrix_csc.col_ptrs_[j + 1];

            // Iterate over the non-zero elements of the cloumn
            for (ITYPE k = multiplicand_start; k < multiplicand_end; ++k) {
                // Get the row index and value of the element
                ITYPE multiplicand_row = multiplicand.matrix_csc.row_indices_[k];
                VTYPE multiplicand_value = multiplicand.matrix_csc.values_[k];

                ITYPE start = matrix_csc.col_ptrs_[multiplicand_row];
                ITYPE end = matrix_csc.col_ptrs_[multiplicand_row + 1];

                // Multiply the element with the corresponding column of the other matrix
                for (ITYPE i = start; i < end; ++i) {
                    ITYPE row = matrix_csc.row_indices_[i];
                    VTYPE value = matrix_csc.values_[i];
                    intermediate[row] += value * multiplicand_value;
                }
            }

            // Add the intermediate results to the result vectors
            for (ITYPE i = 0; i < multiplicand.rows_; ++i) {
                VTYPE result_value = intermediate[i];
                if (std::abs(result_value) >= std::numeric_limits<VTYPE>::min()) {
                    res.matrix_csc.values_.emplace_back(result_value);
                    res.matrix_csc.row_indices_.emplace_back(i);
                }
            }
            res.matrix_csc.col_ptrs_[j+1]=res.matrix_csc.values_.size();
        }

        res.nnz_ = res.matrix_csc.values_.size();

    } else {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix operator*()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }

    return res;
}

// Overload multiplication operator to perform scalar multiplication A*x
template <typename ITYPE, typename VTYPE>
template <typename STYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::operator*(const STYPE &scalar) const {
    static_assert(std::is_convertible<STYPE, VTYPE>::value,
            "MUI Error [matrix_arithmetic.h]: scalar type cannot be converted to matrix element type in scalar multiplication");

    // Create a new sparse matrix object for the result
    sparse_matrix<ITYPE,VTYPE> res(*this);

    if (matrix_format_ == format::COO) {

        for (VTYPE &element : res.matrix_coo.values_) {
            if (static_cast<VTYPE>(scalar) >= std::numeric_limits<VTYPE>::min())
                element *= scalar;
       }

    } else if (matrix_format_ == format::CSR) {

        for (VTYPE &element : res.matrix_csr.values_) {
            if (static_cast<VTYPE>(scalar) >= std::numeric_limits<VTYPE>::min())
                element *= scalar;
       }

    } else if (matrix_format_ == format::CSC) {

        for (VTYPE &element : res.matrix_csc.values_) {
            if (static_cast<VTYPE>(scalar) >= std::numeric_limits<VTYPE>::min())
                element *= scalar;
       }

    } else {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix scalar operator*()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
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
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::hadamard_product(sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    if (rows_ != exist_mat.rows_ || cols_ != exist_mat.cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: matrix size mismatch during matrix Hadamard product" << std::endl;
        std::abort();
    }

    if (exist_mat.matrix_format_ != matrix_format_) {
        exist_mat.format_conversion(this->get_format(), true, true, "overwrite");
    } else {
        if (!exist_mat.is_sorted_unique("matrix_arithmetic.h", "hadamard_product()")){
            if (exist_mat.matrix_format_ == format::COO) {
                exist_mat.sort_coo(true, true, "overwrite");
            } else if (exist_mat.matrix_format_ == format::CSR) {
                exist_mat.sort_csr(true, "overwrite");
            } else if (exist_mat.matrix_format_ == format::CSC) {
                exist_mat.sort_csc(true, "overwrite");
            } else {
                std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised exist_mat matrix format for matrix hadamard_product()" << std::endl;
                std::cerr << "    Please set the exist_mat matrix_format_ as:" << std::endl;
                std::cerr << "    format::COO: COOrdinate format" << std::endl;
                std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                std::abort();
            }
        }
    }

    if (!this->is_sorted_unique("matrix_arithmetic.h", "hadamard_product()")){
        if (matrix_format_ == format::COO) {
            this->sort_coo(true, true, "overwrite");
        } else if (matrix_format_ == format::CSR) {
            this->sort_csr(true, "overwrite");
        } else if (matrix_format_ == format::CSC) {
            this->sort_csc(true, "overwrite");
        } else {
            std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix hadamard_product()" << std::endl;
            std::cerr << "    Please set the matrix_format_ as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }
    }

    // Create a new sparse matrix object for the result
    sparse_matrix<ITYPE,VTYPE> res(rows_, cols_, this->get_format());

    if (matrix_format_ == format::COO) {

        // Perform hadamard product of the COO vectors
        res.matrix_coo.values_.reserve(matrix_coo.values_.size() + exist_mat.matrix_coo.values_.size());
        res.matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size() + exist_mat.matrix_coo.row_indices_.size());
        res.matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size() + exist_mat.matrix_coo.col_indices_.size());

        // Insert the COO vectors of the initial sparse matrix to the result sparse matrix
        res.matrix_coo.values_ = std::vector<VTYPE>(matrix_coo.values_.begin(), matrix_coo.values_.end());
        res.matrix_coo.row_indices_ = std::vector<ITYPE>(matrix_coo.row_indices_.begin(), matrix_coo.row_indices_.end());
        res.matrix_coo.col_indices_ = std::vector<ITYPE>(matrix_coo.col_indices_.begin(), matrix_coo.col_indices_.end());

        // Append the exist_mat COO vectors to the result sparse matrix
        res.matrix_coo.values_.insert(res.matrix_coo.values_.end(), exist_mat.matrix_coo.values_.begin(), exist_mat.matrix_coo.values_.end());
        res.matrix_coo.row_indices_.insert(res.matrix_coo.row_indices_.end(), exist_mat.matrix_coo.row_indices_.begin(), exist_mat.matrix_coo.row_indices_.end());
        res.matrix_coo.col_indices_.insert(res.matrix_coo.col_indices_.end(), exist_mat.matrix_coo.col_indices_.begin(), exist_mat.matrix_coo.col_indices_.end());

        // Sort and deduplicate the result
        res.sort_coo(true, true, "multiply");

    } else if (matrix_format_ == format::CSR) {

        // Perform element-wise hadamard product of the CSR vectors
        res.matrix_csr.values_.reserve(matrix_csr.values_.size() + exist_mat.matrix_csr.values_.size());
        res.matrix_csr.row_ptrs_.resize(rows_ + 1);
        res.matrix_csr.col_indices_.reserve(matrix_csr.col_indices_.size() + exist_mat.matrix_csr.col_indices_.size());

        res.matrix_csr.row_ptrs_[0] = 0;

        ITYPE row = 0;
        while (row < rows_) {
            ITYPE start = matrix_csr.row_ptrs_[row];
            ITYPE end = matrix_csr.row_ptrs_[row + 1];

            ITYPE exist_mat_start = exist_mat.matrix_csr.row_ptrs_[row];
            ITYPE exist_mat_end = exist_mat.matrix_csr.row_ptrs_[row + 1];

            // Merge the values and column indices of the two rows
            ITYPE i = start;
            ITYPE j = exist_mat_start;
            while (i < end && j < exist_mat_end) {
                ITYPE col = matrix_csr.col_indices_[i];
                ITYPE exist_mat_col = exist_mat.matrix_csr.col_indices_[j];

                if (col == exist_mat_col) {
                    // Add the corresponding values if the columns match
                    if (std::abs(matrix_csr.values_[i] * exist_mat.matrix_csr.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                        res.matrix_csr.values_.emplace_back(matrix_csr.values_[i] * exist_mat.matrix_csr.values_[j]);
                        res.matrix_csr.col_indices_.emplace_back(col);
                    }
                    i++;
                    j++;
                } else if (col < exist_mat_col) {
                    i++;
                } else {
                    j++;
                }
            }

            // Update the row pointer
            res.nnz_ = res.matrix_csr.col_indices_.size();
            res.matrix_csr.row_ptrs_[row + 1] = res.nnz_;

            row++;
        }

    } else if (matrix_format_ == format::CSC) {

        // Perform element-wise hadamard product of the CSC vectors
        res.matrix_csc.values_.reserve(matrix_csc.values_.size() + exist_mat.matrix_csc.values_.size());
        res.matrix_csc.row_indices_.reserve(matrix_csc.row_indices_.size() + exist_mat.matrix_csc.row_indices_.size());
        res.matrix_csc.col_ptrs_.resize(cols_ + 1);

        res.matrix_csc.col_ptrs_[0] = 0;

        ITYPE column = 0;
        while (column < cols_) {
            ITYPE start = matrix_csc.col_ptrs_[column];
            ITYPE end = matrix_csc.col_ptrs_[column + 1];

            ITYPE exist_mat_start = exist_mat.matrix_csc.col_ptrs_[column];
            ITYPE exist_mat_end = exist_mat.matrix_csc.col_ptrs_[column + 1];

            // Merge the values and row indices of the two columns
            ITYPE i = start;
            ITYPE j = exist_mat_start;
            while (i < end && j < exist_mat_end) {
                ITYPE row = matrix_csc.row_indices_[i];
                ITYPE exist_mat_row = exist_mat.matrix_csc.row_indices_[j];

                if ((row == exist_mat_row) && std::abs(matrix_csc.values_[i] * exist_mat.matrix_csc.values_[j]) >= std::numeric_limits<VTYPE>::min()) {
                    // Add the corresponding values if the columns match
                    res.matrix_csc.values_.emplace_back(matrix_csc.values_[i] * exist_mat.matrix_csc.values_[j]);
                    res.matrix_csc.row_indices_.emplace_back(row);
                    i++;
                    j++;
                } else if (row < exist_mat_row) {
                    i++;
                } else {
                    j++;
                }
            }

            // Update the column pointer
            res.nnz_ = res.matrix_csc.row_indices_.size();
            res.matrix_csc.col_ptrs_[column + 1] = res.nnz_;

            column++;
        }

    } else {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix hadamard_product()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }

    return res;
}

// Member function to get transpose of matrix
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::transpose(bool performSortAndUniqueCheck) const {

    sparse_matrix<ITYPE,VTYPE> res(*this);

    if (matrix_format_ == format::COO) {

        if (performSortAndUniqueCheck){
            if (!res.is_sorted_unique("matrix_arithmetic.h", "transpose()")){
                res.sort_coo(true, true, "overwrite");
            }
        }

        res.index_reinterpretation();

    } else if (matrix_format_ == format::CSR) {

        res.format_conversion("CSC", performSortAndUniqueCheck, performSortAndUniqueCheck, "overwrite");

        res.format_reinterpretation();

    } else if (matrix_format_ == format::CSC) {

        res.format_conversion("CSR", performSortAndUniqueCheck, performSortAndUniqueCheck, "overwrite");

        res.format_reinterpretation();

    } else {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Unrecognised matrix format for matrix transpose()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }

    return res;

}

// Member function to perform LU decomposition
template <typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::lu_decomposition(sparse_matrix<ITYPE,VTYPE> &L, sparse_matrix<ITYPE,VTYPE> &U) const {
    if (((L.get_rows() != 0) && (L.get_rows() != rows_)) ||
        ((U.get_rows() != 0) && (U.get_rows() != rows_)) ||
        ((L.get_cols() != 0) && (L.get_cols() != cols_)) ||
        ((U.get_cols() != 0) && (U.get_cols() != cols_))) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: L & U Matrices must be null or same size of initial matrix in LU decomposition" << std::endl;
        std::abort();
    }

    if ((!L.empty()) || (!U.empty())) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: L & U Matrices must be empty in LU decomposition" << std::endl;
        std::abort();
    }

    if (rows_ != cols_) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Only square matrix can perform LU decomposition" << std::endl;
        std::abort();
    }

    if ((L.get_rows() != rows_) || (L.get_cols() != cols_)) {
        // Resize the lower triangular matrix
        L.resize(rows_, cols_);
    }

    if ((U.get_rows() != rows_) || (U.get_cols() != cols_)) {
        // Resize the upper triangular matrix
        U.resize(rows_, cols_);
    }

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
    if (((Q.get_rows() != 0) && (Q.get_rows() != rows_)) ||
        ((R.get_rows() != 0) && (R.get_rows() != rows_)) ||
        ((Q.get_cols() != 0) && (Q.get_cols() != cols_)) ||
        ((R.get_cols() != 0) && (R.get_cols() != cols_))) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Q & R Matrices must be null in QR decomposition" << std::endl;
        std::abort();
    }
    if ((!Q.empty()) || (!R.empty())) {
        std::cerr << "MUI Error [matrix_arithmetic.h]: Q & R Matrices must be empty in QR decomposition" << std::endl;
        std::abort();
    }
    assert((rows_ >= cols_) &&
          "MUI Error [matrix_arithmetic.h]: number of rows of matrix should larger or equals to number of columns in QR decomposition");

    if ((Q.get_rows() != rows_) || (Q.get_cols() != cols_)) {
        // Resize the orthogonal matrix
        Q.resize(rows_, cols_);
    }
    if ((R.get_rows() != rows_) || (R.get_cols() != cols_)) {
        // Resize the upper triangular matrix
        R.resize(rows_, cols_);
    }

    // Get a copy of the matrix
    sparse_matrix<ITYPE,VTYPE> mat_copy (*this);
    // Diagonal elements
    std::vector<VTYPE> r_diag (cols_);

    // Calculate the diagonal element values
    for (ITYPE c = 0; c <cols_; ++c)  {
        VTYPE  nrm (0.0);

       // Compute 2-norm of k-th column without under/overflow.
        for (ITYPE r = c; r < rows_; ++r)
            nrm = std::sqrt((nrm * nrm) + (mat_copy.get_value(r, c) * mat_copy.get_value(r, c)));

        if (nrm != static_cast<VTYPE>(0.0))  {

           // Form k-th Householder vector.
            if (mat_copy.get_value(c, c) < static_cast<VTYPE>(0.0))
                nrm = -nrm;

            for (ITYPE r = c; r < rows_; ++r)
                mat_copy.set_value(r, c, (mat_copy.get_value(r, c)/nrm));

            mat_copy.set_value(c, c, (mat_copy.get_value(c, c) + static_cast<VTYPE>(1.0)));

           // Apply transformation to remaining columns.
            for (ITYPE j = c + 1; j < cols_; ++j)  {
                VTYPE  s = 0.0;

                for (ITYPE r = c; r < rows_; ++r)
                    s += mat_copy.get_value(r, c) * mat_copy.get_value(r, j);

                s /= -mat_copy.get_value(c, c);
                for (ITYPE r = c; r < rows_; ++r)
                    mat_copy.set_value(r, j, (mat_copy.get_value(r, j) + s * mat_copy.get_value(r, c)));
            }
        }
        r_diag[c] = -nrm;
    }

    // Calculate the orthogonal matrix
    for (ITYPE c = cols_ - 1; c >= 0; --c)  {
        Q.set_value(c, c, static_cast<VTYPE>(1.0));

        for (ITYPE cc = c; cc < cols_; ++cc)
            if (mat_copy.get_value(c, c) != static_cast<VTYPE>(0.0)) {
                VTYPE s=0.0;

                for (ITYPE r = c; r < rows_; ++r)
                    s += mat_copy.get_value(r, c) * Q.get_value(r, cc);

                s /= -mat_copy.get_value(c, c);
                for (ITYPE r = c; r < rows_; ++r)
                    Q.set_value(r, cc, (Q.get_value(r, cc) + s * mat_copy.get_value(r, c)));
            }
    }

    // Calculate the upper triangular matrix
    for (ITYPE c = 0; c < cols_; ++c)
        for (ITYPE r = 0; r < rows_; ++r)
            if (c < r)
                R.set_value(c, r, mat_copy.get_value(c, r));
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
    sparse_matrix<ITYPE,VTYPE> inverse_mat (rows_,"identity", this->get_format());

    for (ITYPE r = 0; r < rows_; ++r)  {

        ITYPE max_row = r;
        VTYPE max_value= static_cast<VTYPE>(-1.0);

        // Partial pivoting for Gaussian elimination
        ITYPE ppivot;
        for (ITYPE rb = r; rb < rows_; ++rb)  {
            const VTYPE tmp = std::abs(mat_copy.get_value(rb, r));

            if ((tmp > max_value) && (std::abs(tmp) >= std::numeric_limits<VTYPE>::min()))  {
                max_value = tmp;
                max_row = rb;
            }
        }

        assert(std::abs(mat_copy.get_value(max_row, r)) >= std::numeric_limits<VTYPE>::min() &&
                          "MUI Error [matrix_arithmetic.h]: Divide by zero assert for mat_copy.get_value(max_row, r). Cannot perform matrix invert due to singular matrix.");

        if (max_row != r)  {
            for (ITYPE c = 0; c < cols_; ++c)
                mat_copy.swap_elements(r, c, max_row, c);
            ppivot = max_row;
        } else {
            ppivot = 0;
        }

        const ITYPE indx = ppivot;

        if (indx != 0)
            for (ITYPE c = 0; c < cols_; ++c)
                inverse_mat.swap_elements(r, c, indx, c);

        const VTYPE diag = mat_copy.get_value(r, r);

        for (ITYPE c = 0; c < cols_; ++c)  {
            mat_copy.set_value(r, c, (mat_copy.get_value(r, c) / diag));
            inverse_mat.set_value(r, c, (inverse_mat.get_value(r, c) / diag));
        }

        for (ITYPE rr = 0; rr < rows_; ++rr)
            if (rr != r)  {
                const VTYPE off_diag = mat_copy.get_value(rr, r);

                for (ITYPE c = 0; c < cols_; ++c)  {
                    mat_copy.set_value(rr, c, (mat_copy.get_value(rr, c) - off_diag * mat_copy.get_value(r, c)));
                    inverse_mat.set_value(rr, c, (inverse_mat.get_value(rr, c) - off_diag * inverse_mat.get_value(r, c)));
                }
            }
    }
    return inverse_mat;
}

// **************************************************
// ********** Protected member functions ************
// **************************************************

// Protected member function to reinterpret the row and column indexes for sparse matrix with COO format - helper function on matrix transpose
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::index_reinterpretation() {
    assert((matrix_format_ == format::COO) &&
              "MUI Error [matrix_arithmetic.h]: index_reinterpretation() is for COO format only.");

    std::swap(matrix_coo.row_indices_, matrix_coo.col_indices_);
    ITYPE temp_index = rows_;
    rows_ = cols_;
    cols_ = temp_index;

}

// Protected member function to reinterpret the format of sparse matrix between CSR format and CSC format - helper function on matrix transpose
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::format_reinterpretation() {
    assert(((matrix_format_ == format::CSR) || (matrix_format_ == format::CSC)) &&
              "MUI Error [matrix_arithmetic.h]: format_reinterpretation() is for CSR or CSC format.");

    if (matrix_format_ == format::CSR) {

        matrix_csc.col_ptrs_.swap(matrix_csr.row_ptrs_);
        matrix_csc.row_indices_.swap(matrix_csr.col_indices_);
        matrix_csc.values_.swap(matrix_csr.values_);

        ITYPE temp_index = rows_;
        rows_ = cols_;
        cols_ = temp_index;

        matrix_format_ = format::CSC;

        matrix_csr.row_ptrs_.clear();
        matrix_csr.col_indices_.clear();
        matrix_csr.values_.clear();

    } else if (matrix_format_ == format::CSC) {

        matrix_csr.row_ptrs_.swap(matrix_csc.col_ptrs_);
        matrix_csr.col_indices_.swap(matrix_csc.row_indices_);
        matrix_csr.values_.swap(matrix_csc.values_);

        ITYPE temp_index = rows_;
        rows_ = cols_;
        cols_ = temp_index;

        matrix_format_ = format::CSR;

        matrix_csc.col_ptrs_.clear();
        matrix_csc.row_indices_.clear();
        matrix_csc.values_.clear();

    }

}

} // linalg
} // mui

#endif /* MUI_MATRIX_ARITHMETIC_H_ */
