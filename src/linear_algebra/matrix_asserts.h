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
 * @file matrix_asserts.h
 * @author W. Liu
 * @date 15 May 2023
 * @brief Implementation of sparse matrix asserts functions.
 */

#ifndef MUI_MATRIX_ASSERTS_H_
#define MUI_MATRIX_ASSERTS_H_

#include <cassert>

namespace mui {
namespace linalg {

// **************************************************
// ********** Protected member functions ************
// **************************************************

// Member function to assert the matrix vector sizes
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::assert_valid_vector_size(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name.assign("matrix_asserts.h");
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name.assign("assert_valid_vector_size()");
    } else {
        function_name = function_name_input;
    }

    if (matrix_format_ == format::COO) {

    	if (nnz_ < 0) {
            std::cerr << "MUI Error [" << file_name << "]: The number of non-zeros (nnz_=" << nnz_ << ") should be non-negative integer in " << function_name << std::endl;
            std::abort();
        }

        if (rows_ < 0) {
            std::cerr << "MUI Error [" << file_name << "]: The number of rows (rows_=" << rows_ << ") should be non-negative integer in " << function_name << std::endl;
            std::abort();
        }

        if (cols_ < 0) {
            std::cerr << "MUI Error [" << file_name << "]: The number of columns (cols_=" << cols_ << ") should be non-negative integer in " << function_name << std::endl;
            std::abort();
        }

        if ((rows_*cols_) < nnz_) {
            std::cerr << "MUI Warning [" << file_name << "]: Matrix size (" << (rows_*cols_) << ") smaller than the number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << ". Possible duplicated elements occur. " << std::endl;
        }

        if (static_cast<ITYPE>(matrix_coo.values_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: COO values_ matrix size (" << matrix_coo.values_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (static_cast<ITYPE>(matrix_coo.row_indices_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: COO row_indices_ matrix size (" << matrix_coo.row_indices_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (static_cast<ITYPE>(matrix_coo.col_indices_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: COO col_indices_ matrix size (" << matrix_coo.col_indices_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (!matrix_csr.values_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSR values_ matrix (" << matrix_csr.values_.size() << ") under COO matrix format in " << function_name << std::endl;
        }

        if (!matrix_csr.row_ptrs_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSR row_ptrs_ matrix (" << matrix_csr.row_ptrs_.size() << ") under COO matrix format in " << function_name << std::endl;
        }

        if (!matrix_csr.col_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSR col_indices_ matrix (" << matrix_csr.col_indices_.size() << ") under COO matrix format in " << function_name << std::endl;
        }

        if (!matrix_csc.values_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSC values_ matrix (" << matrix_csc.values_.size() << ") under COO matrix format in " << function_name << std::endl;
        }

        if (!matrix_csc.row_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSC row_ptrs_ matrix (" << matrix_csc.row_indices_.size() << ") under COO matrix format in " << function_name << std::endl;
        }

        if (!matrix_csc.col_ptrs_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSC col_indices_ matrix (" << matrix_csc.col_ptrs_.size() << ") under COO matrix format in " << function_name << std::endl;
        }

    } else if (matrix_format_ == format::CSR) {

        if ((rows_*cols_) < nnz_) {
            std::cerr << "MUI Warning [" << file_name << "]: Matrix size (" << (rows_*cols_) << ") smaller than the number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << ". Possible duplicated elements occur. " << std::endl;
        }

        if (!matrix_coo.values_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty COO values_ matrix (" << matrix_coo.values_.size() << ") under CSR matrix format in " << function_name << std::endl;
        }

        if (!matrix_coo.row_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty COO row_indices_ matrix (" << matrix_coo.row_indices_.size() << ") under CSR matrix format in " << function_name << std::endl;
        }

        if (!matrix_coo.col_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty COO col_indices_ matrix (" << matrix_coo.col_indices_.size() << ") under CSR matrix format in " << function_name << std::endl;
        }

        if (static_cast<ITYPE>(matrix_csr.values_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: CSR values_ matrix size (" << matrix_csr.values_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (static_cast<ITYPE>(matrix_csr.row_ptrs_.size()) != (rows_+1)) {
            std::cerr << "MUI Error [" << file_name << "]: CSR row_ptrs_ matrix size (" << matrix_csr.row_ptrs_.size() << ") does not equals to number of rows + 1 (rows_+1=" << (rows_+1) << ") in " << function_name << std::endl;
            std::abort();
        }

        if (static_cast<ITYPE>(matrix_csr.col_indices_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: CSR col_indices_ matrix size (" << matrix_csr.col_indices_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (!matrix_csc.values_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSC values_ matrix (" << matrix_csc.values_.size() << ") under CSR matrix format in " << function_name << std::endl;
        }

        if (!matrix_csc.row_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSC row_ptrs_ matrix (" << matrix_csc.row_indices_.size() << ") under CSR matrix format in " << function_name << std::endl;
        }

        if (!matrix_csc.col_ptrs_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSC col_indices_ matrix (" << matrix_csc.col_ptrs_.size() << ") under CSR matrix format in " << function_name << std::endl;
        }

    } else if (matrix_format_ == format::CSC) {

        if ((rows_*cols_) < nnz_) {
            std::cerr << "MUI Warning [" << file_name << "]: Matrix size (" << (rows_*cols_) << ") smaller than the number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << ". Possible duplicated elements occur. " << std::endl;
        }

        if (!matrix_coo.values_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty COO values_ matrix (" << matrix_coo.values_.size() << ") under CSC matrix format in " << function_name << std::endl;
        }

        if (!matrix_coo.row_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty COO row_indices_ matrix (" << matrix_coo.row_indices_.size() << ") under CSC matrix format in " << function_name << std::endl;
        }

        if (!matrix_coo.col_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty COO col_indices_ matrix (" << matrix_coo.col_indices_.size() << ") under CSC matrix format in " << function_name << std::endl;
        }

        if (!matrix_csr.values_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSR values_ matrix size (" << matrix_csr.values_.size() << ") under CSC matrix format in " << function_name << std::endl;
        }

        if (!matrix_csr.row_ptrs_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSR row_ptrs_ matrix (" << matrix_csr.row_ptrs_.size() << ") under CSC matrix format in " << function_name << std::endl;
        }

        if (!matrix_csr.col_indices_.empty()) {
            std::cerr << "MUI Warning [" << file_name << "]: Non-empty CSR col_indices_ matrix (" << matrix_csr.col_indices_.size() << ") under CSC matrix format in " << function_name << std::endl;
        }

        if (static_cast<ITYPE>(matrix_csc.values_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: CSC values_ matrix size (" << matrix_csc.values_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (static_cast<ITYPE>(matrix_csc.row_indices_.size()) != nnz_) {
            std::cerr << "MUI Error [" << file_name << "]: CSC row_indices_ matrix size (" << matrix_csc.row_indices_.size() << ") does not equals to number of non-zeros (nnz_=" << nnz_ << ") in " << function_name << std::endl;
            std::abort();
        }

        if (static_cast<ITYPE>(matrix_csc.col_ptrs_.size()) != (cols_+1)) {
            std::cerr << "MUI Error [" << file_name << "]: CSC col_indices_ matrix size (" << matrix_csc.col_ptrs_.size() << ") does not equals to number of cols + 1 (cols_+1=" << (cols_+1) << ") in " << function_name << std::endl;
            std::abort();
        }

    } else {
        std::cerr << "MUI Error [" << file_name << "]: unknown matrix format in " << function_name <<  std::endl;
        std::abort();
    }

}

// Member function to assert if the COO matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::assert_coo_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty())
        file_name.assign("matrix_asserts.h");
    else
        file_name = file_name_input;
    if (function_name_input.empty())
        function_name.assign("assert_coo_sorted_unique()");
    else
        function_name = function_name_input;

    ITYPE numEntries = matrix_coo.values_.size();

    if (numEntries > 1) {
        for (ITYPE i = 1; i < numEntries; ++i) {
            // Compare the current entry with the previous one
            if (matrix_coo.row_indices_[i] < matrix_coo.row_indices_[i - 1]) {
                // Row index is not sorted
                std::cerr << "MUI Error [" << file_name << "]: The COO type matrix is not sorted (sorted row index check failed) in " << function_name <<  std::endl;
                std::abort();
            } else if (matrix_coo.row_indices_[i] == matrix_coo.row_indices_[i - 1]) {
                // Row index is the same, check column index
                if (matrix_coo.col_indices_[i] < matrix_coo.col_indices_[i - 1]) {
                    // Column index is not sorted
                    std::cerr << "MUI Error [" << file_name << "]: The COO type matrix is not sorted (sorted column index check failed) in " << function_name <<  std::endl;
                    std::abort();
                } else if (matrix_coo.col_indices_[i] == matrix_coo.col_indices_[i - 1]) {
                    // Column index has duplicate elements
                    std::cerr << "MUI Error [" << file_name << "]: The COO type matrix exists duplicated elements (unique column index check failed) in " << function_name <<  std::endl;
                    std::abort();
                }
            }
        }
    }
}

// Member function to assert if the CSR matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::assert_csr_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty())
        file_name.assign("matrix_asserts.h");
    else
        file_name = file_name_input;
    if (function_name_input.empty())
        function_name.assign("assert_csr_sorted_unique()");
    else
        function_name = function_name_input;

    ITYPE numEntries = matrix_csr.values_.size();

    if (numEntries > 1) {
        for(ITYPE i = 0; i < rows_; ++i){
            if (matrix_csr.row_ptrs_[i] > matrix_csr.row_ptrs_[i+1]) {
                // Row pointers is not sorted
                std::cerr << "MUI Error [" << file_name << "]: The CSR type matrix is not sorted (sorted row pointers check failed) in " << function_name <<  std::endl;
                std::abort();
            }
            for(ITYPE j = matrix_csr.row_ptrs_[i] + 1; j < matrix_csr.row_ptrs_[i+1]; ++j){
                if(matrix_csr.col_indices_[j-1] > matrix_csr.col_indices_[j]){
                    // Column indices is not sorted
                    std::cerr << "MUI Error [" << file_name << "]: The CSR type matrix is not sorted (sorted column index check failed) in " << function_name <<  std::endl;
                    std::abort();
                } else if (matrix_csr.col_indices_[j-1] == matrix_csr.col_indices_[j]) {
                    // Column indices is not unique
                    std::cerr << "MUI Error [" << file_name << "]: The CSR type matrix is not unique (deduplicated column index check failed) in " << function_name <<  std::endl;
                    std::abort();
                }
            }
        }
    }
}

// Member function to assert if the CSC matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::assert_csc_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty())
        file_name.assign("matrix_asserts.h");
    else
      file_name = file_name_input;
    if (function_name_input.empty())
        function_name.assign("assert_csc_sorted_unique()");
    else
      function_name = function_name_input;

    ITYPE numEntries = matrix_csc.values_.size();

    if (numEntries > 1) {
        for(ITYPE i = 0; i < cols_; ++i){
            if (matrix_csc.col_ptrs_[i] > matrix_csc.col_ptrs_[i+1]) {
                // Column pointers is not sorted
                std::cerr << "MUI Error [" << file_name << "]: The CSC type matrix is not sorted (sorted column pointers check failed) in " << function_name <<  std::endl;
                std::abort();
            }
            for(ITYPE j = matrix_csc.col_ptrs_[i] + 1; j < matrix_csc.col_ptrs_[i+1]; ++j){
                if(matrix_csc.row_indices_[j-1] > matrix_csc.row_indices_[j]){
                    // Row indices is not sorted
                    std::cerr << "MUI Error [" << file_name << "]: The CSC type matrix is not sorted (sorted row index check failed) in " << function_name <<  std::endl;
                    std::abort();
                } else if (matrix_csc.row_indices_[j-1] == matrix_csc.row_indices_[j]) {
                    // Row indices is not unique
                    std::cerr << "MUI Error [" << file_name << "]: The CSC type matrix is not unique (deduplicated row index check failed) in " << function_name <<  std::endl;
                    std::abort();
                }
            }
        }
    }
}

} // linalg
} // mui

#endif /* MUI_MATRIX_ASSERTS_H_ */
