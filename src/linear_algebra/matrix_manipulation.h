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
 * @file matrix_manipulation.h
 * @author W. Liu
 * @date 01 February 2023
 * @brief Implementation of sparse matrix manipulation functions.
 */

#ifndef MUI_MATRIX_MANIPULATION_H_
#define MUI_MATRIX_MANIPULATION_H_

#include <cassert>

namespace mui {
namespace linalg {

// Member function to resize a null matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::resize_null(ITYPE r, ITYPE c) {
    assert(((rows_ == 0) && (cols_ == 0)) &&
            "MUI Error [matrix_manipulation.h]: resize_null function only works for null matrix");
    rows_ = r;
    cols_ = c;
}

// Member function to resize an all-zero matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::resize(ITYPE r, ITYPE c) {
    assert(((this->non_zero_elements_count()) == 0) &&
            "MUI Error [matrix_manipulation.h]: resize function only works for all-zero matrix");
    rows_ = r;
    cols_ = c;
}

// Member function to copy a sparse_matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::copy(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
      // Copy the data from the existing matrix
      assert(matrix_.empty() &&
                "MUI Error [matrix_manipulation.h]: copy function only works for empty (all zero elements) matrix");
      assert(((rows_ == exist_mat.rows_) && (cols_ == exist_mat.cols_)) &&
                "MUI Error [matrix_manipulation.h]: matrix size mismatch in copy function ");
      std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
      vec_temp = exist_mat.get_non_zero_elements();
      for (auto element : vec_temp) {
          if (std::abs(exist_mat.get_value(element.first, element.second)) >= std::numeric_limits<VTYPE>::min())
              matrix_[std::make_pair(element.first, element.second)] = exist_mat.get_value(element.first, element.second);
      }
}

// Member function to get a segment of a sparse_matrix
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::segment(ITYPE r_start, ITYPE r_end, ITYPE c_start, ITYPE c_end) {
      // get segment data from the existing matrix
      assert((r_end >= r_start) &&
              "MUI Error [matrix_manipulation.h]: segment function r_end has to be larger or equals to r_start");
      assert((c_end >= c_start) &&
              "MUI Error [matrix_manipulation.h]: segment function c_end has to be larger or equals to c_start");
      assert(((r_end < rows_) && (r_start >= 0) && (c_end < cols_) && (c_start >= 0)) &&
          "MUI Error [matrix_manipulation.h]: Matrix index out of range in segment function");
      sparse_matrix<ITYPE,VTYPE> res((r_end-r_start+1), (c_end-c_start+1));
      for (auto element : matrix_)
          if ((element.first.first >=r_start)  &&
              (element.first.first <=r_end)    &&
              (element.first.second >=c_start) &&
              (element.first.second <=c_end))
              res.set_value((element.first.first-r_start), (element.first.second-c_start), element.second);
      return res;
}

// Member function to insert an element
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_value(ITYPE r, ITYPE c, VTYPE val) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in set_value function");
    if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
        matrix_[std::make_pair(r, c)] = val;
    } else {
        if (matrix_.find(std::make_pair(r, c)) != matrix_.end()) {
            matrix_.erase(std::make_pair(r, c));
        }
    }
}

// Member function to insert the same value to all elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_value(VTYPE val) {
    if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
        for (ITYPE i = 0; i < rows_; ++i) {
            for (ITYPE j = 0; j < cols_; ++j) {
                matrix_[std::make_pair(i, j)] = val;
            }
        }
    } else {
        matrix_.clear();
    }
}

// Member function to swap two elements in a sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::swap_elements(ITYPE r1, ITYPE c1, ITYPE r2, ITYPE c2) {
    assert(((r1 < rows_) && (r1 >= 0) && (c1 < cols_) && (c1 >= 0) &&
            (r2 < rows_) && (r2 >= 0) && (c2 < cols_) && (c2 >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in set_value function");
    VTYPE temp = (*this).get_value(r1, c1);
    (*this).set_value(r1, c1, (*this).get_value(r2, c2));
    (*this).set_value(r2, c2, temp);
}

// Member function to set all elements to zero and empty the sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_zero() {
    matrix_.clear();
}

// Member function to add scalar to a specific elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::add_scalar(ITYPE r, ITYPE c, VTYPE value) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in add_scalar function");
    // check if the element exists
    if (matrix_.find(std::make_pair(r, c)) != matrix_.end()) {
        matrix_[std::make_pair(r, c)] += value;
    } else {
        matrix_[std::make_pair(r, c)] = value;
    }
}

// Member function to subtract a scalar from a specific elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::subtract_scalar(ITYPE r, ITYPE c, VTYPE value) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in subtract_scalar function");
    // check if the element exists
    if (matrix_.find(std::make_pair(r, c)) != matrix_.end()) {
        matrix_[std::make_pair(r, c)] -= value;
    } else {
        matrix_[std::make_pair(r, c)] = -value;
    }
}

// Overloaded assignment operator
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>& sparse_matrix<ITYPE,VTYPE>::operator=(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    if (this != &exist_mat) { // check for self-assignment
        // copy the values from the other matrix to this matrix
        assert(matrix_.empty() &&
                  "MUI Error [matrix_manipulation.h]: assignment operator '=' only works for empty (all zero elements) matrix");

        if((rows_ == 0) && (cols_ == 0)){
            rows_ = exist_mat.rows_;
            cols_ = exist_mat.cols_;
        }

        assert(((rows_ == exist_mat.rows_) && (cols_ == exist_mat.cols_)) &&
                  "MUI Error [matrix_manipulation.h]: matrix size mismatch in assignment operator '='");
        (*this).copy(exist_mat);
    }
    return *this;
}

// Member function to sort the entries by row and column for sparse matrix with COO format
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::sort_coo(bool is_row_major, bool deduplication, const std::string &deduplication_mode) {

    assert((matrix_format_ == format::COO) &&
      "MUI Error [matrix_manipulation.h]: sort_coo_row_major() can only applied to sparse matrix with COO format. Please consider to convert the format into COO before use this function.");

    if (matrix_coo.values_.size() <= 1) return;

	std::string deduplication_mode_trim = string_to_lower(trim(deduplication_mode));

    if ((deduplication_mode_trim != "sum") or (deduplication_mode_trim != "overwrite")) {
    	std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised deduplication mode: " << deduplication_mode << " for sort_coo_row_major() function" << std::endl;
		std::cerr << "    Please set the deduplication mode as:" << std::endl;
		std::cerr << "    'sum': Sum up values for all duplicated elements" << std::endl;
		std::cerr << "    'overwrite' (default): Keeps only the value of the last duplicated element" << std::endl;
		std::abort();
    }

    // Create a vector of indices to hold the original positions
    std::vector<ITYPE> sorted_indices(matrix_coo.values_.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);

    // Sort the indices based on row and column indices
    if (is_row_major) {
        std::sort(sorted_indices.begin(), sorted_indices.end(),
                  [&matrix_coo](ITYPE i, ITYPE j) {
                      return (matrix_coo.row_indices_[i] == matrix_coo.row_indices_[j])
                             ? matrix_coo.col_indices_[i] < matrix_coo.col_indices_[j]
                             : matrix_coo.row_indices_[i] < matrix_coo.row_indices_[j];
                  });
    } else {
        std::sort(sorted_indices.begin(), sorted_indices.end(),
                  [&matrix_coo](ITYPE i, ITYPE j) {
                      return (matrix_coo.col_indices_[i] == matrix_coo.col_indices_[j])
                             ? matrix_coo.row_indices_[i] < matrix_coo.row_indices_[j]
                             : matrix_coo.col_indices_[i] < matrix_coo.col_indices_[j];
                  });
    }

	// Reorder the COO entries based on the sorted indices
	std::vector<VTYPE> sorted_values;
	std::vector<ITYPE> sorted_row_indices;
	std::vector<ITYPE> sorted_column_indices;

	sorted_values.reserve(sorted_indices.size());
	sorted_row_indices.reserve(sorted_indices.size());
	sorted_column_indices.reserve(sorted_indices.size());

	for (ITYPE i = 0; i < sorted_indices.size(); ++i) {
		if (std::abs(matrix_coo.values_[sorted_indices[i]]) >= std::numeric_limits<VTYPE>::min()) {
			sorted_values.emplace_back(matrix_coo.values_[sorted_indices[i]]);
			sorted_row_indices.emplace_back(matrix_coo.row_indices_[sorted_indices[i]]);
			sorted_column_indices.emplace_back(matrix_coo.col_indices_[sorted_indices[i]]);
		}
	}

	if (deduplication) {

		// Deduplicate the COO entries based on the sorted matrix
		std::vector<VTYPE> deduplicated_values;
		std::vector<ITYPE> deduplicated_row_indices;
		std::vector<ITYPE> deduplicated_column_indices;

		deduplicated_values.reserve(sorted_indices.size());
		deduplicated_row_indices.reserve(sorted_indices.size());
		deduplicated_column_indices.reserve(sorted_indices.size());

		ITYPE prev_row = sorted_row_indices[0];
		ITYPE prev_col = sorted_column_indices[0];
		VTYPE sum_value = sorted_values[0];
		VTYPE last_value = sorted_values[0];

		for (ITYPE i = 1; i < sorted_indices.size(); ++i) {

			ITYPE curr_row = sorted_row_indices[i];
			ITYPE curr_col = sorted_column_indices[i];

			if ((curr_row == prev_row) && (curr_col == prev_col)) {
				sum_value += sorted_values[i];
				last_value = sorted_values[i];
			} else {

				prev_row = curr_row;
				prev_col = curr_col;
				sum_value = sorted_values[i];
				last_value = sorted_values[i];

				if ((deduplication_mode_trim == "sum") && (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
					deduplicated_values.emplace_back(sum_value);
					deduplicated_row_indices.emplace_back(prev_row);
					deduplicated_column_indices.emplace_back(prev_col);
				} else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
					deduplicated_values.emplace_back(last_value);
					deduplicated_row_indices.emplace_back(prev_row);
					deduplicated_column_indices.emplace_back(prev_col);
				}
			}
		}

		if ((deduplication_mode_trim == "sum") && (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
			deduplicated_values.emplace_back(sum_value);
			deduplicated_row_indices.emplace_back(prev_row);
			deduplicated_column_indices.emplace_back(prev_col);
		} else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
			deduplicated_values.emplace_back(last_value);
			deduplicated_row_indices.emplace_back(prev_row);
			deduplicated_column_indices.emplace_back(prev_col);
		}

		// Update the COOData struct with the deduplicated entries
		std::swap(matrix_coo.values_, deduplicated_values);
		std::swap(matrix_coo.row_indices_, deduplicated_row_indices);
		std::swap(matrix_coo.col_indices_, deduplicated_column_indices);

		nnz_ = matrix_coo.values_.size();

	} else {

		// Update the COOData struct with the sorted entries
		std::swap(matrix_coo.values_, sorted_values);
		std::swap(matrix_coo.row_indices_, sorted_row_indices);
		std::swap(matrix_coo.col_indices_, sorted_column_indices);

		nnz_ = matrix_coo.values_.size();

	}
}

template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::coo_to_csr() {

    assert((matrix_format_ == format::COO) &&
      "MUI Error [matrix_manipulation.h]: coo_to_csr() can only applied to sparse matrix with COO format.");

    if (matrix_coo.values_.size() == 0) return;

	// Determine the number of non-zero entries in each row
    std::vector<ITYPE> rowPtr(rows_+1, 0);
    for (ITYPE i = 0; i < matrix_coo.row_indices_.size(); ++i) {
    	rowPtr[matrix_coo.row_indices_[i]]++;
    }

	for(ITYPE i = 0, accumulator = 0; i < rows_; ++i){
		ITYPE temp = rowPtr[i];
		rowPtr[i] = accumulator;
		accumulator += temp;
	}
	rowPtr[rows_] = nnz_;

	// Fill the CSRData struct
	matrix_csr.values_.reserve(nnz_);
	matrix_csr.row_ptrs_.reserve(rows_+1);
	matrix_csr.col_indices_.reserve(nnz_);

	std::swap(matrix_csr.values_, matrix_coo.values_);
	std::swap(matrix_csr.row_ptrs_, rowPtr);
	std::swap(matrix_csr.col_indices_, matrix_coo.values_);

	// Deallocate the memory
    matrix_coo.values_.clear();
    matrix_coo.row_indices_.clear();
    matrix_coo.col_indices_.clear();
    rowPtr.clear();

    // Reset the matrix format
    matrix_format_ = format::CSR
}

template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::coo_to_csc() {

    assert((matrix_format_ == format::COO) &&
      "MUI Error [matrix_manipulation.h]: coo_to_csc() can only applied to sparse matrix with COO format.");

    if (matrix_coo.values_.size() == 0) return;

	// Determine the number of non-zero entries in each column
    std::vector<ITYPE> colPtr(cols_+1, 0);
    for (ITYPE i = 0; i < matrix_coo.col_indices_.size(); ++i) {
    	colPtr[matrix_coo.col_indices_[i]]++;
    }

	for(ITYPE i = 0, accumulator = 0; i < cols_; ++i){
		ITYPE temp = colPtr[i];
		colPtr[i] = accumulator;
		accumulator += temp;
	}
	colPtr[cols_] = nnz_;

	// Fill the CSCData struct
	matrix_csc.values_.reserve(nnz_);
	matrix_csc.row_indices_.reserve(nnz_);
	matrix_csc.col_ptrs_.reserve(cols_+1);

	std::swap(matrix_csc.values_, matrix_coo.values_);
	std::swap(matrix_csc.row_indices_, matrix_coo.row_indices_);
	std::swap(matrix_csc.col_ptrs_, colPtr);

	// Deallocate the memory
    matrix_coo.values_.clear();
    matrix_coo.row_indices_.clear();
    matrix_coo.col_indices_.clear();
    colPtr.clear();

    // Reset the matrix format
    matrix_format_ = format::CSC
}

template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csr_to_coo() {

    assert((matrix_format_ == format::CSR) &&
      "MUI Error [matrix_manipulation.h]: csr_to_coo() can only applied to sparse matrix with CSR format.");

    if (matrix_csr.values_.size() == 0) return;

    matrix_coo.values_.reserve(nnz_);
    matrix_coo.row_indices_.reserve(nnz_);
    matrix_coo.col_indices_.reserve(nnz_);

    // Populate the matrix_coo vectors from the matrix_csr vectors
    for (ITYPE row = 0; row < rows_; ++row) {
        for (ITYPE row_nnz = matrix_csr.row_ptrs_[row]; row_nnz < matrix_csr.row_ptrs_[row + 1]; ++row_nnz) {
        	matrix_coo.values_.emplace_back(matrix_csr.values_[row_nnz]);
        	matrix_coo.row_indices_.emplace_back(row);
        	matrix_coo.col_indices_.emplace_back(matrix_csr.col_indices_[row_nnz]);
        }
    }

	// Deallocate the memory
    matrix_csr.values_.clear();
    matrix_csr.row_ptrs_.clear();
    matrix_csr.col_indices_.clear();

    // Reset the matrix format
    matrix_format_ = format::COO
}

template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csr_to_csc() {

    assert((matrix_format_ == format::CSR) &&
      "MUI Error [matrix_manipulation.h]: csr_to_csc() can only applied to sparse matrix with CSR format.");

    if (matrix_csr.values_.size() == 0) return;

    matrix_csc.values_.reserve(nnz_);
    matrix_csc.row_indices_.reserve(nnz_);
    matrix_csc.col_ptrs_.reserve(cols_+1);

	// Determine the number of non-zero entries in each column
    for (ITYPE i = 0; i < matrix_csr.col_indices_.size(); ++i) {
    	matrix_csc.col_ptrs_[matrix_csr.col_indices_[i]]++;
    }

	for(ITYPE i = 0, accumulator = 0; i < cols_; ++i){
		ITYPE temp = matrix_csc.col_ptrs_[i];
		matrix_csc.col_ptrs_[i] = accumulator;
		accumulator += temp;
	}

	matrix_csc.col_ptrs_[cols_] = nnz_;

    for (ITYPE row = 0; row < rows_; ++row) {
        for (ITYPE row_nnz = matrix_csr.row_ptrs_[row]; row_nnz < matrix_csr.row_ptrs_[row + 1]; ++row_nnz) {
        	ITYPE col = matrix_csr.col_indices_[row_nnz];
        	ITYPE dest = matrix_csc.col_ptrs_[col];
        	matrix_csc.row_indices_[dest] = row;
        	matrix_csc.values_[dest] = matrix_csr.values_[row_nnz];
        	matrix_csc.col_ptrs_[col]++;
        }
    }

    for (ITYPE col = 0, last = 0; col <= cols_; ++col) {
    	ITYPE temp = matrix_csc.col_ptrs_[col];
    	matrix_csc.col_ptrs_[col] = last;
    	last = temp;
    }

	// Deallocate the memory
    matrix_csr.values_.clear();
    matrix_csr.row_ptrs_.clear();
    matrix_csr.col_indices_.clear();

    // Reset the matrix format
    matrix_format_ = format::CSC

}

template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csc_to_coo() {

    assert((matrix_format_ == format::CSC) &&
      "MUI Error [matrix_manipulation.h]: csc_to_coo() can only applied to sparse matrix with CSC format.");

    if (matrix_csc.values_.size() == 0) return;

    matrix_coo.values_.reserve(nnz_);
    matrix_coo.row_indices_.reserve(nnz_);
    matrix_coo.col_indices_.reserve(nnz_);

    // Populate the matrix_coo vectors from the matrix_csc vectors
    for (ITYPE col = 0; col < cols_; ++col) {
        for (ITYPE col_nnz = matrix_csc.col_ptrs_[col]; col_nnz < matrix_csc.col_ptrs_[col + 1]; ++col_nnz) {
        	matrix_coo.values_.emplace_back(matrix_csc.values_[col_nnz]);
        	matrix_coo.row_indices_.emplace_back(matrix_csc.row_indices_[col_nnz]);
        	matrix_coo.col_indices_.emplace_back(col);
        }
    }

	// Deallocate the memory
    matrix_csc.values_.clear();
    matrix_csc.row_indices_.clear();
    matrix_csc.col_ptrs_.clear();

    // Reset the matrix format
    matrix_format_ = format::COO
}

template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csc_to_csr() {

    assert((matrix_format_ == format::CSC) &&
      "MUI Error [matrix_manipulation.h]: csc_to_csr() can only applied to sparse matrix with CSC format.");

    if (matrix_csc.values_.size() == 0) return;

    matrix_csr.values_.reserve(nnz_);
    matrix_csr.row_ptrs_.reserve(rows_+1);
    matrix_csr.col_indices_.reserve(nnz_);

	// Determine the number of non-zero entries in each row
    for (ITYPE i = 0; i < matrix_csc.row_indices_.size(); ++i) {
    	matrix_csr.row_ptrs_[matrix_csc.row_indices_[i]]++;
    }

	for(ITYPE i = 0, accumulator = 0; i < rows_; ++i){
		ITYPE temp = matrix_csr.row_ptrs_[i];
		matrix_csr.row_ptrs_[i] = accumulator;
		accumulator += temp;
	}

	matrix_csr.row_ptrs_[rows_] = nnz_;

    for (ITYPE col = 0; col < cols_; ++col) {
        for (ITYPE col_nnz = matrix_csc.col_ptrs_[col]; col_nnz < matrix_csc.col_ptrs_[col + 1]; ++col_nnz) {
        	ITYPE row = matrix_csc.row_indices_[col_nnz];
        	ITYPE dest = matrix_csr.row_ptrs_[row];
        	matrix_csr.col_indices_[dest] = col;
        	matrix_csr.values_[dest] = matrix_csc.values_[col_nnz];
        	matrix_csr.row_ptrs_[row]++;
        }
    }

    for (ITYPE row = 0, last = 0; row <= rows_; ++row) {
    	ITYPE temp = matrix_csr.row_ptrs_[row];
    	matrix_csr.row_ptrs_[row] = last;
    	last = temp;
    }

	// Deallocate the memory
    matrix_csc.values_.clear();
    matrix_csc.row_indices_.clear();
    matrix_csc.col_ptrs_.clear();

    // Reset the matrix format
    matrix_format_ = format::CSR

}

// Member function to convert the format of the sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::format_conversion(const std::string &format, bool sort, bool deduplication, const std::string &deduplication_mode) {

	std::string matrix_format = string_to_upper(trim(format));

    if (matrix_format_ == format::COO) {

    	if (matrix_format == "COO") {

    		std::out << "MUI [matrix_manipulation.h]: Convert matrix format from COO to COO (do nothing)." << std::endl;

    	} else if (matrix_format == "CSR") {

    		if (sort) {
    			this->sparse_matrix<ITYPE,VTYPE>::sort_coo(true, deduplication, deduplication_mode);
    		}

    		this->sparse_matrix<ITYPE,VTYPE>::coo_to_csr();

    	} else if (matrix_format == "CSC") {

    		if (sort) {
    			this->sparse_matrix<ITYPE,VTYPE>::sort_coo(false, deduplication, deduplication_mode);
    		}

    		this->sparse_matrix<ITYPE,VTYPE>::coo_to_csc();

    	} else {

            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised format type: " << format << " for matrix constructor" << std::endl;
            std::cerr << "    Please set the format string as:" << std::endl;
            std::cerr << "    'COO': COOrdinate format" << std::endl;
            std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
            std::abort();

    	}

    } else if (matrix_format_ == format::CSR) {

    	if (matrix_format == "COO") {

    		this->sparse_matrix<ITYPE,VTYPE>::csr_to_coo();

    		if (sort) {
    			this->sparse_matrix<ITYPE,VTYPE>::sort_coo(false, deduplication, deduplication_mode);
    		}

    	} else if (matrix_format == "CSR") {

    		std::out << "MUI [matrix_manipulation.h]: Convert matrix format from CSR to CSR (do nothing)." << std::endl;

    	} else if (matrix_format == "CSC") {

    		this->sparse_matrix<ITYPE,VTYPE>::csr_to_csc();

    	} else {

            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised format type: " << format << " for matrix constructor" << std::endl;
            std::cerr << "    Please set the format string as:" << std::endl;
            std::cerr << "    'COO': COOrdinate format" << std::endl;
            std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
            std::abort();

    	}

    } else if (matrix_format_ == format::CSC) {

    	if (matrix_format == "COO") {

    		this->sparse_matrix<ITYPE,VTYPE>::csc_to_coo();

    		if (sort) {
    			this->sparse_matrix<ITYPE,VTYPE>::sort_coo(true, deduplication, deduplication_mode);
    		}

    	} else if (matrix_format == "CSR") {

    		this->sparse_matrix<ITYPE,VTYPE>::csc_to_csr();

    	} else if (matrix_format == "CSC") {

    		std::out << "MUI [matrix_manipulation.h]: Convert matrix format from CSC to CSC (do nothing)." << std::endl;

    	} else {

            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised format type: " << format << " for matrix constructor" << std::endl;
            std::cerr << "    Please set the format string as:" << std::endl;
            std::cerr << "    'COO': COOrdinate format" << std::endl;
            std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
            std::abort();

    	}

    } else {

        std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format: " << matrix_format_ << " for matrix constructor" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();

    }
}

} // linalg
} // mui

#endif /* MUI_MATRIX_MANIPULATION_H_ */
