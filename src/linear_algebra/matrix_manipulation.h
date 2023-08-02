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
#include <numeric>

namespace mui {
namespace linalg {

// **************************************************
// ************ Public member functions *************
// **************************************************

// Member function to resize an all-zero or null matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::resize(ITYPE r, ITYPE c) {
    assert(((this->non_zero_elements_count()) == 0) &&
            "MUI Error [matrix_manipulation.h]: resize function only works for all-zero matrix");
    rows_ = r;
    cols_ = c;

    if (matrix_format_ == format::CSR) {
        matrix_csr.row_ptrs_.clear();
        matrix_csr.row_ptrs_.resize((r+1), 0);
    } else if (matrix_format_ == format::CSC) {
        matrix_csc.col_ptrs_.clear();
        matrix_csc.col_ptrs_.resize((c+1), 0);
    }
}

// Member function to copy a sparse_matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::copy(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {

    // Copy the data from the existing matrix
    assert(this->empty() &&
            "MUI Error [matrix_manipulation.h]: copy function only works for empty (all zero elements) matrix");
    assert((((rows_ == exist_mat.rows_) && (cols_ == exist_mat.cols_)) || ((rows_ == 0) && (cols_ == 0))) &&
            "MUI Error [matrix_manipulation.h]: matrix size mismatch in copy function ");

    exist_mat.assert_valid_vector_size("matrix_manipulation.h", "copy()");

    std::string format_store = this->get_format();

    if (exist_mat.nnz_>0) {

        if ((rows_ == 0) && (cols_ == 0))
          this->resize(exist_mat.rows_, exist_mat.cols_);

        nnz_ = exist_mat.nnz_;

        if (exist_mat.matrix_format_ == format::COO) {

            matrix_coo.values_.reserve(exist_mat.matrix_coo.values_.size());
            matrix_coo.row_indices_.reserve(exist_mat.matrix_coo.row_indices_.size());
            matrix_coo.col_indices_.reserve(exist_mat.matrix_coo.col_indices_.size());

            matrix_coo.values_ = std::vector<VTYPE>(exist_mat.matrix_coo.values_.begin(), exist_mat.matrix_coo.values_.end());
            matrix_coo.row_indices_ = std::vector<ITYPE>(exist_mat.matrix_coo.row_indices_.begin(), exist_mat.matrix_coo.row_indices_.end());
            matrix_coo.col_indices_ = std::vector<ITYPE>(exist_mat.matrix_coo.col_indices_.begin(), exist_mat.matrix_coo.col_indices_.end());

            matrix_format_ = exist_mat.matrix_format_;

        } else if (exist_mat.matrix_format_ == format::CSR) {

            matrix_csr.values_.reserve(exist_mat.matrix_csr.values_.size());
            matrix_csr.row_ptrs_.reserve(exist_mat.matrix_csr.row_ptrs_.size());
            matrix_csr.col_indices_.reserve(exist_mat.matrix_csr.col_indices_.size());

            matrix_csr.values_ = std::vector<VTYPE>(exist_mat.matrix_csr.values_.begin(), exist_mat.matrix_csr.values_.end());
            matrix_csr.row_ptrs_ = std::vector<ITYPE>(exist_mat.matrix_csr.row_ptrs_.begin(), exist_mat.matrix_csr.row_ptrs_.end());
            matrix_csr.col_indices_ = std::vector<ITYPE>(exist_mat.matrix_csr.col_indices_.begin(), exist_mat.matrix_csr.col_indices_.end());

            matrix_format_ = exist_mat.matrix_format_;

        } else if (exist_mat.matrix_format_ == format::CSC) {

            matrix_csc.values_.reserve(exist_mat.matrix_csc.values_.size());
            matrix_csc.row_indices_.reserve(exist_mat.matrix_csc.row_indices_.size());
            matrix_csc.col_ptrs_.reserve(exist_mat.matrix_csc.col_ptrs_.size());

            matrix_csc.values_ = std::vector<VTYPE>(exist_mat.matrix_csc.values_.begin(), exist_mat.matrix_csc.values_.end());
            matrix_csc.row_indices_ = std::vector<ITYPE>(exist_mat.matrix_csc.row_indices_.begin(), exist_mat.matrix_csc.row_indices_.end());
            matrix_csc.col_ptrs_ = std::vector<ITYPE>(exist_mat.matrix_csc.col_ptrs_.begin(), exist_mat.matrix_csc.col_ptrs_.end());

            matrix_format_ = exist_mat.matrix_format_;

        } else {
              std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised exist_mat format" << std::endl;
              std::cerr << "    Please set the matrix_format_ as:" << std::endl;
              std::cerr << "    format::COO: COOrdinate format" << std::endl;
              std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
              std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
              std::abort();
          }
    }

    if (this->get_format() != format_store)
        this->format_conversion(format_store, true, true, "overwrite");

    this->assert_valid_vector_size("matrix_manipulation.h", "copy()");

}

// Member function to get a segment of a sparse_matrix
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> sparse_matrix<ITYPE,VTYPE>::segment(ITYPE r_start, ITYPE r_end, ITYPE c_start, ITYPE c_end, bool performSortAndUniqueCheck) {
      // get segment data from the existing matrix
      assert((r_end >= r_start) &&
              "MUI Error [matrix_manipulation.h]: segment function r_end has to be larger or equals to r_start");
      assert((c_end >= c_start) &&
              "MUI Error [matrix_manipulation.h]: segment function c_end has to be larger or equals to c_start");
      assert(((r_end < rows_) && (r_start >= 0) && (c_end < cols_) && (c_start >= 0)) &&
          "MUI Error [matrix_manipulation.h]: Matrix index out of range in segment function");

      if (performSortAndUniqueCheck) {
          if (!(this->is_sorted_unique("matrix_manipulation.h", "segment()"))){
              if(matrix_format_ == format::COO) {
                  this->sort_coo(true, true, "overwrite");
              } else if (matrix_format_ == format::CSR) {
                  this->sort_csr(true, "overwrite");
              } else if (matrix_format_ == format::CSC) {
                  this->sort_csc(true, "overwrite");
              } else {
                    std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
                    std::cerr << "    Please set the matrix_format_ as:" << std::endl;
                    std::cerr << "    format::COO: COOrdinate format" << std::endl;
                    std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                    std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                    std::abort();
              }
          }
      }

    sparse_matrix<ITYPE,VTYPE> res((r_end-r_start+1), (c_end-c_start+1), this->get_format());

    if(matrix_format_ == format::COO) {

      // Iterate over the existing non-zero elements
      for (ITYPE i = 0; i < nnz_; ++i) {
          ITYPE row = matrix_coo.row_indices_[i];
          ITYPE col = matrix_coo.col_indices_[i];

          // Check if the current element is within the specified ranges
          if (row >= r_start && row <= r_end && col >= c_start && col <= c_end) {
              // Calculate the indices for the segment
              ITYPE subRow = row - r_start;
              ITYPE subCol = col - c_start;

              res.matrix_coo.values_.reserve(res.matrix_coo.values_.size()+1);
              res.matrix_coo.row_indices_.reserve(res.matrix_coo.row_indices_.size()+1);
              res.matrix_coo.col_indices_.reserve(res.matrix_coo.col_indices_.size()+1);

              // Add the element to the segment matrix_coo struct
              res.matrix_coo.values_.emplace_back(matrix_coo.values_[i]);
              res.matrix_coo.row_indices_.emplace_back(subRow);
              res.matrix_coo.col_indices_.emplace_back(subCol);
              res.nnz_++;

          }
      }

    } else if (matrix_format_ == format::CSR) {

        res.matrix_csr.row_ptrs_.reserve(r_end-r_start+2);

        // Iterate over the row pointers and column indices of the existing non-zero elements
        for (ITYPE row = r_start; row <= r_end; ++row) {
            // Get the starting and ending indices for the current row
            ITYPE start = matrix_csr.row_ptrs_[row];
            ITYPE end = matrix_csr.row_ptrs_[row + 1];

            // Iterate over the non-zero elements in the current row
            for (ITYPE j = start; j < end; ++j) {
                // Get the column index of the current element
                ITYPE col = matrix_csr.col_indices_[j];

                // Check if the current element is within the specified column range
                if (col >= c_start && col <= c_end) {
                    // Calculate the indices for the sub-block
                    ITYPE subCol = col - c_start;

                    res.matrix_csr.values_.reserve(res.matrix_csr.values_.size()+1);
                    res.matrix_csr.col_indices_.reserve(res.matrix_csr.col_indices_.size()+1);

                    // Add the element to the segment matrix_csr struct
                    res.matrix_csr.values_.emplace_back(matrix_csr.values_[j]);
                    res.matrix_csr.col_indices_.emplace_back(subCol);
                    res.nnz_++;
                }
            }
            // Update the row pointer for the segment
            res.matrix_csr.row_ptrs_[row - r_start + 1] = res.matrix_csr.col_indices_.size();
        }

    } else if (matrix_format_ == format::CSC) {

        res.matrix_csc.col_ptrs_.reserve(c_end-c_start+2);

        // Iterate over the column pointers and row indices of the existing non-zero elements
        for (ITYPE col = c_start; col <= c_end; ++col) {
            // Get the starting and ending indices for the current column
            ITYPE start = matrix_csc.col_ptrs_[col];
            ITYPE end = matrix_csc.col_ptrs_[col + 1];

            // Iterate over the non-zero elements in the current column
            for (ITYPE j = start; j < end; ++j) {
                // Get the row index of the current element
                ITYPE row = matrix_csc.row_indices_[j];

                // Check if the current element is within the specified row range
                if (row >= r_start && row <= r_end) {
                    // Calculate the indices for the sub-block
                    ITYPE subRow = row - r_start;

                    res.matrix_csc.values_.reserve(res.matrix_csc.values_.size()+1);
                    res.matrix_csc.row_indices_.reserve(res.matrix_csc.row_indices_.size()+1);

                    // Add the element to the segment matrix_csc struct
                    res.matrix_csc.values_.emplace_back(matrix_csc.values_[j]);
                    res.matrix_csc.row_indices_.emplace_back(subRow);
                    res.nnz_++;
                }
            }

            // Update the column pointer for the segment
            res.matrix_csc.col_ptrs_[col - c_start + 1] = res.matrix_csc.row_indices_.size();
        }

    } else {
      std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
      std::cerr << "    Please set the matrix_format_ as:" << std::endl;
      std::cerr << "    format::COO: COOrdinate format" << std::endl;
      std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
      std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
      std::abort();
    }

    res.assert_valid_vector_size("matrix_manipulation.h", "segment()");

    return res;
}

// Member function to insert an element
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_value(ITYPE r, ITYPE c, VTYPE val, bool performSortAndUniqueCheck) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in set_value function");

    if (performSortAndUniqueCheck) {
        if (!(this->is_sorted_unique("matrix_manipulation.h", "set_value()"))){
            if(matrix_format_ == format::COO) {
                this->sort_coo(true, true, "overwrite");
            } else if (matrix_format_ == format::CSR) {
                this->sort_csr(true, "overwrite");
            } else if (matrix_format_ == format::CSC) {
                this->sort_csc(true, "overwrite");
            } else {
                  std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
                  std::cerr << "    Please set the matrix_format_ as:" << std::endl;
                  std::cerr << "    format::COO: COOrdinate format" << std::endl;
                  std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                  std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                  std::abort();
            }
        }
    }

    if(matrix_format_ == format::COO) {

        if (performSortAndUniqueCheck) {
            this->coo_element_operation(r, c, val, "overwrite", "matrix_manipulation.h", "set_value()");
        } else {
            this->coo_element_operation(r, c, val, "nonsort", "matrix_manipulation.h", "set_value()");
        }

    } else if (matrix_format_ == format::CSR) {

        this->csr_element_operation(r, c, val, "overwrite", "matrix_manipulation.h", "set_value()");

    } else if (matrix_format_ == format::CSC) {

        this->csc_element_operation(r, c, val, "overwrite", "matrix_manipulation.h", "set_value()");

    } else {
          std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
          std::cerr << "    Please set the matrix_format_ as:" << std::endl;
          std::cerr << "    format::COO: COOrdinate format" << std::endl;
          std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
          std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
          std::abort();
    }
}

// Member function to insert the same value to all elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_value(VTYPE val) {

    this->set_zero();

    if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {

        if(matrix_format_ == format::COO) {
            // Resize the vectors to hold all possible elements
            matrix_coo.values_.reserve(rows_*cols_);
            matrix_coo.row_indices_.reserve(rows_*cols_);
            matrix_coo.col_indices_.reserve(rows_*cols_);
            matrix_coo.values_.resize(rows_*cols_);
            matrix_coo.row_indices_.resize(rows_*cols_);
            matrix_coo.col_indices_.resize(rows_*cols_);
            // Fill all elements with the given value
            for (ITYPE i = 0; i < rows_; ++i) {
                std::fill(matrix_coo.row_indices_.begin()+(i*cols_), matrix_coo.row_indices_.end()+((i+1)*cols_), i);
                std::iota(matrix_coo.col_indices_.begin()+(i*cols_), matrix_coo.col_indices_.end()+((i+1)*cols_), 0);
                std::fill(matrix_coo.values_.begin()+(i*cols_), matrix_coo.values_.end()+((i+1)*cols_), val);
            }
        } else if (matrix_format_ == format::CSR) {
            // Resize the vectors to hold all possible elements
            matrix_csr.values_.reserve(rows_*cols_);
            matrix_csr.row_ptrs_.reserve(rows_+1);
            matrix_csr.col_indices_.reserve(rows_*cols_);
            matrix_csr.values_.resize(rows_*cols_);
            matrix_csr.row_ptrs_.resize(rows_+1);
            matrix_csr.col_indices_.resize(rows_*cols_);
            // Fill all elements with the given value
            matrix_csr.row_ptrs_[0] = 0;
            for (ITYPE i = 0; i < rows_; ++i) {
                std::iota(matrix_csr.col_indices_.begin()+(i*cols_), matrix_csr.col_indices_.end()+((i+1)*cols_), 0);
                std::fill(matrix_csr.values_.begin()+(i*cols_), matrix_csr.values_.end()+((i+1)*cols_), val);
                matrix_csr.row_ptrs_[i + 1] = (i+1)*cols_;
            }
        } else if (matrix_format_ == format::CSC) {
            // Resize the vectors to hold all possible elements
            matrix_csc.values_.reserve(rows_*cols_);
            matrix_csc.row_indices_.reserve(rows_*cols_);
            matrix_csc.col_ptrs_.reserve(cols_+1);
            matrix_csc.values_.resize(rows_*cols_);
            matrix_csc.row_indices_.resize(rows_*cols_);
            matrix_csc.col_ptrs_.resize(cols_+1);
            // Fill all elements with the given value
            matrix_csc.col_ptrs_[0] = 0;
            for (ITYPE i = 0; i < cols_; ++i) {
                std::iota(matrix_csc.row_indices_.begin()+(i*rows_), matrix_csc.row_indices_.end()+((i+1)*rows_), 0);
                std::fill(matrix_csc.values_.begin()+(i*rows_), matrix_csc.values_.end()+((i+1)*rows_), val);
                matrix_csc.col_ptrs_[i + 1] = (i+1)*rows_;
            }
        } else {
            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
            std::cerr << "    Please set the matrix_format_ as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }
        nnz_ = rows_*cols_;
    }
}

// Member function to swap two elements in a sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::swap_elements(ITYPE r1, ITYPE c1, ITYPE r2, ITYPE c2) {
    assert(((r1 < rows_) && (r1 >= 0) && (c1 < cols_) && (c1 >= 0) &&
            (r2 < rows_) && (r2 >= 0) && (c2 < cols_) && (c2 >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in swap_elements function");
    VTYPE temp = this->get_value(r1, c1);
    this->set_value(r1, c1, this->get_value(r2, c2));
    this->set_value(r2, c2, temp);
}

// Member function to set all elements to zero and empty the sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_zero() {
    // Clear the existing elements
    if(matrix_format_ == format::COO) {
        matrix_coo.values_.clear();
        matrix_coo.row_indices_.clear();
        matrix_coo.col_indices_.clear();
    } else if (matrix_format_ == format::CSR) {
        matrix_csr.values_.clear();
        matrix_csr.row_ptrs_.clear();
        matrix_csr.col_indices_.clear();
        matrix_csr.row_ptrs_.resize((rows_+1), 0);
    } else if (matrix_format_ == format::CSC) {
        matrix_csc.values_.clear();
        matrix_csc.row_indices_.clear();
        matrix_csc.col_ptrs_.clear();
        matrix_csc.col_ptrs_.resize((cols_+1), 0);
    } else {
        std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }
    nnz_ = 0;
}

// Member function to add scalar to a specific elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::add_scalar(ITYPE r, ITYPE c, VTYPE val, bool performSortAndUniqueCheck) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in add_scalar function");

    if (performSortAndUniqueCheck) {
        if (!(this->is_sorted_unique("matrix_manipulation.h", "add_scalar()"))){
            if(matrix_format_ == format::COO) {
                this->sort_coo(true, true, "overwrite");
            } else if (matrix_format_ == format::CSR) {
                this->sort_csr(true, "overwrite");
            } else if (matrix_format_ == format::CSC) {
                this->sort_csc(true, "overwrite");
            } else {
                  std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
                  std::cerr << "    Please set the matrix_format_ as:" << std::endl;
                  std::cerr << "    format::COO: COOrdinate format" << std::endl;
                  std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                  std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                  std::abort();
            }
        }
    }

    if(matrix_format_ == format::COO) {

        this->coo_element_operation(r, c, val, "plus", "matrix_manipulation.h", "add_scalar()");

    } else if (matrix_format_ == format::CSR) {

        this->csr_element_operation(r, c, val, "plus", "matrix_manipulation.h", "add_scalar()");

    } else if (matrix_format_ == format::CSC) {

        this->csc_element_operation(r, c, val, "plus", "matrix_manipulation.h", "add_scalar()");

    } else {
          std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
          std::cerr << "    Please set the matrix_format_ as:" << std::endl;
          std::cerr << "    format::COO: COOrdinate format" << std::endl;
          std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
          std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
          std::abort();
    }
}

// Member function to subtract a scalar from a specific elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::subtract_scalar(ITYPE r, ITYPE c, VTYPE val, bool performSortAndUniqueCheck) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in subtract_scalar function");
    if (performSortAndUniqueCheck) {
        if (!(this->is_sorted_unique("matrix_manipulation.h", "subtract_scalar()"))){
            if(matrix_format_ == format::COO) {
                this->sort_coo(true, true, "overwrite");
            } else if (matrix_format_ == format::CSR) {
                this->sort_csr(true, "overwrite");
            } else if (matrix_format_ == format::CSC) {
                this->sort_csc(true, "overwrite");
            } else {
                  std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
                  std::cerr << "    Please set the matrix_format_ as:" << std::endl;
                  std::cerr << "    format::COO: COOrdinate format" << std::endl;
                  std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                  std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                  std::abort();
            }
        }
    }

    if(matrix_format_ == format::COO) {

        this->coo_element_operation(r, c, val, "minus", "matrix_manipulation.h", "subtract_scalar()");

    } else if (matrix_format_ == format::CSR) {

        this->csr_element_operation(r, c, val, "minus", "matrix_manipulation.h", "subtract_scalar()");

    } else if (matrix_format_ == format::CSC) {

        this->csc_element_operation(r, c, val, "minus", "matrix_manipulation.h", "subtract_scalar()");

    } else {
          std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
          std::cerr << "    Please set the matrix_format_ as:" << std::endl;
          std::cerr << "    format::COO: COOrdinate format" << std::endl;
          std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
          std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
          std::abort();
    }
}

// Member function to multiply a scalar from a specific elements
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::multiply_scalar(ITYPE r, ITYPE c, VTYPE val, bool performSortAndUniqueCheck) {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_manipulation.h]: Matrix index out of range in multiply_scalar function");

    if (performSortAndUniqueCheck) {
        if (!(this->is_sorted_unique("matrix_manipulation.h", "multiply_scalar()"))){
            if(matrix_format_ == format::COO) {
                this->sort_coo(true, true, "overwrite");
            } else if (matrix_format_ == format::CSR) {
                this->sort_csr(true, "overwrite");
            } else if (matrix_format_ == format::CSC) {
                this->sort_csc(true, "overwrite");
            } else {
                  std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
                  std::cerr << "    Please set the matrix_format_ as:" << std::endl;
                  std::cerr << "    format::COO: COOrdinate format" << std::endl;
                  std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                  std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                  std::abort();
            }
        }
    }

    if(matrix_format_ == format::COO) {

        this->coo_element_operation(r, c, val, "multiply", "matrix_manipulation.h", "multiply_scalar()");

    } else if (matrix_format_ == format::CSR) {

        this->csr_element_operation(r, c, val, "multiply", "matrix_manipulation.h", "multiply_scalar()");

    } else if (matrix_format_ == format::CSC) {

        this->csc_element_operation(r, c, val, "multiply", "matrix_manipulation.h", "multiply_scalar()");

    } else {
          std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
          std::cerr << "    Please set the matrix_format_ as:" << std::endl;
          std::cerr << "    format::COO: COOrdinate format" << std::endl;
          std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
          std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
          std::abort();
    }
}

// Overloaded assignment operator
template <typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>& sparse_matrix<ITYPE,VTYPE>::operator=(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    if (this != &exist_mat) { // check for self-assignment
        // copy the values from the other matrix to this matrix
        assert(this->empty() &&
                  "MUI Error [matrix_manipulation.h]: assignment operator '=' only works for empty (all zero elements) matrix");

        if((rows_ == 0) && (cols_ == 0)){
            rows_ = exist_mat.rows_;
            cols_ = exist_mat.cols_;
        }

        assert(((rows_ == exist_mat.rows_) && (cols_ == exist_mat.cols_)) &&
                  "MUI Error [matrix_manipulation.h]: matrix size mismatch in assignment operator '='");

        std::string format_store = this->get_format();
        if (this->get_format() != exist_mat.get_format()) {
            this->format_conversion(exist_mat.get_format(), true, true, "overwrite");
        }

        // Copy the data from the existing matrix
        if (matrix_format_ == format::COO) {
            matrix_coo.values_ = std::vector<VTYPE>(exist_mat.matrix_coo.values_.begin(), exist_mat.matrix_coo.values_.end());
            matrix_coo.row_indices_ = std::vector<ITYPE>(exist_mat.matrix_coo.row_indices_.begin(), exist_mat.matrix_coo.row_indices_.end());
            matrix_coo.col_indices_ = std::vector<ITYPE>(exist_mat.matrix_coo.col_indices_.begin(), exist_mat.matrix_coo.col_indices_.end());
        } else if (matrix_format_ == format::CSR) {
            matrix_csr.values_ = std::vector<VTYPE>(exist_mat.matrix_csr.values_.begin(), exist_mat.matrix_csr.values_.end());
            matrix_csr.row_ptrs_ = std::vector<ITYPE>(exist_mat.matrix_csr.row_ptrs_.begin(), exist_mat.matrix_csr.row_ptrs_.end());
            matrix_csr.col_indices_ = std::vector<ITYPE>(exist_mat.matrix_csr.col_indices_.begin(), exist_mat.matrix_csr.col_indices_.end());
        } else if (matrix_format_ == format::CSC) {
            matrix_csc.values_ = std::vector<VTYPE>(exist_mat.matrix_csc.values_.begin(), exist_mat.matrix_csc.values_.end());
            matrix_csc.row_indices_ = std::vector<ITYPE>(exist_mat.matrix_csc.row_indices_.begin(), exist_mat.matrix_csc.row_indices_.end());
            matrix_csc.col_ptrs_ = std::vector<ITYPE>(exist_mat.matrix_csc.col_ptrs_.begin(), exist_mat.matrix_csc.col_ptrs_.end());
        } else {
            std::cerr << "MUI Error [matrix_ctor_dtor.h]: Unrecognised matrix format for matrix constructor" << std::endl;
            std::cerr << "    Please set the matrix_format_ as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }

        nnz_ = exist_mat.nnz_;

        if (this->get_format() != format_store) {
            this->format_conversion(format_store, true, true, "overwrite");
        }
    }
    return *this;
}

// Member function to convert the format of the sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::format_conversion(const std::string &format, bool performSortAndUniqueCheck, bool deduplication, const std::string &deduplication_mode) {

    std::string matrix_format = string_to_upper(trim(format));

    if (matrix_format_ == format::COO) {

        if (matrix_format == "COO") {

            std::cout << "MUI [matrix_manipulation.h]: Convert matrix format from COO to COO (do nothing)." << std::endl;

        } else if (matrix_format == "CSR") {

            if (performSortAndUniqueCheck) {
                if (!(this->is_sorted_unique("matrix_manipulation.h", "format_conversion()"))){
                    this->sparse_matrix<ITYPE,VTYPE>::sort_coo(true, deduplication, deduplication_mode);
                }
            }

            this->sparse_matrix<ITYPE,VTYPE>::coo_to_csr();

        } else if (matrix_format == "CSC") {

            this->sparse_matrix<ITYPE,VTYPE>::sort_coo(false, deduplication, deduplication_mode);

            this->sparse_matrix<ITYPE,VTYPE>::coo_to_csc();

        } else {

            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised format type: " << format << " for matrix format conversion" << std::endl;
            std::cerr << "    Please set the format string as:" << std::endl;
            std::cerr << "    'COO': COOrdinate format" << std::endl;
            std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
            std::abort();

        }

    } else if (matrix_format_ == format::CSR) {

        if (matrix_format == "COO") {

            if (performSortAndUniqueCheck) {
                if (!(this->is_sorted_unique("matrix_manipulation.h", "format_conversion()"))){
                    this->sparse_matrix<ITYPE,VTYPE>::sort_csr(deduplication, deduplication_mode);
                }
            }

            this->sparse_matrix<ITYPE,VTYPE>::csr_to_coo();

        } else if (matrix_format == "CSR") {

            std::cout << "MUI [matrix_manipulation.h]: Convert matrix format from CSR to CSR (do nothing)." << std::endl;

        } else if (matrix_format == "CSC") {

            if (performSortAndUniqueCheck) {
                if (!(this->is_sorted_unique("matrix_manipulation.h", "format_conversion()"))){
                    this->sparse_matrix<ITYPE,VTYPE>::sort_csr(deduplication, deduplication_mode);
                }
            }

            this->sparse_matrix<ITYPE,VTYPE>::csr_to_csc();

        } else {

            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised format type: " << format << " for matrix format conversion" << std::endl;
            std::cerr << "    Please set the format string as:" << std::endl;
            std::cerr << "    'COO': COOrdinate format" << std::endl;
            std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
            std::abort();

        }

    } else if (matrix_format_ == format::CSC) {

        if (matrix_format == "COO") {

            if (performSortAndUniqueCheck) {
                if (!(this->is_sorted_unique("matrix_manipulation.h", "format_conversion()"))){
                    this->sparse_matrix<ITYPE,VTYPE>::sort_csc(deduplication, deduplication_mode);
                }
            }

            this->sparse_matrix<ITYPE,VTYPE>::csc_to_coo();

        } else if (matrix_format == "CSR") {

            if (performSortAndUniqueCheck) {
                if (!(this->is_sorted_unique("matrix_manipulation.h", "format_conversion()"))){
                    this->sparse_matrix<ITYPE,VTYPE>::sort_csc(deduplication, deduplication_mode);
                }
            }

            this->sparse_matrix<ITYPE,VTYPE>::csc_to_csr();

        } else if (matrix_format == "CSC") {

            std::cout << "MUI [matrix_manipulation.h]: Convert matrix format from CSC to CSC (do nothing)." << std::endl;

        } else {

            std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised format type: " << format << " for matrix format conversion" << std::endl;
            std::cerr << "    Please set the format string as:" << std::endl;
            std::cerr << "    'COO': COOrdinate format" << std::endl;
            std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
            std::abort();

        }

    } else {

        std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format for matrix format conversion" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();

    }
}

// Member function to sort the entries for the sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::sort_deduplication(bool check_sorted_unique, bool deduplication, const std::string &deduplication_mode, bool is_row_major) {

    if (check_sorted_unique) {
        if (!(this->is_sorted_unique("matrix_manipulation.h", "sort_deduplication()"))){
            if(matrix_format_ == format::COO) {
                this->sort_coo(is_row_major, deduplication, deduplication_mode);
            } else if (matrix_format_ == format::CSR) {
                this->sort_csr(deduplication, deduplication_mode);
            } else if (matrix_format_ == format::CSC) {
                this->sort_csc(deduplication, deduplication_mode);
            } else {
                  std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
                  std::cerr << "    Please set the matrix_format_ as:" << std::endl;
                  std::cerr << "    format::COO: COOrdinate format" << std::endl;
                  std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
                  std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
                  std::abort();
            }
        }
    } else {
        if(matrix_format_ == format::COO) {
            this->sort_coo(is_row_major, deduplication, deduplication_mode);
        } else if (matrix_format_ == format::CSR) {
            this->sort_csr(deduplication, deduplication_mode);
        } else if (matrix_format_ == format::CSC) {
            this->sort_csc(deduplication, deduplication_mode);
        } else {
              std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised matrix format" << std::endl;
              std::cerr << "    Please set the matrix_format_ as:" << std::endl;
              std::cerr << "    format::COO: COOrdinate format" << std::endl;
              std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
              std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
              std::abort();
        }
    }
}


// **************************************************
// ********** Protected member functions ************
// **************************************************

// Protected member function to sort the entries by row and column for sparse matrix with COO format
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::sort_coo(bool is_row_major, bool deduplication, const std::string &deduplication_mode) {

    assert((matrix_format_ == format::COO) &&
      "MUI Error [matrix_manipulation.h]: sort_coo() can only applied to sparse matrix with COO format. Please consider to convert the format into COO before use this function.");

    if (matrix_coo.values_.size() <= 1) return;

    std::string deduplication_mode_trim = string_to_lower(trim(deduplication_mode));

    if ((deduplication_mode_trim != "plus") && (deduplication_mode_trim != "minus") && (deduplication_mode_trim != "multiply") && (deduplication_mode_trim != "overwrite")) {
        std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised deduplication mode: " << deduplication_mode << " for sort_coo() function" << std::endl;
        std::cerr << "    Please set the deduplication mode as:" << std::endl;
        std::cerr << "    'plus': Sum up values for all duplicated elements" << std::endl;
        std::cerr << "    'minus': Take the difference among duplicated elements according to their position (former element minus later element)" << std::endl;
        std::cerr << "    'multiply': Take the product among all duplicated elements" << std::endl;
        std::cerr << "    'overwrite' (default): Keeps only the value of the last duplicated element" << std::endl;
        std::abort();
    }

    // Create a vector of indices to hold the original positions
    std::vector<ITYPE> sorted_indices(matrix_coo.row_indices_.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);

    // Sort the indices based on row and column indices
    if (is_row_major) {
        std::sort(sorted_indices.begin(), sorted_indices.end(),
                  [&](ITYPE i, ITYPE j) {
                      return (matrix_coo.row_indices_[i] == matrix_coo.row_indices_[j])
                             ? matrix_coo.col_indices_[i] < matrix_coo.col_indices_[j]
                             : matrix_coo.row_indices_[i] < matrix_coo.row_indices_[j];
                  });
    } else {
        std::sort(sorted_indices.begin(), sorted_indices.end(),
                  [&](ITYPE i, ITYPE j) {
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

    for (ITYPE i = 0; i < static_cast<ITYPE>(sorted_indices.size()); ++i) {
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
        VTYPE difference_value = sorted_values[0];
        VTYPE product_value = sorted_values[0];
        VTYPE last_value = sorted_values[0];
        bool is_multiply_by_zero = true;

        for (ITYPE i = 1; i < static_cast<ITYPE>(sorted_indices.size()); ++i) {
           ITYPE curr_row = sorted_row_indices[i];
           ITYPE curr_col = sorted_column_indices[i];

           if (i != static_cast<ITYPE>((sorted_indices.size()-1))) {
               if (i == 1) {
                  if ((curr_row == prev_row) && (curr_col == prev_col)) {

                      sum_value += sorted_values[i];
                      difference_value -= sorted_values[i];
                      product_value *= sorted_values[i];
                      last_value = sorted_values[i];
                      is_multiply_by_zero = false;

                  } else {
                      if ((std::abs(sorted_values[0]) >= std::numeric_limits<VTYPE>::min()) && (deduplication_mode_trim != "multiply")) {
                          deduplicated_values.emplace_back(sorted_values[0]);
                          deduplicated_row_indices.emplace_back(prev_row);
                          deduplicated_column_indices.emplace_back(prev_col);
                      }

                      prev_row = curr_row;
                      prev_col = curr_col;
                      sum_value = sorted_values[i];
                      difference_value = sorted_values[i];
                      product_value = sorted_values[i];
                      last_value = sorted_values[i];
                      is_multiply_by_zero = true;

                  }
              } else {
                if ((curr_row == prev_row) && (curr_col == prev_col)) {

                    sum_value += sorted_values[i];
                    difference_value -= sorted_values[i];
                    product_value *= sorted_values[i];
                    last_value = sorted_values[i];
                    is_multiply_by_zero = false;

                } else {

                    if ((deduplication_mode_trim == "plus") && (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(sum_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(difference_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    } else if ((deduplication_mode_trim == "multiply") && (!is_multiply_by_zero) && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(product_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(last_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    }

                    prev_row = curr_row;
                    prev_col = curr_col;
                    sum_value = sorted_values[i];
                    difference_value = sorted_values[i];
                    product_value = sorted_values[i];
                    last_value = sorted_values[i];
                    is_multiply_by_zero = true;

                }
              }
            } else {
                if ((curr_row == prev_row) && (curr_col == prev_col)) {

                    sum_value += sorted_values[i];
                    difference_value -= sorted_values[i];
                    product_value *= sorted_values[i];
                    last_value = sorted_values[i];
                    is_multiply_by_zero = false;

                    if ((deduplication_mode_trim == "plus") && (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(sum_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(difference_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    } else if ((deduplication_mode_trim == "multiply") && (!is_multiply_by_zero) && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(product_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(last_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                    }

                } else {

                    if ((deduplication_mode_trim == "plus") && (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(sum_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                        deduplicated_values.emplace_back(sorted_values[i]);
                        deduplicated_row_indices.emplace_back(curr_row);
                        deduplicated_column_indices.emplace_back(curr_col);
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(difference_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                        deduplicated_values.emplace_back(sorted_values[i]);
                        deduplicated_row_indices.emplace_back(curr_row);
                        deduplicated_column_indices.emplace_back(curr_col);
                    } else if ((deduplication_mode_trim == "multiply") && (!is_multiply_by_zero) && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(product_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                        deduplicated_values.emplace_back(sorted_values[i]);
                        deduplicated_row_indices.emplace_back(curr_row);
                        deduplicated_column_indices.emplace_back(curr_col);
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        deduplicated_values.emplace_back(last_value);
                        deduplicated_row_indices.emplace_back(prev_row);
                        deduplicated_column_indices.emplace_back(prev_col);
                        deduplicated_values.emplace_back(sorted_values[i]);
                        deduplicated_row_indices.emplace_back(curr_row);
                        deduplicated_column_indices.emplace_back(curr_col);
                    }

                }
            }
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

// Protected member function to sort the entries for sparse matrix with CSR format
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::sort_csr(bool deduplication, const std::string &deduplication_mode) {

    assert((matrix_format_ == format::CSR) &&
      "MUI Error [matrix_manipulation.h]: sort_csr() can only applied to sparse matrix with CSR format. Please consider to convert the format into CSR before use this function.");

    if (matrix_csr.values_.size() <= 1) return;

    std::string deduplication_mode_trim = string_to_lower(trim(deduplication_mode));

    if ((deduplication_mode_trim != "plus") && (deduplication_mode_trim != "minus") && (deduplication_mode_trim != "multiply") && (deduplication_mode_trim != "overwrite")) {
        std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised deduplication mode: " << deduplication_mode << " for sort_csr() function" << std::endl;
        std::cerr << "    Please set the deduplication mode as:" << std::endl;
        std::cerr << "    'plus': Sum up values for all duplicated elements" << std::endl;
        std::cerr << "    'minus': Take the difference among duplicated elements according to their position (former element minus later element)" << std::endl;
        std::cerr << "    'multiply': Take the product among all duplicated elements" << std::endl;
        std::cerr << "    'overwrite' (default): Keeps only the value of the last duplicated element" << std::endl;
        std::abort();
    }

    // Temporary vectors for sorted and deduplicated data
    std::vector<ITYPE> sorted_col_indices;
    std::vector<VTYPE> sorted_values;

    // Iterate over each row
    for (ITYPE row = 0; row < rows_; ++row) {

        assert((matrix_csr.row_ptrs_[row] <= matrix_csr.row_ptrs_[row+1]) &&
              "MUI Error [matrix_manipulation.h]: sort_csr() unable to sort and deduplication the sparse matrix with CSR format as the row_ptrs_ vector is not in correct order.");

        // Sort CSR column indices in this row
        std::vector<std::pair<ITYPE,VTYPE>> sorted_indices_pair;
        sorted_indices_pair.reserve(matrix_csr.row_ptrs_[row+1] - matrix_csr.row_ptrs_[row]);

        for(ITYPE j = matrix_csr.row_ptrs_[row]; j < matrix_csr.row_ptrs_[row+1]; ++j){
            sorted_indices_pair.emplace_back(std::make_pair(matrix_csr.col_indices_[j],matrix_csr.values_[j]));
        }

        std::sort(sorted_indices_pair.begin(),sorted_indices_pair.end(),[](const std::pair<ITYPE,VTYPE>& x, const std::pair<ITYPE,VTYPE>& y) {
            return x.first <= y.first;
        });

        for(ITYPE j = matrix_csr.row_ptrs_[row], n = 0; j < matrix_csr.row_ptrs_[row+1]; ++j, ++n){
            matrix_csr.col_indices_[j] = sorted_indices_pair[n].first;
            matrix_csr.values_[j] = sorted_indices_pair[n].second;
        }

        // Clear the temporary vector
        sorted_indices_pair.clear();

        if (deduplication) {
            // Vector to track unique column indices
            std::vector<ITYPE> uniqueColumns;
            uniqueColumns.reserve(matrix_csr.row_ptrs_[row + 1]-matrix_csr.row_ptrs_[row]);

            // Vector to track values corresponding to unique column indices
            std::vector<VTYPE> uniqueValues;
            uniqueValues.reserve(matrix_csr.row_ptrs_[row + 1]-matrix_csr.row_ptrs_[row]);

            // Iterate over the elements in the row
            for (ITYPE i = matrix_csr.row_ptrs_[row]; i < matrix_csr.row_ptrs_[row + 1]; ++i) {
                ITYPE col = matrix_csr.col_indices_[i];
                VTYPE val = matrix_csr.values_[i];

                // Check if the column index already exists
                auto it = std::lower_bound(uniqueColumns.begin(), uniqueColumns.end(), col);
                ITYPE index = std::distance(uniqueColumns.begin(), it);
                if (it != uniqueColumns.end() && *it == col) {
                    // Column index already exists, overwrite the value
                    VTYPE sum_value = uniqueValues[index] + val;
                    VTYPE difference_value = uniqueValues[index] - val;
                    VTYPE product_value = uniqueValues[index] * val;
                    VTYPE last_value = val;
                    if ((deduplication_mode_trim == "plus")&& (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = sum_value;
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = difference_value;
                    } else if ((deduplication_mode_trim == "multiply") && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = product_value;
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = last_value;
                    }
                } else {
                    // Column index does not exist, insert it
                    VTYPE sum_value = val;
                    VTYPE difference_value = val;
                    VTYPE product_value = val;
                    VTYPE last_value = val;
                    if ((deduplication_mode_trim == "plus")&& (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueColumns.emplace_back(col);
                        uniqueValues.emplace_back(sum_value);
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueColumns.emplace_back(col);
                        uniqueValues.emplace_back(difference_value);
                    } else if ((deduplication_mode_trim == "multiply") && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueColumns.emplace_back(col);
                        uniqueValues.emplace_back(product_value);
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueColumns.emplace_back(col);
                        uniqueValues.emplace_back(last_value);
                    }
                }
            }

            // Insert the unique column indices and values into the temporary vector
            sorted_col_indices.insert(sorted_col_indices.end(), uniqueColumns.begin(), uniqueColumns.end());
            sorted_values.insert(sorted_values.end(), uniqueValues.begin(), uniqueValues.end());

            // Update the row_ptr vector
            matrix_csr.row_ptrs_[row + 1] = matrix_csr.row_ptrs_[row] + uniqueColumns.size();
        }
    }

    if (deduplication) {
        // Assign the sorted and deduplicated data back to the col_index and value vectors
        matrix_csr.col_indices_.clear();
        matrix_csr.col_indices_.shrink_to_fit();
        matrix_csr.col_indices_.swap(sorted_col_indices);
        matrix_csr.values_.clear();
        matrix_csr.values_.shrink_to_fit();
        matrix_csr.values_.swap(sorted_values);
        nnz_ = matrix_csr.values_.size();
    }

    // Clear the temporary vector
    sorted_col_indices.clear();
    sorted_values.clear();

}

// Protected member function to sort the entries for sparse matrix with CSC format
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::sort_csc(bool deduplication, const std::string &deduplication_mode) {

    assert((matrix_format_ == format::CSC) &&
      "MUI Error [matrix_manipulation.h]: sort_csc() can only applied to sparse matrix with CSC format. Please consider to convert the format into CSC before use this function.");

    if (matrix_csc.values_.size() <= 1) return;

    std::string deduplication_mode_trim = string_to_lower(trim(deduplication_mode));

    if ((deduplication_mode_trim != "plus") && (deduplication_mode_trim != "minus") && (deduplication_mode_trim != "multiply") && (deduplication_mode_trim != "overwrite")) {
        std::cerr << "MUI Error [matrix_manipulation.h]: Unrecognised deduplication mode: " << deduplication_mode << " for sort_csc() function" << std::endl;
        std::cerr << "    Please set the deduplication mode as:" << std::endl;
        std::cerr << "    'plus': Sum up values for all duplicated elements" << std::endl;
        std::cerr << "    'minus': Take the difference among duplicated elements according to their position (former element minus later element)" << std::endl;
        std::cerr << "    'multiply': Take the product among all duplicated elements" << std::endl;
        std::cerr << "    'overwrite' (default): Keeps only the value of the last duplicated element" << std::endl;
        std::abort();
    }

    // Temporary vectors for sorted and deduplicated data
    std::vector<ITYPE> sorted_row_indices;
    std::vector<VTYPE> sorted_values;

    // Iterate over each column
    for (ITYPE col = 0; col < cols_; ++col) {

        assert((matrix_csc.col_ptrs_[col] <= matrix_csc.col_ptrs_[col+1]) &&
              "MUI Error [matrix_manipulation.h]: sort_csc() unable to sort and deduplication the sparse matrix with CSC format as the col_ptrs_ vector is not in correct order.");

        // Sort CSC row indices in this column
        std::vector<std::pair<ITYPE,VTYPE>> sorted_indices_pair;
        sorted_indices_pair.reserve(matrix_csc.col_ptrs_[col+1] - matrix_csc.col_ptrs_[col]);

        for(ITYPE j = matrix_csc.col_ptrs_[col]; j < matrix_csc.col_ptrs_[col+1]; ++j){
            sorted_indices_pair.emplace_back(std::make_pair(matrix_csc.row_indices_[j],matrix_csc.values_[j]));
        }

        std::sort(sorted_indices_pair.begin(),sorted_indices_pair.end(),[](const std::pair<ITYPE,VTYPE>& x, const std::pair<ITYPE,VTYPE>& y) {
            return x.first <= y.first;
        });

        for(ITYPE j = matrix_csc.col_ptrs_[col], n = 0; j < matrix_csc.col_ptrs_[col+1]; ++j, ++n){
            matrix_csc.row_indices_[j] = sorted_indices_pair[n].first;
            matrix_csc.values_[j] = sorted_indices_pair[n].second;
        }

        // Clear the temporary vector
        sorted_indices_pair.clear();

        if (deduplication) {
            // Vector to track unique row indices
            std::vector<ITYPE> uniqueRows;
            uniqueRows.reserve(matrix_csc.col_ptrs_[col + 1]-matrix_csc.col_ptrs_[col]);

            // Vector to track values corresponding to unique row indices
            std::vector<VTYPE> uniqueValues;
            uniqueValues.reserve(matrix_csc.col_ptrs_[col + 1]-matrix_csc.col_ptrs_[col]);

            // Iterate over the elements in the column
            for (ITYPE i = matrix_csc.col_ptrs_[col]; i < matrix_csc.col_ptrs_[col + 1]; ++i) {
                ITYPE row = matrix_csc.row_indices_[i];
                VTYPE val = matrix_csc.values_[i];

                // Check if the row index already exists
                auto it = std::lower_bound(uniqueRows.begin(), uniqueRows.end(), row);
                ITYPE index = std::distance(uniqueRows.begin(), it);
                if (it != uniqueRows.end() && *it == row) {
                    // Row index already exists, overwrite the value
                    VTYPE sum_value = uniqueValues[index] + val;
                    VTYPE difference_value = uniqueValues[index] - val;
                    VTYPE product_value = uniqueValues[index] * val;
                    VTYPE last_value = val;
                    if ((deduplication_mode_trim == "plus")&& (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = sum_value;
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                         uniqueValues[index] = difference_value;
                    } else if ((deduplication_mode_trim == "multiply") && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = product_value;
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueValues[index] = last_value;
                    }
                } else {
                    // Column index does not exist, insert it
                    VTYPE sum_value = val;
                    VTYPE difference_value = val;
                    VTYPE product_value = val;
                    VTYPE last_value = val;
                    if ((deduplication_mode_trim == "plus")&& (std::abs(sum_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueRows.emplace_back(row);
                        uniqueValues.emplace_back(sum_value);
                    } else if ((deduplication_mode_trim == "minus") && (std::abs(difference_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueRows.emplace_back(row);
                        uniqueValues.emplace_back(difference_value);
                    } else if ((deduplication_mode_trim == "multiply") && (std::abs(product_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueRows.emplace_back(row);
                        uniqueValues.emplace_back(product_value);
                    } else if ((deduplication_mode_trim == "overwrite") && (std::abs(last_value) >= std::numeric_limits<VTYPE>::min())) {
                        uniqueRows.emplace_back(row);
                        uniqueValues.emplace_back(last_value);
                    }
                }
            }

            // Insert the unique row indices and values into the temporary vector
            sorted_row_indices.insert(sorted_row_indices.end(), uniqueRows.begin(), uniqueRows.end());
            sorted_values.insert(sorted_values.end(), uniqueValues.begin(), uniqueValues.end());

            // Update the col_ptrs_ vector
            matrix_csc.col_ptrs_[col + 1] = matrix_csc.col_ptrs_[col] + uniqueRows.size();
        }
    }

    if (deduplication) {
        // Assign the sorted and deduplicated data back to the row_indices and value vectors
        matrix_csc.row_indices_.clear();
        matrix_csc.row_indices_.shrink_to_fit();
        matrix_csc.row_indices_.swap(sorted_row_indices);
        matrix_csc.values_.clear();
        matrix_csc.values_.shrink_to_fit();
        matrix_csc.values_.swap(sorted_values);
        nnz_ = matrix_csc.values_.size();
    }

    // Clear the temporary vector
    sorted_row_indices.clear();
    sorted_values.clear();

}

// Protected member function for element operation of COO matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::coo_element_operation(ITYPE r, ITYPE c, VTYPE val, const std::string &operation_mode, const std::string &file_name_input, const std::string &function_name_input) {

    std::string operation_mode_trim = string_to_lower(trim(operation_mode));

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_manipulation.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "coo_element_operation()";
    } else {
        function_name = function_name_input;
    }

    if ((operation_mode_trim != "plus") and (operation_mode_trim != "minus") and (operation_mode_trim != "multiply") and (operation_mode_trim != "overwrite") and (operation_mode_trim != "nonsort")) {
        std::cerr << "MUI Error [" << file_name << "]: Unrecognised COO element operation mode: " << operation_mode_trim << " for " << function_name << " function" << std::endl;
        std::cerr << "    Please set the COO element operation mode as:" << std::endl;
        std::cerr << "    'plus': Sum up elemental values" << std::endl;
        std::cerr << "    'minus': Take the difference among elemental values" << std::endl;
        std::cerr << "    'multiply': Take the product among elemental values" << std::endl;
        std::cerr << "    'overwrite' (default): Keeps only the last elemental value" << std::endl;
        std::cerr << "    'nonsort': Append the element value without sort or deduplication" << std::endl;
        std::abort();
    }

    if (operation_mode_trim == "nonsort") {
        if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
            matrix_coo.values_.reserve(matrix_coo.values_.size()+1);
            matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size()+1);
            matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size()+1);
            matrix_coo.values_.emplace_back(val);
            matrix_coo.row_indices_.emplace_back(r);
            matrix_coo.col_indices_.emplace_back(c);
            nnz_++;
        }
        return;
    }

    bool isElementAdded = false;

    // Find the range of column indices for the given row using lower_bound and upper_bound
    auto lower = std::lower_bound(matrix_coo.row_indices_.begin(), matrix_coo.row_indices_.end(), r);
    auto upper = std::upper_bound(matrix_coo.row_indices_.begin(), matrix_coo.row_indices_.end(), r);

    ITYPE insert_position = std::distance(matrix_coo.row_indices_.begin(), lower);

    // Iterate over the range of column indices for the given row
    for (auto it = lower; it != upper; ++it) {
        insert_position = std::distance(matrix_coo.row_indices_.begin(), it);
        if (matrix_coo.col_indices_[insert_position] < c) {
            // Do nothing
        } else if (matrix_coo.col_indices_[insert_position] == c) {
            // Found an existing entry with the same row and column, update the value
            if (operation_mode_trim == "plus") {
                // Check if the sum is zero, erase the element if so, overwrite the element if not
                if (std::abs(matrix_coo.values_[insert_position] + val) >= std::numeric_limits<VTYPE>::min()) {
                    matrix_coo.values_[insert_position] += val;
                } else {
                    matrix_coo.row_indices_.erase(matrix_coo.row_indices_.begin() + insert_position);
                    matrix_coo.col_indices_.erase(matrix_coo.col_indices_.begin() + insert_position);
                    matrix_coo.values_.erase(matrix_coo.values_.begin() + insert_position);
                    nnz_--;
                }
                isElementAdded = true;
                break;
            } else if (operation_mode_trim == "minus") {
                // Check if the difference is zero, erase the element if so, overwrite the element if not
                if (std::abs(matrix_coo.values_[insert_position] - val) >= std::numeric_limits<VTYPE>::min()) {
                    matrix_coo.values_[insert_position] -= val;
                } else {
                    matrix_coo.row_indices_.erase(matrix_coo.row_indices_.begin() + insert_position);
                    matrix_coo.col_indices_.erase(matrix_coo.col_indices_.begin() + insert_position);
                    matrix_coo.values_.erase(matrix_coo.values_.begin() + insert_position);
                    nnz_--;
                }
                isElementAdded = true;
                break;
            } else if (operation_mode_trim == "multiply") {
                // Check if the product is zero, erase the element if so, overwrite the element if not
                if (std::abs(matrix_coo.values_[insert_position] * val) >= std::numeric_limits<VTYPE>::min()) {
                    matrix_coo.values_[insert_position] *= val;
                } else {
                    matrix_coo.row_indices_.erase(matrix_coo.row_indices_.begin() + insert_position);
                    matrix_coo.col_indices_.erase(matrix_coo.col_indices_.begin() + insert_position);
                    matrix_coo.values_.erase(matrix_coo.values_.begin() + insert_position);
                    nnz_--;
                }
                isElementAdded = true;
                break;
            } else { // overwrite
                // Check if the value is zero, erase the element if so, overwrite the element if not
                if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                    matrix_coo.values_[insert_position] = val;
                } else {
                    matrix_coo.row_indices_.erase(matrix_coo.row_indices_.begin() + insert_position);
                    matrix_coo.col_indices_.erase(matrix_coo.col_indices_.begin() + insert_position);
                    matrix_coo.values_.erase(matrix_coo.values_.begin() + insert_position);
                    nnz_--;
                }
                isElementAdded = true;
                break;
            }
        } else {
            if (operation_mode_trim != "multiply") {
                if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                    // Insert a new entry at the correct sorted position
                    matrix_coo.values_.reserve(matrix_coo.values_.size()+1);
                    matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size()+1);
                    matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size()+1);

                    matrix_coo.row_indices_.insert(matrix_coo.row_indices_.begin() + insert_position, r);
                    matrix_coo.col_indices_.insert(matrix_coo.col_indices_.begin() + insert_position, c);
                    if (operation_mode_trim == "minus") {
                        matrix_coo.values_.insert(matrix_coo.values_.begin() + insert_position, -val);
                    } else {
                        matrix_coo.values_.insert(matrix_coo.values_.begin() + insert_position, val);
                    }
                    nnz_++;
                }
            }
            isElementAdded = true;
            break;
        }
    }

    if (operation_mode_trim != "multiply") {
        if ((!isElementAdded) && (std::abs(val) >= std::numeric_limits<VTYPE>::min())) {
            insert_position = std::distance(matrix_coo.row_indices_.begin(), upper);
            // Insert a new entry at the correct sorted position
            matrix_coo.values_.reserve(matrix_coo.values_.size()+1);
            matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size()+1);
            matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size()+1);

            matrix_coo.row_indices_.insert(matrix_coo.row_indices_.begin() + insert_position, r);
            matrix_coo.col_indices_.insert(matrix_coo.col_indices_.begin() + insert_position, c);
            if (operation_mode_trim == "minus") {
                matrix_coo.values_.insert(matrix_coo.values_.begin() + insert_position, -val);
            } else {
                matrix_coo.values_.insert(matrix_coo.values_.begin() + insert_position, val);
            }
            nnz_++;
        }
    }
}

// Protected member function for element operation of CSR matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csr_element_operation(ITYPE r, ITYPE c, VTYPE val, const std::string &operation_mode, const std::string &file_name_input, const std::string &function_name_input) {

    std::string operation_mode_trim = string_to_lower(trim(operation_mode));

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_manipulation.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "csr_element_operation()";
    } else {
        function_name = function_name_input;
    }

    if ((operation_mode_trim != "plus") and (operation_mode_trim != "minus") and (operation_mode_trim != "multiply") and (operation_mode_trim != "overwrite")) {
        std::cerr << "MUI Error [" << file_name << "]: Unrecognised CSR element operation mode: " << operation_mode_trim << " for " << function_name << " function" << std::endl;
        std::cerr << "    Please set the CSR element operation mode as:" << std::endl;
        std::cerr << "    'plus': Sum up elemental values" << std::endl;
        std::cerr << "    'minus': Take the difference among elemental values" << std::endl;
        std::cerr << "    'multiply': Take the product among elemental values" << std::endl;
        std::cerr << "    'overwrite' (default): Keeps only the last elemental value" << std::endl;
        std::abort();
    }

    // Find the range of indices for the given row
    ITYPE start = matrix_csr.row_ptrs_[r];
    ITYPE end = matrix_csr.row_ptrs_[r + 1];

    // Find the column index position within the range
    auto it = std::lower_bound(matrix_csr.col_indices_.begin()+start, matrix_csr.col_indices_.begin()+end, c);

    // Check if the column index already exists in the row
    if (it != matrix_csr.col_indices_.begin()+end && *it == c) {
        // Get the index of the found column index
        ITYPE insert_position = std::distance(matrix_csr.col_indices_.begin(), it);
        if (operation_mode_trim == "plus") {
            // Check if the sum value is zero
            if (std::abs(matrix_csr.values_[insert_position] + val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csr.values_[insert_position] += val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csr.col_indices_.erase(matrix_csr.col_indices_.begin()+insert_position);
                matrix_csr.values_.erase(matrix_csr.values_.begin()+insert_position);

                // Adjust the row pointers after the erased element
                for (ITYPE i = r + 1; i < static_cast<ITYPE>(matrix_csr.row_ptrs_.size()); ++i) {
                    matrix_csr.row_ptrs_[i]--;
                }
                nnz_--;
            }
        } else if (operation_mode_trim == "minus") {
            // Check if the difference value is zero
            if (std::abs(matrix_csr.values_[insert_position] - val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csr.values_[insert_position] -= val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csr.col_indices_.erase(matrix_csr.col_indices_.begin()+insert_position);
                matrix_csr.values_.erase(matrix_csr.values_.begin()+insert_position);

                // Adjust the row pointers after the erased element
                for (ITYPE i = r + 1; i < static_cast<ITYPE>(matrix_csr.row_ptrs_.size()); ++i) {
                    matrix_csr.row_ptrs_[i]--;
                }
                nnz_--;
            }
        } else if (operation_mode_trim == "multiply") {
            // Check if the product value is zero
            if (std::abs(matrix_csr.values_[insert_position] * val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csr.values_[insert_position] *= val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csr.col_indices_.erase(matrix_csr.col_indices_.begin()+insert_position);
                matrix_csr.values_.erase(matrix_csr.values_.begin()+insert_position);

                // Adjust the row pointers after the erased element
                for (ITYPE i = r + 1; i < static_cast<ITYPE>(matrix_csr.row_ptrs_.size()); ++i) {
                    matrix_csr.row_ptrs_[i]--;
                }
                nnz_--;
            }
        } else { // overwrite
            // Check if the new value is zero
            if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csr.values_[insert_position] = val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csr.col_indices_.erase(matrix_csr.col_indices_.begin()+insert_position);
                matrix_csr.values_.erase(matrix_csr.values_.begin()+insert_position);

                // Adjust the row pointers after the erased element
                for (ITYPE i = r + 1; i < static_cast<ITYPE>(matrix_csr.row_ptrs_.size()); ++i) {
                    matrix_csr.row_ptrs_[i]--;
                }
                nnz_--;
            }
        }
    } else {
        if (operation_mode_trim != "multiply") {
            // Check if the new value is zero
            if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                // Insert the new value and column index at the determined position
                ITYPE insert_position = std::distance(matrix_csr.col_indices_.begin(), it);
                matrix_csr.col_indices_.reserve(matrix_csr.col_indices_.size()+1);
                matrix_csr.values_.reserve(matrix_csr.values_.size()+1);
                matrix_csr.col_indices_.insert(matrix_csr.col_indices_.begin()+insert_position, c);
                if (operation_mode_trim == "minus") {
                    matrix_csr.values_.insert(matrix_csr.values_.begin()+insert_position, -val);
                } else {
                    matrix_csr.values_.insert(matrix_csr.values_.begin()+insert_position, val);
                }

                // Adjust the row pointers after the inserted element
                for (ITYPE i = r + 1; i < static_cast<ITYPE>(matrix_csr.row_ptrs_.size()); ++i) {
                    matrix_csr.row_ptrs_[i]++;
                }
                nnz_++;
            } else {
                // No existing entry to update, do nothing
            }
        }
    }
}

// Protected member function for element operation of CSC matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csc_element_operation(ITYPE r, ITYPE c, VTYPE val, const std::string &operation_mode, const std::string &file_name_input, const std::string &function_name_input) {

    std::string operation_mode_trim = string_to_lower(trim(operation_mode));

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_manipulation.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "csc_element_operation()";
    } else {
        function_name = function_name_input;
    }

    if ((operation_mode_trim != "plus") and (operation_mode_trim != "minus") and (operation_mode_trim != "multiply") and (operation_mode_trim != "overwrite")) {
        std::cerr << "MUI Error [" << file_name << "]: Unrecognised CSC element operation mode: " << operation_mode_trim << " for " << function_name << " function" << std::endl;
        std::cerr << "    Please set the CSC element operation mode as:" << std::endl;
        std::cerr << "    'plus': Sum up elemental values" << std::endl;
        std::cerr << "    'minus': Take the difference among elemental values" << std::endl;
        std::cerr << "    'multiply': Take the product among elemental values" << std::endl;
        std::cerr << "    'overwrite' (default): Keeps only the last elemental value" << std::endl;
        std::abort();
    }


    // Find the range of column indices for the given column
    ITYPE start = matrix_csc.col_ptrs_[c];
    ITYPE end = matrix_csc.col_ptrs_[c + 1];

    // Find the position of the row index within the range
    auto it = std::lower_bound(matrix_csc.row_indices_.begin()+start, matrix_csc.row_indices_.begin()+end, r);

    // Check if the row index already exists
    if (it != matrix_csc.row_indices_.begin()+end && *it == r) {
        // Get the index of the found column index
        ITYPE insert_position = std::distance(matrix_csc.row_indices_.begin(), it);
        if (operation_mode_trim == "plus") {
            // Check if the sum value is zero
            if (std::abs(matrix_csc.values_[insert_position] + val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csc.values_[insert_position] += val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csc.row_indices_.erase(matrix_csc.row_indices_.begin()+insert_position);
                matrix_csc.values_.erase(matrix_csc.values_.begin()+insert_position);

                // Adjust the column pointers after the erased element
                for (ITYPE i = c + 1; i < static_cast<ITYPE>(matrix_csc.col_ptrs_.size()); ++i) {
                    matrix_csc.col_ptrs_[i]--;
                }
                nnz_--;
            }
        } else if (operation_mode_trim == "minus") {
            // Check if the difference value is zero
            if (std::abs(matrix_csc.values_[insert_position] - val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csc.values_[insert_position] -= val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csc.row_indices_.erase(matrix_csc.row_indices_.begin()+insert_position);
                matrix_csc.values_.erase(matrix_csc.values_.begin()+insert_position);

                // Adjust the column pointers after the erased element
                for (ITYPE i = c + 1; i < static_cast<ITYPE>(matrix_csc.col_ptrs_.size()); ++i) {
                    matrix_csc.col_ptrs_[i]--;
                }
                nnz_--;
            }
        } else if (operation_mode_trim == "multiply") {
            // Check if the product value is zero
            if (std::abs(matrix_csc.values_[insert_position] * val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csc.values_[insert_position] *= val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csc.row_indices_.erase(matrix_csc.row_indices_.begin()+insert_position);
                matrix_csc.values_.erase(matrix_csc.values_.begin()+insert_position);

                // Adjust the column pointers after the erased element
                for (ITYPE i = c + 1; i < static_cast<ITYPE>(matrix_csc.col_ptrs_.size()); ++i) {
                    matrix_csc.col_ptrs_[i]--;
                }
                nnz_--;
            }
        } else { // overwrite
            // Check if the new value is zero
            if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                // Update the existing value with the new value
                matrix_csc.values_[insert_position] = val;
            } else {
                // Erase the existing entry from the vectors
                matrix_csc.row_indices_.erase(matrix_csc.row_indices_.begin()+insert_position);
                matrix_csc.values_.erase(matrix_csc.values_.begin()+insert_position);

                // Adjust the column pointers after the erased element
                for (ITYPE i = c + 1; i < static_cast<ITYPE>(matrix_csc.col_ptrs_.size()); ++i) {
                    matrix_csc.col_ptrs_[i]--;
                }
                nnz_--;
            }
        }
    } else {
        if (operation_mode_trim != "multiply") {
            // Check if the new value is zero
            if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                // Insert the new value and row index at the determined position
                ITYPE insert_position = std::distance(matrix_csc.row_indices_.begin(), it);
                matrix_csc.row_indices_.reserve(matrix_csc.row_indices_.size()+1);
                matrix_csc.values_.reserve(matrix_csc.values_.size()+1);
                matrix_csc.row_indices_.insert(matrix_csc.row_indices_.begin()+insert_position, r);
                if (operation_mode_trim == "minus") {
                    matrix_csc.values_.insert(matrix_csc.values_.begin()+insert_position, -val);
                } else {
                    matrix_csc.values_.insert(matrix_csc.values_.begin()+insert_position, val);
                }
                // Adjust the row pointers after the inserted element
                for (ITYPE i = c + 1; i < static_cast<ITYPE>(matrix_csc.col_ptrs_.size()); ++i) {
                    matrix_csc.col_ptrs_[i]++;
                }
                nnz_++;
            } else {
                // No existing entry to update, do nothing
            }
        }
    }

}

// Protected member function to convert COO matrix into CSR matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::coo_to_csr() {

    assert((matrix_format_ == format::COO) &&
      "MUI Error [matrix_manipulation.h]: coo_to_csr() can only applied to sparse matrix with COO format.");

    if (matrix_coo.values_.size() == 0) {
        matrix_format_ = format::CSR;
        matrix_csr.row_ptrs_.resize((rows_+1), 0);
        return;
    }

    // Determine the number of non-zero entries in each row
    std::vector<ITYPE> rowPtr(rows_+1, 0);

    for (ITYPE i = 0; i < static_cast<ITYPE>(matrix_coo.row_indices_.size()); ++i) {
        if (matrix_coo.row_indices_[i] >= (rows_+1)) {
            std::cerr << "MUI Error [matrix_manipulation.h]: row index: " << matrix_coo.row_indices_[i] << " at the " << i << "th matrix_coo.row_indices_ is out of range (" << (rows_+1) << ") in coo_to_csr()" << std::endl;
            std::abort();
        } else {
            rowPtr[matrix_coo.row_indices_[i]]++;
        }

    }

    for(ITYPE i = 0, accumulator = 0; i < rows_; ++i){
        ITYPE temp = rowPtr[i];
        rowPtr[i] = accumulator;
        accumulator += temp;
    }
    rowPtr[rows_] = nnz_;

    // Fill the CSRData struct
    matrix_csr.row_ptrs_.reserve(rows_+1);
    matrix_csr.values_.reserve(nnz_);
    matrix_csr.col_indices_.reserve(nnz_);

    std::swap(matrix_csr.values_, matrix_coo.values_);
    std::swap(matrix_csr.row_ptrs_, rowPtr);
    std::swap(matrix_csr.col_indices_, matrix_coo.col_indices_);

    // Deallocate the memory
    matrix_coo.values_.clear();
    matrix_coo.row_indices_.clear();
    matrix_coo.col_indices_.clear();
    rowPtr.clear();

    // Reset the matrix format
    matrix_format_ = format::CSR;
}

// Protected member function to convert COO matrix into CSC matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::coo_to_csc() {

    assert((matrix_format_ == format::COO) &&
      "MUI Error [matrix_manipulation.h]: coo_to_csc() can only applied to sparse matrix with COO format.");

    if (matrix_coo.values_.size() == 0) {
        matrix_format_ = format::CSC;
        matrix_csc.col_ptrs_.resize((cols_+1), 0);
        return;
    }

    // Determine the number of non-zero entries in each column
    std::vector<ITYPE> colPtr(cols_+1, 0);

    for (ITYPE i = 0; i < static_cast<ITYPE>(matrix_coo.col_indices_.size()); ++i) {
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
    matrix_format_ = format::CSC;
}

// Protected member function to convert CSR matrix into COO matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csr_to_coo() {

    assert((matrix_format_ == format::CSR) &&
      "MUI Error [matrix_manipulation.h]: csr_to_coo() can only applied to sparse matrix with CSR format.");

    if (matrix_csr.values_.size() == 0) {
        matrix_format_ = format::COO;
        return;
    }

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
    matrix_format_ = format::COO;
}

// Protected member function to convert CSR matrix into CSC matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csr_to_csc() {

    assert((matrix_format_ == format::CSR) &&
      "MUI Error [matrix_manipulation.h]: csr_to_csc() can only applied to sparse matrix with CSR format.");

    if (matrix_csr.values_.size() == 0) {
        matrix_format_ = format::CSC;
        matrix_csc.col_ptrs_.resize((cols_+1), 0);
        return;
    }

    matrix_csc.col_ptrs_.resize((cols_+1), 0);
    matrix_csc.values_.resize(nnz_, 0);
    matrix_csc.row_indices_.resize(nnz_, 0);

    // Determine the number of non-zero entries in each column
    for (ITYPE i = 0; i < cols_; ++i) {
        for (ITYPE j = 0; j < static_cast<ITYPE>(matrix_csr.col_indices_.size()); ++j) {
            if (matrix_csr.col_indices_[j] == i) {
                ++matrix_csc.col_ptrs_[i];
            }
        }
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
    matrix_format_ = format::CSC;

}

// Protected member function to convert CSC matrix into COO matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csc_to_coo() {

    assert((matrix_format_ == format::CSC) &&
      "MUI Error [matrix_manipulation.h]: csc_to_coo() can only applied to sparse matrix with CSC format.");

    if (matrix_csc.values_.size() == 0) {
        matrix_format_ = format::COO;
        return;
    }

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
    matrix_format_ = format::COO;
}

// Protected member function to convert CSC matrix into CSR matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::csc_to_csr() {

    assert((matrix_format_ == format::CSC) &&
      "MUI Error [matrix_manipulation.h]: csc_to_csr() can only applied to sparse matrix with CSC format.");

    if (matrix_csc.values_.size() == 0) {
        matrix_format_ = format::CSR;
        matrix_csr.row_ptrs_.resize((rows_+1), 0);
        return;
    }

    matrix_csr.values_.resize(nnz_, 0);
    matrix_csr.row_ptrs_.resize((rows_+1), 0);
    matrix_csr.col_indices_.resize(nnz_, 0);

    // Determine the number of non-zero entries in each row
    for (ITYPE i = 0; i < static_cast<ITYPE>(matrix_csc.row_indices_.size()); ++i) {
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
    matrix_format_ = format::CSR;

}

// Protected member function to clear all vectors of the sparse matrix
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::clear_vectors() {
    // Clear the existing elements
    matrix_coo.values_.clear();
    matrix_coo.row_indices_.clear();
    matrix_coo.col_indices_.clear();
    matrix_csr.values_.clear();
    matrix_csr.row_ptrs_.clear();
    matrix_csr.col_indices_.clear();
    matrix_csc.values_.clear();
    matrix_csc.row_indices_.clear();
    matrix_csc.col_ptrs_.clear();
    nnz_ = 0;
}

} // linalg
} // mui

#endif /* MUI_MATRIX_MANIPULATION_H_ */
