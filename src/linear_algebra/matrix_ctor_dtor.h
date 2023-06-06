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
 * @file matrix_ctor_dtor.h
 * @author W. Liu
 * @date 01 February 2023
 * @brief Implementation of sparse matrix constructors and destructor.
 */

#ifndef MUI_MATRIX_CTOR_DTOR_H_
#define MUI_MATRIX_CTOR_DTOR_H_

namespace mui {
namespace linalg {

// **************************************************
// ************ Public member functions *************
// **************************************************

// Constructor - takes in size of row and column to generate an matrix, with default arguments of format vectors.
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(ITYPE r, ITYPE c, const std::string &format, const std::vector<VTYPE> &value_vector, const std::vector<ITYPE> &row_vector, const std::vector<ITYPE> &column_vector)
    : rows_(r), cols_(c) {

    this->sparse_matrix<ITYPE,VTYPE>::set_matrix_format(format);


    if (!value_vector.empty()) {

        nnz_ = value_vector.size();

        if (matrix_format_ == format::COO) {

            matrix_coo.values_.reserve(value_vector.size());
            matrix_coo.row_indices_.reserve(row_vector.size());
            matrix_coo.col_indices_.reserve(column_vector.size());

            matrix_coo.values_ = std::vector<VTYPE>(value_vector.begin(), value_vector.end());
            matrix_coo.row_indices_ = std::vector<ITYPE>(row_vector.begin(), row_vector.end());
            matrix_coo.col_indices_ = std::vector<ITYPE>(column_vector.begin(), column_vector.end());

        } else if (matrix_format_ == format::CSR) {

            matrix_csr.values_.reserve(value_vector.size());
            matrix_csr.row_ptrs_.reserve(row_vector.size());
            matrix_csr.col_indices_.reserve(column_vector.size());

            matrix_csr.values_ = std::vector<VTYPE>(value_vector.begin(), value_vector.end());
            matrix_csr.row_ptrs_ = std::vector<ITYPE>(row_vector.begin(), row_vector.end());
            matrix_csr.col_indices_ = std::vector<ITYPE>(column_vector.begin(), column_vector.end());

        } else if (matrix_format_ == format::CSC) {

            matrix_csc.values_.reserve(value_vector.size());
            matrix_csc.row_indices_.reserve(row_vector.size());
            matrix_csc.col_ptrs_.reserve(column_vector.size());

            matrix_csc.values_ = std::vector<VTYPE>(value_vector.begin(), value_vector.end());
            matrix_csc.row_indices_ = std::vector<ITYPE>(row_vector.begin(), row_vector.end());
            matrix_csc.col_ptrs_ = std::vector<ITYPE>(column_vector.begin(), column_vector.end());

        } else {
              std::cerr << "MUI Error [matrix_ctor_dtor.h]: Unrecognised matrix format for matrix constructor" << std::endl;
              std::cerr << "    Please set the matrix_format_ as:" << std::endl;
              std::cerr << "    format::COO: COOrdinate format" << std::endl;
              std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
              std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
              std::abort();
          }

        this->assert_valid_vector_size("matrix_ctor_dtor.h", "sparse_matrix constructor function");
    } else {

        if (matrix_format_ == format::CSR) {
            matrix_csr.row_ptrs_.resize((rows_+1), 0);
        } else if (matrix_format_ == format::CSC) {
            matrix_csc.col_ptrs_.resize((cols_+1), 0);
        }
    }

}

// Constructor - null matrix
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(const std::string &format)
    : rows_(0), cols_(0) {

    this->sparse_matrix<ITYPE,VTYPE>::set_matrix_format(format);

    if (matrix_format_ == format::CSR) {
        matrix_csr.row_ptrs_.resize((rows_+1), 0);
    } else if (matrix_format_ == format::CSC) {
        matrix_csc.col_ptrs_.resize((cols_+1), 0);
    }

}

// Constructor - takes in another sparse_matrix object as an argument
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(const sparse_matrix<ITYPE,VTYPE> &exist_mat)
    : rows_(exist_mat.rows_),
      cols_(exist_mat.cols_),
      nnz_(exist_mat.nnz_),
      matrix_format_(exist_mat.matrix_format_) {

      // Copy the data from the existing matrix
      if (matrix_format_ == format::COO) {
          matrix_coo.values_ = std::move(exist_mat.matrix_coo.values_);
          matrix_coo.row_indices_ = std::move(exist_mat.matrix_coo.row_indices_);
          matrix_coo.col_indices_ = std::move(exist_mat.matrix_coo.col_indices_);
      } else if (matrix_format_ == format::CSR) {
          matrix_csr.values_ = std::move(exist_mat.matrix_csr.values_);
          matrix_csr.row_ptrs_ = std::move(exist_mat.matrix_csr.row_ptrs_);
          matrix_csr.col_indices_ = std::move(exist_mat.matrix_csr.col_indices_);
      } else if (matrix_format_ == format::CSC) {
          matrix_csc.values_ = std::move(exist_mat.matrix_csc.values_);
          matrix_csc.row_indices_ = std::move(exist_mat.matrix_csc.row_indices_);
          matrix_csc.col_ptrs_ = std::move(exist_mat.matrix_csc.col_ptrs_);
      } else {
          std::cerr << "MUI Error [matrix_ctor_dtor.h]: Unrecognised matrix format for matrix constructor" << std::endl;
          std::cerr << "    Please set the matrix_format_ as:" << std::endl;
          std::cerr << "    format::COO: COOrdinate format" << std::endl;
          std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
          std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
          std::abort();
      }

}

// Constructor - takes in a std::vector with row major dense matrix format as an argument
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(const std::vector<std::vector<VTYPE>>& denseVector, const std::string &format) {

    this->sparse_matrix<ITYPE,VTYPE>::set_matrix_format(format);

    // Set the dimensions of the sparse matrix based on the dimensions of the dense matrix vector
    rows_ = denseVector.size();
    cols_ = denseVector[0].size();

    // Calculate the total number of non-zero elements
    for (const auto& row : denseVector) {
        for (ITYPE i = 0; i < row.size(); ++i) {
            if (std::abs(row[i]) >= std::numeric_limits<VTYPE>::min()) {
                nnz_++;
            }
        }
    }

    // Reserve space for the sparse entries
    matrix_coo.values_.reserve(nnz_);
    matrix_coo.row_indices_.reserve(nnz_);
    matrix_coo.col_indices_.reserve(nnz_);

    // Loop over the dense matrix vector and add non-zero entries to COO data
    for (ITYPE i = 0; i < rows_; ++i) {
        for (ITYPE j = 0; j < cols_; ++j) {
            const VTYPE val = denseVector[i][j];
            if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                matrix_coo.values_.emplace_back(val);
                matrix_coo.row_indices_.emplace_back(i);
                matrix_coo.col_indices_.emplace_back(j);
            }
        }
    }

    if (matrix_format_ == format::COO) {

        // Do nothing

    } else if (matrix_format_ == format::CSR) {

        matrix_format_ = format::COO;
        this->sparse_matrix<ITYPE,VTYPE>::coo_to_csr();

    } else if (matrix_format_ == format::CSC) {

        matrix_format_ = format::COO;
        this->sparse_matrix<ITYPE,VTYPE>::sort_coo(false);
        this->sparse_matrix<ITYPE,VTYPE>::coo_to_csc();

    } else {

        std::cerr << "MUI Error [matrix_ctor_dtor.h]: Unrecognised matrix format for matrix constructor" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();

    }

}

// Constructor - generate various square matrices
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(ITYPE n, const std::string &token, const std::string &format)
    : rows_(n), cols_(n) {

    this->set_matrix_format(format);

    if(token.empty()) {

        // empty (all-zero) square matrix (Do nothing from the code perspective)
        if (matrix_format_ == format::CSR) {
            matrix_csr.row_ptrs_.resize((rows_+1), 0);
        } else if (matrix_format_ == format::CSC) {
            matrix_csc.col_ptrs_.resize((cols_+1), 0);
        }

    } else if(string_to_lower(trim(token))=="identity") {

        // set number of non-zero elements as n
        nnz_ = n;

        // identity square matrix
        if (matrix_format_ == format::COO) {

            matrix_coo.values_.reserve(n);
            matrix_coo.row_indices_.reserve(n);
            matrix_coo.col_indices_.reserve(n);

            for (ITYPE i = 0; i < n; ++i) {
                matrix_coo.values_.emplace_back(1.0);
                matrix_coo.row_indices_.emplace_back(i);
                matrix_coo.col_indices_.emplace_back(i);
            }

        } else if (matrix_format_ == format::CSR) {

            matrix_csr.values_.reserve(n);
            matrix_csr.row_ptrs_.reserve(n+1);
            matrix_csr.col_indices_.reserve(n);

            for (ITYPE i = 0; i < n; i++) {
                matrix_csr.values_.emplace_back(1.0);
                matrix_csr.row_ptrs_.emplace_back(i);
                matrix_csr.col_indices_.emplace_back(i);
            }

            matrix_csr.row_ptrs_.emplace_back(n);

        } else if (matrix_format_ == format::CSC) {

            matrix_csc.values_.reserve(n);
            matrix_csc.row_indices_.reserve(n);
            matrix_csc.col_ptrs_.reserve(n+1);

            for (ITYPE i = 0; i < n; ++i) {
                matrix_csc.values_.emplace_back(1.0);
                matrix_csc.row_indices_.emplace_back(i);
                matrix_csc.col_ptrs_.emplace_back(i);
            }

            matrix_csc.col_ptrs_.emplace_back(n);

        } else {
            std::cerr << "MUI Error [matrix_ctor_dtor.h]: Unrecognised matrix format for matrix constructor" << std::endl;
            std::cerr << "    Please set the matrix_format_ as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }

    } else {
        std::cerr << "MUI Error [matrix_ctor_dtor.h]: unidentified token string for square matrix constructor" << std::endl;
        std::cerr << "    Please set the token string as:" << std::endl;
        std::cerr << "    empty string (default): Empty (all-zero) square matrix" << std::endl;
        std::cerr << "    'identity': identity square matrix" << std::endl;
        std::abort();
    }
}

// Destructor
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::~sparse_matrix<ITYPE,VTYPE>() {

    // deallocate the memory for non-zero elements
    matrix_coo.values_.clear();
    matrix_coo.row_indices_.clear();
    matrix_coo.col_indices_.clear();

    matrix_csr.values_.clear();
    matrix_csr.row_ptrs_.clear();
    matrix_csr.col_indices_.clear();

    matrix_csc.values_.clear();
    matrix_csc.row_indices_.clear();
    matrix_csc.col_ptrs_.clear();

    // Reset sparse matrix properties
    rows_ = 0;
    cols_ = 0;
    nnz_ = 0;
    matrix_format_ = format::CSR;
    matrix_csr.row_ptrs_.resize(1, 0);

}

// **************************************************
// ********** Protected member functions ************
// **************************************************

// Protected member function to set matrix format - helper function on matrix constructors
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::set_matrix_format(const std::string &format) {

    std::string matrix_format = string_to_upper(trim(format));

    if (matrix_format == "COO") {
        matrix_format_ = format::COO;
    } else if (matrix_format == "CSR") {
        matrix_format_ = format::CSR;
    } else if (matrix_format == "CSC") {
        matrix_format_ = format::CSC;
    } else {
        std::cerr << "MUI Error [matrix_ctor_dtor.h]: Unrecognised format type: " << format << " for matrix constructor" << std::endl;
        std::cerr << "    Please set the format string as:" << std::endl;
        std::cerr << "    'COO': COOrdinate format" << std::endl;
        std::cerr << "    'CSR' (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    'CSC': Compressed Sparse Column format" << std::endl;
        std::abort();
    }

}

} // linalg
} // mui

#endif /* MUI_MATRIX_CTOR_DTOR_H_ */
