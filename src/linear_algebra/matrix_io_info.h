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
 * @file matrix_io_info.h
 * @author W. Liu
 * @date 01 February 2023
 * @brief Implementation of sparse matrix I/O functions.
 */

#ifndef MUI_MATRIX_IO_INFO_H_
#define MUI_MATRIX_IO_INFO_H_

#include <sstream>
#include <cassert>
#include <limits>
#include <algorithm>
#include <cctype>

namespace mui {
namespace linalg {

// **************************************************
// ************ Public member functions *************
// **************************************************

// Member function to print matrix elements to the console
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::print() const {
    for (ITYPE i = 0; i < rows_; ++i) {
        std::cout << "      ";
       for (ITYPE j = 0; j < cols_; ++j){
           auto it = matrix_.find(std::make_pair(i, j));
           if (it != matrix_.end()) {
               std::cout << it->second << " ";
           } else {
               std::cout << 0 << " ";
           }
       }
       std::cout << std::endl;
    }
}

// Overloading << operator to output matrix in CSV format
template<typename ITYPE, typename VTYPE>
std::ostream& operator << (std::ostream &ofile, const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    for (ITYPE i = 0; i < exist_mat.get_rows(); ++i) {
        for (ITYPE j = 0; j < exist_mat.get_cols(); ++j) {
            if (j == (exist_mat.get_cols() - 1)) {
                ofile << exist_mat.get_value(i,j);
            } else {
                ofile << exist_mat.get_value(i,j) << ",";
            }
        }
        if (i != (exist_mat.get_rows() - 1)) {
            ofile << std::endl;
        }
    }
    return ofile;
}

// Overloading >> operator to read matrix from a file in CSV format with lines start with "//" as comment lines
template<typename ITYPE, typename VTYPE>
std::istream& operator >> (std::istream &ifile, sparse_matrix<ITYPE,VTYPE> &exist_mat) {
    assert((exist_mat.empty()) &&
      "MUI Error [matrix_io_info.h]: Overloading >> operator can only takes in null matrix or empty (all-zero) matrix");
    std::string rawLine;
    std::map<std::pair<ITYPE, ITYPE>, VTYPE> tempMatrix;
    ITYPE row = 0;
    ITYPE col = 0;
    while (std::getline(ifile, rawLine)) {
        std::string line = trim(rawLine);
        // Skips the line if the first two characters are '//'
        if (( line[0] == '/' && line[1] == '/' ) || (line.empty())) continue;
        std::stringstream ss(line);
        std::string value;
        ITYPE colCount = 0;
        while (std::getline(ss, value, ',')) {
            VTYPE val = static_cast<VTYPE>(std::stod(value));
            if (std::abs(val) >= std::numeric_limits<VTYPE>::min()) {
                tempMatrix[std::make_pair(row, colCount)] = val;
            }
            ++colCount;
        }
        if (col == 0){
            col = colCount;
        } else {
            if (col != colCount) {
                std::cout << "MUI Warning [matrix_io_info.h]: The number of columns of the matrix read in at row " <<
                        row << " is " << colCount << ", which is different from previous row (i.e. " <<
                        col << " columns!" << std::endl;
                col = colCount;
            }
        }
        ++row;
    }
    if ((exist_mat.get_rows() == 0) && (exist_mat.get_cols() == 0)) {
        exist_mat.resize_null(row, col);
        for (auto element : tempMatrix) {
            exist_mat.set_value(element.first.first, element.first.second, element.second);
        }
    } else {
        assert(((exist_mat.get_rows() == row) && (exist_mat.get_cols() == col)) &&
          "MUI Error [matrix_io_info.h]: Matrix size mismatching between existing matrix and read in matrix in overloading >> operator ");
        for (auto element : tempMatrix) {
            exist_mat.set_value(element.first.first, element.first.second, element.second);
        }
    }
    return ifile;
}

// Member function to get the value at a given position
template<typename ITYPE, typename VTYPE>
VTYPE sparse_matrix<ITYPE,VTYPE>::get_value(ITYPE r, ITYPE c) const {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_io_info.h]: Matrix index out of range in get_value function");
    auto it = matrix_.find(std::make_pair(r, c));
    if (it != matrix_.end()) {
        return it->second;
    } else {
        return 0;
    }
}

// Member function to get the number of rows
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::get_rows() const {
    return rows_;
}

// Member function to get the number of cols
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::get_cols() const {
    return cols_;
}

// Member function to get non-zero elements
template<typename ITYPE, typename VTYPE>
std::vector<std::pair<ITYPE, ITYPE>> sparse_matrix<ITYPE,VTYPE>::get_non_zero_elements() const {
    std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
    for (auto const &nn_element : matrix_) {
        vec_temp.push_back(std::make_pair(nn_element.first.first, nn_element.first.second));
    }
    return vec_temp;
}

// Member function to get number of non-zero elements
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::non_zero_elements_count() const {
    return matrix_.size();
}

// Member function to check whether the matrix contains all zero elements
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::empty() const {
    return matrix_.empty();
}

// Member function to get the format of the matrix
template<typename ITYPE, typename VTYPE>
std::string sparse_matrix<ITYPE,VTYPE>::getFormat() const
{
	std::string matrix_format;

	if (matrix_format_ == format::COO) {
		matrix_format = "COO";
	} else if (matrix_format_ == format::CSR) {
		matrix_format = "CSR";
	} else if (matrix_format_ == format::CSC) {
		matrix_format = "CSC";
	} else {
        std::cerr << "MUI Error [matrix_io_info.h]: unknown matrix format " << matrix_format_ << std::endl;
        std::abort();
	}

	return matrix_format;
}


// Member function to check if the sparse matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::isSortedUnique(const std::string &file_name, const std::string &function_name) const {

    if (file_name.empty())
        file_name = "matrix_io_info.h";
    if (function_name.empty())
        function_name = "isSortedUnique()";

	if (matrix_format_ == format::COO) {
		return this->isCOOSortedUnique(file_name, function_name);
	} else if (matrix_format_ == format::CSR) {
		return this->isCSRSortedUnique(file_name, function_name);
	} else if (matrix_format_ == format::CSC) {
		return this->isCSCSortedUnique(file_name, function_name);
	} else {
        std::cerr << "MUI Error [matrix_io_info.h]: unknown matrix format " << matrix_format_ << std::endl;
        std::abort();
	}
}

// **************************************************
// ********** Protected member functions ************
// **************************************************

// Protected member function to check if the COO matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::isCOOSortedUnique(const std::string &file_name, const std::string &function_name) const {

    if (file_name.empty())
        file_name = "matrix_io_info.h";
    if (function_name.empty())
        function_name = "isCOOSortedUnique()";

	ITYPE numEntries = matrix_coo.values_.size();

    if (numEntries > 1) {
        for (ITYPE i = 1; i < numEntries; ++i) {
            // Compare the current entry with the previous one
            if (matrix_coo.row_indices_[i] < matrix_coo.row_indices_[i - 1]) {
                // Row index is not sorted
                std::cout << "MUI [" << file_name << "]: The COO type matrix is not sorted (sorted row index check failed) in " << function_name <<  std::endl;
                return false;
            } else if (matrix_coo.row_indices_[i] == matrix_coo.row_indices_[i - 1]) {
                // Row index is the same, check column index
                if (matrix_coo.col_indices_[i] < matrix_coo.col_indices_[i - 1]) {
                    // Column index is not sorted
                    std::cout << "MUI [" << file_name << "]: The COO type matrix is not sorted (sorted column index check failed) in " << function_name <<  std::endl;
                    return false;
                } else if (matrix_coo.col_indices_[i] == matrix_coo.col_indices_[i - 1]) {
                    // Column index has duplicate elements
                    std::cout << "MUI [" << file_name << "]: The COO type matrix exists duplicated elements (unique column index check failed) in " << function_name <<  std::endl;
                    return false;
                }
            }
        }
    }

	return true;

}

// Protected member function to check if the CSR matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::isCSRSortedUnique(const std::string &file_name, const std::string &function_name) const {

    if (file_name.empty())
        file_name = "matrix_io_info.h";
    if (function_name.empty())
        function_name = "isCSRSortedUnique()";

	ITYPE numEntries = matrix_csr.values_.size();

    if (numEntries > 1) {
        for(ITYPE i = 0; i < rows_; ++i){
            if (matrix_csr.row_ptrs_[i] > matrix_csr.row_ptrs_[i+1]) {
                // Row pointers is not sorted
                std::cout << "MUI [" << file_name << "]: The CSR type matrix is not sorted (sorted row pointers check failed) in " << function_name <<  std::endl;
                return false;
            }
            for(ITYPE j = matrix_csr.row_ptrs_[i] + 1; j < matrix_csr.row_ptrs_[i+1]; ++j){
                if(matrix_csr.col_indices_[j-1] > matrix_csr.col_indices_[j]){
                    // Column indices is not sorted
                    std::cout << "MUI [" << file_name << "]: The CSR type matrix is not sorted (sorted column index check failed) in " << function_name <<  std::endl;
                    return false;
                } else if (matrix_csr.col_indices_[j-1] == matrix_csr.col_indices_[j]) {
                    // Column indices is not unique
                    std::cout << "MUI [" << file_name << "]: The CSR type matrix is not unique (deduplicated column index check failed) in " << function_name <<  std::endl;
                    return false;
                }
            }
        }
    }

	return true;

}

// Protected member function to check if the CSC matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::isCSCSortedUnique(const std::string &file_name, const std::string &function_name) const {

    if (file_name.empty())
        file_name = "matrix_io_info.h";
    if (function_name.empty())
        function_name = "isCSCSortedUnique()";

	ITYPE numEntries = matrix_csc.values_.size();

    if (numEntries > 1) {
        for(ITYPE i = 0; i < cols_; ++i){
            if (matrix_csc.col_ptrs_[i] > matrix_csc.col_ptrs_[i+1]) {
                // Column pointers is not sorted
                std::cout << "MUI [" << file_name << "]: The CSC type matrix is not sorted (sorted column pointers check failed) in " << function_name <<  std::endl;
                return false;
            }
            for(ITYPE j = matrix_csc.col_ptrs_[i] + 1; j < matrix_csc.col_ptrs_[i+1]; ++j){
                if(matrix_csc.row_indices_[j-1] > matrix_csc.row_indices_[j]){
                    // Row indices is not sorted
                    std::cout << "MUI [" << file_name << "]: The CSC type matrix is not sorted (sorted row index check failed) in " << function_name <<  std::endl;
                    return false;
                } else if (matrix_csc.row_indices_[j-1] == matrix_csc.row_indices_[j]) {
                    // Row indices is not unique
                    std::cout << "MUI [" << file_name << "]: The CSC type matrix is not unique (deduplicated row index check failed) in " << function_name <<  std::endl;
                    return false;
                }
            }
        }
    }

	return true;

}


} // linalg
} // mui

#endif /* MUI_MATRIX_IO_INFO_H_ */
