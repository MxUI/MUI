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
#include <fstream>
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
           std::cout << this->get_value(i,j) << " ";
       }
       std::cout << std::endl;
    }
}

// Member function to print matrix vectors to the console
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::print_vectors() const {
    if (matrix_format_ == format::COO) {
        std::cout << "The Value vector of matrix in COO format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_coo.values_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        std::cout << "The Row Index vector of matrix in COO format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_coo.row_indices_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        std::cout << "The Column Index vector of matrix in COO format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_coo.col_indices_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    } else if (matrix_format_ == format::CSR) {
        std::cout << "The Value vector of matrix in CSR format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_csr.values_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        std::cout << "The Row Pointers vector of matrix in CSR format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_csr.row_ptrs_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        std::cout << "The Column Index vector of matrix in CSR format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_csr.col_indices_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    } else if (matrix_format_ == format::CSC) {
        std::cout << "The Value vector of matrix in CSC format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_csc.values_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        std::cout << "The Row Index vector of matrix in CSC format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_csc.row_indices_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        std::cout << "The Column Pointers vector of matrix in CSC format: " << std::endl;
        std::cout << "      ";
        for (const auto& element : matrix_csc.col_ptrs_) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format for matrix print_vectors >>" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }
}

// Member function to write matrix vectors to the file
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::write_vectors_to_file(const std::string &format_file_name, const std::string &val_file_name, const std::string &row_file_name, const std::string &col_file_name) const {

    std::string ffn;
    std::string vfn;
    std::string rfn;
    std::string cfn;

    // Open three separate output files
    if ((val_file_name.empty()) ||  (row_file_name.empty()) ||  (col_file_name.empty())) {
        ffn = format_file_name + "_format.dat";
        vfn = format_file_name + "_value.dat";
        rfn = format_file_name + "_row.dat";
        cfn = format_file_name + "_column.dat";
    } else {
        ffn = format_file_name;
        vfn = val_file_name;
        rfn = row_file_name;
        cfn = col_file_name;
    }

    std::ofstream formatFile(ffn);
    std::ofstream valueFile(vfn);
    std::ofstream rowFile(rfn);
    std::ofstream columnFile(cfn);

    // Check if all files were opened successfully
    if (!formatFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening format output files in write_vectors_to_file()." << std::endl;
        std::abort();
    }
    if (!valueFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening value vector output files in write_vectors_to_file()." << std::endl;
        std::abort();
    }
    if (!rowFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening row vector output files in write_vectors_to_file()." << std::endl;
        std::abort();
    }
    if (!columnFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening column vector output files in write_vectors_to_file()." << std::endl;
        std::abort();
    }



    if (matrix_format_ == format::COO) {
        formatFile << "COO" << "\n";
        formatFile << rows_ << "\n";
        formatFile << cols_ << "\n";
        // Write the contents of the vectors to the respective files
        for (ITYPE i = 0; i < matrix_coo.values_.size(); ++i) {
            valueFile << matrix_coo.values_[i] << "\n";
            rowFile << matrix_coo.row_indices_[i] << "\n";
            columnFile << matrix_coo.col_indices_[i] << "\n";
        }
    } else if (matrix_format_ == format::CSR) {
        formatFile << "CSR" << "\n";
        formatFile << rows_ << "\n";
        formatFile << cols_ << "\n";
        // Write the contents of the vectors to the respective files
        for (ITYPE i = 0; i < matrix_csr.values_.size(); ++i) {
            valueFile << matrix_csr.values_[i] << "\n";
            columnFile << matrix_csr.col_indices_[i] << "\n";
        }
        for (ITYPE i = 0; i < (rows_+1); ++i) {
            rowFile << matrix_csr.row_ptrs_[i] << "\n";
        }
    } else if (matrix_format_ == format::CSC) {
        formatFile << "CSC" << "\n";
        formatFile << rows_ << "\n";
        formatFile << cols_ << "\n";
        // Write the contents of the vectors to the respective files
        for (ITYPE i = 0; i < matrix_csc.values_.size(); ++i) {
            valueFile << matrix_csc.values_[i] << "\n";
            rowFile << matrix_csc.row_indices_[i] << "\n";
        }
        for (ITYPE i = 0; i < (cols_+1); ++i) {
            columnFile << matrix_csc.col_ptrs_[i] << "\n";
        }
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format for matrix write_vectors_to_file() >>" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }

    // Close the files
    formatFile.close();
    valueFile.close();
    rowFile.close();
    columnFile.close();

}

// Member function to read matrix vectors to the file
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::read_vectors_from_file(const std::string &format_file_name, const std::string &val_file_name, const std::string &row_file_name, const std::string &col_file_name) {

    std::string ffn;
    std::string vfn;
    std::string rfn;
    std::string cfn;

    // Open three separate input files
    if ((val_file_name.empty()) ||  (row_file_name.empty()) ||  (col_file_name.empty())) {
        ffn = format_file_name + "_format.dat";
        vfn = format_file_name + "_value.dat";
        rfn = format_file_name + "_row.dat";
        cfn = format_file_name + "_column.dat";
    } else {
        ffn = format_file_name;
        vfn = val_file_name;
        rfn = row_file_name;
        cfn = col_file_name;
    }

    std::ifstream formatFile(ffn);
    std::ifstream valueFile(vfn);
    std::ifstream rowFile(rfn);
    std::ifstream columnFile(cfn);

    // Check if all files were opened successfully
    if (!formatFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening format input files in read_vectors_from_file()." << std::endl;
        std::abort();
    }
    if (!valueFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening value vector input files in read_vectors_from_file()." << std::endl;
        std::abort();
    }
    if (!rowFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening row vector input files in read_vectors_from_file()." << std::endl;
        std::abort();
    }
    if (!columnFile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening column vector input files in read_vectors_from_file()." << std::endl;
        std::abort();
    }

    assert((this->empty()) &&
      "MUI Error [matrix_io_info.h]: read_vectors_from_file() can only takes in null matrix or empty (all-zero) matrix");

    format format_store = this->matrix_format_;

    std::string file_matrix_format;

    formatFile >> file_matrix_format;
    formatFile >> rows_;
    formatFile >> cols_;

    std::string file_matrix_format_trim = string_to_upper(trim(file_matrix_format));


    if (file_matrix_format_trim == "COO") {
        this->clear_vectors();
        // Read the contents of the vectors from the respective files
        int val;
        while (valueFile >> val) {
            matrix_coo.values_.reserve(matrix_coo.values_.size()+1);
            matrix_coo.values_.emplace_back(val);
        }
        int row_idx;
        while (rowFile >> row_idx) {
            matrix_coo.row_indices_.reserve(matrix_coo.row_indices_.size()+1);
            matrix_coo.row_indices_.emplace_back(row_idx);
        }
        int col_idx;
        while (columnFile >> col_idx) {
            matrix_coo.col_indices_.reserve(matrix_coo.col_indices_.size()+1);
            matrix_coo.col_indices_.emplace_back(col_idx);
        }
        matrix_format_ = format::COO;
        nnz_ = matrix_coo.values_.size();
    } else if (file_matrix_format_trim == "CSR") {
        this->clear_vectors();
        // Read the contents of the vectors from the respective files
        int val;
        while (valueFile >> val) {
            matrix_csr.values_.reserve(matrix_csr.values_.size()+1);
            matrix_csr.values_.emplace_back(val);
        }

        int row_ptr;
        while (rowFile >> row_ptr) {
            matrix_csr.row_ptrs_.reserve(matrix_csr.row_ptrs_.size()+1);
            matrix_csr.row_ptrs_.emplace_back(row_ptr);
        }

        int col_idx;
        while (columnFile >> col_idx) {
            matrix_csr.col_indices_.reserve(matrix_csr.col_indices_.size()+1);
            matrix_csr.col_indices_.emplace_back(col_idx);
        }

        matrix_format_ = format::CSR;
        nnz_ = matrix_csr.values_.size();

    } else if (file_matrix_format_trim == "CSC") {
        this->clear_vectors();
        // Read the contents of the vectors from the respective files
        int val;
        while (valueFile >> val) {
            matrix_csc.values_.reserve(matrix_csc.values_.size()+1);
            matrix_csc.values_.emplace_back(val);
        }
        int row_idx;
        while (rowFile >> row_idx) {
            matrix_csc.row_indices_.reserve(matrix_csc.row_indices_.size()+1);
            matrix_csc.row_indices_.emplace_back(row_idx);
        }
        int col_ptr;
        while (columnFile >> col_ptr) {
            matrix_csc.col_ptrs_.reserve(matrix_csc.col_ptrs_.size()+1);
            matrix_csc.col_ptrs_.emplace_back(col_ptr);
        }
        matrix_format_ = format::CSC;
        nnz_ = matrix_csc.values_.size();
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format: " << file_matrix_format_trim << " for matrix read_vectors_from_file() >>" << std::endl;
        std::cerr << "    Please set the matrix format as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }

    // Close the files
    formatFile.close();
    valueFile.close();
    rowFile.close();
    columnFile.close();

    this->assert_valid_vector_size("matrix_io_info.h", "read_vectors_from_file()");

    if (format_store != matrix_format_) {
        if (format_store == format::COO) {
            this->format_conversion("COO", true, true, "overwrite");
        } else if (format_store == format::CSR) {
            this->format_conversion("CSR", true, true, "overwrite");
        } else if (format_store == format::CSC) {
            this->format_conversion("CSC", true, true, "overwrite");
        } else {
            std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format for matrix write_vectors_to_file()" << std::endl;
            std::cerr << "    Please set the format_store as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }
    }
}

// Overloading << operator to output matrix in CSV format
template<typename ITYPE, typename VTYPE>
std::ostream& operator << (std::ostream &ofile, const sparse_matrix<ITYPE,VTYPE> &exist_mat) {

    // Check if file was opened successfully
    if (!ofile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening output files in overloaded << operator." << std::endl;
        std::abort();
    }

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

// Function to overloading >> operator to read matrix from a file in CSV format with lines start with "//" as comment lines
template<typename ITYPE, typename VTYPE>
std::istream& operator>>(std::istream &ifile, sparse_matrix<ITYPE,VTYPE> &exist_mat) {

    // Check if file was opened successfully
    if (!ifile) {
        std::cerr << "MUI Error [matrix_io_info.h]: Error opening output files in overloaded >> operator." << std::endl;
        std::abort();
    }

    assert((exist_mat.empty()) &&
      "MUI Error [matrix_io_info.h]: Overloading >> operator can only takes in null matrix or empty (all-zero) matrix");

    std::string format_store = exist_mat.get_format();

    std::string rawLine;

    std::vector<ITYPE> tempRowIndex;
    std::vector<ITYPE> tempColIndex;
    std::vector<VTYPE> tempValue;

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
                tempRowIndex.reserve(tempRowIndex.size()+1);
                tempColIndex.reserve(tempColIndex.size()+1);
                tempValue.reserve(tempValue.size()+1);
                tempRowIndex.emplace_back(row);
                tempColIndex.emplace_back(colCount);
                tempValue.emplace_back(val);
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
        sparse_matrix<ITYPE,VTYPE> temp_matrix(row, col, "COO", tempValue, tempRowIndex, tempColIndex);
        exist_mat.resize(row, col);
        exist_mat.format_conversion("COO", false, false);
        exist_mat.copy(temp_matrix);
        temp_matrix.set_zero();
    } else {
        assert(((exist_mat.get_rows() == row) && (exist_mat.get_cols() == col)) &&
          "MUI Error [matrix_io_info.h]: Matrix size mismatching between existing matrix and read in matrix in overloading >> operator ");
        sparse_matrix<ITYPE,VTYPE> temp_matrix(row, col, "COO", tempValue, tempRowIndex, tempColIndex);
        exist_mat.format_conversion("COO", false, false);
        exist_mat.copy(temp_matrix);
        temp_matrix.set_zero();
    }

    if (format_store != exist_mat.get_format()) {
        if (format_store == "COO") {
            exist_mat.format_conversion("COO", true, true, "overwrite");
        } else if (format_store == "CSR") {
            exist_mat.format_conversion("CSR", true, true, "overwrite");
        } else if (format_store == "CSC") {
            exist_mat.format_conversion("CSC", true, true, "overwrite");
        } else {
            std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format: " << format_store << " for matrix operator >>" << std::endl;
            std::cerr << "    Please set the format_store as:" << std::endl;
            std::cerr << "    format::COO: COOrdinate format" << std::endl;
            std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
            std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
            std::abort();
        }
    }
    return ifile;
}

// Member function to get the value at a given position
template<typename ITYPE, typename VTYPE>
VTYPE sparse_matrix<ITYPE,VTYPE>::get_value(ITYPE r, ITYPE c) const {
    assert(((r < rows_) && (r >= 0) && (c < cols_) && (c >= 0)) &&
        "MUI Error [matrix_io_info.h]: Matrix index out of range in get_value function");

    if (matrix_format_ == format::COO) {
    	for (ITYPE i = 0; i < static_cast<ITYPE>(matrix_coo.row_indices_.size()); ++i) {
            if (matrix_coo.row_indices_[i] == r && matrix_coo.col_indices_[i] == c) {
                return matrix_coo.values_[i];
            }
        }
        // Return 0 if the element is not found
        return static_cast<VTYPE>(0);
    } else if (matrix_format_ == format::CSR) {
        // Find the row range in the row_ptrs_ vector
        ITYPE row_start = matrix_csr.row_ptrs_[r];
        ITYPE row_end = matrix_csr.row_ptrs_[r + 1];

        // Search for the column index within the row range
        auto it = std::lower_bound(matrix_csr.col_indices_.begin()+row_start, matrix_csr.col_indices_.begin()+row_end, c);

        // Check if the element exists in the CSR vector
        if (it != matrix_csr.col_indices_.begin()+row_end && *it == c) {
            ITYPE index = std::distance(matrix_csr.col_indices_.begin(), it);
            return matrix_csr.values_[index];
        }
        // Return 0 if the element is not found
        return static_cast<VTYPE>(0);
    } else if (matrix_format_ == format::CSC) {
        // Find the column range in the col_ptrs_ vector
        ITYPE col_start = matrix_csc.col_ptrs_[c];
        ITYPE col_end = matrix_csc.col_ptrs_[c + 1];

        // Search for the row index within the column range
        auto it = std::lower_bound(matrix_csc.row_indices_.begin()+col_start, matrix_csc.row_indices_.begin()+col_end, r);

        // Check if the element exists in the CSC vector
        if (it != matrix_csc.row_indices_.begin()+col_end && *it == r) {
            ITYPE index = std::distance(matrix_csc.row_indices_.begin(), it);
            return matrix_csc.values_[index];
        }
        // Return 0 if the element is not found
        return static_cast<VTYPE>(0);
    } else {

        std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format for matrix get_value" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();

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
    if (matrix_format_ == format::COO) {
        vec_temp.reserve(matrix_coo.values_.size());
        for (ITYPE i = 0; i < matrix_coo.values_.size(); ++i) {
            vec_temp.emplace_back(std::make_pair(matrix_coo.row_indices_[i], matrix_coo.col_indices_[i]));
        }
    } else if (matrix_format_ == format::CSR) {
        vec_temp.reserve(matrix_csr.values_.size());
        for (ITYPE row = 0; row < matrix_csr.row_ptrs_.size()-1; ++row) {
            ITYPE row_start = matrix_csr.row_ptrs_[row];
            ITYPE row_end = matrix_csr.row_ptrs_[row + 1];

            // Iterate over the non-zero elements in the row
            for (ITYPE i = row_start; i < row_end; ++i) {
                ITYPE column = matrix_csr.col_indices_[i];
                vec_temp.emplace_back(std::make_pair(row, column));
            }
        }
    } else if (matrix_format_ == format::CSC) {
        vec_temp.reserve(matrix_csc.values_.size());
        for (ITYPE column = 0; column < matrix_csc.col_ptrs_.size()-1; ++column) {
            ITYPE col_start = matrix_csc.col_ptrs_[column];
            ITYPE col_end = matrix_csc.col_ptrs_[column + 1];

            // Iterate over the non-zero elements in the column
            for (ITYPE i = col_start; i < col_end; ++i) {
                ITYPE row = matrix_csc.row_indices_[i];
                vec_temp.emplace_back(std::make_pair(row, column));
            }
        }
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format for matrix operator >>" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }
    return vec_temp;
}

// Member function to get number of non-zero elements
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::non_zero_elements_count() const {
    if (matrix_format_ == format::COO) {
        if (nnz_ != static_cast<ITYPE>(matrix_coo.values_.size())) {
            std::cerr << "MUI Error [matrix_io_info.h]: Mismatch matrix_coo.values_ size (" << matrix_coo.values_.size() << ") with number of non-zeros (" << nnz_ <<") in non_zero_elements_count()" << std::endl;
            std::abort();
        }
    } else if (matrix_format_ == format::CSR) {
        if (nnz_ != static_cast<ITYPE>(matrix_csr.values_.size())) {
            std::cerr << "MUI Error [matrix_io_info.h]: Mismatch matrix_csr.values_ size (" << matrix_csr.values_.size() << ") with number of non-zeros (" << nnz_ <<") in non_zero_elements_count()" << std::endl;
            std::abort();
        }
    } else if (matrix_format_ == format::CSC) {
        if (nnz_ != static_cast<ITYPE>(matrix_csc.values_.size())) {
            std::cerr << "MUI Error [matrix_io_info.h]: Mismatch matrix_csc.values_ size (" << matrix_csc.values_.size() << ") with number of non-zeros (" << nnz_ <<") in non_zero_elements_count()" << std::endl;
            std::abort();
        }
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: Unrecognised matrix format for matrix non_zero_elements_count()" << std::endl;
        std::cerr << "    Please set the matrix_format_ as:" << std::endl;
        std::cerr << "    format::COO: COOrdinate format" << std::endl;
        std::cerr << "    format::CSR (default): Compressed Sparse Row format" << std::endl;
        std::cerr << "    format::CSC: Compressed Sparse Column format" << std::endl;
        std::abort();
    }

    return nnz_;
}

// Member function to check whether the matrix contains all zero elements
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::empty() const {
    if ((matrix_coo.values_.empty()) &&
        (matrix_coo.row_indices_.empty()) &&
        (matrix_coo.col_indices_.empty()) &&
        (matrix_csr.values_.empty()) &&
        (matrix_csr.col_indices_.empty()) &&
        (matrix_csc.values_.empty()) &&
        (matrix_csc.row_indices_.empty())) {
        return true;
    } else {
        return false;
    }
}

// Member function to get the format of the matrix
template<typename ITYPE, typename VTYPE>
std::string sparse_matrix<ITYPE,VTYPE>::get_format() const
{
    std::string matrix_format;

    if (matrix_format_ == format::COO) {
        matrix_format = "COO";
    } else if (matrix_format_ == format::CSR) {
        matrix_format = "CSR";
    } else if (matrix_format_ == format::CSC) {
        matrix_format = "CSC";
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: unknown matrix format" << std::endl;
        std::abort();
    }

    return matrix_format;
}

// Member function to check if the sparse matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::is_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_io_info.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "is_sorted_unique()";
    } else {
        function_name = function_name_input;
    }

    if (matrix_format_ == format::COO) {
        return this->is_coo_sorted_unique(file_name, function_name);
    } else if (matrix_format_ == format::CSR) {
        return this->is_csr_sorted_unique(file_name, function_name);
    } else if (matrix_format_ == format::CSC) {
        return this->is_csc_sorted_unique(file_name, function_name);
    } else {
        std::cerr << "MUI Error [matrix_io_info.h]: unknown matrix format" << std::endl;
        std::abort();
    }
}

// **************************************************
// ********** Protected member functions ************
// **************************************************

// Protected member function to check if the COO matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::is_coo_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_io_info.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "is_coo_sorted_unique()";
    } else {
        function_name = function_name_input;
    }

    ITYPE numEntries = matrix_coo.values_.size();

    if (numEntries > 1) {
        for (ITYPE i = 1; i < numEntries; ++i) {
            // Compare the current entry with the previous one
            if (matrix_coo.row_indices_[i] < matrix_coo.row_indices_[i - 1]) {
                // Row index is not sorted
                if (DEBUG) {
                    std::cout << "MUI [" << file_name << "]: The COO type matrix is not sorted (sorted row index check failed) in " << function_name <<  std::endl;
                }
                return false;
            } else if (matrix_coo.row_indices_[i] == matrix_coo.row_indices_[i - 1]) {
                // Row index is the same, check column index
                if (matrix_coo.col_indices_[i] < matrix_coo.col_indices_[i - 1]) {
                    // Column index is not sorted
                    if (DEBUG) {
                        std::cout << "MUI [" << file_name << "]: The COO type matrix is not sorted (sorted column index check failed) in " << function_name <<  std::endl;
                    }
                    return false;
                } else if (matrix_coo.col_indices_[i] == matrix_coo.col_indices_[i - 1]) {
                    // Column index has duplicate elements
                    if (DEBUG) {
                        std::cout << "MUI [" << file_name << "]: The COO type matrix exists duplicated elements (unique column index check failed) in " << function_name <<  std::endl;
                    }
                    return false;
                }
            }
        }
    }
    return true;
}

// Protected member function to check if the CSR matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::is_csr_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_io_info.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "is_csr_sorted_unique()";
    } else {
        function_name = function_name_input;
    }

    ITYPE numEntries = matrix_csr.values_.size();

    if (numEntries > 1) {
        for(ITYPE i = 0; i < rows_; ++i){
            if (matrix_csr.row_ptrs_[i] > matrix_csr.row_ptrs_[i+1]) {
                // Row pointers is not sorted
                if (DEBUG) {
                    std::cout << "MUI [" << file_name << "]: The CSR type matrix is not sorted (sorted row pointers check failed) in " << function_name <<  std::endl;
                }
                return false;
            }
            for(ITYPE j = matrix_csr.row_ptrs_[i] + 1; j < matrix_csr.row_ptrs_[i+1]; ++j){
                if(matrix_csr.col_indices_[j-1] > matrix_csr.col_indices_[j]){
                    // Column indices is not sorted
                    if (DEBUG) {
                        std::cout << "MUI [" << file_name << "]: The CSR type matrix is not sorted (sorted column index check failed) in " << function_name <<  std::endl;
                    }
                    return false;
                } else if (matrix_csr.col_indices_[j-1] == matrix_csr.col_indices_[j]) {
                    // Column indices is not unique
                    if (DEBUG) {
                        std::cout << "MUI [" << file_name << "]: The CSR type matrix is not unique (deduplicated column index check failed) in " << function_name <<  std::endl;
                    }
                    return false;
                }
            }
        }
    }
    return true;
}

// Protected member function to check if the CSC matrix is sorted and deduplicated
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::is_csc_sorted_unique(const std::string &file_name_input, const std::string &function_name_input) const {

    std::string file_name;
    std::string function_name;

    if (file_name_input.empty()) {
        file_name = "matrix_io_info.h";
    } else {
        file_name = file_name_input;
    }

    if (function_name_input.empty()) {
        function_name = "is_csc_sorted_unique()";
    } else {
        function_name = function_name_input;
    }

    ITYPE numEntries = matrix_csc.values_.size();

    if (numEntries > 1) {
        for(ITYPE i = 0; i < cols_; ++i){
            if (matrix_csc.col_ptrs_[i] > matrix_csc.col_ptrs_[i+1]) {
                // Column pointers is not sorted
                if (DEBUG) {
                    std::cout << "MUI [" << file_name << "]: The CSC type matrix is not sorted (sorted column pointers check failed) in " << function_name <<  std::endl;
                }
                return false;
            }
            for(ITYPE j = matrix_csc.col_ptrs_[i] + 1; j < matrix_csc.col_ptrs_[i+1]; ++j){
                if(matrix_csc.row_indices_[j-1] > matrix_csc.row_indices_[j]){
                    // Row indices is not sorted
                    if (DEBUG) {
                        std::cout << "MUI [" << file_name << "]: The CSC type matrix is not sorted (sorted row index check failed) in " << function_name <<  std::endl;
                    }
                    return false;
                } else if (matrix_csc.row_indices_[j-1] == matrix_csc.row_indices_[j]) {
                    // Row indices is not unique
                    if (DEBUG) {
                        std::cout << "MUI [" << file_name << "]: The CSC type matrix is not unique (deduplicated row index check failed) in " << function_name <<  std::endl;
                    }
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
