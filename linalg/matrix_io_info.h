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
 * @file matrix_io_info.h
 * @author W. Liu
 * @date 01 February 2023
 * @brief Implemantation of sparse matrix I/O functions.
 */

#ifndef MUI_MATRIX_IO_INFO_H_
#define MUI_MATRIX_IO_INFO_H_

#include <sstream>
#include <cassert>
#include <limits>

namespace mui {
namespace linalg {

// Function to left trim a string - helper function on matrix file I/O
std::string ltrim(const std::string &s) {
    const std::string WHITESPACE = " \n\r\t\f\v";
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

// Function to right trim a string - helper function on matrix file I/O
std::string rtrim(const std::string &s) {
    const std::string WHITESPACE = " \n\r\t\f\v";
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

// Function to trim a string on both sides - helper function on matrix file I/O
std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

// Member function to print matrix elements to the console
template<typename ITYPE, typename VTYPE>
void sparse_matrix<ITYPE,VTYPE>::print() {
    for (ITYPE i = 0; i < rows; ++i) {
        std::cout << "      ";
       for (ITYPE j = 0; j < cols; ++j){
           auto it = matrix.find(std::make_pair(i, j));
           if (it != matrix.end()) {
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
    assert(((r < rows) && (r >= 0) && (c < cols) && (c >= 0)) &&
        "MUI Error [matrix_io_info.h]: Matrix index out of range in get_value function");
    auto it = matrix.find(std::make_pair(r, c));
    if (it != matrix.end()) {
        return it->second;
    } else {
        return 0;
    }
}

// Member function to get the number of rows
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::get_rows() const {
    return rows;
}

// Member function to get the number of cols
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::get_cols() const {
    return cols;
}

// Member function to get non-zero elements
template<typename ITYPE, typename VTYPE>
std::vector<std::pair<ITYPE, ITYPE>> sparse_matrix<ITYPE,VTYPE>::get_non_zero_elements() const {
    std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
    for (auto const &nn_element : matrix) {
        vec_temp.push_back(std::make_pair(nn_element.first.first, nn_element.first.second));
    }
    return vec_temp;
}

// Member function to get number of non-zero elements
template<typename ITYPE, typename VTYPE>
ITYPE sparse_matrix<ITYPE,VTYPE>::non_zero_elements_count() const {
    return matrix.size();
}

// Member function to check whether the matrix contains all zero elements
template<typename ITYPE, typename VTYPE>
bool sparse_matrix<ITYPE,VTYPE>::empty() {
    return matrix.empty();
}

} // linalg
} // mui

#endif /* MUI_MATRIX_IO_INFO_H_ */
