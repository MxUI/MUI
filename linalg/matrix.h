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
 * @file matrix.h
 * @author W. Liu
 * @date 27 January 2023
 * @brief Base class for sparse matrix includes basic arithmetic operations 
 * such as addition, subtraction, and multiplication.
 */

#ifndef MUI_SPARSE_MATRIX_H_
#define MUI_SPARSE_MATRIX_H_

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cassert>
#include <limits>

namespace mui {
namespace linalg {

template<typename ITYPE, typename VTYPE>
class sparse_matrix {
    private:
        std::map<std::pair<ITYPE, ITYPE>, VTYPE> matrix;
        ITYPE rows;
        ITYPE cols;

    public:
        // Constructor
        sparse_matrix<ITYPE,VTYPE>(ITYPE r, ITYPE c)
            : rows(r), cols(c) {}

        // Overload constructor - null
        sparse_matrix<ITYPE,VTYPE>()
            : rows(0), cols(0) {}

        // Overload constructor - taken another sparse_matrix object as an argument
        sparse_matrix<ITYPE,VTYPE>(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
              // Copy the data from the existing matrix
              rows = exist_mat.rows;
              cols = exist_mat.cols;
              std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
              vec_temp = exist_mat.get_non_zero_elements();
              for (auto elememt : vec_temp) {
                  if (std::abs(exist_mat.get_value(elememt.first, elememt.second)) > std::numeric_limits<VTYPE>::min())
                      matrix[std::make_pair(elememt.first, elememt.second)] = exist_mat.get_value(elememt.first, elememt.second);
              }
          }

        // Function to print matrix elements to the console
        void print() {
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

        // Overloading the operator to output matrix in CSV format
        friend std::ostream& operator << (std::ostream& ofile, const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
            for (ITYPE i = 0; i < exist_mat.get_rows(); ++i) {
                for (ITYPE j = 0; j < exist_mat.get_cols(); ++j) {
                    auto it = exist_mat.matrix.find(std::make_pair(i, j));
                    if (it != exist_mat.matrix.end()) {
                        if (j == (exist_mat.get_cols() - 1)) {
                            ofile << it->second;
                        } else {
                            ofile << it->second << ",";
                        }
                    } else {
                        if (j == (exist_mat.get_cols() - 1)) {
                            ofile << 0;
                        } else {
                            ofile << 0 << ",";
                        }
                    }
                }
                ofile << std::endl;
            }
            return ofile;
        }

        // Function to resize a null matrix
        void resize_null(ITYPE r, ITYPE c) {
            assert(((rows == 0) && (cols == 0)) &&
                    "MUI Error [matrix.h]: resize_null function only works for null matrix");
            rows = r;
            cols = c;
        }

        // Function to copy a sparse_matrix
        void copy(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
              // Copy the data from the existing matrix
              assert(matrix.empty() &&
                        "MUI Error [matrix.h]: copy function only works for empty (all zero elements) matrix");
              assert(((rows == exist_mat.rows) && (cols == exist_mat.cols)) &&
                        "MUI Error [matrix.h]: matrix size mismatch in copy function ");
              std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
              vec_temp = exist_mat.get_non_zero_elements();
              for (auto elememt : vec_temp) {
                  if (std::abs(exist_mat.get_value(elememt.first, elememt.second)) > std::numeric_limits<VTYPE>::min())
                      matrix[std::make_pair(elememt.first, elememt.second)] = exist_mat.get_value(elememt.first, elememt.second);
              }
          }

        // Function to get a segment of a sparse_matrix
        sparse_matrix<ITYPE,VTYPE> segment(ITYPE r_start, ITYPE r_end, ITYPE c_start, ITYPE c_end) {
              // get segment data from the existing matrix
              assert((r_end >= r_start) &&
                      "MUI Error [matrix.h]: segment function r_end has to be larger or equals to r_start");
              assert((c_end >= c_start) &&
                      "MUI Error [matrix.h]: segment function c_end has to be larger or equals to c_start");
              sparse_matrix<ITYPE,VTYPE> res((r_end-r_start+1), (c_end-c_start+1));
              for (auto elememt : matrix)
                  if ((elememt.first.first >=r_start)  &&
                      (elememt.first.first <=r_end)    &&
                      (elememt.first.second >=c_start) &&
                      (elememt.first.second <=c_end))
                      res.set_value(elememt.first.first, elememt.first.second, elememt.second);
              return res;
          }

        // Function to insert an element
        void set_value(ITYPE r, ITYPE c, VTYPE val) {
            if (std::abs(val) > std::numeric_limits<VTYPE>::min())
                matrix[std::make_pair(r, c)] = val;
        }

        // Function to get the value at a given position
        VTYPE get_value(ITYPE r, ITYPE c) const {
            auto it = matrix.find(std::make_pair(r, c));
            if (it != matrix.end()) {
                return it->second;
            } else {
                return 0;
            }
        }

        // Function to get the number of rows
        ITYPE get_rows() const {
            return rows;
        }

        // Function to get the number of cols
        ITYPE get_cols() const {
            return cols;
        }

        // Function to get non-zero elements
        std::vector<std::pair<ITYPE, ITYPE>> get_non_zero_elements() const {
            std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
            for (auto const &nn_element : matrix)
            {
                vec_temp.push_back(std::make_pair(nn_element.first.first, nn_element.first.second));
            }
            return vec_temp;
        }

        // Function to set all elements to zero and empty the sparse matrix
        void set_zero() {
            matrix.clear();
        }

        // Function to check whether the matrix contains all zero elements
        bool empty() {
            return matrix.empty();
        }

        // Function to add scalar to a specific elements
        void add_scalar(ITYPE i, ITYPE j, VTYPE value) {
            // check if the element exists
            if (matrix.find(std::make_pair(i, j)) != matrix.end()) {
                matrix[std::make_pair(i, j)] += value;
            } else {
                matrix[std::make_pair(i, j)] = value;
            }
        }

        // Function to subtract a scalar from a specific elements
        void subtract_scalar(ITYPE i, ITYPE j, VTYPE value) {
            // check if the element exists
            if (matrix.find(std::make_pair(i, j)) != matrix.end()) {
                matrix[std::make_pair(i, j)] -= value;
            } else {
                matrix[std::make_pair(i, j)] = -value;
            }
        }

        // Function to get transpose of matrix
        sparse_matrix<ITYPE,VTYPE> transpose() {
            sparse_matrix<ITYPE,VTYPE> res(cols, rows);
            for (auto elememt : matrix)
                res.set_value(elememt.first.second, elememt.first.first, elememt.second);
            return res;
        }

        // Overload addition operator to perform sparse matrix addition
        sparse_matrix<ITYPE,VTYPE> operator+(sparse_matrix<ITYPE,VTYPE> &addend) {

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
        sparse_matrix<ITYPE,VTYPE> operator-(sparse_matrix<ITYPE,VTYPE> &subtrahend) {
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
        sparse_matrix<ITYPE,VTYPE> operator*(sparse_matrix<ITYPE,VTYPE> &multiplicand) {
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

        // Overloaded assignment operator
        sparse_matrix<ITYPE,VTYPE>& operator=(const sparse_matrix<ITYPE,VTYPE> &exist_mat) {
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

        // Function of dot product
        VTYPE dot_product(sparse_matrix<ITYPE,VTYPE> &b_) const {
            assert(((cols == 1)&&(b_.cols == 1)) &&
                "MUI Error [matrix.h]: dot_product function only works for column vectors");
            sparse_matrix<ITYPE,VTYPE> tempThis(*this);
            sparse_matrix<ITYPE,VTYPE> thisT(tempThis.transpose());
            sparse_matrix<ITYPE,VTYPE> tempMat(thisT * b_);
            assert(((tempMat.get_rows() == 1)&&(tempMat.get_cols() == 1)) &&
                            "MUI Error [matrix.h]: result of dot_product function should be a scalar");
            return (tempMat.get_value(0,0));
        }

};

} // linalg
} // mui

#endif /* MUI_SPARSE_MATRIX_H_ */
