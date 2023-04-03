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
                  "MUI Error [matrix_arithmetic.h]: assignment operator '=' only works for empty (all zero elements) matrix");

        if((rows_ == 0) && (cols_ == 0)){
            rows_ = exist_mat.rows_;
            cols_ = exist_mat.cols_;
        }

        assert(((rows_ == exist_mat.rows_) && (cols_ == exist_mat.cols_)) &&
                  "MUI Error [matrix_arithmetic.h]: matrix size mismatch in assignment operator '='");
        (*this).copy(exist_mat);
    }
    return *this;
}

} // linalg
} // mui

#endif /* MUI_MATRIX_MANIPULATION_H_ */
