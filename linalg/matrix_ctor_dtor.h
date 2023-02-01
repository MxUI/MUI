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
 * @file matrix_ctor_dtor.h
 * @author W. Liu
 * @date 01 February 2023
 * @brief Implemantation of sparse matrix constructors and destructor.
 */

#ifndef MUI_MATRIX_CTOR_DTOR_H_
#define MUI_MATRIX_CTOR_DTOR_H_

#include "matrix_io_info.h"

namespace mui {
namespace linalg {

// Constructor - takes in size of row and column to generate an empty matrix
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(ITYPE r, ITYPE c)
    : rows(r), cols(c), dummy_(0) {
}

// Constructor - null matrix
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix()
    : rows(0), cols(0), dummy_(0) {
}

// Constructor - takes in another sparse_matrix object as an argument
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(const sparse_matrix<ITYPE,VTYPE> &exist_mat)
    : rows(exist_mat.rows), cols(exist_mat.cols), dummy_(0) {
      // Copy the data from the existing matrix
      std::vector<std::pair<ITYPE, ITYPE>> vec_temp;
      vec_temp = exist_mat.get_non_zero_elements();
      for (auto element : vec_temp) {
          if (std::abs(exist_mat.get_value(element.first, element.second)) >= std::numeric_limits<VTYPE>::min())
              matrix[std::make_pair(element.first, element.second)] = exist_mat.get_value(element.first, element.second);
      }
}

// Constructor - generate various square matrices
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE>::sparse_matrix(ITYPE n, const std::string &token)
    : rows(n), cols(n), dummy_(0) {
    if(token.empty()) {
        // empty (all-zero) square matrix (Do nothing from the code perspective)
    } else if(trim(token)=="identity") {
        // identity square matrix
        for (ITYPE i = 0; i < n; ++i) {
            matrix[std::make_pair(i, i)] = static_cast<VTYPE>(1);
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
    matrix.clear();
    // set matrix properties to null
    rows = 0;
    cols = 0;
    dummy_ = 0;
}

} // linalg
} // mui

#endif /* MUI_MATRIX_CTOR_DTOR_H_ */
