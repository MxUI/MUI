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
 * @file solver_ge.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Implementation to solve problem A.x = b using the Gaussian Elimination method.
 *
 */

#ifndef MUI_GAUSSIAN_ELINATION_H_
#define MUI_GAUSSIAN_ELINATION_H_

namespace mui {
namespace linalg {

// Constructor for one-dimensional Gaussian Elimination solver
template<typename ITYPE, typename VTYPE>
gaussian_elimination_1d<ITYPE, VTYPE>::gaussian_elimination_1d(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b)
    : A_(A),
      b_(b){
        assert(A_.get_rows() == b_.get_rows() &&
                "MUI Error [solver_ge.h]: Number of rows of A matrix must be the same as the number of rows of b matrix");
        assert(b_.get_cols() == 1 &&
                "MUI Error [solver_ge.h]: Number of column of b matrix must be 1");
        x_.resize(b_.get_rows(),b_.get_cols());
}

// Constructor for multidimensional Gaussian Elimination solver
template<typename ITYPE, typename VTYPE>
gaussian_elimination<ITYPE, VTYPE>::gaussian_elimination(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b)
    : A_(A),
      b_(b){
        assert(A_.get_rows() == b_.get_rows() &&
                "MUI Error [solver_ge.h]: Number of rows of A matrix must be the same as the number of rows of b matrix");
        x_.resize(b_.get_rows(),b_.get_cols());
        b_column_.resize(b_.get_rows(), 1);
}

// Destructor for one-dimensional Gaussian Elimination solver
template<typename ITYPE, typename VTYPE>
gaussian_elimination_1d<ITYPE, VTYPE>::~gaussian_elimination_1d() {
    // Deallocate the memory for matrices
    A_.set_zero();
    x_.set_zero();
    b_.set_zero();
}

// Destructor for multidimensional Gaussian Elimination solver
template<typename ITYPE, typename VTYPE>
gaussian_elimination<ITYPE, VTYPE>::~gaussian_elimination() {
    // Deallocate the memory for matrices
    A_.set_zero();
    x_.set_zero();
    b_.set_zero();
    b_column_.set_zero();
}

// Member function for one-dimensional Gaussian Elimination solver to solve
template<typename ITYPE, typename VTYPE>
std::pair<ITYPE, VTYPE> gaussian_elimination_1d<ITYPE, VTYPE>::solve(sparse_matrix<ITYPE,VTYPE> x_init) {
    // x_init does not in used in Gaussian Elimination

    for (ITYPE i = 0; i < A_.get_rows(); ++i) {
        // Find the pivot
        ITYPE pivot = i;
        for (ITYPE j = i+1; j < A_.get_cols(); ++j) {
            if (std::abs(A_.get_value(j, i)) > std::abs(A_.get_value(pivot, i))) {
                pivot = j;
            }
        }

        // Swap the pivot row with the current row
        for (ITYPE j = 0; j < A_.get_cols(); ++j) {
            A_.swap_elements(i, j, pivot, j);
        }
        b_.swap_elements(i, 0, pivot, 0);

      // Elimination
      for (ITYPE j = i+1; j < A_.get_cols(); ++j) {
          assert(std::abs(A_.get_value(i, i)) >= std::numeric_limits<VTYPE>::min() &&
                  "MUI Error [solver_ge.h]: Divide by zero assert for A_.get_value(i, i)");
          VTYPE factor = A_.get_value(j, i) / A_.get_value(i, i);
          for (ITYPE k = i; k < A_.get_cols(); ++k) {
              A_.set_value(j, k, (A_.get_value(j, k) - (factor * A_.get_value(i, k))));
          }
          b_.set_value(j, 0, (b_.get_value(j, 0) - (factor * b_.get_value(i, 0))));
      }
    }

    // Back substitution
    for (ITYPE i = A_.get_cols()-1; i >= 0; --i) {
        x_.set_value(i, 0, b_.get_value(i, 0));
        for (ITYPE j = i+1; j < A_.get_cols(); ++j) {
            x_.set_value(i, 0, (x_.get_value(i, 0) - (A_.get_value(i, j) * x_.get_value(j, 0))));
        }
        assert(std::abs(A_.get_value(i, i)) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_ge.h]: Divide by zero assert for A_.get_value(i, i)");
        x_.set_value(i, 0, (x_.get_value(i, 0) / A_.get_value(i, i)));
    }

    return std::make_pair(static_cast<ITYPE>(0),static_cast<VTYPE>(0));
}

// Member function for multidimensional Gaussian Elimination solver to solve
template<typename ITYPE, typename VTYPE>
std::pair<ITYPE, VTYPE> gaussian_elimination<ITYPE, VTYPE>::solve(sparse_matrix<ITYPE,VTYPE> x_init) {
    // x_init does not in used in Gaussian Elimination

    for (ITYPE j = 0; j < b_.get_cols(); ++j) {
        b_column_.set_zero();
        b_column_ = b_.segment(0,(b_.get_rows()-1),j,j);
        gaussian_elimination_1d<ITYPE, VTYPE> ge(A_, b_column_);
        std::pair<ITYPE, VTYPE> geReturnTemp = ge.solve();
        sparse_matrix<ITYPE,VTYPE> x_column(b_.get_rows(),1);
        x_column = ge.getSolution();
        for (ITYPE i = 0; i < x_column.get_rows(); ++i) {
            x_.set_value(i, j, x_column.get_value(i,0));
        }
    }

    return std::make_pair(static_cast<ITYPE>(0),static_cast<VTYPE>(0));
}

// Member function for one-dimensional Gaussian Elimination solver to get the solution
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> gaussian_elimination_1d<ITYPE, VTYPE>::getSolution() {
    return x_;
}

// Member function for multidimensional Gaussian Elimination solver to get the solution
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> gaussian_elimination<ITYPE, VTYPE>::getSolution() {
    return x_;
}

} // linalg
} // mui

#endif /* MUI_GAUSSIAN_ELINATION_H_ */
