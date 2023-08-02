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
 * @file preconditioner_ic.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Implementation of Incomplete Cholesky preconditioner.
 */

#ifndef MUI_PRECONDITIONER_IC_H_
#define MUI_PRECONDITIONER_IC_H_

#include <math.h>
#include <limits>

namespace mui {
namespace linalg {

// Constructor
template<typename ITYPE, typename VTYPE>
incomplete_cholesky_preconditioner<ITYPE,VTYPE>::incomplete_cholesky_preconditioner(const sparse_matrix<ITYPE,VTYPE>& A) {
    // Initialise the lower triangular matrix
    L_.resize(A.get_rows(), A.get_cols());
    // Construct the lower triangular matrix
    for (ITYPE i = 0; i < A.get_rows(); ++i) {
        for (ITYPE j = 0; j <= i; ++j) {
            if (i == j) {
                VTYPE sum = 0;
                for (ITYPE k = 0; k < j; ++k) {
                 sum += std::pow(L_.get_value(j, k), 2);
                }
                L_.set_value(j, j, (std::sqrt(A.get_value(j, j) - sum)));
            } else {
                VTYPE sum = 0;
                for (ITYPE k = 0; k < j; ++k) {
                 sum += L_.get_value(i, k) * L_.get_value(j, k);
                }
                assert(std::abs(L_.get_value(j, j)) >= std::numeric_limits<VTYPE>::min() &&
                        "MUI Error [preconditioner_ic.h]: Divide by zero assert for L_.get_value(j, j)");
                L_.set_value(i, j, ((A.get_value(i, j) - sum) / L_.get_value(j, j)));
            }
        }
    }
 }

// Destructor
template<typename ITYPE, typename VTYPE>
incomplete_cholesky_preconditioner<ITYPE,VTYPE>::~incomplete_cholesky_preconditioner() {
    // Deallocate the memory for the lower triangular matrix
    L_.set_zero();
}

// Member function on preconditioner apply
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> incomplete_cholesky_preconditioner<ITYPE,VTYPE>::apply(const sparse_matrix<ITYPE,VTYPE>& x) {
    assert((x.get_cols()==1) &&
        "MUI Error [preconditioner_ic.h]: apply only works for column vectors");
    sparse_matrix<ITYPE,VTYPE> y(x.get_rows(), x.get_cols());
    sparse_matrix<ITYPE,VTYPE> z(x.get_rows(), x.get_cols());

    // Forward substitution
    for (ITYPE i = 0; i < x.get_rows(); ++i) {
        VTYPE sum = 0;
        for (ITYPE j = 0; j < i; ++j) {
            sum += L_.get_value(i, j) * y.get_value(j,0);
        }
        assert(std::abs(L_.get_value(i, i)) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [preconditioner_ic.h]: Divide by zero assert for L_.get_value(i, i)");
        y.set_value(i, 0, ((x.get_value(i, 0) - sum) / L_.get_value(i, i)));
    }

    // Backward substitution
    for (ITYPE i = x.get_rows() - 1; i >= 0; --i) {
        VTYPE sum = 0;
        for (ITYPE j = i + 1; j < x.get_rows(); ++j) {
            sum += L_.get_value(j, i) * z.get_value(j,0);
        }
        assert(std::abs(L_.get_value(i, i)) >= std::numeric_limits<VTYPE>::min() &&
          "MUI Error [preconditioner_ic.h]: Divide by zero assert for L_.get_value(i, i)");
        z.set_value(i, 0, ((y.get_value(i, 0) - sum) / L_.get_value(i, i)));
    }
    return z;
}

} // linalg
} // mui

#endif /* MUI_PRECONDITIONER_IC_H_ */
