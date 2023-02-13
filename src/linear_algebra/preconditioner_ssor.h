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
 * @file preconditioner_ssor.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Class of Symmetric Successive Over-relaxation preconditioner.
 */

#ifndef MUI_PRECONDITIONER_SSOR_H_
#define MUI_PRECONDITIONER_SSOR_H_

#include <math.h>
#include <limits>

namespace mui {
namespace linalg {

// Constructor
template<typename ITYPE, typename VTYPE>
symmetric_successive_over_relaxation_preconditioner<ITYPE,VTYPE>::symmetric_successive_over_relaxation_preconditioner(const sparse_matrix<ITYPE,VTYPE>& A, VTYPE omega)
 : A_(A), omega_(omega) {}

// Destructor
template<typename ITYPE, typename VTYPE>
symmetric_successive_over_relaxation_preconditioner<ITYPE,VTYPE>::~symmetric_successive_over_relaxation_preconditioner() {
    // Deallocate the memory for the coefficient matrix
    A_.set_zero();
    // Set the relaxation parameter to null
    omega_ = 0;
}

// Member function on preconditioner apply
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> symmetric_successive_over_relaxation_preconditioner<ITYPE,VTYPE>::apply(const sparse_matrix<ITYPE,VTYPE>& x) {
    assert((x.get_cols()==1) &&
        "MUI Error [preconditioner_ssor.h]: apply only works for column vectors");
    sparse_matrix<ITYPE,VTYPE> y(x.get_rows(), x.get_cols());
    sparse_matrix<ITYPE,VTYPE> z(x.get_rows(), x.get_cols());

    // Forward substitution
    for (ITYPE i = 0; i < A_.get_rows(); ++i) {
        VTYPE sum = 0;
        for (ITYPE j = 0; j < A_.get_cols(); ++j) {
            if (i < j) {
                sum += A_.get_value(i,j) * y.get_value(j,0);
            } else if (i == j) {
                assert(std::abs(omega_ * A_.get_value(i,j)) >= std::numeric_limits<VTYPE>::min() &&
                  "MUI Error [preconditioner_ssor.h]: Divide by zero assert for omega_ * A_.get_value(i,j)");
                y.set_value(i, 0, ((x.get_value(i,0) - sum) / (omega_ * A_.get_value(i,j))));
            }
        }
    }

    // Back substitution
    for (ITYPE i = A_.get_rows() - 1; i >= 0; i--) {
        VTYPE sum = 0;
        for (ITYPE j = 0; j < A_.get_cols(); ++j) {
            if (i > j) {
                sum += A_.get_value(i,j) * z.get_value(j,0);
            } else if (i == j) {
                assert(std::abs(A_.get_value(i,j)) >= std::numeric_limits<VTYPE>::min() &&
                  "MUI Error [preconditioner_ssor.h]: Divide by zero assert for A_.get_value(i,j)");
                z.set_value(i, 0, (omega_ * (y.get_value(i,0) - sum) / A_.get_value(i,j)));
            }
        }
    }
    return z;
}

} // linalg
} // mui

#endif /* MUI_PRECONDITIONER_SSOR_H_ */
