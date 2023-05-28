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
 * @file preconditioner_diagonal.h
 * @author W. Liu
 * @date 28 March 2023
 * @brief Diagonal (Jacobi) preconditioner.
 */

#ifndef MUI_PRECONDITIONER_DIAGONAL_H_
#define MUI_PRECONDITIONER_DIAGONAL_H_

#include <math.h>
#include <limits>

namespace mui {
namespace linalg {

// Constructor
template<typename ITYPE, typename VTYPE>
diagonal_preconditioner<ITYPE,VTYPE>::diagonal_preconditioner(const sparse_matrix<ITYPE,VTYPE>& A) {
    // Initialise the lower triangular matrix
    inv_diag_.resize(A.get_rows(), A.get_cols());
    // Construct the inverse diagonal matrix
    for (int i = 0; i < A.get_rows(); i++) {
        if (std::abs(A.get_value(i,i)) >= std::numeric_limits<VTYPE>::min()) {
            inv_diag_.set_value(i, i, 1.0 / A.get_value(i,i));
        } else {
            inv_diag_.set_value(i, i, 1.0);
        }
    }
 }

// Destructor
template<typename ITYPE, typename VTYPE>
diagonal_preconditioner<ITYPE,VTYPE>::~diagonal_preconditioner() {
    // Deallocate the memory for the inverse diagonal matrix
    inv_diag_.set_zero();
}

// Member function on preconditioner apply
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> diagonal_preconditioner<ITYPE,VTYPE>::apply(const sparse_matrix<ITYPE,VTYPE>& x) {
    assert((x.get_cols()==1) &&
        "MUI Error [preconditioner_diagonal.h]: apply only works for column vectors");
    sparse_matrix<ITYPE,VTYPE> z(x.get_rows(), x.get_cols());

    for (int i = 0; i < x.get_rows(); i++) {
        if (std::abs(inv_diag_.get_value(i,i)) >= std::numeric_limits<VTYPE>::min()) {
            z.set_value(i, 0, inv_diag_.get_value(i,i)*x.get_value(i,0));
        } else {
            z.set_value(i, 0, 0.0);
        }
    }

    return z;
}

} // linalg
} // mui

#endif /* MUI_PRECONDITIONER_DIAGONAL_H_ */
