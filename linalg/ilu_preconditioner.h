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
 * @file ilu_preconditioner.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Class of ilu preconditioner.
 */

#ifndef MUI_ILU_PRECONDITIONER_H_
#define MUI_ILU_PRECONDITIONER_H_

#include "preconditioner_base.h"

namespace mui {
namespace linalg {

template<typename ITYPE, typename VTYPE>
class ilu_preconditioner : public preconditioner<ITYPE,VTYPE> {
private:
	sparse_matrix<ITYPE,VTYPE> L, U;

public:
    ilu_preconditioner(const sparse_matrix<ITYPE,VTYPE>& A) {
        // Initialize L and U matrices
        L.resize_null(A.get_rows(), A.get_cols());
        U.resize_null(A.get_rows(), A.get_cols());

        // Perform the ILU factorization
        for (ITYPE i = 0; i < A.get_rows(); ++i) {
            // Copy the diagonal element
            L.set_value(i, i, A.get_value(i, i));
            U.set_value(i, i, static_cast<VTYPE>(1.0));

            for (ITYPE j = 0; j < A.get_cols(); ++j) {
                if (j < i) {
                    // Copy the lower triangular elements
                	L.set_value(i, j, A.get_value(i, j));
                } else if (j > i) {
                	// Copy the upper triangular elements
                	U.set_value(i, j, A.get_value(i, j));
                }
            }

            // Perform the forward substitution step
            for (ITYPE j1 = 0; j1 < L.get_cols(); ++j1) {
                if (j1 < i) {
                    for (ITYPE j2 = 0; j2 < U.get_cols(); ++j2) {
                        if (j2 > j1) {
                            L.set_value(i, j2, (L.get_value(i, j2) - (L.get_value(i, j1) * U.get_value(j1, j2))));
                        }
                    }
                }
            }

            // Perform the backward substitution step
            for (ITYPE j1 = 0; j1 < U.get_cols(); ++j1) {
                if (j1 > i) {
                	for (ITYPE j2 = 0; j2 < L.get_cols(); ++j2) {
                        if (j2 < j1) {
                        	U.set_value(i, j2, (U.get_value(i, j2) / L.get_value(j1, j1)));
                        }
                	}

                }
            }
        }
    }

    sparse_matrix<ITYPE,VTYPE> apply(const sparse_matrix<ITYPE,VTYPE>& x) {
        assert((x.get_cols()==1) &&
            "MUI Error [ilu_preconditioner.h]: apply only works for column vectors");
    	sparse_matrix<ITYPE,VTYPE> y(x.get_rows(), x.get_cols());
    	sparse_matrix<ITYPE,VTYPE> z(x.get_rows(), x.get_cols());

        // Perform the forward substitution step
        for (ITYPE i = 0; i < x.get_rows(); ++i) {
        	VTYPE sum = 0;
            for (ITYPE j = 0; j < L.get_cols(); ++j) {
                if (j < i) {
                    sum += L.get_value(i, j) * y.get_value(j,0);
                }
            }
            y.set_value(i,0,(x.get_value(i,0)-sum));
        }

        // Perform the backward substitution step
        for (ITYPE i = x.get_rows() - 1; i >= 0; i--) {
        	VTYPE sum = 0;
            for (ITYPE j = 0; j < U.get_cols(); ++j) {
                if (j > i) {
                    sum += U.get_value(i, j) * z.get_value(j,0);
                }
            }
            z.set_value(i,0,((y.get_value(i,0) - sum) / L.get_value(i, i)));
        }
        return z;
    }
};

} // linalg
} // mui

#endif /* MUI_ILU_PRECONDITIONER_H_ */
