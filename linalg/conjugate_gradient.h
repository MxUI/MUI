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
 * @file conjugate_gradient.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Class to solve problem A.x = b using the Conjugate Gradient method.
 */

#ifndef MUI_CONJUGATE_GRADIENT_H_
#define MUI_CONJUGATE_GRADIENT_H_

#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include "matrix.h"

namespace mui {
namespace linalg {

template<typename ITYPE, typename VTYPE>
class conjugate_gradient {
    private:
        sparse_matrix<ITYPE,VTYPE> A_;
        sparse_matrix<ITYPE,VTYPE> x_;
        sparse_matrix<ITYPE,VTYPE> b_;
        sparse_matrix<ITYPE,VTYPE> r_;
        sparse_matrix<ITYPE,VTYPE> p_;
        VTYPE cg_solve_tol_;
        ITYPE cg_max_iter_;

    public:
        conjugate_gradient(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b, VTYPE cg_solve_tol= 1e-6, ITYPE cg_max_iter= 0) 
            : A_(A), 
              b_(b),
              cg_solve_tol_(cg_solve_tol),
              cg_max_iter_(cg_max_iter){
                assert(b_.get_cols() == 1 && 
                        "MUI Error [conjugate_gradient.h]: Number of column of b matrix must be 1");
                x_.resize_null(A_.get_rows(),1);
                r_.resize_null(A_.get_rows(),1);
                p_.resize_null(A_.get_rows(),1);
        }

        std::pair<ITYPE, VTYPE> solve(sparse_matrix<ITYPE,VTYPE> x_init = sparse_matrix<ITYPE,VTYPE>()) {
            if (!x_init.empty()){
                assert(((x_init.get_rows() == x_.get_rows()) && (x_init.get_cols() == x_.get_cols())) && 
                        "MUI Error [conjugate_gradient.h]: Size of x_init matrix mismatch with size of x_ matrix");
                // Initialize x_ with x_init
                x_.copy(x_init);
                // Initialise r_ with b-Ax0
                sparse_matrix<ITYPE,VTYPE> Ax0 = A_* x_init;
                r_.copy(b_-Ax0);
            } else {
                // Initialise r_ with b
                r_.copy(b_);
            }

            // Initialise p_ with r_
            p_.copy(r_);

            VTYPE r_norm0 = r_.dot_product(r_);
            assert(r_norm0 >= cg_solve_tol_ && 
                    "MUI Error [conjugate_gradient.h]: Divide by zero assert for r_norm0");
            VTYPE r_norm = r_norm0;
            VTYPE r_norm_rel = std::sqrt(r_norm/r_norm0);

            ITYPE kIter;
            if(cg_max_iter_ == 0) {
                kIter = std::numeric_limits<ITYPE>::max();
            } else {
                kIter = cg_max_iter_;
            }

            ITYPE acturalKIterCount = 0;
            
            for (ITYPE k = 0; k < kIter; ++k) {
                ++acturalKIterCount;
                sparse_matrix<ITYPE,VTYPE> Ap = A_*p_;
                VTYPE p_dot_Ap = p_.dot_product(Ap);
                assert(p_dot_Ap >= cg_solve_tol_ && 
                        "MUI Error [conjugate_gradient.h]: Divide by zero assert for p_dot_Ap");
                VTYPE alpha = r_norm / p_dot_Ap;
                for (ITYPE j = 0; j < A_.get_rows(); ++j) {
                    x_.add_scalar(j, 0, (alpha * (p_.get_value(j,0))));
                    r_.subtract_scalar(j, 0, (alpha * (Ap.get_value(j,0))));
                }
                VTYPE updated_r_norm = r_.dot_product(r_);
                assert(r_norm >= cg_solve_tol_ && 
                        "MUI Error [conjugate_gradient.h]: Divide by zero assert for r_norm");
                VTYPE beta = updated_r_norm / r_norm;
                r_norm = updated_r_norm;
                for (ITYPE j = 0; j < A_.get_rows(); ++j) {
                    p_.set_value(j, 0, (r_.get_value(j,0)+(beta*p_.get_value(j,0))));
                }

                r_norm_rel = std::sqrt(r_norm/r_norm0);
                if (r_norm_rel <= cg_solve_tol_) {
                    break;
                }
            }
            return std::make_pair(acturalKIterCount,r_norm_rel);
        }

        sparse_matrix<ITYPE,VTYPE> getSolution() {
            return x_;
        }

    private:

};

} // linalg
} // mui

#endif /* MUI_CONJUGATE_GRADIENT_H_ */