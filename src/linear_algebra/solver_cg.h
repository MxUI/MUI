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
 * @file solver_cg.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Implementation to solve problem A.x = b using the Conjugate Gradient method.
 *        Based on Baratta, Igor, Chris Richardson, and Garth Wells. "Performance 
 *        analysis of matrix-free conjugate gradient kernels using SYCL." 
 *        In International Workshop on OpenCL, pp. 1-10. 2022.
 */

#ifndef MUI_CONJUGATE_GRADIENT_H_
#define MUI_CONJUGATE_GRADIENT_H_

#include <cmath>

namespace mui {
namespace linalg {

// Constructor for one-dimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
conjugate_gradient_1d<ITYPE, VTYPE>::conjugate_gradient_1d(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b, VTYPE cg_solve_tol, ITYPE cg_max_iter, preconditioner<ITYPE,VTYPE>* M)
    : A_(A),
      b_(b),
      cg_solve_tol_(cg_solve_tol),
      cg_max_iter_(cg_max_iter),
      M_(M){
        assert(b_.get_cols() == 1 &&
                "MUI Error [solver_cg.h]: Number of column of b matrix must be 1");
        x_.resize(A_.get_rows(),1);
        r_.resize(A_.get_rows(),1);
        z_.resize(A_.get_rows(),1);
        p_.resize(A_.get_rows(),1);
}

// Constructor for multidimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
conjugate_gradient<ITYPE, VTYPE>::conjugate_gradient(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b, VTYPE cg_solve_tol, ITYPE cg_max_iter, preconditioner<ITYPE,VTYPE>* M)
    : A_(A),
      b_(b),
      cg_solve_tol_(cg_solve_tol),
      cg_max_iter_(cg_max_iter),
      M_(M){
        assert(A_.get_rows() == b_.get_rows() &&
                "MUI Error [solver_cg.h]: Number of rows of A matrix must be the same as the number of rows of b matrix");
        b_column_.resize(b_.get_rows(),1);
        x_.resize(b_.get_rows(),b_.get_cols());
        x_init_column_.resize(b_.get_rows(),1);
}

// Destructor for one-dimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
conjugate_gradient_1d<ITYPE, VTYPE>::~conjugate_gradient_1d() {
    // Deallocate the memory for matrices
    A_.set_zero();
    x_.set_zero();
    b_.set_zero();
    r_.set_zero();
    z_.set_zero();
    p_.set_zero();
    // Set properties to null
    cg_solve_tol_ = 0;
    cg_max_iter_ = 0;
    // Deallocate the memory for preconditioner pointer
    if(M_!=nullptr) {
        M_ = nullptr;
        delete[] M_;
    }
}

// Destructor for multidimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
conjugate_gradient<ITYPE, VTYPE>::~conjugate_gradient() {
    // Deallocate the memory for matrices
    A_.set_zero();
    x_.set_zero();
    b_.set_zero();
    b_column_.set_zero();
    x_init_column_.set_zero();
    // Set properties to null
    cg_solve_tol_ = 0;
    cg_max_iter_ = 0;
    // Deallocate the memory for preconditioner pointer
    if(M_!=nullptr) {
        M_ = nullptr;
        delete[] M_;
    }
}

// Member function for one-dimensional Conjugate Gradient solver to solve
template<typename ITYPE, typename VTYPE>
std::pair<ITYPE, VTYPE> conjugate_gradient_1d<ITYPE, VTYPE>::solve(sparse_matrix<ITYPE,VTYPE> x_init) {
    if (!x_init.empty()){
        assert(((x_init.get_rows() == x_.get_rows()) && (x_init.get_cols() == x_.get_cols())) &&
                "MUI Error [solver_cg.h]: Size of x_init matrix mismatch with size of x_ matrix");
        // Initialize x_ with x_init
        x_.copy(x_init);
        // Initialise r_ with b-Ax0
        sparse_matrix<ITYPE,VTYPE> Ax0 = A_* x_init;
        r_.copy(b_-Ax0);
    } else {
        // Initialise r_ with b
        r_.copy(b_);
    }

    // Initialise z_ with r_
    z_.copy(r_);

    if (M_) {
        sparse_matrix<ITYPE,VTYPE> tempZ(z_.get_rows(), z_.get_cols());
        tempZ = M_->apply(z_);
        z_.set_zero();
        z_.copy(tempZ);
    }

    // Initialise p_ with z_
    p_.copy(z_);

    VTYPE r_norm0 = r_.dot_product(z_);
    assert(std::abs(r_norm0) >= std::numeric_limits<VTYPE>::min() &&
            "MUI Error [solver_cg.h]: Divide by zero assert for r_norm0");
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
        assert(std::abs(p_dot_Ap) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_cg.h]: Divide by zero assert for p_dot_Ap");
        VTYPE alpha = r_norm / p_dot_Ap;
        for (ITYPE j = 0; j < A_.get_rows(); ++j) {
            x_.add_scalar(j, 0, (alpha * (p_.get_value(j,0))));
            r_.subtract_scalar(j, 0, (alpha * (Ap.get_value(j,0))));
        }

        z_.set_zero();
        z_.copy(r_);

        if (M_) {
            sparse_matrix<ITYPE,VTYPE> tempZ(z_.get_rows(), z_.get_cols());
            tempZ = M_->apply(z_);
            z_.set_zero();
            z_.copy(tempZ);
        }

        VTYPE updated_r_norm = r_.dot_product(z_);
        assert(std::abs(r_norm) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_cg.h]: Divide by zero assert for r_norm");
        VTYPE beta = updated_r_norm / r_norm;
        r_norm = updated_r_norm;
        for (ITYPE j = 0; j < A_.get_rows(); ++j) {
            p_.set_value(j, 0, (z_.get_value(j,0)+(beta*p_.get_value(j,0))));
        }

        r_norm_rel = std::sqrt(r_norm/r_norm0);
        if (r_norm_rel <= cg_solve_tol_) {
            break;
        }
    }
    return std::make_pair(acturalKIterCount,r_norm_rel);
}

// Member function for multidimensional Conjugate Gradient solver to solve
template<typename ITYPE, typename VTYPE>
std::pair<ITYPE, VTYPE> conjugate_gradient<ITYPE, VTYPE>::solve(sparse_matrix<ITYPE,VTYPE> x_init) {
    if (!x_init.empty()){
        assert(((x_init.get_rows() == b_.get_rows()) && (x_init.get_cols() == b_.get_cols())) &&
                "MUI Error [solver_cg.h]: Size of x_init matrix mismatch with size of b_ matrix");
    }

    std::pair<ITYPE, VTYPE> cgReturn;
    for (ITYPE j = 0; j < b_.get_cols(); ++j) {
        b_column_.set_zero();
        b_column_ = b_.segment(0,(b_.get_rows()-1),j,j);
        conjugate_gradient_1d<ITYPE, VTYPE> cg(A_, b_column_, cg_solve_tol_, cg_max_iter_, M_);
        if (!x_init.empty()) {
            x_init_column_.set_zero();
            x_init_column_ = x_init.segment(0,(x_init.get_rows()-1),j,j);
        }
        std::pair<ITYPE, VTYPE> cgReturnTemp = cg.solve(x_init_column_);
        if (cgReturn.first < cgReturnTemp.first)
            cgReturn.first = cgReturnTemp.first;
        cgReturn.second += cgReturnTemp.second;
        sparse_matrix<ITYPE,VTYPE> x_column(b_.get_rows(),1);
        x_column = cg.getSolution();
        for (ITYPE i = 0; i < x_column.get_rows(); ++i) {
            x_.set_value(i, j, x_column.get_value(i,0));
        }
    }
    cgReturn.second /= b_.get_cols();

    return cgReturn;
}

// Member function for one-dimensional Conjugate Gradient solver to get the solution
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> conjugate_gradient_1d<ITYPE, VTYPE>::getSolution() {
    return x_;
}

// Member function for multidimensional Conjugate Gradient solver to get the solution
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> conjugate_gradient<ITYPE, VTYPE>::getSolution() {
    return x_;
}

} // linalg
} // mui

#endif /* MUI_CONJUGATE_GRADIENT_H_ */
