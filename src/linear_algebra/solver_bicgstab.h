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
 * @file solver_bicgstab.h
 * @author W. Liu
 * @date 13 Arpil 2023
 * @brief Implementation to solve problem A.x = b using the Biconjugate Gradient Stabilized method.
 *
 */

#ifndef MUI_BICONJUGATE_GRADIENT_STABILIZED_H_
#define MUI_BICONJUGATE_GRADIENT_STABILIZED_H_

#include <cmath>

namespace mui {
namespace linalg {

// Constructor for one-dimensional Biconjugate Gradient Stabilized solver
template<typename ITYPE, typename VTYPE>
biconjugate_gradient_stabilized_1d<ITYPE, VTYPE>::biconjugate_gradient_stabilized_1d(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b, VTYPE bicgstab_solve_tol, ITYPE bicgstab_max_iter, preconditioner<ITYPE,VTYPE>* M)
    : A_(A),
      b_(b),
      bicgstab_solve_tol_(bicgstab_solve_tol),
      bicgstab_max_iter_(bicgstab_max_iter),
      M_(M){
        assert(b_.get_cols() == 1 &&
                "MUI Error [solver_bicgstab.h]: Number of column of b matrix must be 1");
        x_.resize(A_.get_rows(),1);
        r_.resize(A_.get_rows(),1);
        rTilde_.resize(A_.get_rows(),1);
        v_.resize(A_.get_rows(),1);
        t_.resize(A_.get_rows(),1);
        p_.resize(A_.get_rows(),1);
        s_.resize(A_.get_rows(),1);
        h_.resize(A_.get_rows(),1);
        y_.resize(A_.get_rows(),1);
        z_.resize(A_.get_rows(),1);
        alpha_ = 1.0;
        beta_ = 0.0;
        omega_ = 1.0;
        rho_ = 1.0;
        rhoTilde_ = 1.0;
}

// Constructor for multidimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
biconjugate_gradient_stabilized<ITYPE, VTYPE>::biconjugate_gradient_stabilized(sparse_matrix<ITYPE,VTYPE> A, sparse_matrix<ITYPE,VTYPE> b, VTYPE bicgstab_solve_tol, ITYPE bicgstab_max_iter, preconditioner<ITYPE,VTYPE>* M)
    : A_(A),
      b_(b),
      bicgstab_solve_tol_(bicgstab_solve_tol),
      bicgstab_max_iter_(bicgstab_max_iter),
      M_(M){
        assert(A_.get_rows() == b_.get_rows() &&
                "MUI Error [solver_bicgstab.h]: Number of rows of A matrix must be the same as the number of rows of b matrix");
        b_column_.resize(b_.get_rows(),1);
        x_.resize(b_.get_rows(),b_.get_cols());
        x_init_column_.resize(b_.get_rows(),1);
}

// Destructor for one-dimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
biconjugate_gradient_stabilized_1d<ITYPE, VTYPE>::~biconjugate_gradient_stabilized_1d() {
    // Deallocate the memory for matrices
    A_.set_zero();
    x_.set_zero();
    b_.set_zero();
    r_.set_zero();
    rTilde_.set_zero();
    v_.set_zero();
    t_.set_zero();
    p_.set_zero();
    s_.set_zero();
    h_.set_zero();
    y_.set_zero();
    z_.set_zero();
    // Set properties to null
    alpha_ = 0.0;
    beta_ = 0.0;
    omega_ = 0.0;
    rho_ = 0.0;
    rhoTilde_ = 0.0;
    bicgstab_solve_tol_ = 0;
    bicgstab_max_iter_ = 0;
    // Deallocate the memory for preconditioner pointer
    if(M_!=nullptr) {
        M_ = nullptr;
        delete[] M_;
    }
}

// Destructor for multidimensional Conjugate Gradient solver
template<typename ITYPE, typename VTYPE>
biconjugate_gradient_stabilized<ITYPE, VTYPE>::~biconjugate_gradient_stabilized() {
    // Deallocate the memory for matrices
    A_.set_zero();
    x_.set_zero();
    b_.set_zero();
    b_column_.set_zero();
    x_init_column_.set_zero();
    // Set properties to null
    bicgstab_solve_tol_ = 0;
    bicgstab_max_iter_ = 0;
    // Deallocate the memory for preconditioner pointer
    if(M_!=nullptr) {
        M_ = nullptr;
        delete[] M_;
    }
}

// Member function for one-dimensional Conjugate Gradient solver to solve
template<typename ITYPE, typename VTYPE>
std::pair<ITYPE, VTYPE> biconjugate_gradient_stabilized_1d<ITYPE, VTYPE>::solve(sparse_matrix<ITYPE,VTYPE> x_init) {
    if (!x_init.empty()){
        assert(((x_init.get_rows() == x_.get_rows()) && (x_init.get_cols() == x_.get_cols())) &&
                "MUI Error [solver_bicgstab.h]: Size of x_init matrix mismatch with size of x_ matrix");
        // Initialize x_ with x_init
        x_.copy(x_init);
        // Initialise r_ with b-Ax0
        sparse_matrix<ITYPE,VTYPE> Ax0 = A_* x_init;
        r_.copy(b_-Ax0);
    } else {
        // Initialise r_ with b
        r_.copy(b_);
    }

    // Initialise rTilde_ with r_
    rTilde_.copy(r_);
    // Initialise p_ with r_
    p_.copy(r_);

    bool debug_switch = true;

    if (M_) {
        std::cout << "MUI Warning [solver_bicgstab.h]: Preconditioner is not yet supported by BiCGStab yet. "
                << "The preconditioner is ignored."<< std::endl;
        debug_switch = false;
    }

    if (M_ && debug_switch) {
        sparse_matrix<ITYPE,VTYPE> tempR(r_.get_rows(), r_.get_cols());
        tempR = M_->apply(r_);
        r_.set_zero();
        r_.copy(tempR);
        rTilde_.set_zero();
        rTilde_.copy(r_);
    }

    VTYPE r_norm0 = r_.dot_product(r_);
    assert(std::abs(r_norm0) >= std::numeric_limits<VTYPE>::min() &&
            "MUI Error [solver_bicgstab.h]: Divide by zero assert for r_norm0");
    VTYPE r_norm = r_norm0;
    VTYPE r_norm_rel = std::sqrt(r_norm/r_norm0);

    ITYPE kIter;
    if(bicgstab_max_iter_ == 0) {
        kIter = std::numeric_limits<ITYPE>::max();
    } else {
        kIter = bicgstab_max_iter_;
    }

    ITYPE acturalKIterCount = 0;

    for (ITYPE k = 0; k < kIter; ++k) {
        ++acturalKIterCount;

        rho_ = rTilde_.dot_product(r_);
        assert(std::abs(rhoTilde_) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_bicgstab.h]: Divide by zero assert for rhoTilde_");
        assert(std::abs(omega_) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_bicgstab.h]: Divide by zero assert for omega_");
        if (k>0) {
            beta_ = (rho_ / rhoTilde_) * (alpha_ / omega_);
            sparse_matrix<ITYPE,VTYPE> omega_dot_v = omega_ * v_;
            sparse_matrix<ITYPE,VTYPE> beta_p_omega_dot_v = beta_ * (p_ - omega_dot_v);
            p_.set_zero();
            p_.copy(r_ + beta_p_omega_dot_v);
        } else {
            p_.set_zero();
            p_.copy(r_);
        }

        y_.set_zero();
        y_.copy(p_);
        if (M_ && debug_switch) {
            sparse_matrix<ITYPE,VTYPE> tempY(y_.get_rows(), y_.get_cols());
            tempY = M_->apply(y_);
            y_.set_zero();
            y_.copy(tempY);
        }

        v_.set_zero();
        v_.copy(A_*y_);
        VTYPE rTilde_dot_v = rTilde_.dot_product(v_);
        assert(std::abs(rTilde_dot_v) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_bicgstab.h]: Divide by zero assert for rTilde_dot_v");
        alpha_ = rho_ / rTilde_dot_v;
        sparse_matrix<ITYPE,VTYPE> alpha_dot_p = alpha_ * p_;
        h_.set_zero();
        h_.copy(x_ + alpha_dot_p);
        sparse_matrix<ITYPE,VTYPE> alpha_dot_v = alpha_ * v_;
        s_.set_zero();
        s_.copy(r_ - alpha_dot_v);

        z_.set_zero();
        z_.copy(s_);
        if (M_ && debug_switch) {
            sparse_matrix<ITYPE,VTYPE> tempZ(z_.get_rows(), z_.get_cols());
            tempZ = M_->apply(z_);
            z_.set_zero();
            z_.copy(tempZ);
        }

        t_.set_zero();
        t_.copy(A_ * z_);
        VTYPE t_dot_t = t_.dot_product(t_);
        assert(std::abs(t_dot_t) >= std::numeric_limits<VTYPE>::min() &&
                "MUI Error [solver_bicgstab.h]: Divide by zero assert for t_dot_t");
        omega_ = t_.dot_product(s_) / t_dot_t;
        sparse_matrix<ITYPE,VTYPE> omega_dot_s = omega_ * s_;
        x_.set_zero();
        x_.copy(h_ + omega_dot_s);
        sparse_matrix<ITYPE,VTYPE> omega_dot_t = omega_ * t_;
        r_.set_zero();
        r_.copy(s_ - omega_dot_t);

        if (M_ && debug_switch) {
            sparse_matrix<ITYPE,VTYPE> tempR(r_.get_rows(), r_.get_cols());
            tempR = M_->apply(r_);
            r_.set_zero();
            r_.copy(tempR);
            rTilde_.set_zero();
            rTilde_.copy(r_);
        }

        r_norm = r_.dot_product(r_);
        r_norm_rel = std::sqrt(r_norm/r_norm0);
        if (r_norm_rel <= bicgstab_solve_tol_) {
            break;
        }

        rhoTilde_ = rho_;
    }
    return std::make_pair(acturalKIterCount,r_norm_rel);
}

// Member function for multidimensional Conjugate Gradient solver to solve
template<typename ITYPE, typename VTYPE>
std::pair<ITYPE, VTYPE> biconjugate_gradient_stabilized<ITYPE, VTYPE>::solve(sparse_matrix<ITYPE,VTYPE> x_init) {
    if (!x_init.empty()){
        assert(((x_init.get_rows() == b_.get_rows()) && (x_init.get_cols() == b_.get_cols())) &&
                "MUI Error [solver_bicgstab.h]: Size of x_init matrix mismatch with size of b_ matrix");
    }

    std::pair<ITYPE, VTYPE> bicgstabReturn;
    for (ITYPE j = 0; j < b_.get_cols(); ++j) {
        b_column_.set_zero();
        b_column_ = b_.segment(0,(b_.get_rows()-1),j,j);
        biconjugate_gradient_stabilized_1d<ITYPE, VTYPE> bicgstab(A_, b_column_, bicgstab_solve_tol_, bicgstab_max_iter_, M_);
        if (!x_init.empty()) {
            x_init_column_.set_zero();
            x_init_column_ = x_init.segment(0,(x_init.get_rows()-1),j,j);
        }
        std::pair<ITYPE, VTYPE> bicgstabReturnTemp = bicgstab.solve(x_init_column_);
        if (bicgstabReturn.first < bicgstabReturnTemp.first)
            bicgstabReturn.first = bicgstabReturnTemp.first;
        bicgstabReturn.second += bicgstabReturnTemp.second;
        sparse_matrix<ITYPE,VTYPE> x_column(b_.get_rows(),1);
        x_column = bicgstab.getSolution();
        for (ITYPE i = 0; i < x_column.get_rows(); ++i) {
            x_.set_value(i, j, x_column.get_value(i,0));
        }
    }
    bicgstabReturn.second /= b_.get_cols();

    return bicgstabReturn;
}

// Member function for one-dimensional Conjugate Gradient solver to get the solution
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> biconjugate_gradient_stabilized_1d<ITYPE, VTYPE>::getSolution() {
    return x_;
}

// Member function for multidimensional Conjugate Gradient solver to get the solution
template<typename ITYPE, typename VTYPE>
sparse_matrix<ITYPE,VTYPE> biconjugate_gradient_stabilized<ITYPE, VTYPE>::getSolution() {
    return x_;
}

} // linalg
} // mui

#endif /* MUI_BICONJUGATE_GRADIENT_STABILIZED_H_ */
