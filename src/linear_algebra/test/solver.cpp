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
 * @file solver.cpp
 * @author W. Liu
 * @date 28 January 2023
 * @brief Unit test on Conjugate Gradient solver with preconditioners.
 */

#include <iostream>
#include "../solver.h"

void test00 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 00: 2-by-2 Matrix A ==============" << std::endl;
    std::cout << "================== MUI Conjugate Gradient Solver ===========" << std::endl;
    std::cout << "======================= w/o Preconditioner =================" << std::endl;
    std::cout << "======================= w/o Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> cgReturn;

    int Css_size = 2;
    int Aas_size = 2;
    double cg_solve_tol = 1e-6;
    int cg_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 1);
    H_.resize(Aas_size, 1);
    H_ref.resize(Aas_size, 1);
    H_diff.resize(Aas_size, 1);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);

    Css.set_value(0, 0, 2.5409);
    Css.set_value(0, 1, -0.0113);
    Css.set_value(1, 0, -0.0113);
    Css.set_value(1, 1, 0.5287);

    Aas.set_value(0, 0, 1.3864);
    Aas.set_value(1, 0, 0.3719);

    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter);
    cgReturn = cg.solve();
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    std::cout << std::endl;
}

void test01 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 01: 2-by-2 Matrix A ==============" << std::endl;
    std::cout << "================== MUI Conjugate Gradient Solver ===========" << std::endl;
    std::cout << "======================= ILU Preconditioner =================" << std::endl;
    std::cout << "======================= w/o Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> cgReturn;

    int Css_size = 2;
    int Aas_size = 2;
    double cg_solve_tol = 1e-6;
    int cg_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 1);
    H_.resize(Aas_size, 1);
    H_ref.resize(Aas_size, 1);
    H_diff.resize(Aas_size, 1);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);

    Css.set_value(0, 0, 2.5409);
    Css.set_value(0, 1, -0.0113);
    Css.set_value(1, 0, -0.0113);
    Css.set_value(1, 1, 0.5287);

    Aas.set_value(0, 0, 1.3864);
    Aas.set_value(1, 0, 0.3719);

    mui::linalg::incomplete_lu_preconditioner<int,double> M(Css);
    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);
    cgReturn = cg.solve();
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    std::cout << std::endl;

}

void test02 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 02: 3-by-3 Matrix A ==============" << std::endl;
    std::cout << "================== MUI Conjugate Gradient Solver ===========" << std::endl;
    std::cout << "====================== SSOR Preconditioner =================" << std::endl;
    std::cout << "======================= w/o Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> cgReturn;

    int Css_size = 3;
    int Aas_size = 3;
    double cg_solve_tol = 1e-6;
    int cg_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 1);
    H_.resize(Aas_size, 1);
    H_ref.resize(Aas_size, 1);
    H_diff.resize(Aas_size, 1);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);

    Css.set_value(0, 0, 0.7444);
    Css.set_value(0, 1, -0.5055);
    Css.set_value(0, 2, -0.0851);
    Css.set_value(1, 0, -0.5055);
    Css.set_value(1, 1, 3.4858);
    Css.set_value(1, 2, 0.0572);
    Css.set_value(2, 0, -0.0851);
    Css.set_value(2, 1, 0.0572);
    Css.set_value(2, 2, 0.4738);

    Aas.set_value(0, 0, -0.0043);
    Aas.set_value(1, 0, 2.2501);
    Aas.set_value(2, 0, 0.2798);

    mui::linalg::symmetric_successive_over_relaxation_preconditioner<int,double> M(Css, 1.2);
    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);
    cgReturn = cg.solve();
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    std::cout << std::endl;

}

void test03 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 03: 6-by-6 Matrix A ==============" << std::endl;
    std::cout << "================== MUI Conjugate Gradient Solver ===========" << std::endl;
    std::cout << "======================= IC Preconditioner ==================" << std::endl;
    std::cout << "====================== with Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> cgReturn;

    int Css_size = 6;
    int Aas_size = 6;
    double cg_solve_tol = 1e-6;
    int cg_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_init; //< Initial Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 1);
    H_init.resize(Aas_size, 1);
    H_.resize(Aas_size, 1);
    H_ref.resize(Aas_size, 1);
    H_diff.resize(Aas_size, 1);

    H_init.set_value(0, 0, 9);
    H_init.set_value(1, 0, 0);
    H_init.set_value(2, 0, -2);
    H_init.set_value(3, 0, 3);
    H_init.set_value(4, 0, -2);
    H_init.set_value(5, 0, 5);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);
    H_ref.set_value(3, 0, 0.5449);
    H_ref.set_value(4, 0, 0.4237);
    H_ref.set_value(5, 0, 0.6459);

    Css.set_value(0, 0, 3.4430);
    Css.set_value(0, 1, -0.3963);
    Css.set_value(0, 2, 2.5012);
    Css.set_value(0, 3, 0.9525);
    Css.set_value(0, 4, 0.6084);
    Css.set_value(0, 5, -1.2728);

    Css.set_value(1, 0, -0.3963);
    Css.set_value(1, 1, 0.6015);
    Css.set_value(1, 2, -0.4108);
    Css.set_value(1, 3, -0.1359);
    Css.set_value(1, 4, -0.0295);
    Css.set_value(1, 5, 0.2630);

    Css.set_value(2, 0, 2.5012);
    Css.set_value(2, 1, -0.4108);
    Css.set_value(2, 2, 2.5927);
    Css.set_value(2, 3, 0.7072);
    Css.set_value(2, 4, 0.5587);
    Css.set_value(2, 5, -1.0613);

    Css.set_value(3, 0, 0.9525);
    Css.set_value(3, 1, -0.1359);
    Css.set_value(3, 2, 0.7072);
    Css.set_value(3, 3, 1.1634);
    Css.set_value(3, 4, 0.1920);
    Css.set_value(3, 5, -0.4344);

    Css.set_value(4, 0, 0.6084);
    Css.set_value(4, 1, -0.0295);
    Css.set_value(4, 2, 0.5587);
    Css.set_value(4, 3, 0.1920);
    Css.set_value(4, 4, 0.7636);
    Css.set_value(4, 5, -0.3261);

    Css.set_value(5, 0, -1.2728);
    Css.set_value(5, 1, 0.2630);
    Css.set_value(5, 2, -1.0613);
    Css.set_value(5, 3, -0.4344);
    Css.set_value(5, 4, -0.3261);
    Css.set_value(5, 5, 1.0869);

    Aas.set_value(0, 0, 3.0685);
    Aas.set_value(1, 0, 0.0484);
    Aas.set_value(2, 0, 2.5783);
    Aas.set_value(3, 0, 1.2865);
    Aas.set_value(4, 0, 0.8671);
    Aas.set_value(5, 0, -0.8230);

    mui::linalg::incomplete_cholesky_preconditioner<int,double> M(Css);
    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);
    cgReturn = cg.solve(H_init);
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    std::cout << std::endl;

}

void test04 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "= TEST 04: 6-by-6 Matrix A with multidimensional Matrix b ==" << std::endl;
    std::cout << "================== MUI Conjugate Gradient Solver ===========" << std::endl;
    std::cout << "======================= IC Preconditioner ==================" << std::endl;
    std::cout << "====================== with Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> cgReturn;

    int Css_size = 6;
    int Aas_size = 6;
    double cg_solve_tol = 1e-6;
    int cg_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_init; //< Initial Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 4);
    H_init.resize(Aas_size, 4);
    H_.resize(Aas_size, 4);
    H_ref.resize(Aas_size, 4);
    H_diff.resize(Aas_size, 4);

    H_init.set_value(0, 0, 9);
    H_init.set_value(1, 0, 0);
    H_init.set_value(2, 0, -2);
    H_init.set_value(3, 0, 3);
    H_init.set_value(4, 0, -2);
    H_init.set_value(5, 0, 5);
    H_init.set_value(0, 1, 9);
    H_init.set_value(1, 1, 0);
    H_init.set_value(2, 1, -2);
    H_init.set_value(3, 1, 3);
    H_init.set_value(4, 1, -2);
    H_init.set_value(5, 1, 5);
    H_init.set_value(0, 2, 9);
    H_init.set_value(1, 2, 0);
    H_init.set_value(2, 2, -2);
    H_init.set_value(3, 2, 3);
    H_init.set_value(4, 2, -2);
    H_init.set_value(5, 2, 5);
    H_init.set_value(0, 3, 9);
    H_init.set_value(1, 3, 0);
    H_init.set_value(2, 3, -2);
    H_init.set_value(3, 3, 3);
    H_init.set_value(4, 3, -2);
    H_init.set_value(5, 3, 5);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);
    H_ref.set_value(3, 0, 0.5449);
    H_ref.set_value(4, 0, 0.4237);
    H_ref.set_value(5, 0, 0.6459);
    H_ref.set_value(0, 1, 0.5488);
    H_ref.set_value(1, 1, 0.7152);
    H_ref.set_value(2, 1, 0.6028);
    H_ref.set_value(3, 1, 0.5449);
    H_ref.set_value(4, 1, 0.4237);
    H_ref.set_value(5, 1, 0.6459);
    H_ref.set_value(0, 2, 0.5488);
    H_ref.set_value(1, 2, 0.7152);
    H_ref.set_value(2, 2, 0.6028);
    H_ref.set_value(3, 2, 0.5449);
    H_ref.set_value(4, 2, 0.4237);
    H_ref.set_value(5, 2, 0.6459);
    H_ref.set_value(0, 3, 0.5488);
    H_ref.set_value(1, 3, 0.7152);
    H_ref.set_value(2, 3, 0.6028);
    H_ref.set_value(3, 3, 0.5449);
    H_ref.set_value(4, 3, 0.4237);
    H_ref.set_value(5, 3, 0.6459);

    Css.set_value(0, 0, 3.4430);
    Css.set_value(0, 1, -0.3963);
    Css.set_value(0, 2, 2.5012);
    Css.set_value(0, 3, 0.9525);
    Css.set_value(0, 4, 0.6084);
    Css.set_value(0, 5, -1.2728);

    Css.set_value(1, 0, -0.3963);
    Css.set_value(1, 1, 0.6015);
    Css.set_value(1, 2, -0.4108);
    Css.set_value(1, 3, -0.1359);
    Css.set_value(1, 4, -0.0295);
    Css.set_value(1, 5, 0.2630);

    Css.set_value(2, 0, 2.5012);
    Css.set_value(2, 1, -0.4108);
    Css.set_value(2, 2, 2.5927);
    Css.set_value(2, 3, 0.7072);
    Css.set_value(2, 4, 0.5587);
    Css.set_value(2, 5, -1.0613);

    Css.set_value(3, 0, 0.9525);
    Css.set_value(3, 1, -0.1359);
    Css.set_value(3, 2, 0.7072);
    Css.set_value(3, 3, 1.1634);
    Css.set_value(3, 4, 0.1920);
    Css.set_value(3, 5, -0.4344);

    Css.set_value(4, 0, 0.6084);
    Css.set_value(4, 1, -0.0295);
    Css.set_value(4, 2, 0.5587);
    Css.set_value(4, 3, 0.1920);
    Css.set_value(4, 4, 0.7636);
    Css.set_value(4, 5, -0.3261);

    Css.set_value(5, 0, -1.2728);
    Css.set_value(5, 1, 0.2630);
    Css.set_value(5, 2, -1.0613);
    Css.set_value(5, 3, -0.4344);
    Css.set_value(5, 4, -0.3261);
    Css.set_value(5, 5, 1.0869);

    Aas.set_value(0, 0, 3.0685);
    Aas.set_value(1, 0, 0.0484);
    Aas.set_value(2, 0, 2.5783);
    Aas.set_value(3, 0, 1.2865);
    Aas.set_value(4, 0, 0.8671);
    Aas.set_value(5, 0, -0.8230);
    Aas.set_value(0, 1, 3.0685);
    Aas.set_value(1, 1, 0.0484);
    Aas.set_value(2, 1, 2.5783);
    Aas.set_value(3, 1, 1.2865);
    Aas.set_value(4, 1, 0.8671);
    Aas.set_value(5, 1, -0.8230);
    Aas.set_value(0, 2, 3.0685);
    Aas.set_value(1, 2, 0.0484);
    Aas.set_value(2, 2, 2.5783);
    Aas.set_value(3, 2, 1.2865);
    Aas.set_value(4, 2, 0.8671);
    Aas.set_value(5, 2, -0.8230);
    Aas.set_value(0, 3, 3.0685);
    Aas.set_value(1, 3, 0.0484);
    Aas.set_value(2, 3, 2.5783);
    Aas.set_value(3, 3, 1.2865);
    Aas.set_value(4, 3, 0.8671);
    Aas.set_value(5, 3, -0.8230);

    mui::linalg::incomplete_cholesky_preconditioner<int,double> M(Css);
    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);
    cgReturn = cg.solve(H_init);
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    std::cout << std::endl;

}

void test05 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 05: 6-by-6 Matrix A ==============" << std::endl;
    std::cout << "================ MUI Gaussian Elimination Solver ===========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> geReturn;

    int Css_size = 6;
    int Aas_size = 6;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 1);
    H_.resize(Aas_size, 1);
    H_ref.resize(Aas_size, 1);
    H_diff.resize(Aas_size, 1);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);
    H_ref.set_value(3, 0, 0.5449);
    H_ref.set_value(4, 0, 0.4237);
    H_ref.set_value(5, 0, 0.6459);

    Css.set_value(0, 0, 3.4430);
    Css.set_value(0, 1, -0.3963);
    Css.set_value(0, 2, 2.5012);
    Css.set_value(0, 3, 0.9525);
    Css.set_value(0, 4, 0.6084);
    Css.set_value(0, 5, -1.2728);

    Css.set_value(1, 0, -0.3963);
    Css.set_value(1, 1, 0.6015);
    Css.set_value(1, 2, -0.4108);
    Css.set_value(1, 3, -0.1359);
    Css.set_value(1, 4, -0.0295);
    Css.set_value(1, 5, 0.2630);

    Css.set_value(2, 0, 2.5012);
    Css.set_value(2, 1, -0.4108);
    Css.set_value(2, 2, 2.5927);
    Css.set_value(2, 3, 0.7072);
    Css.set_value(2, 4, 0.5587);
    Css.set_value(2, 5, -1.0613);

    Css.set_value(3, 0, 0.9525);
    Css.set_value(3, 1, -0.1359);
    Css.set_value(3, 2, 0.7072);
    Css.set_value(3, 3, 1.1634);
    Css.set_value(3, 4, 0.1920);
    Css.set_value(3, 5, -0.4344);

    Css.set_value(4, 0, 0.6084);
    Css.set_value(4, 1, -0.0295);
    Css.set_value(4, 2, 0.5587);
    Css.set_value(4, 3, 0.1920);
    Css.set_value(4, 4, 0.7636);
    Css.set_value(4, 5, -0.3261);

    Css.set_value(5, 0, -1.2728);
    Css.set_value(5, 1, 0.2630);
    Css.set_value(5, 2, -1.0613);
    Css.set_value(5, 3, -0.4344);
    Css.set_value(5, 4, -0.3261);
    Css.set_value(5, 5, 1.0869);

    Aas.set_value(0, 0, 3.0685);
    Aas.set_value(1, 0, 0.0484);
    Aas.set_value(2, 0, 2.5783);
    Aas.set_value(3, 0, 1.2865);
    Aas.set_value(4, 0, 0.8671);
    Aas.set_value(5, 0, -0.8230);

    mui::linalg::gaussian_elimination_1d<int,double> ge(Css, Aas);
    geReturn = ge.solve();
    H_ = ge.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << std::endl;

}

void test06 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "= TEST 06: 6-by-6 Matrix A with multidimensional Matrix b ==" << std::endl;
    std::cout << "================ MUI Gaussian Elimination Solver ===========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> geReturn;

    int Css_size = 6;
    int Aas_size = 6;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 4);
    H_.resize(Aas_size, 4);
    H_ref.resize(Aas_size, 4);
    H_diff.resize(Aas_size, 4);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);
    H_ref.set_value(3, 0, 0.5449);
    H_ref.set_value(4, 0, 0.4237);
    H_ref.set_value(5, 0, 0.6459);
    H_ref.set_value(0, 1, 0.5488);
    H_ref.set_value(1, 1, 0.7152);
    H_ref.set_value(2, 1, 0.6028);
    H_ref.set_value(3, 1, 0.5449);
    H_ref.set_value(4, 1, 0.4237);
    H_ref.set_value(5, 1, 0.6459);
    H_ref.set_value(0, 2, 0.5488);
    H_ref.set_value(1, 2, 0.7152);
    H_ref.set_value(2, 2, 0.6028);
    H_ref.set_value(3, 2, 0.5449);
    H_ref.set_value(4, 2, 0.4237);
    H_ref.set_value(5, 2, 0.6459);
    H_ref.set_value(0, 3, 0.5488);
    H_ref.set_value(1, 3, 0.7152);
    H_ref.set_value(2, 3, 0.6028);
    H_ref.set_value(3, 3, 0.5449);
    H_ref.set_value(4, 3, 0.4237);
    H_ref.set_value(5, 3, 0.6459);

    Css.set_value(0, 0, 3.4430);
    Css.set_value(0, 1, -0.3963);
    Css.set_value(0, 2, 2.5012);
    Css.set_value(0, 3, 0.9525);
    Css.set_value(0, 4, 0.6084);
    Css.set_value(0, 5, -1.2728);

    Css.set_value(1, 0, -0.3963);
    Css.set_value(1, 1, 0.6015);
    Css.set_value(1, 2, -0.4108);
    Css.set_value(1, 3, -0.1359);
    Css.set_value(1, 4, -0.0295);
    Css.set_value(1, 5, 0.2630);

    Css.set_value(2, 0, 2.5012);
    Css.set_value(2, 1, -0.4108);
    Css.set_value(2, 2, 2.5927);
    Css.set_value(2, 3, 0.7072);
    Css.set_value(2, 4, 0.5587);
    Css.set_value(2, 5, -1.0613);

    Css.set_value(3, 0, 0.9525);
    Css.set_value(3, 1, -0.1359);
    Css.set_value(3, 2, 0.7072);
    Css.set_value(3, 3, 1.1634);
    Css.set_value(3, 4, 0.1920);
    Css.set_value(3, 5, -0.4344);

    Css.set_value(4, 0, 0.6084);
    Css.set_value(4, 1, -0.0295);
    Css.set_value(4, 2, 0.5587);
    Css.set_value(4, 3, 0.1920);
    Css.set_value(4, 4, 0.7636);
    Css.set_value(4, 5, -0.3261);

    Css.set_value(5, 0, -1.2728);
    Css.set_value(5, 1, 0.2630);
    Css.set_value(5, 2, -1.0613);
    Css.set_value(5, 3, -0.4344);
    Css.set_value(5, 4, -0.3261);
    Css.set_value(5, 5, 1.0869);

    Aas.set_value(0, 0, 3.0685);
    Aas.set_value(1, 0, 0.0484);
    Aas.set_value(2, 0, 2.5783);
    Aas.set_value(3, 0, 1.2865);
    Aas.set_value(4, 0, 0.8671);
    Aas.set_value(5, 0, -0.8230);
    Aas.set_value(0, 1, 3.0685);
    Aas.set_value(1, 1, 0.0484);
    Aas.set_value(2, 1, 2.5783);
    Aas.set_value(3, 1, 1.2865);
    Aas.set_value(4, 1, 0.8671);
    Aas.set_value(5, 1, -0.8230);
    Aas.set_value(0, 2, 3.0685);
    Aas.set_value(1, 2, 0.0484);
    Aas.set_value(2, 2, 2.5783);
    Aas.set_value(3, 2, 1.2865);
    Aas.set_value(4, 2, 0.8671);
    Aas.set_value(5, 2, -0.8230);
    Aas.set_value(0, 3, 3.0685);
    Aas.set_value(1, 3, 0.0484);
    Aas.set_value(2, 3, 2.5783);
    Aas.set_value(3, 3, 1.2865);
    Aas.set_value(4, 3, 0.8671);
    Aas.set_value(5, 3, -0.8230);

    mui::linalg::gaussian_elimination<int,double> ge(Css, Aas);
    geReturn = ge.solve();
    H_ = ge.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << std::endl;

}


void test07 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "= TEST 07: 6-by-6 Matrix A with multidimensional Matrix b ==" << std::endl;
    std::cout << "================== MUI Conjugate Gradient Solver ===========" << std::endl;
    std::cout << "=================== Diagonal Preconditioner ================" << std::endl;
    std::cout << "====================== with Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> cgReturn;

    int Css_size = 6;
    int Aas_size = 6;
    double cg_solve_tol = 1e-6;
    int cg_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_init; //< Initial Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 4);
    H_init.resize(Aas_size, 4);
    H_.resize(Aas_size, 4);
    H_ref.resize(Aas_size, 4);
    H_diff.resize(Aas_size, 4);

    H_init.set_value(0, 0, 9);
    H_init.set_value(1, 0, 0);
    H_init.set_value(2, 0, -2);
    H_init.set_value(3, 0, 3);
    H_init.set_value(4, 0, -2);
    H_init.set_value(5, 0, 5);
    H_init.set_value(0, 1, 9);
    H_init.set_value(1, 1, 0);
    H_init.set_value(2, 1, -2);
    H_init.set_value(3, 1, 3);
    H_init.set_value(4, 1, -2);
    H_init.set_value(5, 1, 5);
    H_init.set_value(0, 2, 9);
    H_init.set_value(1, 2, 0);
    H_init.set_value(2, 2, -2);
    H_init.set_value(3, 2, 3);
    H_init.set_value(4, 2, -2);
    H_init.set_value(5, 2, 5);
    H_init.set_value(0, 3, 9);
    H_init.set_value(1, 3, 0);
    H_init.set_value(2, 3, -2);
    H_init.set_value(3, 3, 3);
    H_init.set_value(4, 3, -2);
    H_init.set_value(5, 3, 5);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);
    H_ref.set_value(3, 0, 0.5449);
    H_ref.set_value(4, 0, 0.4237);
    H_ref.set_value(5, 0, 0.6459);
    H_ref.set_value(0, 1, 0.5488);
    H_ref.set_value(1, 1, 0.7152);
    H_ref.set_value(2, 1, 0.6028);
    H_ref.set_value(3, 1, 0.5449);
    H_ref.set_value(4, 1, 0.4237);
    H_ref.set_value(5, 1, 0.6459);
    H_ref.set_value(0, 2, 0.5488);
    H_ref.set_value(1, 2, 0.7152);
    H_ref.set_value(2, 2, 0.6028);
    H_ref.set_value(3, 2, 0.5449);
    H_ref.set_value(4, 2, 0.4237);
    H_ref.set_value(5, 2, 0.6459);
    H_ref.set_value(0, 3, 0.5488);
    H_ref.set_value(1, 3, 0.7152);
    H_ref.set_value(2, 3, 0.6028);
    H_ref.set_value(3, 3, 0.5449);
    H_ref.set_value(4, 3, 0.4237);
    H_ref.set_value(5, 3, 0.6459);

    Css.set_value(0, 0, 3.4430);
    Css.set_value(0, 1, -0.3963);
    Css.set_value(0, 2, 2.5012);
    Css.set_value(0, 3, 0.9525);
    Css.set_value(0, 4, 0.6084);
    Css.set_value(0, 5, -1.2728);

    Css.set_value(1, 0, -0.3963);
    Css.set_value(1, 1, 0.6015);
    Css.set_value(1, 2, -0.4108);
    Css.set_value(1, 3, -0.1359);
    Css.set_value(1, 4, -0.0295);
    Css.set_value(1, 5, 0.2630);

    Css.set_value(2, 0, 2.5012);
    Css.set_value(2, 1, -0.4108);
    Css.set_value(2, 2, 2.5927);
    Css.set_value(2, 3, 0.7072);
    Css.set_value(2, 4, 0.5587);
    Css.set_value(2, 5, -1.0613);

    Css.set_value(3, 0, 0.9525);
    Css.set_value(3, 1, -0.1359);
    Css.set_value(3, 2, 0.7072);
    Css.set_value(3, 3, 1.1634);
    Css.set_value(3, 4, 0.1920);
    Css.set_value(3, 5, -0.4344);

    Css.set_value(4, 0, 0.6084);
    Css.set_value(4, 1, -0.0295);
    Css.set_value(4, 2, 0.5587);
    Css.set_value(4, 3, 0.1920);
    Css.set_value(4, 4, 0.7636);
    Css.set_value(4, 5, -0.3261);

    Css.set_value(5, 0, -1.2728);
    Css.set_value(5, 1, 0.2630);
    Css.set_value(5, 2, -1.0613);
    Css.set_value(5, 3, -0.4344);
    Css.set_value(5, 4, -0.3261);
    Css.set_value(5, 5, 1.0869);

    Aas.set_value(0, 0, 3.0685);
    Aas.set_value(1, 0, 0.0484);
    Aas.set_value(2, 0, 2.5783);
    Aas.set_value(3, 0, 1.2865);
    Aas.set_value(4, 0, 0.8671);
    Aas.set_value(5, 0, -0.8230);
    Aas.set_value(0, 1, 3.0685);
    Aas.set_value(1, 1, 0.0484);
    Aas.set_value(2, 1, 2.5783);
    Aas.set_value(3, 1, 1.2865);
    Aas.set_value(4, 1, 0.8671);
    Aas.set_value(5, 1, -0.8230);
    Aas.set_value(0, 2, 3.0685);
    Aas.set_value(1, 2, 0.0484);
    Aas.set_value(2, 2, 2.5783);
    Aas.set_value(3, 2, 1.2865);
    Aas.set_value(4, 2, 0.8671);
    Aas.set_value(5, 2, -0.8230);
    Aas.set_value(0, 3, 3.0685);
    Aas.set_value(1, 3, 0.0484);
    Aas.set_value(2, 3, 2.5783);
    Aas.set_value(3, 3, 1.2865);
    Aas.set_value(4, 3, 0.8671);
    Aas.set_value(5, 3, -0.8230);

    mui::linalg::diagonal_preconditioner<int,double> M(Css);
    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);
    cgReturn = cg.solve(H_init);
    H_ = cg.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << cgReturn.first <<" with final r_norm_rel: " << cgReturn.second << std::endl;

    std::cout << std::endl;

}

void test08 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 08: 6-by-6 Matrix A ==============" << std::endl;
    std::cout << "========= MUI Biconjugate Gradient Stabilized Solver =======" << std::endl;
    std::cout << "======================= w/o Preconditioner =================" << std::endl;
    std::cout << "====================== with Initial Guess ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::pair<int, double> bicgstabReturn;

    int Css_size = 6;
    int Aas_size = 6;
    double bicgstab_solve_tol = 1e-6;
    int bicgstab_max_iter = 1000;

    mui::linalg::sparse_matrix<int,double> Css; //< Matrix of radial basis function evaluations between prescribed points
    mui::linalg::sparse_matrix<int,double> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points
    mui::linalg::sparse_matrix<int,double> H_init; //< Initial Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_; //< Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    Css.resize(Css_size, Css_size);
    Aas.resize(Aas_size, 4);
    H_init.resize(Aas_size, 4);
    H_.resize(Aas_size, 4);
    H_ref.resize(Aas_size, 4);
    H_diff.resize(Aas_size, 4);

    H_init.set_value(0, 0, 9);
    H_init.set_value(1, 0, 0);
    H_init.set_value(2, 0, -2);
    H_init.set_value(3, 0, 3);
    H_init.set_value(4, 0, -2);
    H_init.set_value(5, 0, 5);
    H_init.set_value(0, 1, 9);
    H_init.set_value(1, 1, 0);
    H_init.set_value(2, 1, -2);
    H_init.set_value(3, 1, 3);
    H_init.set_value(4, 1, -2);
    H_init.set_value(5, 1, 5);
    H_init.set_value(0, 2, 9);
    H_init.set_value(1, 2, 0);
    H_init.set_value(2, 2, -2);
    H_init.set_value(3, 2, 3);
    H_init.set_value(4, 2, -2);
    H_init.set_value(5, 2, 5);
    H_init.set_value(0, 3, 9);
    H_init.set_value(1, 3, 0);
    H_init.set_value(2, 3, -2);
    H_init.set_value(3, 3, 3);
    H_init.set_value(4, 3, -2);
    H_init.set_value(5, 3, 5);

    H_ref.set_value(0, 0, 0.5488);
    H_ref.set_value(1, 0, 0.7152);
    H_ref.set_value(2, 0, 0.6028);
    H_ref.set_value(3, 0, 0.5449);
    H_ref.set_value(4, 0, 0.4237);
    H_ref.set_value(5, 0, 0.6459);
    H_ref.set_value(0, 1, 0.5488);
    H_ref.set_value(1, 1, 0.7152);
    H_ref.set_value(2, 1, 0.6028);
    H_ref.set_value(3, 1, 0.5449);
    H_ref.set_value(4, 1, 0.4237);
    H_ref.set_value(5, 1, 0.6459);
    H_ref.set_value(0, 2, 0.5488);
    H_ref.set_value(1, 2, 0.7152);
    H_ref.set_value(2, 2, 0.6028);
    H_ref.set_value(3, 2, 0.5449);
    H_ref.set_value(4, 2, 0.4237);
    H_ref.set_value(5, 2, 0.6459);
    H_ref.set_value(0, 3, 0.5488);
    H_ref.set_value(1, 3, 0.7152);
    H_ref.set_value(2, 3, 0.6028);
    H_ref.set_value(3, 3, 0.5449);
    H_ref.set_value(4, 3, 0.4237);
    H_ref.set_value(5, 3, 0.6459);

    Css.set_value(0, 0, 3.4430);
    Css.set_value(0, 1, -0.3963);
    Css.set_value(0, 2, 2.5012);
    Css.set_value(0, 3, 0.9525);
    Css.set_value(0, 4, 0.6084);
    Css.set_value(0, 5, -1.2728);

    Css.set_value(1, 0, -0.3963);
    Css.set_value(1, 1, 0.6015);
    Css.set_value(1, 2, -0.4108);
    Css.set_value(1, 3, -0.1359);
    Css.set_value(1, 4, -0.0295);
    Css.set_value(1, 5, 0.2630);

    Css.set_value(2, 0, 2.5012);
    Css.set_value(2, 1, -0.4108);
    Css.set_value(2, 2, 2.5927);
    Css.set_value(2, 3, 0.7072);
    Css.set_value(2, 4, 0.5587);
    Css.set_value(2, 5, -1.0613);

    Css.set_value(3, 0, 0.9525);
    Css.set_value(3, 1, -0.1359);
    Css.set_value(3, 2, 0.7072);
    Css.set_value(3, 3, 1.1634);
    Css.set_value(3, 4, 0.1920);
    Css.set_value(3, 5, -0.4344);

    Css.set_value(4, 0, 0.6084);
    Css.set_value(4, 1, -0.0295);
    Css.set_value(4, 2, 0.5587);
    Css.set_value(4, 3, 0.1920);
    Css.set_value(4, 4, 0.7636);
    Css.set_value(4, 5, -0.3261);

    Css.set_value(5, 0, -1.2728);
    Css.set_value(5, 1, 0.2630);
    Css.set_value(5, 2, -1.0613);
    Css.set_value(5, 3, -0.4344);
    Css.set_value(5, 4, -0.3261);
    Css.set_value(5, 5, 1.0869);

    Aas.set_value(0, 0, 3.0685);
    Aas.set_value(1, 0, 0.0484);
    Aas.set_value(2, 0, 2.5783);
    Aas.set_value(3, 0, 1.2865);
    Aas.set_value(4, 0, 0.8671);
    Aas.set_value(5, 0, -0.8230);
    Aas.set_value(0, 1, 3.0685);
    Aas.set_value(1, 1, 0.0484);
    Aas.set_value(2, 1, 2.5783);
    Aas.set_value(3, 1, 1.2865);
    Aas.set_value(4, 1, 0.8671);
    Aas.set_value(5, 1, -0.8230);
    Aas.set_value(0, 2, 3.0685);
    Aas.set_value(1, 2, 0.0484);
    Aas.set_value(2, 2, 2.5783);
    Aas.set_value(3, 2, 1.2865);
    Aas.set_value(4, 2, 0.8671);
    Aas.set_value(5, 2, -0.8230);
    Aas.set_value(0, 3, 3.0685);
    Aas.set_value(1, 3, 0.0484);
    Aas.set_value(2, 3, 2.5783);
    Aas.set_value(3, 3, 1.2865);
    Aas.set_value(4, 3, 0.8671);
    Aas.set_value(5, 3, -0.8230);

    mui::linalg::biconjugate_gradient_stabilized<int,double> bicgstab(Css, Aas, bicgstab_solve_tol, bicgstab_max_iter);

    bicgstabReturn = bicgstab.solve(H_init);
    H_ = bicgstab.getSolution();

    std::cout << "Matrix H_: " << std::endl;
    H_.print();

    std::cout << "Reference value of Matrix H_: " << std::endl;
    H_ref.print();

    H_diff = H_ - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    std::cout << "Total CG iteration number: " << bicgstabReturn.first <<" with final r_norm_rel: " << bicgstabReturn.second << std::endl;

    std::cout << std::endl;
}

int main() {

    // Perform test 00
    test00();
    // Perform test 01
    test01();
    // Perform test 02
    test02();
    // Perform test 03
    test03();
    // Perform test 04
    test04();
    // Perform test 05
    test05();
    // Perform test 06
    test06();
    // Perform test 07
    test07();
    // Perform test 08
    test08();

    return 0;
}
