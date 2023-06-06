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
 * @file matrix_arithmetics.cpp
 * @author W. Liu
 * @date 28 January 2023
 * @brief Unit test on Sparse Matrix arithmetic operations.
 */

#include <iostream>
#include <fstream>
#include "../matrix.h"

void test02 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "============= TEST 02: Matrix Arithmetic Operations ========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> A(3, 3, "CSR");
    A.set_value(0, 0, 1);
    A.set_value(0, 2, 2);
    A.set_value(1, 1, 3);
    std::cout << "Matrix A: " << std::endl;
    A.print();

    // Create matrix B
    std::vector<double> row0_vector{4,5,0};
    std::vector<double> row1_vector{0,6,0};
    std::vector<double> row2_vector{0,0,0};
    std::vector<std::vector<double>> dense_vector;
    dense_vector.push_back(row0_vector);
    dense_vector.push_back(row1_vector);
    dense_vector.push_back(row2_vector);
    mui::linalg::sparse_matrix<int,double> B(dense_vector, "CSR");
    std::cout << "Matrix B: " << std::endl;
    B.print();


    mui::linalg::sparse_matrix<int,double> C(3, 3, "CSR");
    C = A + B;
    std::cout << "Addition of matrices (A + B): " << std::endl;
    C.print();

    C.set_zero();
    C = A - B;
    std::cout << "Subtraction of matrices (A - B): " << std::endl;
    C.print();

    C.set_zero();
    C = A * B;
    std::cout << "Multiplication of matrices (A * B): " << std::endl;
    C.print();

    C.set_zero();
    C = 8 * A;
    std::cout << "Scalar multiplication (8 * A): " << std::endl;
    C.print();

    C.set_zero();
    C = A * 8;
    std::cout << "Scalar multiplication (A * 8): " << std::endl;
    C.print();

    mui::linalg::sparse_matrix<int,double> D(A.get_rows(), 1, "CSR");
    mui::linalg::sparse_matrix<int,double> E(B.get_rows(), 1, "CSR");
    double dot_result = 0;
    D = A.segment(0,(A.get_rows()-1),1,1);
    E = B.segment(0,(B.get_rows()-1),1,1);
    dot_result = D.dot_product(E);
    std::cout << "The 2nd row of matrix A dot product with the 2nd row of matrix B: " << std::endl;
    std::cout << "      " << dot_result << std::endl;

    C.set_zero();
    C = A.hadamard_product(B);
    std::cout << "Hadamard product (A {*} B): " << std::endl;
    C.print();

    C.set_zero();
    C = A.transpose();
    std::cout << "Transpose of matrix A (A^T): " << std::endl;
    C.print();

    mui::linalg::sparse_matrix<int,double> F(2, 2, "CSR");
    F.set_value(0, 0, -1);
    F.set_value(0, 1, 1.5);
    F.set_value(1, 0, 1);
    F.set_value(1, 1, -1);

    std::cout << "Matrix F: " << std::endl;
    F.print();

    mui::linalg::sparse_matrix<int,double> G = F.inverse();
    std::cout << "Inverse of matrix F^(-1) : " << std::endl;
    G.print();

    mui::linalg::sparse_matrix<int,double> H = F * G;
    std::cout << "Proof inverse of matrix by F * F^(-1) (should be an identity matrix): " << std::endl;
    H.print();
}

void test03 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "=============== TEST 03: Matrix LU Decomposition ===========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> G(2,2, "CSR");
    G.set_value(0, 0, 4);
    G.set_value(0, 1, 3);
    G.set_value(1, 0, 6);
    G.set_value(1, 1, 3);

    std::cout << "Matrix G: " << std::endl;
    G.print();

    mui::linalg::sparse_matrix<int,double> L;
    mui::linalg::sparse_matrix<int,double> U;
    G.lu_decomposition(L,U);
    std::cout << "L matrix of LU decomposition of matrix G: " << std::endl;
    L.print();
    std::cout << "U matrix of LU decomposition of matrix G: " << std::endl;
    U.print();

    mui::linalg::sparse_matrix<int,double> H = L * U;
    std::cout << "Proof LU decomposition by L * U (should equals to matrix G): " << std::endl;
    H.print();
}

void test04 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "=============== TEST 04: Matrix QR Decomposition ===========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> H(3,3, "CSR");
    H.set_value(0, 0, 12);
    H.set_value(0, 1, -51);
    H.set_value(0, 2, 4);
    H.set_value(1, 0, 6);
    H.set_value(1, 1, 167);
    H.set_value(1, 2, -68);
    H.set_value(2, 0, -4);
    H.set_value(2, 1, 24);
    H.set_value(2, 2, -41);

    std::cout << "Matrix H: " << std::endl;
    H.print();

    mui::linalg::sparse_matrix<int,double> Q;
    mui::linalg::sparse_matrix<int,double> R;
    H.qr_decomposition(Q,R);
    std::cout << "Q matrix of QR decomposition of matrix H: " << std::endl;
    Q.print();
    std::cout << "R matrix of QR decomposition of matrix H: " << std::endl;
    R.print();

    mui::linalg::sparse_matrix<int,double> I = Q * R;
    std::cout << "Proof QR decomposition by Q * R (should equals to matrix H): " << std::endl;
    I.print();
}

int main()
{
    // Perform test 02
    test02();
    // Perform test 03
    test03();
    // Perform test 04
    test04();

    return 0;
}


