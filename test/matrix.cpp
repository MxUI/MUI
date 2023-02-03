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
 * @file matrix.cpp
 * @author W. Liu
 * @date 28 January 2023
 * @brief Unit test on Sparse Matrix.
 */

#include <iostream>
#include <fstream>
#include "../linalg/matrix.h"

void test00 (mui::linalg::sparse_matrix<int,double> A,
            mui::linalg::sparse_matrix<int,double> B,
            mui::linalg::sparse_matrix<int,double> C,
            mui::linalg::sparse_matrix<int,double> D) {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "=============== TEST 00: Basic Matrix Setup and I/O ========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    std::cout << "Matrix A: " << std::endl;
    A.print();

    std::cout << "Matrix B: " << std::endl;
    B.print();

    std::cout << "Empty matrix C: " << std::endl;
    C.print();

    std::cout << "Identity matrix D: " << std::endl;
    D.print();

    // Output A matrix to a file in CSV format
    std::ofstream ofile("matrix.csv");
    ofile<< "// **************";
    ofile << "\n";
    ofile<< "// **** TEST ****";
    ofile << "\n";
    ofile<< "// **************";
    ofile << "\n";
    ofile << "//  ";
    ofile << "\n";
    ofile << A;
    ofile.close();

    mui::linalg::sparse_matrix<int,double> E;
    // Reads matrix from a file
    std::ifstream ifile("matrix.csv");
    ifile >> E;
    ifile.close();

    std::cout << "Matrix File I/O Test in CSV format" << std::endl;
    std::cout << "Read in matrix E (should equals to matrix A): " << std::endl;
    E.print();
}

void test01 (mui::linalg::sparse_matrix<int,double> A) {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==================== TEST 01: Matrix Segments ==============" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> B;
    B.resize(A.get_rows(),1);
    for (int j = 0; j < A.get_cols(); ++j) {
        B.set_zero();
        B = A.segment(0,(A.get_rows()-1),j,j);
        std::cout << "segment of matrix A -- the " << (j+1) << "th column: " << std::endl;
        B.print();
    }
}

void test02 (mui::linalg::sparse_matrix<int,double> A,
        mui::linalg::sparse_matrix<int,double> B) {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "============= TEST 02: Matrix Arithmetic Operations ========" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> C = A + B;
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
    C = A.hadamard_product(B);
    std::cout << "Hadamard product (A {*} B): " << std::endl;
    C.print();

    C.set_zero();
    C = A.transpose();
    std::cout << "Transpose of matrix A (A^T): " << std::endl;
    C.print();

    mui::linalg::sparse_matrix<int,double> F(2,2);
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

    mui::linalg::sparse_matrix<int,double> G(2,2);
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

    mui::linalg::sparse_matrix<int,double> H(3,3);
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
    // Create matrix A
    mui::linalg::sparse_matrix<int,double> A(3, 3);
    A.set_value(0, 0, 1);
    A.set_value(0, 2, 2);
    A.set_value(1, 1, 3);
    // Create matrix B
    mui::linalg::sparse_matrix<int,double> B(3, 3);
    B.set_value(0, 0, 4);
    B.set_value(0, 1, 5);
    B.set_value(1, 1, 6);
    // Create empty matrix C
    mui::linalg::sparse_matrix<int,double> C(3);
    // Create identity matrix D
    mui::linalg::sparse_matrix<int,double> D(3,"identity");

    // Perform test 00
    test00(A, B, C, D);
    // Perform test 01
    test01(A);
    // Perform test 02
    test02(A, B);
    // Perform test 03
    test03();
    // Perform test 04
    test04();

    return 0;
   }
