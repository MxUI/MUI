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

int main()
{
    mui::linalg::sparse_matrix<int,double> a(3, 3);
    a.set_value(0, 0, 1);
    a.set_value(0, 2, 2);
    a.set_value(1, 1, 3);

    mui::linalg::sparse_matrix<int,double> b(3, 3);
    b.set_value(0, 0, 4);
    b.set_value(0, 1, 5);
    b.set_value(1, 1, 6);

    mui::linalg::sparse_matrix<int,double> c(3);

    mui::linalg::sparse_matrix<int,double> d(3,"identity");

    std::cout << "Matrix A: " << std::endl;
    a.print();

    std::cout << "Matrix B: " << std::endl;
    b.print();

    std::cout << "Empty matrix C: " << std::endl;
    c.print();

    std::cout << "Identity matrix D: " << std::endl;
    d.print();

    mui::linalg::sparse_matrix<int,double> e = a + b;
    std::cout << "Addition of matrices (A + B): " << std::endl;
    e.print();

    e.set_zero();
    e = a - b;
    std::cout << "Subtraction of matrices (A - B): " << std::endl;
    e.print();

    e.set_zero();
    e = a * b;
    std::cout << "Multiplication of matrices (A * B): " << std::endl;
    e.print();

    e.set_zero();
    e = 8 * a;
    std::cout << "Scalar multiplication (8 * A): " << std::endl;
    e.print();

    e.set_zero();
    e = a.hadamard_product(b);
    std::cout << "Hadamard product (A {*} B): " << std::endl;
    e.print();

    e.set_zero();
    e = a.transpose();
    std::cout << "Transpose of A matrix (A^T): " << std::endl;
    e.print();

    // Outputs matrix to a file in CSV format
    std::ofstream ofile("matrix.csv");
    ofile<< "// **************";
    ofile << "\n";
    ofile<< "// **** TEST ****";
    ofile << "\n";
    ofile<< "// **************";
    ofile << "\n";
    ofile << "//  ";
    ofile << "\n";
    ofile << a;
    ofile.close();

    mui::linalg::sparse_matrix<int,double> f;
    // Reads matrix from a file
    std::ifstream ifile("matrix.csv");
    ifile >> f;
    ifile.close();

    std::cout << "Matrix File I/O Test in CSV format F = A: " << std::endl;
    f.print();

    mui::linalg::sparse_matrix<int,double> g(2,2);
    g.set_value(0, 0, 4);
    g.set_value(0, 1, 3);
    g.set_value(1, 0, 6);
    g.set_value(1, 1, 3);

    std::cout << "Matrix G: " << std::endl;
    g.print();

    mui::linalg::sparse_matrix<int,double> h;
    mui::linalg::sparse_matrix<int,double> i;
    g.lu_decomposition(h,i);
    std::cout << "L matrix of LU decomposition of matrix G: " << std::endl;
    h.print();
    std::cout << "U matrix of LU decomposition of matrix G: " << std::endl;
    i.print();

    mui::linalg::sparse_matrix<int,double> j = h * i;
    std::cout << "Proof LU decomposition by L * U matrix: " << std::endl;
    j.print();

//    mui::linalg::sparse_matrix<int,double> I = g * h;
//    std::cout << "Proof inverse of matrix by (G * G^(-1)): " << std::endl;
//    I.print();

    return 0;
   }
