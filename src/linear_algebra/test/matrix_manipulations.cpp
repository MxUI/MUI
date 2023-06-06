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
 * @file matrix_manipulations.cpp
 * @author W. Liu
 * @date 28 January 2023
 * @brief Unit test on Sparse Matrix manipulations.
 */

#include <iostream>
#include <fstream>
#include "../matrix.h"

void test00 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "============= TEST 00: Basic Matrix Setup, I/O & Info ======" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    // Create matrix A
    mui::linalg::sparse_matrix<int,double> A(3, 3, "CSR");
    A.set_value(0, 0, 1);
    A.set_value(0, 2, 2);
    A.set_value(1, 1, 3);

    // Create matrix B
    std::vector<double> row0_vector{4,5,0};
    std::vector<double> row1_vector{0,6,0};
    std::vector<double> row2_vector{0,0,0};
    std::vector<std::vector<double>> dense_vector;
    dense_vector.push_back(row0_vector);
    dense_vector.push_back(row1_vector);
    dense_vector.push_back(row2_vector);
    mui::linalg::sparse_matrix<int,double> B(dense_vector, "COO");

    // Create empty matrix C
    mui::linalg::sparse_matrix<int,double> C(3, "", "CSR");
    // Create identity matrix D
    mui::linalg::sparse_matrix<int,double> D(3,"identity", "CSC");

    std::cout << "Matrix A: " << std::endl;
    A.print();
    A.print_vectors();
    std::cout << "Get value at (1, 1) of Matrix A: " << A.get_value(1,1) << std::endl;
    std::cout << "Get format of Matrix A: " << A.get_format() << std::endl;
    std::cout << std::endl;

    std::cout << "Matrix B: " << std::endl;
    B.print();
    B.print_vectors();
    std::cout << "Get format of Matrix B: " << B.get_format() << std::endl;
    std::cout << "Get the number of rows of Matrix B: " << B.get_rows() << std::endl;
    std::cout << "Get the number of columns of Matrix B: " << B.get_cols() << std::endl;
    std::cout << std::endl;

    std::cout << "Empty matrix C: " << std::endl;
    C.print();
    C.print_vectors();
    std::cout << "Check whether Matrix C contains all zero elements: " << C.empty() << std::endl;
    std::cout << std::endl;

    std::cout << "Identity matrix D: " << std::endl;
    D.print();
    D.print_vectors();
    std::vector<std::pair<int, int>> d_non_zero_index;
    d_non_zero_index = D.get_non_zero_elements();
    std::cout << "Check whether Matrix D contains all zero elements: " << D.empty() << std::endl;
    std::cout << "Get the non-zero elements of Matrix D: " << std::endl;
    for (const auto& element : d_non_zero_index) {
        std::cout << element.first << ", " <<element.second << std::endl;
    }
    std::cout << "Get the number of non-zero elements of Matrix D: " << D.non_zero_elements_count() << std::endl;
    std::cout << std::endl;

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

    mui::linalg::sparse_matrix<int,double> E("CSC");
    // Reads matrix from a file
    std::ifstream ifile("matrix.csv");
    ifile >> E;
    ifile.close();

    std::cout << "Matrix File I/O Test in CSV format" << std::endl;
    std::cout << "Read in matrix E (should equals to matrix A): " << std::endl;
    E.print();
    E.print_vectors();
    std::cout << std::endl;

    // Output A matrix vectors to a file
    A.write_vectors_to_file("matrix_vector_format.csv",
      "matrix_vector_value.csv",
      "matrix_vector_row.csv",
      "matrix_vector_column.csv");

    mui::linalg::sparse_matrix<int,double> F("COO");
    // Reads matrix from a file
    F.read_vectors_from_file("matrix_vector_format.csv",
                            "matrix_vector_value.csv",
                            "matrix_vector_row.csv",
                            "matrix_vector_column.csv");

    std::cout << "Matrix vectors file I/O Test" << std::endl;
    std::cout << "Read in matrix F (should equals to matrix A): " << std::endl;
    F.print();
    F.print_vectors();
    std::cout << std::endl;

    // Create matrix G
    std::vector<double> value_vector{1,4,2,3};
    std::vector<int> row_vector{0,1,0,1};
    std::vector<int> col_vector{0,1,2,1};

    mui::linalg::sparse_matrix<int,double> G(3, 3, "COO", value_vector, row_vector, col_vector);
    std::cout << "Matrix G: " << std::endl;
    G.print();
    G.print_vectors();
    std::cout << "Check if the G sparse matrix is sorted and deduplicated: "<< G.is_sorted_unique() << std::endl;
    std::cout << std::endl;

}

void test01 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "================= TEST 01: Matrix manipulations ============" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    // Create matrix A
    mui::linalg::sparse_matrix<int,double> A(3, 3, "CSR");
    A.set_value(0, 0, 1);
    A.set_value(0, 2, 2);
    A.set_value(1, 1, 3);

    mui::linalg::sparse_matrix<int,double> B("CSC");
    B.copy(A);
    std::cout << "Copy of matrix A: " << std::endl;
    B.print();
    B.print_vectors();
    B.swap_elements(1,1,1,0);
    B.swap_elements(2,2,2,1);
    std::cout << "Swap value of elements (1,1) and (1,0) of matrix A: " << std::endl;
    B.print();
    B.print_vectors();
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> C("CSR");
    C.resize(A.get_rows(),1);
    for (int j = 0; j < A.get_cols(); ++j) {
        C.set_zero();
        C = A.segment(0,(A.get_rows()-1),j,j);
        std::cout << "segment of matrix A -- the " << (j+1) << "th column: " << std::endl;
        C.print();
    }
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> D(3,3,"CSC");
    D.set_value(6);
    std::cout << "Set value '6' for all elements of matrix D: " << std::endl;
    D.print();
    D.print_vectors();
    D.set_zero();
    std::cout << "Set zero of matrix D: " << std::endl;
    D.print();
    D.print_vectors();
    std::cout << std::endl;

    mui::linalg::sparse_matrix<int,double> E("CSR");
    E = A;
    std::cout << "Matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.add_scalar(1, 1, 3, true);
    std::cout << "Add '3' at element (1,1) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.subtract_scalar(1, 1, 2, true);
    std::cout << "Subtract '2' at element (1,1) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.multiply_scalar(1, 1, 3, true);
    std::cout << "multiply '3' at element (1,1) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.add_scalar(1, 1, -12, true);
    std::cout << "Add '-12' at element (1,1) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.subtract_scalar(0, 0, 1, true);
    std::cout << "Subtract '1' at element (0,0) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.multiply_scalar(0, 2, 0, true);
    std::cout << "multiply '0' at element (0,2) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.add_scalar(2, 2, 3, true);
    std::cout << "Add '3' at element (2,2) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.subtract_scalar(2, 1, 1, true);
    std::cout << "Subtract '1' at element (2,1) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    E.multiply_scalar(2, 0, 5, true);
    std::cout << "multiply '5' at element (2,0) of matrix E: " << std::endl;
    E.print();
    E.print_vectors();
    std::cout << std::endl;

}

int main()
{
    // Perform test 00
    test00();
    // Perform test 01
    test01();

    return 0;
}


