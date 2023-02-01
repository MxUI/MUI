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

    std::cout << "Matrix A: " << std::endl;
    a.print();

    std::cout << "Matrix B: " << std::endl;
    b.print();

    std::cout << "Identity matrix C: " << std::endl;
    c.print();

    mui::linalg::sparse_matrix<int,double> d = a + b;
    std::cout << "Addition of matrices (A + B): " << std::endl;
    d.print();

    d.set_zero();
    d = a - b;
    std::cout << "Subtraction of matrices (A - B): " << std::endl;
    d.print();

    d.set_zero();
    d = a * b;
    std::cout << "Multiplication of matrices (A * B): " << std::endl;
    d.print();

    d.set_zero();
    d = 8*a;
    std::cout << "Scalar multiplication (8 * A): " << std::endl;
    d.print();

    d.set_zero();
    d = a.hadamard_product(b);
    std::cout << "Hadamard product (A {*} B): " << std::endl;
    d.print();

    d.set_zero();
    d = a.transpose();
    std::cout << "Transpose of A matrix (A^T): " << std::endl;
    d.print();

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

    mui::linalg::sparse_matrix<int,double> e;
    // Reads matrix from a file
    std::ifstream ifile("matrix.csv");
    ifile >> e;
    ifile.close();

    std::cout << "Matrix File I/O Test in CSV format E = A: " << std::endl;
    e.print();

    return 0;
   }
