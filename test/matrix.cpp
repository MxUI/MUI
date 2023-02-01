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

    std::cout << "Matrix A: " << std::endl;
    a.print();

    std::cout << "Matrix B: " << std::endl;
    b.print();

    mui::linalg::sparse_matrix<int,double> c = a + b;
    std::cout << "Addition of matrices (A + B): " << std::endl;
    c.print();

    c.set_zero();
    c = a - b;
    std::cout << "Subtraction of matrices (A - B): " << std::endl;
    c.print();

    c.set_zero();
    c = a * b;
    std::cout << "Multiplication of matrices (A * B): " << std::endl;
    c.print();

    c.set_zero();
    c = 8*a;
    std::cout << "Scalar multiplication (8 * A): " << std::endl;
    c.print();

    c.set_zero();
    c = a.hadamard_product(b);
    std::cout << "Hadamard product (A {*} B): " << std::endl;
    c.print();

    c.set_zero();
    c = a.transpose();
    std::cout << "Transpose of A matrix (A^T): " << std::endl;
    c.print();

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

    mui::linalg::sparse_matrix<int,double> d;
    // Reads matrix from a file
    std::ifstream ifile("matrix.csv");
    ifile >> d;
    ifile.close();

    std::cout << "Matrix File I/O Test in CSV format D = A: " << std::endl;
    d.print();

    return 0;
   }
