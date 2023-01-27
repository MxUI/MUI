#include <iostream>
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
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
	   for (int j = 0; j < 3; j++)
		   std::cout << a.get_value(i, j) << " ";
	   std::cout << std::endl;
	}
    

	std::cout << "Matrix B: " << std::endl;
	for (int i = 0; i < 3; i++) {
       std::cout << "      ";
	   for (int j = 0; j < 3; j++)
		   std::cout << b.get_value(i, j)  << " ";
	   std::cout << std::endl;
	}
	
	mui::linalg::sparse_matrix<int,double> c = a + b;
    std::cout << "Addition of matrices (A+B): " << std::endl;
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
	   for (int j = 0; j < 3; j++)
		   std::cout << c.get_value(i, j)  << " ";
	   std::cout << std::endl;
	}

	c.set_zero();
	c = a - b;
	std::cout << "Subtraction of matrices (A-B): " << std::endl;
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
	   for (int j = 0; j < 3; j++)
		   std::cout << c.get_value(i, j)  << " ";
	   std::cout << std::endl;
	}

	c.set_zero();
	c = a * b;
	std::cout << "Multiplication of matrices (A*B): " << std::endl;
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
	   for (int j = 0; j < 3; j++)
		   std::cout << c.get_value(i, j) << " ";
	   std::cout << std::endl;
	}

	c.set_zero();
	c = a.transpose();
	std::cout << "Transpose of A matrix (A^T): " << std::endl;
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
	   for (int j = 0; j < 3; j++)
		   std::cout << c.get_value(i, j) << " ";
	   std::cout << std::endl;
	}

    return 0;
   }
