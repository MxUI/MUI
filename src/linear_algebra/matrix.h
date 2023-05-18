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
 * @file matrix.h
 * @author W. Liu
 * @date 27 January 2023
 * @brief Base class for sparse matrix based on COO format includes basic
 * arithmetic operations such as addition, subtraction, and multiplication.
 */

#ifndef MUI_SPARSE_MATRIX_H_
#define MUI_SPARSE_MATRIX_H_

#include <map>
#include <vector>
#include "linalg_util.H"

namespace mui {
namespace linalg {

// Class of sparse matrix
template<typename ITYPE, typename VTYPE>
class sparse_matrix {

    public:
        // *****************************************
        // ****** Constructors & Destructor ********
        // *****************************************

        // Constructor - takes in size of row and column to generate an empty matrix
        sparse_matrix<ITYPE,VTYPE>(ITYPE, ITYPE, const std::string & = "CSR");
        // Constructor - null matrix
        sparse_matrix<ITYPE,VTYPE>(const std::string & = "CSR");
        // Constructor - takes in another sparse_matrix object as an argument
        sparse_matrix<ITYPE,VTYPE>(const sparse_matrix<ITYPE,VTYPE> &);
        // Constructor - takes in a std::vector with row major dense matrix format as an argument
        sparse_matrix<ITYPE,VTYPE>(const std::vector<std::vector<VTYPE>> &, const std::string & = "CSR");
        // Constructor - generate various square matrices
        sparse_matrix<ITYPE,VTYPE>(ITYPE, const std::string & = {}, const std::string & = "CSR");
        // Destructor
        ~sparse_matrix<ITYPE,VTYPE>();

        // *****************************************
        // ********** Matrix I/O & info ************
        // *****************************************

        // Member function to print matrix elements to the console
        void print() const;
        // Member function to get the value at a given position
        VTYPE get_value(ITYPE, ITYPE) const;
        // Member function to get the number of rows
        ITYPE get_rows() const;
        // Member function to get the number of cols
        ITYPE get_cols() const;
        // Member function to get non-zero elements
        std::vector<std::pair<ITYPE, ITYPE>> get_non_zero_elements() const;
        // Member function to get number of non-zero elements
        ITYPE non_zero_elements_count() const;
        // Member function to check whether the matrix contains all zero elements
        bool empty() const;
        // Member function to get the format of the matrix
        std::string getFormat() const;

        // *****************************************
        // ********* Matrix manipulations **********
        // *****************************************

        // Member function to resize a null matrix
        void resize_null(ITYPE, ITYPE);
        // Member function to resize an all-zero matrix
        void resize(ITYPE, ITYPE);
        // Member function to copy a sparse_matrix
        void copy(const sparse_matrix<ITYPE,VTYPE> &);
        // Member function to get a segment of a sparse_matrix
        sparse_matrix<ITYPE,VTYPE> segment(ITYPE, ITYPE, ITYPE, ITYPE);
        // Member function to insert an element
        void set_value(ITYPE, ITYPE, VTYPE);
        // Member function to insert the same value to all elements
        void set_value(VTYPE);
        // Member function to swap two elements in a sparse matrix
        void swap_elements(ITYPE, ITYPE, ITYPE, ITYPE);
        // Member function to set all elements to zero and empty the sparse matrix
        void set_zero();
        // Member function to add scalar to a specific elements
        void add_scalar(ITYPE, ITYPE, VTYPE);
        // Member function to subtract a scalar from a specific elements
        void subtract_scalar(ITYPE, ITYPE, VTYPE);
        // Overloaded assignment operator
        sparse_matrix<ITYPE,VTYPE>& operator=(const sparse_matrix<ITYPE,VTYPE> &);
        // Member function to convert the format of the sparse matrix
        void format_conversion(const std::string & = "COO", bool = false, bool = false, const std::string & = "overwrite");

        // *****************************************
        // ********* Arithmetic operations *********
        // *****************************************

        // Overload addition operator to perform sparse matrix addition
        sparse_matrix<ITYPE,VTYPE> operator+(sparse_matrix<ITYPE,VTYPE> &) const;
        // Overload subtraction operator to perform sparse matrix subtraction
        sparse_matrix<ITYPE,VTYPE> operator-(sparse_matrix<ITYPE,VTYPE> &) const;
        // Overload multiplication operator to perform sparse matrix multiplication
        sparse_matrix<ITYPE,VTYPE> operator*(sparse_matrix<ITYPE,VTYPE> &) const;
        // Overload multiplication operator to perform scalar multiplication
        template <typename STYPE>
        sparse_matrix<ITYPE,VTYPE> operator*(const STYPE &) const;
        // Member function of dot product
        VTYPE dot_product(sparse_matrix<ITYPE,VTYPE> &) const;
        // Member function of Hadamard product
        sparse_matrix<ITYPE,VTYPE> hadamard_product(const sparse_matrix<ITYPE,VTYPE> &) const;
        // Member function to get transpose of matrix
        sparse_matrix<ITYPE,VTYPE> transpose() const;
        // Member function to perform LU decomposition
        void lu_decomposition(sparse_matrix<ITYPE,VTYPE> &, sparse_matrix<ITYPE,VTYPE> &) const;
        // Member function to perform QR decomposition
        void qr_decomposition(sparse_matrix<ITYPE,VTYPE> &, sparse_matrix<ITYPE,VTYPE> &) const;
        // Member function to get the inverse of matrix
        sparse_matrix<ITYPE,VTYPE> inverse() const;

        // *****************************************
        // **************** Asserts ****************
        // *****************************************

        // Member function to assert the matrix format
        //void assertValidFormat() const;

    protected:

        // *****************************************
        // ****** Constructors & Destructor ********
        // *****************************************

        // Protected member function to set matrix format - helper function on matrix constructors
        void set_matrix_format(const std::string & = "CSR");

        // *****************************************
        // ********* Matrix manipulations **********
        // *****************************************

        // Protected member function to sort the entries by row and column for sparse matrix with COO format
        void sort_coo(bool = true, bool = false, const std::string & = "overwrite");
        // Protected member function to convert COO matrix into CSR matrix
        void coo_to_csr();
        // Protected member function to convert COO matrix into CSC matrix
        void coo_to_csc();
        // Protected member function to convert CSR matrix into COO matrix
        void csr_to_coo();
        // Protected member function to convert CSR matrix into CSC matrix
        void csr_to_csc();
        // Protected member function to convert CSC matrix into COO matrix
        void csc_to_coo();
        // Protected member function to convert CSC matrix into CSR matrix
        void csc_to_csr();

	private:
		// Format of sparse matrix
		enum class format {
			COO,
			CSR,
			CSC
		};

		// COO format data
		struct matrix_coo {
			// Values of non-zero elements of sparse matrix
			std::vector<VTYPE> values_;
			// Row index of each element in the values_ vector
			std::vector<ITYPE> row_indices_;
			// Column index of each element in the values_ vector
			std::vector<ITYPE> col_indices_;
		}

		// CSR format data
		struct matrix_csr {
			// Values of non-zero elements of sparse matrix
			std::vector<VTYPE> values_;
			// Row pointers of each element in the values_ vector
			std::vector<ITYPE> row_ptrs_;
			// Column index of each element in the values_ vector
			std::vector<ITYPE> col_indices_;
		}

		// CSC format data
		struct matrix_csc {
			// Values of non-zero elements of sparse matrix
			std::vector<VTYPE> values_;
			// Row index of each element in the values_ vector
			std::vector<ITYPE> row_indices_;
			// Column pointers of each element in the values_ vector
			std::vector<ITYPE> col_ptrs_;
		}

		// Number of rows of sparse matrix
		ITYPE rows_ = 0;
		// Number of columns of sparse matrix
		ITYPE cols_ = 0;
		// Number of non-zero elements of sparse matrix
		ITYPE nnz_ = 0;
		// Format indicator with default value of format::CSR
		format matrix_format_ = format::CSR;

		// Dummy member variable for invalid or unassigned elements in sparse matrix
		VTYPE dummy_ = 0;

};

} // linalg
} // mui

// Include implementations
// #include "../linear_algebra/matrix_asserts.h"
#include "../linear_algebra/matrix_ctor_dtor.h"
#include "../linear_algebra/matrix_arithmetic.h"
#include "../linear_algebra/matrix_manipulation.h"
#include "../linear_algebra/matrix_io_info.h"

#endif /* MUI_SPARSE_MATRIX_H_ */
