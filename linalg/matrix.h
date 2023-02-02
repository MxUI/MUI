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

namespace mui {
namespace linalg {

template<typename ITYPE, typename VTYPE>
class sparse_matrix {

    public:
        // *****************************************
        // ****** Constructors & Destructor ********
        // *****************************************

        // Constructor - takes in size of row and column to generate an empty matrix
        sparse_matrix<ITYPE,VTYPE>(ITYPE, ITYPE);
        // Constructor - null matrix
        sparse_matrix<ITYPE,VTYPE>();
        // Constructor - takes in another sparse_matrix object as an argument
        sparse_matrix<ITYPE,VTYPE>(const sparse_matrix<ITYPE,VTYPE> &);
        // Constructor - generate various square matrices
        sparse_matrix<ITYPE,VTYPE>(ITYPE, const std::string & = {});
        // Destructor
        ~sparse_matrix<ITYPE,VTYPE>();

        // *****************************************
        // ********** Matrix I/O & info ************
        // *****************************************

        // Member function to print matrix elements to the console
        void print();
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
        bool empty();

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
        // Member function to set all elements to zero and empty the sparse matrix
        void set_zero();
        // Member function to add scalar to a specific elements
        void add_scalar(ITYPE, ITYPE, VTYPE);
        // Member function to subtract a scalar from a specific elements
        void subtract_scalar(ITYPE, ITYPE, VTYPE);

        // *****************************************
        // ********* Arithmetic operations *********
        // *****************************************

        // Overload addition operator to perform sparse matrix addition
        sparse_matrix<ITYPE,VTYPE> operator+(sparse_matrix<ITYPE,VTYPE> &);
        // Overload subtraction operator to perform sparse matrix subtraction
        sparse_matrix<ITYPE,VTYPE> operator-(sparse_matrix<ITYPE,VTYPE> &);
        // Overload multiplication operator to perform sparse matrix multiplication
        sparse_matrix<ITYPE,VTYPE> operator*(sparse_matrix<ITYPE,VTYPE> &);
        // Overload multiplication operator to perform scalar multiplication
        template <typename STYPE>
        sparse_matrix<ITYPE,VTYPE> operator*(const STYPE &) const;
        // Overloaded assignment operator
        sparse_matrix<ITYPE,VTYPE>& operator=(const sparse_matrix<ITYPE,VTYPE> &);
        // Member function of dot product
        VTYPE dot_product(sparse_matrix<ITYPE,VTYPE> &) const;
        // Member function of Hadamard product
        sparse_matrix<ITYPE,VTYPE> hadamard_product(const sparse_matrix<ITYPE,VTYPE> &);
        // Member function to get transpose of matrix
        sparse_matrix<ITYPE,VTYPE> transpose();
        // Member function to get the inverse of matrix
        sparse_matrix<ITYPE,VTYPE> inverse();

    private:
        // Non-zero sparse matrix elements in COO format
        std::map<std::pair<ITYPE, ITYPE>, VTYPE> matrix_;
        // Number of rows of sparse matrix
        ITYPE rows_;
        // Number of columns of sparse matrix
        ITYPE cols_;
        // Dummy member variable for invalid or unassigned elements in sparse matrix
        VTYPE dummy_;
};

} // linalg
} // mui

// Include implementations
#include "matrix_ctor_dtor.h"
#include "matrix_arithmetic.h"
#include "matrix_manipulation.h"
#include "matrix_io_info.h"

#endif /* MUI_SPARSE_MATRIX_H_ */
