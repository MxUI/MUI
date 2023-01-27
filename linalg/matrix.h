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
 * @brief Base class for sparse matrix includes basic arithmetic operations 
 * such as addition, subtraction, and multiplication.
 */

#ifndef MUI_SPARSE_MATRIX_H_
#define MUI_SPARSE_MATRIX_H_

#include <iostream>
#include <map>
#include <vector>

namespace mui {
namespace linalg {

template<typename ITYPE, typename VTYPE>
class sparse_matrix {
    private:
        std::map<std::pair<ITYPE, ITYPE>, VTYPE> matrix;
        ITYPE rows;
        ITYPE cols;

    public:
        // Constructor
        sparse_matrix<ITYPE,VTYPE>(ITYPE r, ITYPE c)
			: rows(r), cols(c) {}

        // Function to insert an element
        void set_value(ITYPE r, ITYPE c, VTYPE val) {
        	if (val != 0)
        		matrix[std::make_pair(r, c)] = val;
        }

        // Function to get the value at a given position
        VTYPE get_value(ITYPE r, ITYPE c) {
            auto it = matrix.find(std::make_pair(r, c));
            if (it != matrix.end()) {
                return it->second;
            } else {
                return 0;
            }
        }

        // Function to set all elements to zero and empty the sparse matrix
        void set_zero() {
        	matrix.clear();
        }

        // Function to get transpose of matrix
        sparse_matrix<ITYPE,VTYPE> transpose() {
        	sparse_matrix<ITYPE,VTYPE> res(cols, rows);
            for (auto elememt : matrix)
            	res.set_value(elememt.first.second, elememt.first.first, elememt.second);
            return res;
        }

        // Overload addition operator to perform sparse matrix addition
        sparse_matrix<ITYPE,VTYPE> operator+(sparse_matrix<ITYPE,VTYPE> &addend) {

			if (rows != addend.rows || cols != addend.cols) {
				std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix addition" << std::endl;
				std::abort();
			}

            sparse_matrix<ITYPE,VTYPE> res(rows, cols);
            for (auto elememt : matrix) {
                res.set_value(elememt.first.first, elememt.first.second, elememt.second);
            }
            for (auto elememt : addend.matrix) {
                res.set_value(elememt.first.first, elememt.first.second, res.get_value(elememt.first.first, elememt.first.second) + elememt.second);
            }
            return res;
        }

        // Overload subtraction operator to perform sparse matrix subtraction
        sparse_matrix<ITYPE,VTYPE> operator-(sparse_matrix<ITYPE,VTYPE> &subtrahend) {
            if (rows != subtrahend.rows || cols != subtrahend.cols) {
				std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix subtraction" << std::endl;
				std::abort();
            }
            sparse_matrix<ITYPE,VTYPE> res(rows, cols);
            for (auto elememt : matrix) {
                res.set_value(elememt.first.first, elememt.first.second, elememt.second);
            }
            for (auto elememt : subtrahend.matrix) {
                res.set_value(elememt.first.first, elememt.first.second, res.get_value(elememt.first.first, elememt.first.second) - elememt.second);
            }
            return res;
        }

        // Overload multiplication operator to perform sparse matrix multiplication
        sparse_matrix<ITYPE,VTYPE> operator*(sparse_matrix<ITYPE,VTYPE> &multiplicand) {
            if (cols != multiplicand.rows) {
            	std::cerr << "MUI Error [matrix.h]: matrix size mismatch during matrix multiplication" << std::endl;
            	std::abort();
            }
            sparse_matrix<ITYPE,VTYPE> res(rows, multiplicand.cols);
            for (auto elememt1 : matrix) {
                for (auto elememt2 : multiplicand.matrix) {
                    if (elememt1.first.second == elememt2.first.first) {
                        res.set_value(elememt1.first.first, elememt2.first.second, res.get_value(elememt1.first.first, elememt2.first.second) + elememt1.second * elememt2.second);
                    }
                }
            }
            return res;
        }

};

} // linalg
} // mui

#endif /* MUI_SPARSE_MATRIX_H_ */
