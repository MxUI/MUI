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
 * @file preconditioner.h
 * @author W. Liu
 * @date 28 January 2023
 * @brief Preconditioner classes.
 */

#ifndef MUI_PRECONDITIONER_H_
#define MUI_PRECONDITIONER_H_

namespace mui {
namespace linalg {

// Base preconditioner class
template<typename ITYPE, typename VTYPE>
class preconditioner {

    public:
        // Abstract function on preconditioner apply
        virtual sparse_matrix<ITYPE,VTYPE> apply(const sparse_matrix<ITYPE,VTYPE> &) = 0;
        // Destructor
        virtual ~preconditioner(){};
};

// Class of Incomplete LU preconditioner
template<typename ITYPE, typename VTYPE>
class incomplete_lu_preconditioner : public preconditioner<ITYPE,VTYPE> {

    public:
        // Constructor
        incomplete_lu_preconditioner(const sparse_matrix<ITYPE,VTYPE>&);
        // Destructor
        ~incomplete_lu_preconditioner();
        // Member function on preconditioner apply
        sparse_matrix<ITYPE,VTYPE> apply(const sparse_matrix<ITYPE,VTYPE>&);

    private:
        // Lower triangular matrix for Incomplete LU preconditioner
        sparse_matrix<ITYPE,VTYPE> L_;
        // Upper triangular matrix for Incomplete LU preconditioner
        sparse_matrix<ITYPE,VTYPE> U_;

};

// Class of Incomplete Cholesky preconditioner
template<typename ITYPE, typename VTYPE>
class incomplete_cholesky_preconditioner : public preconditioner<ITYPE,VTYPE> {

    public:
        // Constructor
        incomplete_cholesky_preconditioner(const sparse_matrix<ITYPE,VTYPE>&);
        // Destructor
        ~incomplete_cholesky_preconditioner();
        // Member function on preconditioner apply
        sparse_matrix<ITYPE,VTYPE> apply(const sparse_matrix<ITYPE,VTYPE>&);

    private:
        // Lower triangular matrix for Incomplete Cholesky preconditioner
        sparse_matrix<ITYPE,VTYPE> L_;
};

// Class of Symmetric Successive Over-relaxation preconditioner
template<typename ITYPE, typename VTYPE>
class symmetric_successive_over_relaxation_preconditioner : public preconditioner<ITYPE,VTYPE> {

    public:
        // Constructor
        symmetric_successive_over_relaxation_preconditioner(const sparse_matrix<ITYPE,VTYPE>&, VTYPE = 1.0);
        // Destructor
        ~symmetric_successive_over_relaxation_preconditioner();
        // Member function on preconditioner apply
        sparse_matrix<ITYPE,VTYPE> apply(const sparse_matrix<ITYPE,VTYPE>&);

    private:
        // The coefficient matrix of the matrix equation
        sparse_matrix<ITYPE,VTYPE> A_;
        // The relaxation parameter
        VTYPE omega_;
};

// Class of diagonal (Jacobi) preconditioner
template<typename ITYPE, typename VTYPE>
class diagonal_preconditioner : public preconditioner<ITYPE,VTYPE> {

    public:
        // Constructor
        diagonal_preconditioner(const sparse_matrix<ITYPE,VTYPE>&);
        // Destructor
        ~diagonal_preconditioner();
        // Member function on preconditioner apply
        sparse_matrix<ITYPE,VTYPE> apply(const sparse_matrix<ITYPE,VTYPE>&);

    private:
        // The inverse diagonal matrix
        sparse_matrix<ITYPE,VTYPE> inv_diag_;

};

} // linalg
} // mui

// Include implementations
#include "../linear_algebra/preconditioner_ilu.h"
#include "../linear_algebra/preconditioner_ic.h"
#include "../linear_algebra/preconditioner_ssor.h"
#include "../linear_algebra/preconditioner_diagonal.h"

#endif /* MUI_PRECONDITIONER_H_ */
