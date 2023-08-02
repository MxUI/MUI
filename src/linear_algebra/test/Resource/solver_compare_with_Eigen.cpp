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
 * @file solver_compare_with_Eigen.cpp
 * @author W. Liu
 * @date 01 April 2023
 * @brief Conjugate Gradient solver with preconditioners compared with Eigen.
 */

#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include <ctime>
#include "../../solver.h"

int runType = -1;

inline double test00 () {

	std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==== TEST 00: Find Matrix H_i by Matrix Css and Aas ========" << std::endl;
    std::cout << "================== Simplified RBF Matrices =================" << std::endl;
    std::cout << "======================= mui::linalg ========================" << std::endl;
    std::cout << "================ Diagonal Preconditioner ===================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    const clock_t begin_time = clock();

    double cg_solve_tol = 1e-6;
    int cg_max_iter = 100;

    std::pair<int, double> cgReturn;

    // Reads Css matrix from a file
    mui::linalg::sparse_matrix<int,double> Css;
    std::ifstream ifileCss("Css.csv");
    ifileCss >> Css;
    ifileCss.close();

    // Reads Aas matrix from a file
    mui::linalg::sparse_matrix<int,double> Aas;
    std::ifstream ifileAas("Aas.csv");
    ifileAas >> Aas;
    ifileAas.close();

    // Reads H_ref matrix from a file
    mui::linalg::sparse_matrix<int,double> H_ref; //< Reference value of Transformation Matrix
    std::ifstream ifileHRef("Hi.csv");
    ifileHRef >> H_ref;
    ifileHRef.close();

    mui::linalg::sparse_matrix<int,double> H_i;   //< Coupling Matrix
    mui::linalg::sparse_matrix<int,double> H_diff; //< Difference between reference value and calculated value of Transformation Matrix

    H_i.resize(Aas.get_rows(), Aas.get_cols());
    H_diff.resize(Aas.get_rows(), Aas.get_cols());

    mui::linalg::diagonal_preconditioner<int,double> M(Css);

    mui::linalg::conjugate_gradient<int,double> cg(Css, Aas, cg_solve_tol, cg_max_iter, &M);

    cgReturn = cg.solve();
    H_i = cg.getSolution();

    std::cout << "Matrix H_i: " << std::endl;
    H_i.print();

    std::cout << "Reference value of Matrix H_i: " << std::endl;
    H_ref.print();

    H_diff = H_i - H_ref;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    H_diff.print();

    double mui_time = static_cast<double>(clock() - begin_time) / CLOCKS_PER_SEC;

    std::cout << "MUI CG iteration number: " << cgReturn.first <<" with final MUI r_norm_rel: " << cgReturn.second << " solved in " << mui_time << "s " << std::endl;

    std::cout << std::endl;

    return mui_time;
}

inline double test01 () {

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "==== TEST 01: Find Matrix H_i by Matrix Css and Aas ========" << std::endl;
    std::cout << "================== Simplified RBF Matrices =================" << std::endl;
    std::cout << "========================= Eigen ============================" << std::endl;
    std::cout << "================ Diagonal Preconditioner ===================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;

    const clock_t begin_time = clock();

    double cg_solve_tol = 1e-6;
    int cg_max_iter = 100;

    Eigen::SparseMatrix<double> Css;
    Eigen::SparseMatrix<double> Aas;
    Eigen::SparseMatrix<double> H_ref; //< Difference between reference value and calculated value of Transformation Matrix

    std::pair<int, double> cgReturn;

    // Reads Css matrix from a file
    mui::linalg::sparse_matrix<int,double> Css_raw;
    std::ifstream ifileCss("Css.csv");
    ifileCss >> Css_raw;
    ifileCss.close();

    // Reads Aas matrix from a file
    mui::linalg::sparse_matrix<int,double> Aas_raw;
    std::ifstream ifileAas("Aas.csv");
    ifileAas >> Aas_raw;
    ifileAas.close();

    // Reads H_ref matrix from a file
    mui::linalg::sparse_matrix<int,double> H_ref_raw; //< Reference value of Transformation Matrix
    std::ifstream ifileHRef("Hi.csv");
    ifileHRef >> H_ref_raw;
    ifileHRef.close();

	Css.resize(Css_raw.get_rows(), Css_raw.get_cols());
	Aas.resize(Aas_raw.get_rows(), Aas_raw.get_cols());
	H_ref.resize(H_ref_raw.get_rows(), H_ref_raw.get_cols());

	//set Css
    std::vector<Eigen::Triplet<double> > coefsC;

	for ( size_t i = 0; i < Css_raw.get_rows(); i++ ) {
		for ( size_t j = 0; j < Css_raw.get_cols(); j++ ) {
			coefsC.emplace_back(Eigen::Triplet<double> (i, j, Css_raw.get_value(i, j)));
		}
	}

	Css.reserve(coefsC.size());
	Css.setFromTriplets(coefsC.begin(), coefsC.end());

	//set Aas
	std::vector<Eigen::Triplet<double>> coefs;

	for ( size_t i = 0; i < Aas_raw.get_rows(); i++ ) {
		for ( size_t j = 0; j < Aas_raw.get_cols(); j++ ) {
			coefs.emplace_back(Eigen::Triplet<double> (i, j, Aas_raw.get_value(i, j)));
		}
	}

	Aas.reserve(coefs.size());
	Aas.setFromTriplets(coefs.begin(), coefs.end());

	//set H_ref
    std::vector<Eigen::Triplet<double> > coefsRef;

	for ( size_t i = 0; i < H_ref_raw.get_rows(); i++ ) {
		for ( size_t j = 0; j < H_ref_raw.get_cols(); j++ ) {
			coefsRef.emplace_back(Eigen::Triplet<double> (i, j, H_ref_raw.get_value(i, j)));
		}
	}

	H_ref.reserve(coefsRef.size());
	H_ref.setFromTriplets(coefsRef.begin(), coefsRef.end());

	//Eigen solve
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
			Eigen::Lower | Eigen::Upper,
			Eigen::DiagonalPreconditioner<double>> solver(Css);

		solver.setMaxIterations(cg_max_iter);
		solver.setTolerance(cg_solve_tol);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_i = solver.solve(Aas);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_diff = H_i - H_ref;

    std::cout << "Matrix H_i: " << std::endl;
    std::cout << H_i << std::endl;

    std::cout << "Reference value of Matrix H_i: " << std::endl;
    std::cout << H_ref << std::endl;

    std::cout << "Difference between calculated value and reference value: " << std::endl;
    std::cout << H_diff << std::endl;

    double eigen_time = static_cast<double>(clock() - begin_time) / CLOCKS_PER_SEC;

    std::cout << "Eigen CG iteration number: " << solver.iterations() <<" with final Eigen r_norm_rel: " << solver.error() << " solved in " << eigen_time << "s " << std::endl;

    std::cout << std::endl;

    return eigen_time;

}

int main(int argc, char** argv) {

	//Parse input parameters
	if( argc < 2 ) {  //Check the number of parameters
		std::cerr << "Usage: " << argv[0] << " [linear algebra solution (eigen/mui/both)]" << std::endl;
		exit(-1);
	}

	std::string runTypeArg(argv[1]);

	if ( runTypeArg.compare("mui") == 0 )
		runType = 0;
	else if ( runTypeArg.compare("eigen") == 0 )
		runType = 1;
	else if ( runTypeArg.compare("both") == 0 )
		runType = 2;
	else {
		std::cerr << "Usage: " << argv[0] << " [linear algebra solution (eigen/mui/both)]" << std::endl;
		exit(-1);
	}

	double mui_time = 0;
	double eigen_time = 0;

	if( runType == 0 ) { // mui::linalg only
		double mui_time = test00();

		std::cout << "MUI::linalg time: " << mui_time << std::endl;
	}
    if( runType == 1 ) { // eigen only
    	double eigen_time = test01();

    	std::cout << "Eigen time: " << eigen_time << std::endl;
    }
    if( runType == 2 ) { // both tests
    	double mui_time = test00();
    	double eigen_time = test01();

    	std::cout << "MUI::linalg time: " << mui_time << " Eigen time: " << eigen_time << std::endl;
		std::cout << "Eigen is " <<  mui_time/eigen_time << " times faster than MUI::linalg" << std::endl;
    }

    return 0;
}
