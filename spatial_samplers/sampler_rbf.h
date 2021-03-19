/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
*                    A. Skillen                                              *
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
 * @file sampler_rbf.h
 * @author A. Skillen, W. Liu
 * @date November 2018
 * @brief Spatial sampler using Gaussian Radial Basis Function interpolation.
 */

#ifndef MUI_SAMPLER_RBF_H_
#define MUI_SAMPLER_RBF_H_

#include "../util.h"
#include <Eigen/Sparse>
#include "../uniface.h"
#include <iterator>
#include <ctime>

namespace mui {

template<typename CONFIG=default_config, typename O_TP=typename CONFIG::REAL, typename I_TP=O_TP>
class sampler_rbf 
{
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;
	using EXCEPTION  = typename CONFIG::EXCEPTION;

	static const bool QUIET = CONFIG::QUIET;

/**
 * Input parameters:
 * 1.REAL r:
 *    r is the searching radius
 * 2. std::vector<point_type>& pts:
 *    pts is vector of points that pre-set for RBF interpolation
 * 3. INT basisFunc:
 *    basisFunc is a parameter for basis function selection. Implemented functions are as follows:
 *    Gaussian(default): basisFunc_=0
 *    Wendland's C0:     basisFunc_=1
 *    Wendland's C2:     basisFunc_=2
 *    Wendland's C4:     basisFunc_=3
 *    Wendland's C6:     basisFunc_=4
 * 4. bool conservative:
 *    conservative is a boolen switch on the mode of RBF interpolation:
 *    consistent mode(default): conservative=false
 *    conservative mode:        conservative=true
 * 5. bool polynomial:
 *    polynomial is a boolen switch on the polynomial term of the transformation matrix:
 *    without polynomial terms(default): polynomial=false
 *    with polynomial terms:             polynomial=true
 * 6. bool smoothFunc:
 *    smoothFunc is a boolen switch on the smoothing function of the transformation matrix:
 *    without smoothing function(default): smoothFunc=false
 *    with smoothing function:             smoothFunc=true
 * 7. bool readMatrix:
 *    readMatrix is a parameter to select whether to construct transformation matrix or read matrix from file:
 *    construct transformation matrix and write it to file(default): readMatrix=false
 *    read transformation matrix from file:                          readMatrix=true
 * 8. const std::string& fileAddress:
 *    fileAddress is a parameter to set the address that the transformation matrix I/O.
 *    The default value of fileAddress is an empty string.
 * 9. double cutoff:
 *    cutoff is a parameter to set the cut-off of the Gaussian basis function (only valid for basisFunc_=0).
 *    The default value of cutoff is 1e-9
 */

	sampler_rbf( REAL r, std::vector<point_type>& pts, INT basisFunc=0, bool conservative=false, bool polynomial=false, bool smoothFunc=false, bool readMatrix=false, const std::string& fileAddress=std::string(), REAL cutoff=1e-9 ) :
	r_(r), 
	initialised_(false),
	CABrow_(0),
	CABcol_(0),
	Hrow_(0),
	Hcol_(0),
	conservative_(conservative), 
	consistent_(!conservative),
	polynomial_(polynomial),
	smoothFunc_(smoothFunc),
	readMatrix_(readMatrix),
	fileAddress_(fileAddress),
	N_sp_(50),
	M_ap_(50),
	basisFunc_(basisFunc),
	pts_(pts)
	{
		//set s to give rbf(r)=cutoff (default 1e-9)
		s_=std::pow(-std::log(cutoff),0.5) / r_;
	}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) {
		if( !initialised_ )
			computeRBFtransformation( data_points );

		auto p = std::find_if(pts_.begin(),pts_.end(),[focus](point_type b) { return normsq(focus - b) < 1e-10; });

		if( p==std::end(pts_) ) {
			EXCEPTION( std::runtime_error("Point not found. Must pre-set points for RBF interpolation" ));
		}

		auto i = std::distance(pts_.begin(),p);

		REAL sum = 0.0;
		for( INT j=0 ; j < data_points.size() ; j++) {
			sum += H_(i, j) * data_points[j].second;
		}

		return sum;
	}

	inline geometry::any_shape<CONFIG> support( point_type focus, REAL domain_mag ) const {
		return geometry::any_shape<CONFIG>();
	}

	void preSetFetchPoints( std::vector<point_type>& pts ) {
		pts_ = pts;
		initialised_ = false;
	}
	void addFetchPoint( point_type pt ) {
		pts_.push_back( pt );
		initialised_ = false;
	}


private:
	template<template<typename,typename> class CONTAINER>
	void computeRBFtransformation( const CONTAINER<ITYPE,CONFIG> &data_points ) {

		if(conservative_){
			if(pts_.size()<=50){
				N_sp_=pts_.size();
			}
		}else{
			if(data_points.size()<=50){
				N_sp_=data_points.size();
			}
		}
		if(smoothFunc_){
			if(pts_.size()<=50){
				M_ap_=(pts_.size()-1);
			}
		}

        if (readMatrix_) {

            std::ifstream inputFileMatrixSize(fileAddress_+"/matrixSize.dat");

            if (!inputFileMatrixSize) {
                std::cerr << "Could not locate the file address of matrixSize.dat!" << std::endl;
            } else{
				std::string tempS;
				std::vector<INT> tempV;
				while (std::getline(inputFileMatrixSize, tempS)) {
					// Skips the line if the first two characters are '//'
					if (tempS[0] == '/' && tempS[1] == '/') continue;
					std::stringstream lineStream(tempS);
                    std::string tempSS;
					while(std::getline(lineStream, tempSS, ',')) {
						tempV.push_back(std::stoi(tempSS));
					}
				}
				CABrow_ = tempV[0];
				CABcol_ = tempV[1];
				CAArow_ = tempV[2];
				CAAcol_ = tempV[3];
				Hrow_ = tempV[4];
				Hcol_ = tempV[5];
            }

            std::ifstream inputFileCAB(fileAddress_+"/connectivityAB.dat");

            if (!inputFileCAB) {
                std::cerr << "Could not locate the file address on the connectivityAB.dat!" << std::endl;
            } else{

                connectivityAB_.resize(CABrow_);

                for(int i=0; i<CABrow_; ++i){

                    connectivityAB_[i].resize(CABcol_);
                    std::string tempS;
                    while (std::getline(inputFileCAB, tempS)) {
						// Skips the line if the first two characters are '//'
						if (tempS[0] == '/' && tempS[1] == '/') continue;
						std::stringstream lineStream(tempS);
						std::string tempSS;
                        std::vector<INT> tempV;
                        while(std::getline(lineStream, tempSS, ',')) {
                            tempV.push_back(std::stoi(tempSS));
                        }
                        connectivityAB_.push_back(tempV);
                    }
                }
            }
            if (smoothFunc_) {
				std::ifstream inputFileCAA(fileAddress_+"/connectivityAA.dat");

				if (!inputFileCAA) {
					std::cerr << "Could not locate the file address on the connectivityAA.dat!" << std::endl;
				} else{

					if ((CAArow_== 0)||(CAAcol_== 0)) {
						std::cerr << "Error on the size of connectivityAA matrix in matrixSize.dat! number of rows: " <<
						CAArow_ << "number of columns: " << CAAcol_ <<
						". Please make sure the maxtixes were generated with the smoothing function switched on"<< std::endl;
					} else{

						connectivityAA_.resize(CAArow_);

						for(int i=0; i<CAArow_; ++i){

							connectivityAA_[i].resize(CAAcol_);
							std::string tempS;
							while (std::getline(inputFileCAA, tempS)) {
								// Skips the line if the first two characters are '//'
								if (tempS[0] == '/' && tempS[1] == '/') continue;
								std::stringstream lineStream(tempS);
								std::string tempSS;
								std::vector<INT> tempV;
								while(std::getline(lineStream, tempSS, ',')) {
									tempV.push_back(std::stoi(tempSS));
								}
								connectivityAA_.push_back(tempV);
							}
						}
					}
				}
			}

            H_.resize( Hrow_, Hcol_);

            std::ifstream inputFileHMatrix(fileAddress_+"/Hmatrix.dat");

            if (!inputFileHMatrix) {
                std::cerr << "Could not locate the file address on the Hmatrix.dat!" << std::endl;
            } else{
                std::string tempS;
                //std::vector<double> tempV;
                int tempRow = 0;
                int tempPoints = 0;
                while (std::getline(inputFileHMatrix, tempS)){
					// Skips the line if the first two characters are '//'
					if (tempS[0] == '/' && tempS[1] == '/') continue;
                    std::stringstream lineStream(tempS);
                    std::string tempSS;
                    int tempCol = 0;
                    while (std::getline(lineStream, tempSS, ',')) {
                        //tempV.push_back(std::stod(tempSS));
                        H_(tempRow,tempCol) = std::stod(tempSS);
                        ++tempCol;
                        ++tempPoints;
                    }
                    ++tempRow;
                }

                if ((tempRow!=Hrow_) || ((tempPoints/tempRow)!=Hcol_)) {
                    std::cerr << "tempRow (" << tempRow << ") does NOT equals to Hrow_ (" << Hrow_ << "), or" << std::endl
                        <<"(tempPoints/tempRow) (" << (tempPoints/tempRow) << ") does NOT equals to Hcol_ (" << Hcol_ << ")"<< std::endl;
                }
                // else{
                    //H_=Eigen::Map< Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(tempV.data(), tempRow, tempV.size()/tempRow);
               // }
            }


        } else {

            buildConnectivity( data_points, N_sp_ );

            H_.resize( pts_.size(), data_points.size() );

			if (smoothFunc_) {
				buildConnectivityAA( M_ap_ );
				H_toSmooth_.resize( pts_.size(), data_points.size() );
			}

            std::ofstream outputFileMatrixSize(fileAddress_+"/matrixSize.dat");

            if (!outputFileMatrixSize) {
                std::cerr << "Could not locate the file address of matrixSize.dat!" << std::endl;
            } else{
                outputFileMatrixSize << "// *********************************************************************************************************************************************";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// **** This is the 'matrixSize.dat' file of the RBF spatial sampler of the MUI library";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// **** This file contains the size (number of rows and number of columns) of the Point Connectivity Matrix (N) and the Coupling Matrix (H).";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// **** The file uses the Comma-Separated Values (CSV) format and the ASCII format with the meanings as follows: ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// ****			The number of rows of the Point Connectivity Matrix (N), ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// ****			The number of columns of the Point Connectivity Matrix (N),";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// ****			The number of rows of the Point Connectivity Matrix (M) (for smoothing), ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// ****			The number of columns of the Point Connectivity Matrix (M) (for smoothing),";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// ****			The number of rows of the Coupling Matrix (H),";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// ****			The number of columns of the Coupling Matrix (H)";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "// *********************************************************************************************************************************************";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "//  ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << connectivityAB_.size();
                outputFileMatrixSize << ",";
                outputFileMatrixSize << connectivityAB_[0].size();
                outputFileMatrixSize << ",";
                if (smoothFunc_) {
					outputFileMatrixSize << connectivityAA_.size();
					outputFileMatrixSize << ",";
					outputFileMatrixSize << connectivityAA_[0].size();
					outputFileMatrixSize << ",";
				}else{
					outputFileMatrixSize << "0";
					outputFileMatrixSize << ",";
					outputFileMatrixSize << "0";
					outputFileMatrixSize << ",";
				}
                outputFileMatrixSize << H_.rows();
                outputFileMatrixSize << ",";
                outputFileMatrixSize << H_.cols();
                outputFileMatrixSize << "\n";
            }
			const clock_t begin_time = clock();
            if( consistent_ ) {

                if (polynomial_){

                    for( int row=0; row<pts_.size(); row++ ) {

                        Eigen::SparseMatrix<REAL> Css; //< Matrix of radial basis function evaluations between prescribed points
                        Eigen::SparseMatrix<REAL> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points

						Css.resize((1+N_sp_+CONFIG::D), (1+N_sp_+CONFIG::D));
						Aas.resize(1, (1+N_sp_+CONFIG::D));

                        std::vector<Eigen::Triplet<REAL> > coefsC;

                        //set Css
                        for( INT i=0 ; i < N_sp_ ; i++) {
                            for( INT j=i ; j < N_sp_ ; j++) {
                                int glob_i = connectivityAB_[row][i];
                                int glob_j = connectivityAB_[row][j];

                                auto d = norm( data_points[glob_i].first - data_points[glob_j].first );

                                if ( d < r_ ) {
                                    REAL w = rbf(d);
                                    coefsC.push_back( Eigen::Triplet<REAL>(i, j, w ) );

                                    if( i!=j ) {
                                        coefsC.push_back( Eigen::Triplet<REAL>( j, i, w ) );
                                    }
                                }
                            }
                        }

                        for( INT i=0; i < N_sp_; i++) {
                            coefsC.push_back( Eigen::Triplet<REAL>(i, N_sp_, 1 ) );
                            coefsC.push_back( Eigen::Triplet<REAL>(N_sp_, i, 1 ) );

                            int glob_i = connectivityAB_[row][i];

                            for (int tempDim=0; tempDim < CONFIG::D; tempDim++){
                                coefsC.push_back(Eigen::Triplet<REAL>(i, (N_sp_+tempDim+1), data_points[glob_i].first[tempDim]));
                                coefsC.push_back(Eigen::Triplet<REAL>((N_sp_+tempDim+1), i, data_points[glob_i].first[tempDim]));
                            }
                        }

                        Css.reserve(coefsC.size());
                        Css.setFromTriplets(coefsC.begin(), coefsC.end());

                        //set Aas
                        std::vector<Eigen::Triplet<REAL> > coefs;

                        for( INT j=0 ; j < N_sp_ ; j++) {
                            int glob_j = connectivityAB_[row][j];

                            auto d = norm( pts_[row] - data_points[glob_j].first );

                            if ( d < r_ ) {
                                coefs.push_back( Eigen::Triplet<REAL>(0, j, rbf(d) ) );
                            }
                        }

                        coefs.push_back( Eigen::Triplet<REAL>(0, N_sp_, 1 ) );

                        for (int tempDim=0; tempDim < CONFIG::D; tempDim++){
                            coefs.push_back(Eigen::Triplet<REAL>(0, (N_sp_+tempDim+1), pts_[row][tempDim]));
                        }

                        Aas.reserve(coefs.size());
                        Aas.setFromTriplets(coefs.begin(), coefs.end());

                        //invert Css
						Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> I(Css.rows(), Css.cols());
						I.setIdentity();
						Eigen::ConjugateGradient<Eigen::SparseMatrix<REAL> > solver(Css);
						solver.setTolerance(1e-6);
						Eigen::SparseMatrix<REAL> invCss = (solver.solve(I)).sparseView(1e8);

						if( CONFIG::DEBUG ) {
							std::cout << "#iterations of invCss:     " << solver.iterations() << ". Error of invCss: "<< solver.error()<< std::endl;
						}

						Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_i = (Aas * invCss).pruned(1e8);

                        if (smoothFunc_) {
							for( INT j=0 ; j < N_sp_ ; j++) {
								int glob_j = connectivityAB_[row][j];
								H_toSmooth_(row, glob_j) = H_i(0,j);
							}
						}else{
							for( INT j=0 ; j < N_sp_ ; j++) {
								int glob_j = connectivityAB_[row][j];
								H_(row, glob_j) = H_i(0,j);
							}
						}
                    }
					if (smoothFunc_) {
						for( int row=0; row<pts_.size(); row++ ) {
							for( INT j=0 ; j < N_sp_ ; j++) {
								int glob_j = connectivityAB_[row][j];
								REAL h_j_sum = 0.;
								REAL f_sum = 0.;

								for( INT k=0 ; k < M_ap_; k++) {
									INT row_k=connectivityAA_[row][k];
									if (row_k==row) {
										std::cerr << "Invalid row_k value: " << row_k << std::endl;
										} else{
											h_j_sum += std::pow(dist_h_i(row, row_k),(-2));
										}
								}

								for( INT k=0 ; k < M_ap_; k++) {
									INT row_k=connectivityAA_[row][k];
									if (row_k==row) {
										std::cerr << "Invalid row_k value: " << row_k << std::endl;
										} else{
											REAL w_i=((std::pow(dist_h_i(row, row_k),(-2)))/(h_j_sum));
											f_sum += w_i*H_toSmooth_(row_k, glob_j);
										}
								}

								H_(row, glob_j)=0.5*(f_sum+H_toSmooth_(row, glob_j));
							}
						}
					}
                } else{ // Without polynomial terms

                    for( int row=0; row<pts_.size(); row++ ) {
                        Eigen::SparseMatrix<REAL> AA; //< Matrix of radial basis function evaluations between prescribed points
                        Eigen::SparseMatrix<REAL> AB; //< Matrix of RBF evaluations between prescribed and interpolation points

                        AA.resize( N_sp_, N_sp_ );
                        AB.resize( 1, N_sp_ );

                        std::vector<Eigen::Triplet<REAL> > coefs;

                        //set AA
                        for( INT i=0 ; i < N_sp_ ; i++) {
                            for( INT j=i ; j < N_sp_ ; j++) {
                                int glob_i = connectivityAB_[row][i];
                                int glob_j = connectivityAB_[row][j];

                                auto d = norm( data_points[glob_i].first - data_points[glob_j].first );
                                if ( d < r_ ) {
                                    REAL w = rbf(d);
                                    coefs.push_back( Eigen::Triplet<REAL>(i, j , w ) );
                                    if( i!=j ) {
                                        coefs.push_back( Eigen::Triplet<REAL>( j, i, w ) );
                                    }
                                }
                            }
                        }

                        AA.reserve(coefs.size());
                        AA.setFromTriplets(coefs.begin(), coefs.end());

                        //set AB
                        coefs.clear();
                        for( INT j=0 ; j < N_sp_ ; j++) {
                            int glob_j = connectivityAB_[row][j];

                            auto d = norm( pts_[row] - data_points[glob_j].first );
                            if ( d < r_ ) {
                                coefs.push_back( Eigen::Triplet<REAL>(0, j, rbf(d) ) );
                            }
                        }

                        AB.reserve(coefs.size());
                        AB.setFromTriplets(coefs.begin(), coefs.end());

                        //invert AA
                        Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> I(AA.rows(), AA.cols());
                        I.setIdentity();
                        Eigen::ConjugateGradient<Eigen::SparseMatrix<REAL> > solver(AA);
                        solver.setTolerance(1e-6);
                        Eigen::SparseMatrix<REAL> invAA = (solver.solve(I)).sparseView(1e8);

                        if( CONFIG::DEBUG ) {
                        	std::cout << "MUI [sampler_rbf.h]: invCss iteration count: " << solver.iterations()
									  << "                          invCss error: " << solver.error() << std::endl;
                        }

                        Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_i = (AB * invAA).pruned(1e8);

						if (smoothFunc_) {
							for( INT j=0 ; j < N_sp_ ; j++) {
								int glob_j = connectivityAB_[row][j];
								H_toSmooth_(row, glob_j) = H_i(0,j);
							}
						}else{
							for( INT j=0 ; j < N_sp_ ; j++) {
								int glob_j = connectivityAB_[row][j];
								H_(row, glob_j) = H_i(0,j);
							}
						}
                    }
                    if (smoothFunc_) {
						for( int row=0; row<pts_.size(); row++ ) {
							for( INT j=0 ; j < N_sp_ ; j++) {
								int glob_j = connectivityAB_[row][j];
								REAL h_j_sum = 0.;
								REAL f_sum = 0.;

								for( INT k=0 ; k < M_ap_; k++) {
									INT row_k=connectivityAA_[row][k];
									if (row_k==row) {
										std::cerr << "Invalid row_k value: " << row_k << std::endl;
										} else{
											h_j_sum += std::pow(dist_h_i(row, row_k),(-2));
										}
								}

								for( INT k=0 ; k < M_ap_; k++) {
									INT row_k=connectivityAA_[row][k];
									if (row_k==row) {
										std::cerr << "Invalid row_k value: " << row_k << std::endl;
										} else{
											REAL w_i=((std::pow(dist_h_i(row, row_k),(-2)))/(h_j_sum));
											f_sum += w_i*H_toSmooth_(row_k, glob_j);
										}
								}

								H_(row, glob_j)=0.5*(f_sum+H_toSmooth_(row, glob_j));
							}
						}
					}
				}
            }else { //conservative

                if (polynomial_){

                    Eigen::SparseMatrix<REAL> Css; //< Matrix of radial basis function evaluations between prescribed points
                    Eigen::SparseMatrix<REAL> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points

                    Css.resize((1+pts_.size()+CONFIG::D), (1+pts_.size()+CONFIG::D));
                    Aas.resize((data_points.size()), (1+pts_.size()+CONFIG::D));

                    std::vector<Eigen::Triplet<REAL> > coefsC;

                    //set Css
                    for( INT i=0 ; i < pts_.size() ; i++) {
                        for( INT j=i ; j < pts_.size() ; j++) {

                            auto d = norm( pts_[i] - pts_[j] );

                            if ( d < r_ ) {
                                REAL w = rbf(d);
                                coefsC.push_back( Eigen::Triplet<REAL>((i+CONFIG::D+1), (j+CONFIG::D+1), w ) );

                                if( i!=j ) {
                                    coefsC.push_back( Eigen::Triplet<REAL>( (j+CONFIG::D+1), (i+CONFIG::D+1), w ) );
                                }
                            }
                        }
                    }

    /*                 for( INT j=0; j < pts_.size(); j++) {
                        coefsC.push_back( Eigen::Triplet<REAL>(0, (j+CONFIG::D+1), 1 ) );
                        coefsC.push_back( Eigen::Triplet<REAL>((j+CONFIG::D+1), 0, 1 ) );

                        for (int tempDim=0; tempDim < CONFIG::D; tempDim++){
                            coefsC.push_back(Eigen::Triplet<REAL>((tempDim+1), (j+CONFIG::D+1), pts_[j][tempDim]));
                            coefsC.push_back(Eigen::Triplet<REAL>((j+CONFIG::D+1), (tempDim+1), pts_[j][tempDim]));
                        }
                    } */

                    Css.reserve(coefsC.size());
                    Css.setFromTriplets(coefsC.begin(), coefsC.end());

                    //set Aas
                    std::vector<Eigen::Triplet<REAL> > coefs;

                    for( INT i=0 ; i < data_points.size() ; i++) {
                        for( INT j=0 ; j < pts_.size() ; j++) {

                            auto d = norm( data_points[i].first - pts_[j] );

                            if ( d < r_ ) {
                                coefs.push_back( Eigen::Triplet<REAL>(i, (j+CONFIG::D+1), rbf(d) ) );
                            }
                        }
                    }

    /*                  for( INT i=0 ; i < data_points.size() ; i++) {

                        coefs.push_back( Eigen::Triplet<REAL>(i, 0, 1 ) );
                        for (int tempDim=0; tempDim < CONFIG::D; tempDim++){
                            coefs.push_back(Eigen::Triplet<REAL>(i, (tempDim+1), data_points[i].first[tempDim]));
                        }
                    }
     */
                    Aas.reserve(coefs.size());
                    Aas.setFromTriplets(coefs.begin(), coefs.end());

                    //invert Css
                    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> I(Css.rows(), Css.cols());
                    I.setIdentity();
                    Eigen::ConjugateGradient<Eigen::SparseMatrix<REAL> > solver(Css);
                    solver.setTolerance(1e-4);
                    Eigen::SparseMatrix<REAL> invCss = (solver.solve(I)).sparseView(1e8);

                    if( CONFIG::DEBUG ) {
                        std::cout << "#iterations of invCss:     " << solver.iterations() << ". Error of invCss: "<< solver.error()<< std::endl;
                    }


                    Eigen::SparseMatrix<REAL> AasTrans = Aas.transpose();

                    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_more = (invCss * AasTrans).pruned(1e8);
                    int r = H_more.rows();
                    int c = H_more.cols();

                    if (smoothFunc_) {
						for( INT i=0 ; i < pts_.size() ; i++) {
							for( INT j=0 ; j < data_points.size() ; j++) {
								H_toSmooth_(i, j) = H_more((i+CONFIG::D+1),j);
							}
						}
					}else{
						for( INT i=0 ; i < pts_.size() ; i++) {
							for( INT j=0 ; j < data_points.size() ; j++) {
								H_(i, j) = H_more((i+CONFIG::D+1),j);
							}
						}
					}

                    if (smoothFunc_) {
						for( INT row=0 ; row < pts_.size() ; row++) {
							for( INT j=0 ; j < data_points.size() ; j++) {
								REAL h_j_sum = 0.;
								REAL f_sum = 0.;
								for( INT k=0 ; k < M_ap_; k++) {
									INT row_k=connectivityAA_[row][k];
									if (row_k==row) {
										std::cerr << "Invalid row_k value: " << row_k << std::endl;
										} else{
											h_j_sum += std::pow(dist_h_i(row, row_k),(-2));
										}
								}

								for( INT k=0 ; k < M_ap_; k++) {
									INT row_k=connectivityAA_[row][k];
									if (row_k==row) {
										std::cerr << "Invalid row_k value: " << row_k << std::endl;
										} else{
											REAL w_i=((std::pow(dist_h_i(row, row_k),(-2)))/(h_j_sum));
											f_sum += w_i*H_toSmooth_(row_k, j);
										}
								}

								H_(row, j)=0.5*(f_sum+H_toSmooth_(row, j));
							}
						}
					}
                } else{ // Without polynomial terms

					for( int row=0; row<data_points.size(); row++ ) {
						Eigen::SparseMatrix<REAL> AA; //< Matrix of radial basis function evaluations between prescribed points
						Eigen::SparseMatrix<REAL> AB; //< Matrix of RBF evaluations between prescribed and interpolation points

						AA.resize( pts_.size(), pts_.size() );
						AB.resize( 1, pts_.size() );

						std::vector<Eigen::Triplet<REAL> > coefs;

						//set AA
						for( INT i=0 ; i < pts_.size() ; i++) {
							for( INT j=i ; j < pts_.size() ; j++) {

								int glob_i = connectivityAB_[row][i];
								int glob_j = connectivityAB_[row][j];

								auto d = norm( pts_[glob_i] - pts_[glob_j] );
								if ( d < r_ ) {
									REAL w = rbf(d);
									coefs.push_back( Eigen::Triplet<REAL>(i, j , w ) );
									if( i!=j ) {
										coefs.push_back( Eigen::Triplet<REAL>( j, i, w ) );
									}
								}
							}
						}

						AA.reserve(coefs.size());
						AA.setFromTriplets(coefs.begin(), coefs.end());

						//set AB
						coefs.clear();
						for( INT j=0 ; j < pts_.size() ; j++) {
							int glob_j = connectivityAB_[row][j];

							auto d = norm( data_points[row].first - pts_[glob_j] );
							if ( d < r_ ) {
									coefs.push_back( Eigen::Triplet<REAL>(0, j, rbf(d) ) );
								}
						}

						AB.reserve(coefs.size());
						AB.setFromTriplets(coefs.begin(), coefs.end());

						//invert AA
						Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> I(AA.rows(), AA.cols());
						I.setIdentity();
						Eigen::ConjugateGradient<Eigen::SparseMatrix<REAL> > solver(AA);
						solver.setTolerance(1e-6);
						Eigen::SparseMatrix<REAL> invAA = (solver.solve(I)).sparseView(1e8);

						if( CONFIG::DEBUG ) {
							std::cout << "#iterations:     " << solver.iterations() << ". Error: "<< solver.error()<< std::endl;
						}

						Eigen::SparseMatrix<REAL> ABTrans = AB.transpose();

						Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_j = (invAA * ABTrans).pruned(1e8);

						if (smoothFunc_) {
							for( INT i=0 ; i < pts_.size() ; i++) {
								int glob_i = connectivityAB_[row][i];
								H_toSmooth_(glob_i, row) = H_j(i,0);
							}
						}else{
							for( INT i=0 ; i < pts_.size() ; i++) {
								int glob_i = connectivityAB_[row][i];
								H_(glob_i, row) = H_j(i,0);
							}
						}
					}
                    if (smoothFunc_) {
						for( int row=0; row<data_points.size(); row++ ) {
							for( INT i=0 ; i < pts_.size() ; i++) {
								int glob_i = connectivityAB_[row][i];
								REAL h_j_sum = 0.;
								REAL f_sum = 0.;

								for( INT k=0 ; k < M_ap_; k++) {
									INT global_k=connectivityAA_[glob_i][k];
									if (global_k==glob_i) {
										std::cerr << "Invalid global_k value: " << global_k << std::endl;
										} else{
											h_j_sum += std::pow(dist_h_i(glob_i, global_k),(-2));
										}
								}

								for( INT k=0 ; k < M_ap_; k++) {
									INT global_k=connectivityAA_[glob_i][k];
									if (global_k==glob_i) {
										std::cerr << "Invalid global_k value: " << global_k << std::endl;
										} else{
											REAL w_i=((std::pow(dist_h_i(glob_i, global_k),(-2)))/(h_j_sum));
											f_sum += w_i*H_toSmooth_(global_k, row);
										}
								}

								H_(glob_i, row)=0.5*(f_sum+H_toSmooth_(glob_i, row));
							}
						}
					}
                }
            }
			if( CONFIG::DEBUG ) {
				std::cout << "run time [s]: " << float( clock () - begin_time )/  CLOCKS_PER_SEC << std::endl;
            }
			const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

            std::ofstream outputFileHMatrix(fileAddress_+"/Hmatrix.dat");

            if (!outputFileHMatrix) {
                std::cerr << "Could not locate the file address of Hmatrix.dat!" << std::endl;
            } else{
				outputFileHMatrix << "// ************************************************************************************************";
                outputFileHMatrix << "\n";
				outputFileHMatrix << "// **** This is the 'Hmatrix.dat' file of the RBF spatial sampler of the MUI library";
                outputFileHMatrix << "\n";
                outputFileHMatrix << "// **** This file contains the entire matrix of the Coupling Matrix (H).";
                outputFileHMatrix << "\n";
                outputFileHMatrix << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire H matrix";
                outputFileHMatrix << "\n";
                outputFileHMatrix << "// ************************************************************************************************";
                outputFileHMatrix << "\n";
                outputFileHMatrix << "// ";
                outputFileHMatrix << "\n";
                outputFileHMatrix << H_.format(CSVFormat);
            }
        }
		initialised_ = true;
	}

	///Radial basis function
	REAL rbf( point_type x1, point_type x2 ) {
		auto d = norm( x1 - x2 );
		return rbf(d);
	}
	
	REAL rbf(REAL d) {

		REAL w = 0.0;

		switch(basisFunc_) {

		   case 0:
			  //Gaussian
			  w = ( d < r_ ) ? std::exp( -(s_*s_ * d * d) ) : 0.0;
			  return w;
			  break;

		   case 1:
			  //Wendland's C0
			  w = std::pow(( 1. - d ), 2);
			  return w;
			  break;

		   case 2:
			  //Wendland's C2
			  w = (std::pow((1. - d), 4))*((4 * d) + 1);
			  return w;
			  break;

		   case 3:
			  //Wendland's C4
			  w = (std::pow((1. - d), 6))*((35 * d * d)+(18 * d) + 3);
			  return w;
			  break;

		   case 4:
			  //Wendland's C6
			  w = (std::pow((1. - d), 8))*((32 * d * d * d)+(25 * d * d)+(8 * d) + 1);
			  return w;
			  break;

		   default :
			  std::cerr << "MUI Error [sampler_rbf.h]: invalid RBF basis function number (" << basisFunc_ << ")" << std::endl
			  << "Please set the RBF basis function number (basisFunc_) as: "<< std::endl 
			  << "basisFunc_=0 (Gaussian); "<< std::endl
			  << "basisFunc_=1 (Wendland's C0); " << std::endl
			  << "basisFunc_=2 (Wendland's C2); " << std::endl
			  << "basisFunc_=3 (Wendland's C4); " << std::endl
			  << "basisFunc_=4 (Wendland's C6); " << std::endl;
		}

	}
	
	///Distances function
	REAL dist_h_i(INT pts_i, INT pts_j) {

		REAL h_i = 0.0;

		switch (CONFIG::D) {
			case 1:
				h_i = std::sqrt( (std::pow((pts_[pts_i][0] - pts_[pts_j][0]),2)));
				break;
			case 2:
				h_i = std::sqrt( (std::pow((pts_[pts_i][0] - pts_[pts_j][0]),2))+
										(std::pow((pts_[pts_i][1] - pts_[pts_j][1]),2)));
				break;
			case 3:
				h_i = std::sqrt( (std::pow((pts_[pts_i][0] - pts_[pts_j][0]),2))+
										(std::pow((pts_[pts_i][1] - pts_[pts_j][1]),2))+
										(std::pow((pts_[pts_i][2] - pts_[pts_j][2]),2)));;
				break;
			default:
				std::cerr << "Wrong number of CONFIG::D" << std::endl;;
		}

		return h_i;
	}

	template<template<typename,typename> class CONTAINER>
	void buildConnectivity( const CONTAINER<ITYPE,CONFIG> &data_points, const INT NP ) {

        bool warningSent = false;
        int  pointsCount = 0;


		std::ofstream outputFileCAB(fileAddress_+"/connectivityAB.dat");

		if (!outputFileCAB) {
			std::cerr << "Could not locate the file address on the connectivityAB.dat!" << std::endl;
		} else {
			outputFileCAB << "// ************************************************************************************************";
			outputFileCAB << "\n";
			outputFileCAB << "// **** This is the 'connectivityAB.dat' file of the RBF spatial sampler of the MUI library";
			outputFileCAB << "\n";
			outputFileCAB << "// **** This file contains the entire matrix of the Point Connectivity Matrix (N).";
			outputFileCAB << "\n";
			outputFileCAB << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
			outputFileCAB << "\n";
			outputFileCAB << "// ************************************************************************************************";
			outputFileCAB << "\n";
			outputFileCAB << "// ";
			outputFileCAB << "\n";
		}

        if( consistent_ ) {

            connectivityAB_.resize(pts_.size());

            for( INT i=0; i<pts_.size(); i++ ) {

                pointsCount = 0;

                for( INT n=0; n<NP; n++ ) {
                    REAL cur = 1e10;
                    INT bestj = -1;
                    for( INT j=0; j<data_points.size(); j++ ) {

                        bool added = false;
                        for( int k=0; k<connectivityAB_[i].size(); k++ ) {
                            if( connectivityAB_[i][k] == j ) {
                                added=true;
                                break;
                            }
                        }

                        if( added ) {
                            continue;
                        }

                        auto d = normsq( pts_[i] - data_points[j].first );

                        if( d<cur ) {
                            cur=d;
                            bestj = j;
                        }

                        if ((n==0) && ( d < r_*r_ ))
                            ++pointsCount;
                    }

                    connectivityAB_[i].push_back(bestj);

					outputFileCAB << bestj << ",";

                    if ((!warningSent) && (pointsCount > 120) && (n ==0)) {
                        if( !QUIET )
                        	std::cout << "MUI Warning [sampler_rbf.h]: RBF search radius too large (No. points found " << pointsCount << ")" << std::endl;
                        warningSent = true;
                    }

                }
				outputFileCAB << '\n';
            }
        } else{ //conservative

            connectivityAB_.resize(data_points.size());

            for( INT i=0; i<data_points.size(); i++ ) {

                pointsCount = 0;

                for( INT n=0; n<pts_.size(); n++ ) {
                    REAL cur = 1e10;
                    INT bestj = -1;
                    for( INT j=0; j<pts_.size(); j++ ) {

                        bool added = false;
                        for( int k=0; k<connectivityAB_[i].size(); k++ ) {
                            if( connectivityAB_[i][k] == j ) {
                                added=true;
                                break;
                            }
                        }

                        if( added ) {
                            continue;
                        }

                        auto d = normsq( data_points[i].first - pts_[j] );

                        if( d<cur ) {
                            cur=d;
                            bestj = j;
                        }

                        if ((n==0) && ( d < r_*r_ ))
                            ++pointsCount;
                    }

                    connectivityAB_[i].push_back(bestj);

					outputFileCAB << bestj << ",";

                    if ((!warningSent) && (pointsCount > 120) && (n ==0)) {
                    	if( !QUIET )
                    		std::cout << "MUI Warning [sampler_rbf.h]: RBF search radius too large (No. points found " << pointsCount << ")" << std::endl;
                        warningSent = true;
                    }
                }
				outputFileCAB << '\n';
            }
        }
		outputFileCAB.close();
		return;
	}

	void buildConnectivityAA( const INT MP ) {

        bool warningSent = false;

		std::ofstream outputFileCAA(fileAddress_+"/connectivityAA.dat");

		if (!outputFileCAA) {
			std::cerr << "Could not locate the file address on the connectivityAA.dat!" << std::endl;
		} else {
			outputFileCAA << "// ************************************************************************************************";
			outputFileCAA << "\n";
			outputFileCAA << "// **** This is the 'connectivityAA.dat' file of the RBF spatial sampler of the MUI library";
			outputFileCAA << "\n";
			outputFileCAA << "// **** This file contains the entire matrix of the Point Connectivity Matrix (M) (for smoothing).";
			outputFileCAA << "\n";
			outputFileCAA << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
			outputFileCAA << "\n";
			outputFileCAA << "// ************************************************************************************************";
			outputFileCAA << "\n";
			outputFileCAA << "// ";
			outputFileCAA << "\n";
		}

            connectivityAA_.resize(pts_.size());

            for( INT i=0; i<pts_.size(); i++ ) {

                for( INT n=0; n<MP; n++ ) {
                    REAL cur = 1e10;
                    INT bestj = -1;
                    for( INT j=0; j<pts_.size(); j++ ) {

                        bool added = false;
                        for( int k=0; k<connectivityAA_[i].size(); k++ ) {
                            if( connectivityAA_[i][k] == j ) {
                                added=true;
                                break;
                            }
                        }

                        if(( added ) || (i == j)) {
                            continue;
                        }

                        auto d = normsq( pts_[i] - pts_[j] );

                        if( d<cur ) {
                            cur=d;
                            bestj = j;
                        }
                    }

                    connectivityAA_[i].push_back(bestj);

					outputFileCAA << bestj << ",";

                }
				outputFileCAA << '\n';
            }

		outputFileCAA.close();
		return;
	}

protected:
	REAL r_;
	REAL s_;

	bool initialised_;
    INT CABrow_;
	INT CABcol_;
    INT CAArow_;
	INT CAAcol_;
	INT Hrow_;
	INT Hcol_;
	const bool conservative_;
	const bool consistent_;
	const bool polynomial_;
	const bool smoothFunc_;
    const bool readMatrix_;
    const INT basisFunc_;
	const std::string fileAddress_;

	INT N_sp_;
	INT M_ap_;
	
	std::vector<point_type> pts_;
	Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_;  //< Transformation Matrix
	Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_toSmooth_;

	std::vector<std::vector<INT> > connectivityAB_;
	std::vector<std::vector<INT> > connectivityAA_;
};

}

#endif /* MUI_SAMPLER_RBF_H_ */
