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
 * @author A. Skillen
 * @date November 2018
 * @brief Spatial sampler using Gaussian Radial Basis Function interpolation.
 */

#ifndef MUI_SAMPLER_RBF_H_
#define MUI_SAMPLER_RBF_H_

#include "../util.h"
#include <Eigen/Sparse>
#include "../uniface.h"

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

	sampler_rbf( REAL r, std::vector<point_type>& pts, bool conservative=false, double cutoff=1e-9 ) : 
	r_(r), 
	initialised_(false), 
	conservative_(conservative), 
	consistent_(!conservative),
	N_sp_(50),
	pts_(pts)
	{
		if( !CONFIG::FIXEDPOINTS ) {
			 EXCEPTION(std::runtime_error("Not (yet) implemented for dynamic points."));
		}

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

		buildConnectivity( data_points, N_sp_ );
		H_.resize( pts_.size(), data_points.size() );

		for( int row=0; row<pts_.size(); row++ ) {
			Eigen::SparseMatrix<REAL> AA; //< Matrix of radial basis function evaluations between prescribed points
			Eigen::SparseMatrix<REAL> AB; //< Matrix of RBF evaluations between prescribed and interpolation points

			if( consistent_ ) {
				AA.resize( N_sp_, N_sp_ );
				AB.resize( 1, N_sp_ );
			}
			else { //conservative
				//AA.resize( pts_.size(), pts_.size() );
				//AB.resize( data_points.size(), pts_.size() );
			}

			std::vector<Eigen::Triplet<REAL> > coefs;

			//set AA
			if( consistent_ ) {
				for( INT i=0 ; i < N_sp_ ; i++) {
					for( INT j=i ; j < N_sp_ ; j++) {
						int glob_i = connectivityAB_[row][i];
						int glob_j = connectivityAB_[row][j];

						auto d = normsq( data_points[glob_i].first - data_points[glob_j].first );
						if ( d < r_*r_ ) {
							REAL w = rbf(d);
							coefs.push_back( Eigen::Triplet<REAL>(i, j , w ) );
							if( i!=j ) {
								coefs.push_back( Eigen::Triplet<REAL>( j, i, w ) );
							}
						}
					}
				}
			}
			else { //conservative
				//for( INT i=0 ; i < pts_.size() ; i++) {
				//	for( INT j=i ; j < pts_.size() ; j++) {
				//		auto d = normsq( pts_[i] - pts_[j] );
				//		if ( d < r_*r_ ) {
				//			REAL w = rbf(d);
				//			coefs.push_back( Eigen::Triplet<REAL>(i, j , w ) );
				//			if( i!=j ) {
				//				coefs.push_back( Eigen::Triplet<REAL>( j, i, w ) );
				//			}
				//		}
				//	}
				//}
			}
			
			AA.reserve(coefs.size());
			AA.setFromTriplets(coefs.begin(), coefs.end());


			//invert AA
			Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> I(AA.rows(), AA.cols());
			I.setIdentity();
			Eigen::ConjugateGradient<Eigen::SparseMatrix<REAL> > solver(AA);
			solver.setTolerance(1e-6);
			Eigen::SparseMatrix<REAL> invAA = (solver.solve(I)).sparseView(1e8);
			
			if( CONFIG::DEBUG ) {
				std::cout << "#iterations:     " << solver.iterations() << ". Error: "<< solver.error()<< std::endl;
			}
			
			
			//set AB
			coefs.clear();
			for( INT j=0 ; j < N_sp_ ; j++) {
				int glob_j = connectivityAB_[row][j];
				
				auto d = normsq( pts_[row] - data_points[glob_j].first );
				if ( d < r_*r_ ) {
					if( consistent_ ) {
						coefs.push_back( Eigen::Triplet<REAL>(0, j, rbf(d) ) );
					}
					else { // conservative
						coefs.push_back( Eigen::Triplet<REAL>(j, 0, rbf(d) ) );
					}
				}
			}

			AB.reserve(coefs.size());
			AB.setFromTriplets(coefs.begin(), coefs.end());
			
			Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_i = (AB * invAA).pruned(1e8);

			for( INT j=0 ; j < N_sp_ ; j++) {
				int glob_j = connectivityAB_[row][j];
				H_(row, glob_j) = H_i(0,j);
			}
		}


		if( conservative_ )
			H_.transposeInPlace();

		initialised_ = true;
	}
	
	
	REAL rbf( point_type x1, point_type x2 ) {
		auto d = normsq( x1 - x2 );
		return rbf(d);
	}
	
	///Gaussian radial basis function 
	REAL rbf(REAL d) {
		REAL w = ( d < r_*r_ ) ? std::exp( -(s_*s_ * d) ) : 0.0;
		return w;
	}
	
	template<template<typename,typename> class CONTAINER>
	void buildConnectivity( const CONTAINER<ITYPE,CONFIG> &data_points, const INT NP ) {
		connectivityAB_.resize(pts_.size());
		
		for( INT i=0; i<pts_.size(); i++ ) {
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
				}
				
				connectivityAB_[i].push_back(bestj);
			}
		}
		
		return;
	}

protected:
	REAL r_;
	REAL s_;

	bool initialised_;
	const bool conservative_;
	const bool consistent_;
	
	INT N_sp_;
	
	std::vector<point_type> pts_;
	Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> H_;  //< Transformation Matrix
	
	std::vector<std::vector<INT> > connectivityAB_;
};

}

#endif /* MUI_SAMPLER_RBF_H_ */
