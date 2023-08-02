/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
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
 * @file sampler_sph_quintic.h
 * @author X. Bian
 * @date 21 July 2014
 * @brief Spatial sampler that provides a value at a point using a Smoothed Particle
 * Hydrodynamics (SPH) derived interpolation method with a quintic spline kernel.
 *
 * M. Liu, & G. Liu, Smoothed particle hydrodynamics (SPH): an overview and recent
 * developments,Archives of computational methods in engineering, 17.1 (2010), pp. 25-76.
 */

#ifndef MUI_SAMPLER_SPH_QUINTIC_H_
#define MUI_SAMPLER_SPH_QUINTIC_H_

#include "../../general/util.h"
#include "../../config.h"
#include "../sampler.h"

namespace mui
{

template<typename CONFIG=default_config, typename O_TP=typename CONFIG::REAL, typename I_TP=O_TP>
class sampler_sph_quintic
{
public:
    using OTYPE      = O_TP;
    using ITYPE      = I_TP;
    using REAL       = typename CONFIG::REAL;
    using INT        = typename CONFIG::INT;
    using point_type = typename CONFIG::point_type;
    const static int D = CONFIG::D;

    sampler_sph_quintic( REAL r_ ) : r( r_ ), hinv( REAL( 3 ) / r_ )
    {
        static_assert( D == 1 || D == 2 || D == 3, "Quintic kernel for dimension other than 1,2,3 not defined." );
        REAL sigma;
        switch( D ) {
        case 1:
            sigma = 1.0 / 120.0;
            break;
        case 2:
            sigma = 7.0 / 478.0 / PI;
            break;
        case 3:
            sigma = 1.0 / 120.0 / PI;
            break;
        }
        norm_factor = sigma * powr<D>( hinv );
    }

    template<template<typename, typename> class CONTAINER>
    inline OTYPE filter( point_type focus, const CONTAINER<ITYPE, CONFIG> &data_points ) const
    {
        OTYPE vsum = 0;
        for( size_t i = 0 ; i < data_points.size() ; i++ ) {
            auto dist2 = normsq( focus - data_points[i].first );
            if( dist2 < r * r ) {
                REAL w = quintic_polynomial( sqrt( dist2 ) );
                vsum += data_points[i].second * w;
            }
        }
        return vsum;
    }

    inline geometry::any_shape<CONFIG> support( point_type focus, REAL domain_mag ) const
    {
        return geometry::sphere<CONFIG>( focus, r );
    }

protected:
    REAL r, hinv, norm_factor;

    inline REAL quintic_polynomial( const REAL dist ) const
    {
        REAL s, w;
        REAL s1, s2, s3;
        REAL s1_5, s2_5, s3_5;
        s = dist * hinv;
        s1 = 1.0 - s;
        s2 = 2.0 - s;
        s3 = 3.0 - s;
        s1_5 = 15.0 * powr<5>( s1 );
        s2_5 = -6.0 * powr<5>( s2 );
        s3_5 = powr<5>( s3 );
        if( s < 1.0 ) {
            w = s3_5 + s2_5 + s1_5;
        } else if( s < 2.0 ) {
            w = s3_5 + s2_5;
        } else if( s < 3.0 ) {
            w = s3_5;
        } else {
            w = 0.0;
        }
        w *= norm_factor;
        return w;
    }
};

}

#endif /* MUI_SAMPLER_SPH_QUINTIC_H_ */
