/*
 * sampler_sph_quintic.h
 *
 *  Created on: July 21, 2014
 *      Author: xbian
 *
 *  SPH interpolation using quintic kernel
 */

#ifndef MUI_SAMPLER_SPH_QUINTIC_H_
#define MUI_SAMPLER_SPH_QUINTIC_H_

#include "util.h"

namespace mui
{

template<typename O_TP, typename I_TP = O_TP, typename CONFIG = default_config>
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
        for( INT i = 0 ; i < data_points.size() ; i++ ) {
            auto dist2 = normsq( focus - data_points[i].first );
            if( dist2 < r * r ) {
                REAL w = quintic_polynomial( sqrt( dist2 ) );
                vsum += data_points[i].second * w;
            }
        }
        return vsum;
    }

    inline geometry::any_shape<CONFIG> support( point_type focus ) const
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
