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
 * @file algo_fixed_relaxation.h
 * @author W. Liu
 * @date 05 October 2022
 * @brief Fixed Relaxation coupling algorithm
 */

#ifndef MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_
#define MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_

#include "../util.h"
#include "../config.h"
#include <iterator>

namespace mui {

template<typename CONFIG=default_config> class algo_fixed_relaxation {
public:
    using REAL       = typename CONFIG::REAL;
    using INT        = typename CONFIG::INT;
    using time_type  = typename CONFIG::time_type;
    using point_type = typename CONFIG::point_type;

    algo_fixed_relaxation( REAL under_relaxation_factor = 1.0,
       std::vector<std::pair<point_type, REAL>> pts_value_init =
        std::vector<std::pair<point_type, REAL>>() ):
       init_under_relaxation_factor_(under_relaxation_factor) {

        under_relaxation_factor_ = init_under_relaxation_factor_;

        if (!pts_value_init.empty()) {
            pts_time_value_.insert(pts_time_value_.begin(),
                std::make_pair(
                    std::make_pair(
                        std::numeric_limits<time_type>::lowest(),
                        std::numeric_limits<time_type>::lowest()
                    ),pts_value_init
                )
            );
        }
    }

    //- relaxation based on single time value
    template<typename OTYPE>
    OTYPE relaxation(std::pair<time_type, time_type> t, point_type focus, OTYPE filtered_value) const {

        OTYPE filtered_old_value = 0.0;

        if (pts_time_value_.empty()) {

            filtered_old_value = 0.0;

            std::vector<std::pair<point_type, REAL>> pts_value_temp{
                std::make_pair(focus,calculate_relaxed_value(filtered_value,filtered_old_value))
                };

            pts_time_value_.insert(pts_time_value_.begin(),std::make_pair(t,pts_value_temp));

            return calculate_relaxed_value(filtered_value,filtered_old_value);

        } else { // pts_time_value_ not empty

            auto present_iter = std::find_if(pts_time_value_.begin(), pts_time_value_.end(),
                [t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
                    return (((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) &&
                            ((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon()));
                });

            auto previous_iter = std::find_if(pts_time_value_.begin(), pts_time_value_.end(),
                [t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
                    return ((((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) &&
                            b.first.second < t.second) || (b.first.first < t.first));
                });

            if ((present_iter == std::end(pts_time_value_)) &&
                (previous_iter == std::end(pts_time_value_)) ) {

                std::cerr << "Non-monotonic time marching does not (yet) supported for the Fixed Relaxation coupling method! " << std::endl;

            } else if ((present_iter != std::end(pts_time_value_)) &&
                (previous_iter == std::end(pts_time_value_)) ) {

                    auto pts_relx_val_iter = std::find_if(present_iter->second.begin(),
                        present_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                    return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                    });

                    if ( pts_relx_val_iter == std::end(present_iter->second) ) {

                        // Interpolate the relaxed value by N2_linear
                        REAL r2min_1st = std::numeric_limits<REAL>::max();
                        REAL r2min_2nd = std::numeric_limits<REAL>::max();
                        OTYPE value_1st = 0, value_2nd = 0;
                        for( size_t i = 0 ; i < present_iter->second.size() ; i++ ) {
                            REAL dr2 = normsq( focus - present_iter->second[i].first );
                            if ( dr2 < r2min_1st ) {
                                r2min_2nd = r2min_1st;
                                value_2nd = value_1st;
                                r2min_1st = dr2;
                                value_1st = present_iter->second[i].second ;
                            } else if ( dr2 < r2min_2nd ) {
                                r2min_2nd = dr2;
                                value_2nd = present_iter->second[i].second ;
                            }
                        }

                        REAL r1 = std::sqrt( r2min_1st );
                        REAL r2 = std::sqrt( r2min_2nd );

                        auto relaxed_value_temp = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

                        present_iter->second.insert(present_iter->second.begin(),std::make_pair(focus, relaxed_value_temp));

                        return relaxed_value_temp;

                    } else {

                        return pts_relx_val_iter->second;

                    }

            } else if ((present_iter == std::end(pts_time_value_)) &&
                (previous_iter != std::end(pts_time_value_)) ) {

                    auto pts_relx_val_iter = std::find_if(previous_iter->second.begin(),
                        previous_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                    return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                    });

                    if ( pts_relx_val_iter == std::end(previous_iter->second) ) {

                        // Interpolate the relaxed value by N2_linear
                        REAL r2min_1st = std::numeric_limits<REAL>::max();
                        REAL r2min_2nd = std::numeric_limits<REAL>::max();
                        OTYPE value_1st = 0, value_2nd = 0;
                        for( size_t i = 0 ; i < previous_iter->second.size() ; i++ ) {
                            REAL dr2 = normsq( focus - previous_iter->second[i].first );
                            if ( dr2 < r2min_1st ) {
                                r2min_2nd = r2min_1st;
                                value_2nd = value_1st;
                                r2min_1st = dr2;
                                value_1st = previous_iter->second[i].second ;
                            } else if ( dr2 < r2min_2nd ) {
                                r2min_2nd = dr2;
                                value_2nd = previous_iter->second[i].second ;
                            }
                        }

                        REAL r1 = std::sqrt( r2min_1st );
                        REAL r2 = std::sqrt( r2min_2nd );

                        filtered_old_value = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

                    } else {

                        filtered_old_value = pts_relx_val_iter->second;

                    }

                    std::vector<std::pair<point_type, REAL>> pts_value_temp{
                        std::make_pair(focus,calculate_relaxed_value(filtered_value,filtered_old_value))
                    };

                    pts_time_value_.insert(pts_time_value_.begin(),std::make_pair(t, pts_value_temp));

                    return calculate_relaxed_value(filtered_value,filtered_old_value);

            } else {

                    auto pts_relx_val_iter = std::find_if(previous_iter->second.begin(),
                        previous_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                    return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                    });

                    if ( pts_relx_val_iter == std::end(previous_iter->second) ) {

                        // Interpolate the relaxed value by N2_linear
                        REAL r2min_1st = std::numeric_limits<REAL>::max();
                        REAL r2min_2nd = std::numeric_limits<REAL>::max();
                        OTYPE value_1st = 0, value_2nd = 0;
                        for( size_t i = 0 ; i < previous_iter->second.size() ; i++ ) {
                            REAL dr2 = normsq( focus - previous_iter->second[i].first );
                            if ( dr2 < r2min_1st ) {
                                r2min_2nd = r2min_1st;
                                value_2nd = value_1st;
                                r2min_1st = dr2;
                                value_1st = previous_iter->second[i].second ;
                            } else if ( dr2 < r2min_2nd ) {
                                r2min_2nd = dr2;
                                value_2nd = previous_iter->second[i].second ;
                            }
                        }

                        REAL r1 = std::sqrt( r2min_1st );
                        REAL r2 = std::sqrt( r2min_2nd );

                        filtered_old_value = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

                    } else {

                        filtered_old_value = pts_relx_val_iter->second;

                    }

                    auto pts_present_relx_val_iter = std::find_if(present_iter->second.begin(),
                        present_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                    return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                    });

                    if ( pts_present_relx_val_iter == std::end(present_iter->second) ) {

                        present_iter->second.insert(present_iter->second.begin(),
                            std::make_pair(focus,calculate_relaxed_value(
                                filtered_value,filtered_old_value)
                            )
                        );

                    } else {

                        pts_present_relx_val_iter->second = calculate_relaxed_value(filtered_value,filtered_old_value);

                    }

                    return calculate_relaxed_value(filtered_value,filtered_old_value);

            }

        }
    }

private:
    template<typename OTYPE>
    OTYPE calculate_relaxed_value(OTYPE filtered_value, OTYPE filtered_old_value) {

        return (under_relaxation_factor_ * filtered_value) + ((1 - under_relaxation_factor_) * filtered_old_value);

    }

protected:
    REAL init_under_relaxation_factor_;
    REAL under_relaxation_factor_;

    mutable std::vector<std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>>> pts_time_value_;

};

}

#endif /* MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_ */
