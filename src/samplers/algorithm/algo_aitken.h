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
 * @file algo_aitken.h
 * @author W. Liu
 * @date 05 October 2022
 * @brief Aitken coupling algorithm
 */

#ifndef MUI_COUPLING_ALGORITHM_AITKEN_H_
#define MUI_COUPLING_ALGORITHM_AITKEN_H_

#include "../../general/util.h"
#include "../../config.h"
#include <mpi.h>
#include <iterator>

namespace mui {

template<typename CONFIG=default_config> class algo_aitken {
public:
    using REAL       = typename CONFIG::REAL;
    using INT        = typename CONFIG::INT;
    using time_type  = typename CONFIG::time_type;
    using iterator_type = typename CONFIG::iterator_type;
    using point_type = typename CONFIG::point_type;

    algo_aitken( REAL under_relaxation_factor = 1.0,
       REAL under_relaxation_factor_max = 1.0,
      MPI_Comm local_comm = MPI_COMM_NULL,
       std::vector<std::pair<point_type, REAL>> pts_vlu_init =
        std::vector<std::pair<point_type, REAL>>(),
       REAL res_l2_norm_nm1 = 0.0):
      init_under_relaxation_factor_(under_relaxation_factor),
      under_relaxation_factor_max_(under_relaxation_factor_max),
      local_mpi_comm_world_(local_comm) {
        minimum_iterator_ = 1;
        under_relaxation_factor_.insert(under_relaxation_factor_.begin(),
            std::make_pair(std::make_pair(
                std::numeric_limits<time_type>::lowest(),
                (minimum_iterator_ -1)), init_under_relaxation_factor_
            )
        );

        if (!pts_vlu_init.empty()) {
            pts_time_value_.insert(pts_time_value_.begin(),
                std::make_pair(std::make_pair(
                    std::numeric_limits<time_type>::lowest(),
                    (minimum_iterator_ -1)),pts_vlu_init
                )
            );
        }

        if(res_l2_norm_nm1 != 0.0){
            residual_l2_norm_.insert(residual_l2_norm_.begin(),
                std::make_pair(std::make_pair(
                    std::numeric_limits<time_type>::lowest(),
                    (minimum_iterator_ -1)),
                    (std::make_pair(static_cast<INT>(0), res_l2_norm_nm1)
                    )
                )
            );
        }
    }

    //- relaxation based on single time value
    template<typename OTYPE>
    OTYPE relaxation(std::pair<time_type,iterator_type> t, point_type focus, OTYPE filtered_value) const {

        OTYPE filtered_old_value = 0.0;

        if (pts_time_value_.empty()) {

            assert(pts_time_res_.empty());

            filtered_old_value = 0.0;

            std::vector<std::pair<point_type, REAL>> pts_value_temp{
                std::make_pair(focus,calculate_relaxed_value(t,filtered_value,filtered_old_value))
                };

            pts_time_value_.insert(pts_time_value_.begin(),std::make_pair(t,pts_value_temp));

            std::vector<std::pair<point_type, REAL>> pts_res_temp{
                std::make_pair(focus,calculate_point_residual(t,filtered_value,filtered_old_value))
                };

            pts_time_res_.insert(pts_time_res_.begin(),std::make_pair(t,pts_res_temp));

            return calculate_relaxed_value(t,filtered_value,filtered_old_value);

        } else {

                auto present_iter = std::find_if(pts_time_value_.begin(), pts_time_value_.end(),
                    [t](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                            (t.second == b.first.second);});

                auto mi = minimum_iterator_;
                auto previous_iter = std::find_if(pts_time_value_.begin(), pts_time_value_.end(),
                    [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((t.second == mi) ?
                                (b.first.first < t.first) ||
                                 (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                  (b.first.second == (mi - 1))) :
                                ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                 (b.first.second < t.second));});

                auto present_res_iter = std::find_if(pts_time_res_.begin(), pts_time_res_.end(),
                    [t](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                            (t.second == b.first.second);});

                mi = minimum_iterator_;
                auto previous_res_iter = std::find_if(pts_time_res_.begin(), pts_time_res_.end(),
                    [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((t.second == mi) ?
                                (b.first.first < t.first) ||
                                 (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                  (b.first.second == (mi - 1))) :
                                ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                 (b.first.second < t.second));});

                if ((present_iter == std::end(pts_time_value_)) &&
                    (previous_iter == std::end(pts_time_value_)) ) {

                    assert((present_res_iter == std::end(pts_time_res_)) || pts_time_res_.empty());

                    std::cerr << "Non-monotonic time marching does not (yet) supported for the Aitken coupling method! " << std::endl;

                } else if ((present_iter != std::end(pts_time_value_)) &&
                    (previous_iter == std::end(pts_time_value_)) ) {

                        assert((((present_iter->first.first - present_res_iter->first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (present_iter->first.second == present_res_iter->first.second)) || pts_time_res_.empty());

                        if (!pts_time_res_.empty()) {
                            assert(!residual_l2_norm_.empty());
                        }

                        auto pts_relx_value_iter = std::find_if(present_iter->second.begin(),
                            present_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                        return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                        });

                        auto pts_relx_res_iter = std::find_if(present_res_iter->second.begin(),
                            present_res_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                        return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                        });

                        if ( pts_relx_value_iter == std::end(present_iter->second) ) {

                            assert((pts_relx_res_iter == std::end(present_res_iter->second)) || pts_time_res_.empty());

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

                            if (pts_time_res_.empty()) {

                                std::vector<std::pair<point_type, REAL>> pts_res_temp{
                                    std::make_pair(focus, (filtered_value-relaxed_value_temp))
                                    };

                                pts_time_res_.insert(pts_time_res_.begin(),std::make_pair(t,pts_res_temp));

                            } else {
                                present_res_iter->second.insert(present_res_iter->second.begin(),std::make_pair(focus, (filtered_value-relaxed_value_temp)));
                            }

                            return relaxed_value_temp;

                        } else {

                            assert((normsq(pts_relx_value_iter->first - pts_relx_res_iter->first) < std::numeric_limits<REAL>::epsilon()) || pts_time_res_.empty());

                            return pts_relx_value_iter->second;

                        }

                } else if ((present_iter == std::end(pts_time_value_)) &&
                    (previous_iter != std::end(pts_time_value_)) ) {

                        assert((present_res_iter == std::end(pts_time_res_)) || pts_time_res_.empty());

                        auto pts_relx_value_iter = std::find_if(previous_iter->second.begin(),
                            previous_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                        return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                        });

                        if ( pts_relx_value_iter == std::end(previous_iter->second) ) {

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

                            filtered_old_value = pts_relx_value_iter->second;

                        }

                        std::vector<std::pair<point_type, REAL>> pts_value_temp{
                            std::make_pair(focus,calculate_relaxed_value(t,filtered_value,filtered_old_value))
                        };

                        pts_time_value_.insert(pts_time_value_.begin(),std::make_pair(t, pts_value_temp));

                        std::vector<std::pair<point_type, REAL>> pts_res_temp{
                            std::make_pair(focus,calculate_point_residual(t,filtered_value,filtered_old_value))
                            };

                        pts_time_res_.insert(pts_time_res_.begin(),std::make_pair(t,pts_res_temp));

                        present_iter = std::find_if(pts_time_value_.begin(), pts_time_value_.end(),
                            [t](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                                return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                        (t.second == b.first.second);});

                        auto mi = minimum_iterator_;
                        previous_iter = std::find_if(pts_time_value_.begin(), pts_time_value_.end(),
                            [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                                return ((t.second == mi) ?
                                        (b.first.first < t.first) ||
                                         (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                          (b.first.second == (mi - 1))) :
                                        ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                         (b.first.second < t.second));});

                        present_res_iter = std::find_if(pts_time_res_.begin(), pts_time_res_.end(),
                            [t](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                                return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                        (t.second == b.first.second);});

                        mi = minimum_iterator_;
                        previous_res_iter = std::find_if(pts_time_res_.begin(), pts_time_res_.end(),
                            [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                            return ((t.second == mi) ?
                                    (b.first.first < t.first) ||
                                     (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                      (b.first.second == (mi - 1))) :
                                    ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                     (b.first.second < t.second));});

                        auto pts_residual_l2_norm_iter = std::find_if(residual_l2_norm_.begin(),
                            residual_l2_norm_.end(), [previous_iter](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
                        return ((previous_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (previous_iter->first.second == b.first.second);});

                        if(pts_residual_l2_norm_iter == std::end(residual_l2_norm_)) {

                            if (previous_res_iter!=std::end(pts_time_res_)) {

                                assert(((previous_iter->first.first - previous_res_iter->first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                        (previous_iter->first.second == previous_res_iter->first.second) );

                                REAL local_residual_mag_sq_sum_temp = 0.0;
                                REAL residual_mag_sq_sum_temp = 0.0;

                                for (auto & element_pair : previous_res_iter->second) {
                                    local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                                }

                                if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                                    residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                                } else {
                                    MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                                }

                                if((residual_mag_sq_sum_temp != 0) || (!residual_l2_norm_.empty())){
                                    residual_l2_norm_.insert(residual_l2_norm_.begin(),
                                        std::make_pair(
                                            previous_res_iter->first, (
                                                std::make_pair(
                                                    static_cast<INT>(
                                                        previous_res_iter->second.size()
                                                    ), std::sqrt(residual_mag_sq_sum_temp)
                                                )
                                            )
                                        )
                                    );

                                    pts_residual_l2_norm_iter = std::find_if(residual_l2_norm_.begin(),
                                        residual_l2_norm_.end(), [previous_iter](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
                                    return ((previous_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                            (previous_iter->first.second == b.first.second);});
                                }
                            }

                        } else {

                            if(pts_residual_l2_norm_iter->second.first != 0) {

                                auto pts_time_res_iter = std::find_if(pts_time_res_.begin(),
                                pts_time_res_.end(), [previous_res_iter](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                                return ((previous_res_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                        (previous_res_iter->first.second == b.first.second);});

                                assert(pts_time_res_iter != std::end(pts_time_res_) );

                                if(pts_time_res_iter->second.size() != static_cast<size_t>(pts_residual_l2_norm_iter->second.first)) {

                                    REAL local_residual_mag_sq_sum_temp = 0.0;
                                    REAL residual_mag_sq_sum_temp = 0.0;

                                    for (auto & element_pair : pts_time_res_iter->second) {
                                        local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                                    }

                                    if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                                        residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                                    } else {
                                        MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                                    }

                                    pts_residual_l2_norm_iter->second.second = std::sqrt(residual_mag_sq_sum_temp);
                                }
                            }
                        }

                        return calculate_relaxed_value(t,filtered_value,filtered_old_value);

                } else {

                        assert((((present_iter->first.first - present_res_iter->first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (present_iter->first.second == present_res_iter->first.second)) || pts_time_res_.empty());

                        auto pts_relx_value_iter = std::find_if(previous_iter->second.begin(),
                            previous_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                        return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                        });

                        if ( pts_relx_value_iter == std::end(previous_iter->second) ) {

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

                            filtered_old_value = pts_relx_value_iter->second;

                        }

                        auto pts_present_relx_value_iter = std::find_if(present_iter->second.begin(),
                            present_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                        return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                        });

                        auto pts_relx_res_iter = std::find_if(present_res_iter->second.begin(),
                            present_res_iter->second.end(), [focus](std::pair<point_type, REAL> b) {
                        return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
                        });

                        if ( pts_present_relx_value_iter == std::end(present_iter->second) ) {

                            assert((pts_relx_res_iter == std::end(present_res_iter->second)) || pts_time_res_.empty());

                            present_iter->second.insert(present_iter->second.begin(),
                                std::make_pair(focus,calculate_relaxed_value(
                                    t,filtered_value,filtered_old_value)
                                )
                            );

                            if (pts_time_res_.empty()) {

                                std::vector<std::pair<point_type, REAL>> pts_res_temp{
                                    std::make_pair(focus,calculate_point_residual(
                                        t,filtered_value,filtered_old_value))
                                    };

                                pts_time_res_.insert(pts_time_res_.begin(),std::make_pair(t,pts_res_temp));

                            } else {

                                present_res_iter->second.insert(present_res_iter->second.begin(),
                                    std::make_pair(focus,calculate_point_residual(
                                        t,filtered_value,filtered_old_value)
                                    )
                                );
                            }
                        } else {

                            assert((normsq(pts_present_relx_value_iter->first - pts_relx_res_iter->first) < std::numeric_limits<REAL>::epsilon()) || pts_time_res_.empty());

                            pts_present_relx_value_iter->second = calculate_relaxed_value(t,filtered_value,filtered_old_value);

                            if (pts_time_res_.empty()) {

                                std::vector<std::pair<point_type, REAL>> pts_res_temp{
                                    std::make_pair(focus,calculate_point_residual(
                                        t,filtered_value,filtered_old_value))
                                    };

                                pts_time_res_.insert(pts_time_res_.begin(),std::make_pair(t,pts_res_temp));

                            } else {
                                pts_relx_res_iter->second = calculate_point_residual(t,filtered_value,filtered_old_value);
                            }
                        }

                        return calculate_relaxed_value(t,filtered_value,filtered_old_value);
                }
                return calculate_relaxed_value(t,filtered_value,filtered_old_value);
        }
    }

    REAL get_under_relaxation_factor (time_type t_single) {

        std::pair<time_type,iterator_type> t = std::make_pair(std::numeric_limits<time_type>::lowest(),t_single);

        update_under_relaxation_factor(t);

        auto under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
            under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
            return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
            (t.second == b.first.second);});

        assert(under_relaxation_present_iter != std::end(under_relaxation_factor_) );

        return under_relaxation_present_iter->second;

    }

    REAL get_under_relaxation_factor (time_type t, iterator_type it) {

        std::pair<time_type,iterator_type> time = std::make_pair(t,it);

        update_under_relaxation_factor(time);

        auto under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
            under_relaxation_factor_.end(), [time](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
            return ((time.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
            (time.second == b.first.second);});

        assert(under_relaxation_present_iter != std::end(under_relaxation_factor_) );

        return under_relaxation_present_iter->second;

    }

    REAL get_residual_L2_Norm (time_type t_single) {

        std::pair<time_type,iterator_type> t = std::make_pair(std::numeric_limits<time_type>::lowest(),t_single);

        update_under_relaxation_factor(t);

        auto mi = minimum_iterator_;
        auto res_l2_norm_nm1_iter = std::find_if(residual_l2_norm_.begin(),
            residual_l2_norm_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
        return ((t.second == mi) ?
                (b.first.first < t.first) ||
                 (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                  (b.first.second == (mi - 1))) :
                ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                 (b.first.second < t.second));});

        if(res_l2_norm_nm1_iter != std::end(residual_l2_norm_) ) {
            return res_l2_norm_nm1_iter->second.second;
        } else {
            return 0.0;
        }
    }

    REAL get_residual_L2_Norm (time_type t, iterator_type it) {

        std::pair<time_type,iterator_type> time = std::make_pair(t,it);

        update_under_relaxation_factor(time);

        auto mi = minimum_iterator_;
        auto res_l2_norm_nm1_iter = std::find_if(residual_l2_norm_.begin(),
            residual_l2_norm_.end(), [time, &mi](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
        return ((time.second == mi) ?
                (b.first.first < time.first) ||
                 (((time.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                  (b.first.second == (mi - 1))) :
                ((time.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                 (b.first.second < time.second));});

        if(res_l2_norm_nm1_iter != std::end(residual_l2_norm_) ) {
            return res_l2_norm_nm1_iter->second.second;
        } else {
            return 0.0;
        }
    }

private:
    template<typename OTYPE>
    OTYPE calculate_relaxed_value(std::pair<time_type,iterator_type> t, OTYPE filtered_value, OTYPE filtered_old_value) const {
 
        update_under_relaxation_factor(t);

        auto under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
            under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
            return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
            (t.second == b.first.second);});

        assert(under_relaxation_present_iter != std::end(under_relaxation_factor_) );

        return (under_relaxation_present_iter->second * filtered_value) +
            ((1 - under_relaxation_present_iter->second) * filtered_old_value);

    }

    void update_under_relaxation_factor(std::pair<time_type,iterator_type> t) const {

        auto under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
            under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
        return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
            (t.second == b.first.second);});

        auto mi = minimum_iterator_;
        auto under_relaxation_prev_iter = std::find_if(under_relaxation_factor_.begin(),
            under_relaxation_factor_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
        return ((t.second == mi) ?
                (b.first.first < t.first) ||
                 (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                  (b.first.second == (mi - 1))) :
                ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                 (b.first.second < t.second));});

        mi = minimum_iterator_;
        auto res_l2_norm_nm1_iter = std::find_if(residual_l2_norm_.begin(),
            residual_l2_norm_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
        return ((t.second == mi) ?
                (b.first.first < t.first) ||
                 (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                  (b.first.second == (mi - 1))) :
                ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                 (b.first.second < t.second));});

        mi = minimum_iterator_;
        auto res_l2_norm_nm2_iter = std::find_if(residual_l2_norm_.begin(),
            residual_l2_norm_.end(), [res_l2_norm_nm1_iter, &mi](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
        return ((res_l2_norm_nm1_iter->first.second == mi) ?
                (b.first.first < res_l2_norm_nm1_iter->first.first) ||
                 (((res_l2_norm_nm1_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                  (b.first.second == (mi - 1))) :
                ((res_l2_norm_nm1_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                 (b.first.second < res_l2_norm_nm1_iter->first.second));});


        if (t.second <= minimum_iterator_) {

            under_relaxation_factor_.insert(under_relaxation_factor_.begin(),std::make_pair(t, init_under_relaxation_factor_));

            under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
                under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
            return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                    (t.second == b.first.second);});

            auto mi = minimum_iterator_;
            under_relaxation_prev_iter = std::find_if(under_relaxation_factor_.begin(),
                under_relaxation_factor_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
            return ((t.second == mi) ?
                    (b.first.first < t.first) ||
                     (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                      (b.first.second == (mi - 1))) :
                    ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                     (b.first.second < t.second));});

            auto present_res_iter = std::find_if(pts_time_res_.begin(), pts_time_res_.end(),
                [t](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                    return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                            (t.second == b.first.second);});

            mi = minimum_iterator_;
            auto res_l2_norm_iter = std::find_if(residual_l2_norm_.begin(),
                residual_l2_norm_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>> b) {
            return (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                     (b.first.second == t.second));});

            if ((present_res_iter != std::end(pts_time_res_)) && (res_l2_norm_iter == std::end(residual_l2_norm_))) {

                REAL local_residual_mag_sq_sum_temp = 0.0;
                REAL residual_mag_sq_sum_temp = 0.0;

                for (auto & element_pair : present_res_iter->second) {
                    local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                }

                if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                    residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                } else {
                    MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                }

                if((residual_mag_sq_sum_temp != 0) || (!residual_l2_norm_.empty())){
                    residual_l2_norm_.insert(residual_l2_norm_.begin(),
                        std::make_pair(
                            present_res_iter->first, (
                                std::make_pair(
                                    static_cast<INT>(
                                        present_res_iter->second.size()
                                    ), std::sqrt(residual_mag_sq_sum_temp)
                                )
                            )
                        )
                    );
                }
            }

        } else {

            if (under_relaxation_present_iter==std::end(under_relaxation_factor_)) {

                if(res_l2_norm_nm2_iter == std::end(residual_l2_norm_) ) {

                    under_relaxation_factor_.insert(under_relaxation_factor_.begin(),std::make_pair(t, init_under_relaxation_factor_));

                    under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
                        under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                    return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                            (t.second == b.first.second);});

                    auto mi = minimum_iterator_;
                    under_relaxation_prev_iter = std::find_if(under_relaxation_factor_.begin(),
                        under_relaxation_factor_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                    return ((t.second == mi) ?
                            (b.first.first < t.first) ||
                             (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                              (b.first.second == (mi - 1))) :
                            ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                             (b.first.second < t.second));});


                } else {

                    if(res_l2_norm_nm2_iter->second.first != 0 ) {
                        auto pts_time_res_nm2_iter = std::find_if(pts_time_res_.begin(),
                        pts_time_res_.end(), [res_l2_norm_nm2_iter](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((res_l2_norm_nm2_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (res_l2_norm_nm2_iter->first.second == b.first.second);});

                        assert(pts_time_res_nm2_iter != std::end(pts_time_res_) );

                        if(pts_time_res_nm2_iter->second.size() != static_cast<size_t>(res_l2_norm_nm2_iter->second.first)) {

                            REAL local_residual_mag_sq_sum_temp = 0.0;
                            REAL residual_mag_sq_sum_temp = 0.0;

                            for (auto & element_pair : pts_time_res_nm2_iter->second) {
                                local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                            }

                            if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                                residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                            } else {
                                MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                            }

                            res_l2_norm_nm2_iter->second.second = std::sqrt(residual_mag_sq_sum_temp);

                            update_under_relaxation_factor(res_l2_norm_nm2_iter->first);
                        }

                        auto pts_time_res_nm1_iter = std::find_if(pts_time_res_.begin(),
                        pts_time_res_.end(), [res_l2_norm_nm1_iter](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((res_l2_norm_nm1_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (res_l2_norm_nm1_iter->first.second == b.first.second);});

                        if(pts_time_res_nm1_iter->second.size() != static_cast<size_t>(res_l2_norm_nm1_iter->second.first)) {

                            REAL local_residual_mag_sq_sum_temp = 0.0;
                            REAL residual_mag_sq_sum_temp = 0.0;

                            for (auto & element_pair : pts_time_res_nm1_iter->second) {
                                local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                            }

                            if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                                residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                            } else {
                                MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                            }

                            res_l2_norm_nm1_iter->second.second = std::sqrt(residual_mag_sq_sum_temp);

                            update_under_relaxation_factor(res_l2_norm_nm1_iter->first);
                        }
                    }

                    double nominator   = res_l2_norm_nm2_iter->second.second;
                    double denominator = res_l2_norm_nm1_iter->second.second - res_l2_norm_nm2_iter->second.second;

                    if (denominator != 0.0 ) {

                        if (under_relaxation_prev_iter==std::end(under_relaxation_factor_)) {
                            under_relaxation_factor_.insert(under_relaxation_factor_.begin(),
                                std::make_pair(
                                    t, calculate_aitken_constraint_pnz_control(init_under_relaxation_factor_)
                                )
                            );

                            under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
                                under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                            return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (t.second == b.first.second);
                            });

                            auto mi = minimum_iterator_;
                            under_relaxation_prev_iter = std::find_if(under_relaxation_factor_.begin(),
                                under_relaxation_factor_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                            return ((t.second == mi) ?
                                    (b.first.first < t.first) ||
                                     (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                      (b.first.second == (mi - 1))) :
                                    ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                     (b.first.second < t.second));});
                        } else {
                            under_relaxation_factor_.insert(under_relaxation_factor_.begin(),
                                std::make_pair(
                                    t, calculate_aitken_constraint_pnz_control(
                                        -under_relaxation_prev_iter->second * (nominator/denominator)
                                    )
                                )
                            );

                            under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
                                under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                            return ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (t.second == b.first.second);});

                            auto mi = minimum_iterator_;
                            under_relaxation_prev_iter = std::find_if(under_relaxation_factor_.begin(),
                                under_relaxation_factor_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                            return ((t.second == mi) ?
                                    (b.first.first < t.first) ||
                                     (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                      (b.first.second == (mi - 1))) :
                                    ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                     (b.first.second < t.second));});
                        }

                    } else {
                        under_relaxation_factor_.insert(under_relaxation_factor_.begin(),std::make_pair(t, init_under_relaxation_factor_));

                        under_relaxation_present_iter = std::find_if(under_relaxation_factor_.begin(),
                            under_relaxation_factor_.end(), [t](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                        return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) &&
                                (t.second == b.first.second);});

                        auto mi = minimum_iterator_;
                        under_relaxation_prev_iter = std::find_if(under_relaxation_factor_.begin(),
                            under_relaxation_factor_.end(), [t, &mi](std::pair<std::pair<time_type,iterator_type>, REAL> b) {
                        return ((t.second == mi) ?
                                (b.first.first < t.first) ||
                                 (((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                  (b.first.second == (mi - 1))) :
                                ((t.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                 (b.first.second < t.second));});
                    }
                }
            } else {

                assert(under_relaxation_present_iter!=std::end(under_relaxation_factor_));

                if(res_l2_norm_nm2_iter == std::end(residual_l2_norm_) ) {
                    if(under_relaxation_present_iter->second != init_under_relaxation_factor_) {
                        std::cout << "change under Relax Factor to its initial value." << std::endl;
                        under_relaxation_present_iter->second = init_under_relaxation_factor_;
                    }
                } else {

                    if(res_l2_norm_nm2_iter->second.first != 0 ) {
                        auto pts_time_res_nm2_iter = std::find_if(pts_time_res_.begin(),
                        pts_time_res_.end(), [res_l2_norm_nm2_iter](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return ((res_l2_norm_nm2_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (res_l2_norm_nm2_iter->first.second == b.first.second);});

                        assert(pts_time_res_nm2_iter != std::end(pts_time_res_) );

                        if(pts_time_res_nm2_iter->second.size() != static_cast<size_t>(res_l2_norm_nm2_iter->second.first)) {

                            REAL local_residual_mag_sq_sum_temp = 0.0;
                            REAL residual_mag_sq_sum_temp = 0.0;

                            for (auto & element_pair : pts_time_res_nm2_iter->second) {
                                local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                            }

                            if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                                residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                            } else {
                                MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                            }

                            res_l2_norm_nm2_iter->second.second = std::sqrt(residual_mag_sq_sum_temp);

                        }

                        auto pts_time_res_nm1_iter = std::find_if(pts_time_res_.begin(),
                        pts_time_res_.end(), [res_l2_norm_nm1_iter](std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>> b) {
                        return (((res_l2_norm_nm1_iter->first.first - b.first.first) < std::numeric_limits<time_type>::epsilon()) &&
                                (res_l2_norm_nm1_iter->first.second == b.first.second));});

                        if(pts_time_res_nm1_iter->second.size() != static_cast<size_t>(res_l2_norm_nm1_iter->second.first)) {

                            REAL local_residual_mag_sq_sum_temp = 0.0;
                            REAL residual_mag_sq_sum_temp = 0.0;

                            for (auto & element_pair : pts_time_res_nm1_iter->second) {
                                local_residual_mag_sq_sum_temp += std::pow(element_pair.second, 2);
                            }

                            if ( local_mpi_comm_world_ == MPI_COMM_NULL ) {
                                residual_mag_sq_sum_temp = local_residual_mag_sq_sum_temp;
                            } else {
                                MPI_Allreduce(&local_residual_mag_sq_sum_temp, &residual_mag_sq_sum_temp, 1, MPI_DOUBLE, MPI_SUM, local_mpi_comm_world_);
                            }

                            res_l2_norm_nm1_iter->second.second = std::sqrt(residual_mag_sq_sum_temp);

                            update_under_relaxation_factor(res_l2_norm_nm1_iter->first);
                        }
                    }

                    double nominator   = res_l2_norm_nm2_iter->second.second;
                    double denominator = res_l2_norm_nm1_iter->second.second - res_l2_norm_nm2_iter->second.second;

                    if (denominator != 0.0 ) {

                        if (under_relaxation_prev_iter==std::end(under_relaxation_factor_)) {
                            if(under_relaxation_present_iter->second != calculate_aitken_constraint_pnz_control(init_under_relaxation_factor_)
                            ) {
                                std::cout << "Update under Relx Factor. Position 01." << std::endl;
                                under_relaxation_present_iter->second = calculate_aitken_constraint_pnz_control(init_under_relaxation_factor_);
                            }
                        } else {
                            if(under_relaxation_present_iter->second != calculate_aitken_constraint_pnz_control(
                                -under_relaxation_prev_iter->second * (nominator/denominator))
                            ) {
                                std::cout << "Update under Relx Factor. Position 02." << std::endl;
                                under_relaxation_present_iter->second = calculate_aitken_constraint_pnz_control(
                                    -under_relaxation_prev_iter->second * (nominator/denominator)
                                );
                            }
                        }
                    } else {
                        if(under_relaxation_present_iter->second != init_under_relaxation_factor_) {
                                std::cout << "Update under Relx Factor. Position 03." << std::endl;
                                under_relaxation_present_iter->second = init_under_relaxation_factor_;
                        }
                    }
                }
            }
        }
    }

    REAL calculate_aitken_constraint(REAL under_relaxation_factor) {

        return (sign(under_relaxation_factor) * std::min(std::abs(under_relaxation_factor),under_relaxation_factor_max_));

    }

    REAL calculate_aitken_constraint_pn_control(REAL under_relaxation_factor) {

        return (std::min(std::abs(under_relaxation_factor),under_relaxation_factor_max_));
    }

    REAL calculate_aitken_constraint_pnz_control(REAL under_relaxation_factor) const {

        return ((std::min(std::abs(under_relaxation_factor),under_relaxation_factor_max_)) < init_under_relaxation_factor_) ? init_under_relaxation_factor_ : (std::min(std::abs(under_relaxation_factor),under_relaxation_factor_max_));
    }

    template<typename OTYPE>
    OTYPE calculate_point_residual(std::pair<time_type,iterator_type> t, OTYPE filtered_value, OTYPE filtered_old_value) const {

        return (filtered_value - calculate_relaxed_value(t, filtered_value, filtered_old_value));

    }

    INT sign(REAL value)
    {
        return (value < 0) ? -1 : ((value > 0) ? 1 : 0);
    }

protected:
    const REAL init_under_relaxation_factor_;

    const REAL under_relaxation_factor_max_;

    MPI_Comm local_mpi_comm_world_;

    iterator_type minimum_iterator_;

    mutable std::vector<std::pair<std::pair<time_type,iterator_type>, REAL>> under_relaxation_factor_;

    mutable std::vector<std::pair<std::pair<time_type,iterator_type>,std::pair<INT, REAL>>> residual_l2_norm_;

    mutable std::vector<std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>>> pts_time_value_;

    mutable std::vector<std::pair<std::pair<time_type,iterator_type>,std::vector<std::pair<point_type, REAL>>>> pts_time_res_;

};

}

#endif /* MUI_COUPLING_ALGORITHM_AITKEN_H_ */
