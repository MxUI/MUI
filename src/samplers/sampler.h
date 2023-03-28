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
 * @file sampler.h
 * @author Y. H. Tang
 * @date 10 February 2014
 * @brief A reference file for making custom samplers. The new sampler does
 * not have to derive from this class, it just needs to implement all the
 * interfaces with the signatures specified.
 */

#ifndef MUI_SAMPLER_H_
#define MUI_SAMPLER_H_

#include "../config.h"
#include "../geometry/geometry.h"
#include "../storage/virtual_container.h"

namespace mui {
/*
class sampler {
public:
    sampler() {
        r = 1.0;
    }

    double filter( point<double,3> focus, const virtual_container<double> &data_points ) const {
        double sum = 0.;
        int n = 0;
        for( size_t i = 0 ; i < data_points.size() ; i++ ) {
            if ( ( focus - data_points[i].first ).norm() < r ) {
                sum += data_points[i].second;
                n++;
            }
        }
        if (n) return sum / n;
        else return 0;
    }
    inline span<> support() const {
        return span<>() || geometry::sphere<>( point<double,3>(0), r );
    }

protected:
    double r;
};
*/
}

#endif /* MUI_SAMPLER_H_ */
