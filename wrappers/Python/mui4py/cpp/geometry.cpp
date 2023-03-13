/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2023 C. Richardson                                           *
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
 * @file geometry.cpp
 * @author C. Richardson
 * @date 20 January 2023
 * @brief Geometries for MUI Python wrapper.
 */

#include <mui.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include "config_name.h"

// Geometry

template <typename Tconfig>
void declare_geometry(py::module &m)
{
    std::string name = "_geometry_Shape" + config_name<Tconfig>();
    py::class_<mui::geometry::shape<Tconfig>>(m, name.c_str());
    name = "_geometry_Any_shape" + config_name<Tconfig>();
    using Tclass = mui::geometry::any_shape<Tconfig>;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<const mui::geometry::shape<Tconfig> &>());
    name = "_geometry_Point" + config_name<Tconfig>();
    py::class_<mui::geometry::point<Tconfig>, mui::geometry::shape<Tconfig>>(m, name.c_str())
        .def(py::init<const mui::point<typename Tconfig::REAL, Tconfig::D> &>());
    name = "_geometry_Box" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    py::class_<mui::geometry::box<Tconfig>, mui::geometry::shape<Tconfig>>(m, name.c_str())
        .def(py::init<const mui::point<Treal, Tconfig::D> &,
                      const mui::point<Treal, Tconfig::D> &>());
    name = "_geometry_Sphere" + config_name<Tconfig>();
    py::class_<mui::geometry::sphere<Tconfig>, mui::geometry::shape<Tconfig>>(m, name.c_str())
        .def(py::init<const mui::point<typename Tconfig::REAL, Tconfig::D> &,
                      typename Tconfig::REAL>());

    py::implicitly_convertible<mui::geometry::shape<Tconfig>,
                               mui::geometry::any_shape<Tconfig>>();
}

void geometry(py::module &m)
{
#ifdef PYTHON_INT_64

    declare_geometry<mui::mui_config_1dx>(m);
    declare_geometry<mui::mui_config_2dx>(m);
    declare_geometry<mui::mui_config_3dx>(m);
    declare_geometry<mui::mui_config_1fx>(m);
    declare_geometry<mui::mui_config_2fx>(m);
    declare_geometry<mui::mui_config_3fx>(m);

#elif defined PYTHON_INT_32

    declare_geometry<mui::mui_config_1d>(m);
    declare_geometry<mui::mui_config_2d>(m);
    declare_geometry<mui::mui_config_3d>(m);
    declare_geometry<mui::mui_config_1f>(m);
    declare_geometry<mui::mui_config_2f>(m);
    declare_geometry<mui::mui_config_3f>(m);

#else
#error PYTHON_INT_[32|64] not defined.
#endif
}
