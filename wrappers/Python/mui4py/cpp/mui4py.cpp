/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2023 C. Richardson, E. R. Fernandez, W. Liu                  *
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
 * @file mui4py.cpp
 * @author C. Richardson, E. R. Fernandez, W. Liu
 * @date 20 January 2023
 * @brief Main c++ file for MUI Python wrapper.
 */

#include <mui.h>
#include <pybind11/chrono.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi4py/mpi4py.h>
#include <string>

#include "compiler_info.h"

// Declaration of other files
void geometry(py::module &m);
void sampler(py::module &m);
void temporal_sampler(py::module &m);
void algorithm(py::module &m);
void uniface1d(py::module &m);
void uniface2d(py::module &m);
void uniface3d(py::module &m);
void uniface1f(py::module &m);
void uniface2f(py::module &m);
void uniface3f(py::module &m);

PYBIND11_MODULE(mui4py_mod, m)
{
  m.doc() = "MUI bindings for Python.";

  // Expose numerical limits from C++
  m.attr("numeric_limits_real") = std::numeric_limits<double>::min();
  m.attr("numeric_limits_int") = std::numeric_limits<int>::min();
  m.attr("numeric_limits_uint") = std::numeric_limits<unsigned int>::lowest();

  geometry(m);
  sampler(m);
  temporal_sampler(m);
  algorithm(m);
  uniface1d(m);
  uniface2d(m);
  uniface3d(m);
  uniface1f(m);
  uniface2f(m);
  uniface3f(m);

  m.def("set_quiet", &mui::set_quiet, "");
  m.def(
      "mpi_split_by_app", []() -> py::handle
      {
    if (import_mpi4py() < 0)
      Py_RETURN_NONE;
    return PyMPIComm_New(mui::mpi_split_by_app()); },
      "");
  m.def("get_mpi_version", &get_mpi_version, "");
  m.def("get_compiler_config", &get_compiler_config, "");
  m.def("get_compiler_version", &get_compiler_version, "");
}
