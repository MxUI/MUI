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
 * @file uniface.cpp
 * @author C. Richardson, E. R. Fernandez, W. Liu
 * @date 11 March 2023
 * @brief Uniface calss for MUI Python wrapper.
 */

#include <mui.h>
#include <pybind11/chrono.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi4py/mpi4py.h>
#include <string>

#include "config_name.h"
#include "sampler_name.h"
#include "temporal_name.h"

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal>
void declare_uniface_fetch(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;

  std::string fetch_name = "fetch_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
              (T(Tclass::*)(
                  const std::string &, const mui::point<Treal, Tconfig::D> &, const Ttime,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &, bool)) &
                  Tclass::fetch,
              "")
      .def(fetch_name.c_str(),
           (T(Tclass::*)(
               const std::string &, const mui::point<Treal, Tconfig::D> &, const Ttime,
               const Titer,
               const Tsampler<Tconfig, T, T> &,
               const Ttemporal<Tconfig> &, bool)) &
               Tclass::fetch,
           "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal>
void declare_uniface_fetch_many(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;

  std::string fetch_many_name = "fetch_many_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_many_name.c_str(),
              (py::array_t<T, py::array::c_style>(Tclass::*)(
                  const std::string &attr,
                  const py::array_t<Treal, py::array::c_style> points, const Ttime t,
                  const Tsampler<Tconfig, T, T> &sampler,
                  const Ttemporal<Tconfig> &t_sampler)) &
                  Tclass::fetch_many,
              "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler>
void declare_uniface_fetch_all_temporal(py::class_<mui::uniface<Tconfig>> &uniface)
{
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);
}

template <typename Tconfig, typename T>
void declare_uniface_funcs(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;

  std::string push_name = "push_" + type_name<T>();
  std::string push_many_name = "push_many_" + type_name<T>();
  std::string fetch_name = "fetch_" + type_name<T>();
  uniface.def(push_name.c_str(),
              (void(Tclass::*)(const std::string &, const mui::point<Treal, Tconfig::D> &,
                               const T &)) &
                  Tclass::push,
              "")
      .def(push_name.c_str(),
           (void(Tclass::*)(const std::string &, const T &)) & Tclass::push,
           "")
      .def(fetch_name.c_str(),
           (T(Tclass::*)(const std::string &)) & Tclass::fetch, "");

  uniface.def(push_many_name.c_str(),
              (void(Tclass::*)(const std::string &attr, const py::array_t<Treal> &points,
                               const py::array_t<T> &values)) &
                  Tclass::push_many,
              "");
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_gauss>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_moving_average>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_nearest_neighbor>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_pseudo_n2_linear>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_pseudo_nearest_neighbor>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_shepard_quintic>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_sum_quintic>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_sph_quintic>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_exact>(uniface);
  declare_uniface_fetch_all_temporal<Tconfig, T, mui::sampler_rbf>(uniface);
}

template <typename Tconfig>
void declare_uniface_class(py::module &m)
{
  std::string name = "_create_uniface" + config_name<Tconfig>();

  m.def(name.c_str(), [](std::string domain, std::vector<std::string> interfaces, pybind11::handle const& world)
        {
			// Import mpi4py if it does not exist.
			if (!PyMPIComm_Get)
			{
				if (import_mpi4py() < 0)
				{
					throw std::runtime_error(
						"MUI Error [uniface.cpp]: mpi4py not loaded correctly\n");
				}
			}

			MPI_Comm comm;
			MPI_Comm ric_mpiComm;

			if (world.is_none()) {
				comm = MPI_COMM_WORLD;
				ric_mpiComm = reinterpret_cast<MPI_Comm>(comm);
			} else {
				PyObject *py_src = world.ptr();
				MPI_Comm *comm_p = PyMPIComm_Get(py_src);
				ric_mpiComm = reinterpret_cast<MPI_Comm>(*comm_p);
			}
            return mui::create_uniface<Tconfig>(domain, interfaces, ric_mpiComm);
        });

  m.def(name.c_str(), [](std::string domain, std::vector<std::string> interfaces)
        {
            return mui::create_uniface<Tconfig>(domain, interfaces, MPI_COMM_WORLD);
        });

  name = "_Uniface" + config_name<Tconfig>();
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;
  py::class_<Tclass> uniface(m, name.c_str());

  uniface
      .def("commit", (int(Tclass::*)(Ttime)) & Tclass::commit, "")
      .def("forecast", (void(Tclass::*)(Ttime)) & Tclass::forecast, "")
      //.def("is_ready", &Tclass::is_ready, "")
      .def("barrier", (void(Tclass::*)(Ttime)) & Tclass::barrier, "")
      .def("barrier", (void(Tclass::*)(Ttime, Titer)) & Tclass::barrier, "")
      .def("forget", (void(Tclass::*)(Ttime, bool)) & Tclass::forget, "")
      .def("forget", (void(Tclass::*)(Ttime, Ttime, bool)) & Tclass::forget, "")
      .def("set_memory", (void(Tclass::*)(Ttime)) & Tclass::set_memory, "")
      .def("announce_send_span",
           (void(Tclass::*)(Ttime, Ttime, mui::geometry::any_shape<Tconfig>,
                            bool synchronised)) &
               Tclass::announce_send_span,
           "")
      .def("announce_recv_span",
           (void(Tclass::*)(Ttime, Ttime, mui::geometry::any_shape<Tconfig>,
                            bool synchronised)) &
               Tclass::announce_recv_span,
           "")
      .def("announce_send_disable",
           (void(Tclass::*)()) & Tclass::announce_send_disable, "")
      .def("announce_recv_disable",
           (void(Tclass::*)()) & Tclass::announce_recv_disable, "")
      //    DEFINE_MUI_UNIFACE_FETCH_5ARGS() DEFINE_MUI_UNIFACE_FETCH_6ARGS()

      .def(py::init<const std::string &>());

  // Special cases for string
  uniface.def("push_string",
              (void(Tclass::*)(const std::string &, const mui::point<Treal, Tconfig::D> &,
                               const std::string &)) &
                  Tclass::push,
              "")
      .def("push_string",
           (void(Tclass::*)(const std::string &, const std::string &)) & Tclass::push,
           "")
      .def("fetch_string",
           (std::string(Tclass::*)(const std::string &)) & Tclass::fetch, "");
  declare_uniface_fetch<Tconfig, std::string, mui::sampler_exact, mui::temporal_sampler_exact>(uniface);

  declare_uniface_funcs<Tconfig, double>(uniface);
  declare_uniface_funcs<Tconfig, float>(uniface);
  //declare_uniface_funcs<Tconfig, std::int32_t>(uniface);
  //declare_uniface_funcs<Tconfig, std::int64_t>(uniface);
}

void uniface(py::module &m)
{
#ifdef PYTHON_INT_64
  declare_uniface_class<mui::mui_config_1dx>(m);
 // declare_uniface_class<mui::mui_config_2dx>(m);
  declare_uniface_class<mui::mui_config_3dx>(m);
//  declare_uniface_class<mui::mui_config_1fx>(m);
  declare_uniface_class<mui::mui_config_2fx>(m);
 // declare_uniface_class<mui::mui_config_3fx>(m);

#elif defined PYTHON_INT_32
  declare_uniface_class<mui::mui_config_1d>(m);
  //declare_uniface_class<mui::mui_config_2d>(m);
  declare_uniface_class<mui::mui_config_3d>(m);
  //declare_uniface_class<mui::mui_config_1f>(m);
  declare_uniface_class<mui::mui_config_2f>(m);
  //declare_uniface_class<mui::mui_config_3f>(m);

#else
#error PYTHON_INT_[32|64] not defined.
#endif
}