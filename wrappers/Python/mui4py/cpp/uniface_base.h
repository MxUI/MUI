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
 * @file uniface_base.h
 * @author C. Richardson, E. R. Fernandez, W. Liu
 * @date 11 March 2023
 * @brief Uniface base calss for MUI Python wrapper.
 */

#include <mui.h>
#include <pybind11/pybind11.h>
#include <mpi4py/mpi4py.h>

#include "config_name.h"
#include "sampler_name.h"
#include "temporal_name.h"
#include "algorithm_name.h"

template <typename Tconfig, typename T>
void declare_uniface_fetch_single(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;

  std::string fetch_name = "fetch_single_" + type_name<T>();
  uniface.def(fetch_name.c_str(), (T(Tclass::*)(const std::string &)) & Tclass::fetch, "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal>
void declare_uniface_fetch(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;

  std::string fetch_name = "fetch_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
              (T(Tclass::*)(
                  const std::string &,
				  const mui::point<Treal, Tconfig::D> &,
				  const Ttime,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &, bool)) &
                  Tclass::fetch,
              "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal>
void declare_uniface_fetch_dual(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;

  std::string fetch_name = "fetch_dual_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
           (T(Tclass::*)(
               const std::string &,
			   const mui::point<Treal, Tconfig::D> &,
			   const Ttime,
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
                  const std::string &,
                  const py::array_t<Treal, py::array::c_style>,
				  const Ttime,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &, bool)) &
                  Tclass::fetch_many,
              "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal>
void declare_uniface_fetch_many_dual(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;

  std::string fetch_many_name = "fetch_many_dual_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_many_name.c_str(),
              (py::array_t<T, py::array::c_style>(Tclass::*)(
                  const std::string &,
                  const py::array_t<Treal, py::array::c_style>,
				  const Ttime,
				  const Titer,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &, bool)) &
                  Tclass::fetch_many,
              "");
}


template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal, template <typename> class Talgorithm>
void declare_uniface_fetch_algo(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;

  std::string fetch_name = "fetch_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>() + "_" + algorithm_name<Tconfig, Talgorithm>();
  uniface.def(fetch_name.c_str(),
              (T(Tclass::*)(
                  const std::string &,
				  const mui::point<Treal, Tconfig::D> &,
				  const Ttime,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &,
                  const Talgorithm<Tconfig> &,
				  bool)) &
                  Tclass::fetch,
              "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal, template <typename> class Talgorithm>
void declare_uniface_fetch_algo_dual(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;

  std::string fetch_name = "fetch_dual_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>() + "_" + algorithm_name<Tconfig, Talgorithm>();
  uniface.def(fetch_name.c_str(),
           (T(Tclass::*)(
               const std::string &,
			   const mui::point<Treal, Tconfig::D> &,
			   const Ttime,
               const Titer,
               const Tsampler<Tconfig, T, T> &,
               const Ttemporal<Tconfig> &,
			   const Talgorithm<Tconfig> &,
			   bool)) &
               Tclass::fetch,
           "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal, template <typename> class Talgorithm>
void declare_uniface_fetch_algo_many(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;

  std::string fetch_many_name = "fetch_many_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>() + "_" + algorithm_name<Tconfig, Talgorithm>();
  uniface.def(fetch_many_name.c_str(),
              (py::array_t<T, py::array::c_style>(Tclass::*)(
                  const std::string &,
                  const py::array_t<Treal, py::array::c_style>,
				  const Ttime,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &,
				  const Talgorithm<Tconfig> &, bool)) &
                  Tclass::fetch_many,
              "");
}

template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal, template <typename> class Talgorithm>
void declare_uniface_fetch_algo_many_dual(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;

  std::string fetch_many_name = "fetch_many_dual_" + type_name<T>() + "_" + sampler_name<Tconfig, T, Tsampler>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>() + "_" + algorithm_name<Tconfig, Talgorithm>();
  uniface.def(fetch_many_name.c_str(),
              (py::array_t<T, py::array::c_style>(Tclass::*)(
                  const std::string &,
                  const py::array_t<Treal, py::array::c_style>,
				  const Ttime,
				  const Titer,
                  const Tsampler<Tconfig, T, T> &,
                  const Ttemporal<Tconfig> &,
				  const Talgorithm<Tconfig> &, bool)) &
                  Tclass::fetch_many,
              "");
}

template <typename Tconfig, typename T, template <typename> class Ttemporal>
void declare_uniface_fetch_points(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Ttime = typename Tconfig::time_type;
  using Treal = typename Tconfig::REAL;

  std::string fetch_name = "fetch_points_" + type_name<T>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
              (py::array_t<Treal, py::array::c_style>(Tclass::*)(
                  const std::string &,
				  const Ttime,
                  const Ttemporal<Tconfig> &,
				  bool,
				  T)) &
                  Tclass::fetch_points_np,
              "");
}

template <typename Tconfig, typename T, template <typename> class Ttemporal>
void declare_uniface_fetch_points_dual(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;
  using Treal = typename Tconfig::REAL;

  std::string fetch_name = "fetch_points_dual_" + type_name<T>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
              (py::array_t<Treal, py::array::c_style>(Tclass::*)(
                  const std::string &,
				  const Ttime,
				  const Titer,
                  const Ttemporal<Tconfig> &,
				  bool,
				  T)) &
                  Tclass::fetch_points_np,
              "");
}

template <typename Tconfig, typename T, template <typename> class Ttemporal>
void declare_uniface_fetch_values(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Ttime = typename Tconfig::time_type;

  std::string fetch_name = "fetch_values_" + type_name<T>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
              (std::vector<T>(Tclass::*)(
                  const std::string &,
				  const Ttime,
                  const Ttemporal<Tconfig> &,
				  bool)) &
                  Tclass::fetch_values,
              "");
}

template <typename Tconfig, typename T, template <typename> class Ttemporal>
void declare_uniface_fetch_values_dual(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;

  std::string fetch_name = "fetch_values_dual_" + type_name<T>() + "_" + temporal_sampler_name<Tconfig, Ttemporal>();
  uniface.def(fetch_name.c_str(),
              (std::vector<T>(Tclass::*)(
                  const std::string &,
				  const Ttime,
				  const Titer,
                  const Ttemporal<Tconfig> &,
				  bool)) &
                  Tclass::fetch_values,
              "");
}


template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler, template <typename> class Ttemporal>
void declare_uniface_fetch_all_algorithm(py::class_<mui::uniface<Tconfig>> &uniface)
{

  declare_uniface_fetch_algo<Tconfig, T, Tsampler, Ttemporal, mui::algo_fixed_relaxation>(uniface);
  declare_uniface_fetch_algo<Tconfig, T, Tsampler, Ttemporal, mui::algo_aitken>(uniface);

  declare_uniface_fetch_algo_dual<Tconfig, T, Tsampler, Ttemporal, mui::algo_fixed_relaxation>(uniface);
  declare_uniface_fetch_algo_dual<Tconfig, T, Tsampler, Ttemporal, mui::algo_aitken>(uniface);

  declare_uniface_fetch_algo_many<Tconfig, T, Tsampler, Ttemporal, mui::algo_fixed_relaxation>(uniface);
  declare_uniface_fetch_algo_many<Tconfig, T, Tsampler, Ttemporal, mui::algo_aitken>(uniface);

  declare_uniface_fetch_algo_many_dual<Tconfig, T, Tsampler, Ttemporal, mui::algo_fixed_relaxation>(uniface);
  declare_uniface_fetch_algo_many_dual<Tconfig, T, Tsampler, Ttemporal, mui::algo_aitken>(uniface);
}


template <typename Tconfig, typename T, template <typename, typename, typename> class Tsampler>
void declare_uniface_fetch_all_temporal(py::class_<mui::uniface<Tconfig>> &uniface)
{

  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_dual<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_dual<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_dual<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_dual<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_many<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_many_dual<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_many_dual<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_many_dual<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_many_dual<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_all_algorithm<Tconfig, T, Tsampler, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_all_algorithm<Tconfig, T, Tsampler, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_all_algorithm<Tconfig, T, Tsampler, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_all_algorithm<Tconfig, T, Tsampler, mui::temporal_sampler_sum>(uniface);
}

template <typename Tconfig, typename T>
void declare_uniface_fetch_points_values_temporal(py::class_<mui::uniface<Tconfig>> &uniface)
{

  declare_uniface_fetch_points<Tconfig, T, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_points<Tconfig, T, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_points<Tconfig, T, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_points<Tconfig, T, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_points_dual<Tconfig, T, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_points_dual<Tconfig, T, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_points_dual<Tconfig, T, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_points_dual<Tconfig, T, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_values<Tconfig, T, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_values<Tconfig, T, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_values<Tconfig, T, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_values<Tconfig, T, mui::temporal_sampler_sum>(uniface);

  declare_uniface_fetch_values_dual<Tconfig, T, mui::temporal_sampler_exact>(uniface);
  declare_uniface_fetch_values_dual<Tconfig, T, mui::temporal_sampler_gauss>(uniface);
  declare_uniface_fetch_values_dual<Tconfig, T, mui::temporal_sampler_mean>(uniface);
  declare_uniface_fetch_values_dual<Tconfig, T, mui::temporal_sampler_sum>(uniface);
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
void declare_uniface_string(py::class_<mui::uniface<Tconfig>> &uniface)
{
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;

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

}

template <typename Tconfig>
void declare_uniface_class (py::module &m)
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
						"MUI Error [wrappers/Python/mui4py/cpp/uniface.cpp]: mpi4py not loaded correctly\n");
				}
			}

			MPI_Comm ric_mpiComm;

			if (world.is_none()) {
				MPI_Comm comm = MPI_COMM_WORLD;
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
  using Ttime = typename Tconfig::time_type;
  using Titer = typename Tconfig::iterator_type;
  py::class_<Tclass> uniface(m, name.c_str());

  uniface
      .def("commit", (int(Tclass::*)(Ttime, Titer)) & Tclass::commit, "")
      .def("forecast", (void(Tclass::*)(Ttime, Titer)) & Tclass::forecast, "")
      .def("is_ready", (bool(Tclass::*)(const std::string &, Ttime) const) &Tclass::is_ready, "")
      .def("is_ready", (bool(Tclass::*)(const std::string &, Ttime, Titer) const) &Tclass::is_ready, "")
      .def("barrier", (void(Tclass::*)(Ttime)) & Tclass::barrier, "")
      .def("barrier", (void(Tclass::*)(Ttime, Titer)) & Tclass::barrier, "")
      .def("forget", (void(Tclass::*)(Ttime, bool)) & Tclass::forget, "")
      .def("forget", (void(Tclass::*)(std::pair<Ttime,Titer>, bool)) & Tclass::forget, "")
      .def("forget", (void(Tclass::*)(Ttime, Ttime, bool)) & Tclass::forget, "")
      .def("forget", (void(Tclass::*)(std::pair<Ttime,Titer>, std::pair<Ttime,Titer>, bool)) & Tclass::forget, "")
      .def("set_memory", (void(Tclass::*)(Ttime)) & Tclass::set_memory, "")
      .def("uri_host", (std::string(Tclass::*)()) & Tclass::uri_host, "")
      .def("uri_path", (std::string(Tclass::*)()) & Tclass::uri_path, "")
      .def("uri_protocol", (std::string(Tclass::*)()) & Tclass::uri_protocol, "")
      .def("announce_send_span",
           (void(Tclass::*)(Ttime, Ttime, mui::geometry::any_shape<Tconfig>,
                            bool)) & Tclass::announce_send_span, "")
      .def("announce_recv_span",
           (void(Tclass::*)(Ttime, Ttime, mui::geometry::any_shape<Tconfig>,
                            bool)) & Tclass::announce_recv_span, "")
      .def("announce_send_disable",
           (void(Tclass::*)(bool)) & Tclass::announce_send_disable, "")
      .def("announce_recv_disable",
           (void(Tclass::*)(bool)) & Tclass::announce_recv_disable, "")
      .def("update_smart_send",
           (void(Tclass::*)(Ttime)) & Tclass::update_smart_send, "")
      .def("barrier_ss_send",
           (void(Tclass::*)()) & Tclass::barrier_ss_send, "")
      .def("barrier_ss_recv",
           (void(Tclass::*)()) & Tclass::barrier_ss_recv, "")
      .def(py::init<const std::string &>());

  declare_uniface_string<Tconfig> (uniface);

  declare_uniface_funcs<Tconfig, double>(uniface);
  declare_uniface_funcs<Tconfig, float>(uniface);
  declare_uniface_funcs<Tconfig, std::int32_t>(uniface);
  declare_uniface_funcs<Tconfig, std::int64_t>(uniface);

  declare_uniface_fetch_single<Tconfig, double>(uniface);
  declare_uniface_fetch_single<Tconfig, float>(uniface);
  declare_uniface_fetch_single<Tconfig, std::int32_t>(uniface);
  declare_uniface_fetch_single<Tconfig, std::int64_t>(uniface);

  declare_uniface_fetch_points_values_temporal<Tconfig, double>(uniface);
  declare_uniface_fetch_points_values_temporal<Tconfig, float>(uniface);
  declare_uniface_fetch_points_values_temporal<Tconfig, std::int32_t>(uniface);
  declare_uniface_fetch_points_values_temporal<Tconfig, std::int64_t>(uniface);
}
