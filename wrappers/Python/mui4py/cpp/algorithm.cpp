/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2023 W. Liu                                                  *
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
 * @file algorithm.cpp
 * @author W. Liu
 * @date 18 March 2023
 * @brief Algorithms for MUI Python wrapper.
 */

#include <mui.h>
#include <pybind11/pybind11.h>
#include <mpi4py/mpi4py.h>
#include <string>
#include "config_name.h"

template <typename Tconfig>
void declare_algorithms(py::module &m)
{
    std::string name = "_Algorithm_fixed_relaxation" + config_name<Tconfig>();
    using TclassFR = mui::algo_fixed_relaxation<Tconfig>;
    using Treal = typename Tconfig::REAL;
    using Tpoint = typename Tconfig::point_type;
    using Ttime = typename Tconfig::time_type;
    using Titer = typename Tconfig::iterator_type;
    py::class_<TclassFR> algo_fixed_relaxation(m, name.c_str());
    algo_fixed_relaxation.def(py::init([](Treal under_relaxation_factor = 1.0,
                                pybind11::handle const& local_comm = py::none(),
                                std::vector<std::pair<Tpoint, Treal>> pts_vlu_init = std::vector<std::pair<Tpoint, Treal>>()) {

                // Import mpi4py if it does not exist.
                if (!PyMPIComm_Get)
                {
                    if (import_mpi4py() < 0)
                    {
                        throw std::runtime_error(
                            "MUI Error [wrappers/Python/mui4py/cpp/algorithm.cpp]: mpi4py not loaded correctly\n");
                    }
                }

                MPI_Comm ric_mpiComm;

                if (local_comm.is_none()) {
                    MPI_Comm comm = MPI_COMM_NULL;
                    ric_mpiComm = reinterpret_cast<MPI_Comm>(comm);
                } else {
                    PyObject *py_src = local_comm.ptr();
                    MPI_Comm *comm_p = PyMPIComm_Get(py_src);
                    ric_mpiComm = reinterpret_cast<MPI_Comm>(*comm_p);
                }

                return new mui::algo_fixed_relaxation<Tconfig>(under_relaxation_factor,
                                                     ric_mpiComm,
                                                     pts_vlu_init);
    }))
               .def("get_under_relaxation_factor", (Treal(TclassFR::*)(Ttime)) &TclassFR::get_under_relaxation_factor, "")
               .def("get_under_relaxation_factor", (Treal(TclassFR::*)(Ttime, Titer)) &TclassFR::get_under_relaxation_factor, "")
               .def("get_residual_L2_Norm", (Treal(TclassFR::*)(Ttime)) &TclassFR::get_residual_L2_Norm, "")
               .def("get_residual_L2_Norm", (Treal(TclassFR::*)(Ttime, Titer)) &TclassFR::get_residual_L2_Norm, "");

    name = "_Algorithm_aitken" + config_name<Tconfig>();
    using TclassAitken = mui::algo_aitken<Tconfig>;
    using Treal = typename Tconfig::REAL;
    using Tpoint = typename Tconfig::point_type;
    using Ttime = typename Tconfig::time_type;
    using Titer = typename Tconfig::iterator_type;
    py::class_<TclassAitken> algo_aitken(m, name.c_str());
    algo_aitken.def(py::init([](Treal under_relaxation_factor = 1.0,
                                Treal under_relaxation_factor_max = 1.0,
                                pybind11::handle const& local_comm = py::none(),
                                std::vector<std::pair<Tpoint, Treal>> pts_vlu_init = std::vector<std::pair<Tpoint, Treal>>(),
                                Treal res_l2_norm_nm1 = 0.0) {

                // Import mpi4py if it does not exist.
                if (!PyMPIComm_Get)
                {
                    if (import_mpi4py() < 0)
                    {
                        throw std::runtime_error(
                            "MUI Error [wrappers/Python/mui4py/cpp/algorithm.cpp]: mpi4py not loaded correctly\n");
                    }
                }

                MPI_Comm ric_mpiComm;

                if (local_comm.is_none()) {
                    MPI_Comm comm = MPI_COMM_NULL;
                    ric_mpiComm = reinterpret_cast<MPI_Comm>(comm);
                } else {
                    PyObject *py_src = local_comm.ptr();
                    MPI_Comm *comm_p = PyMPIComm_Get(py_src);
                    ric_mpiComm = reinterpret_cast<MPI_Comm>(*comm_p);
                }

                return new mui::algo_aitken<Tconfig>(under_relaxation_factor,
                                                     under_relaxation_factor_max,
                                                     ric_mpiComm,
                                                     pts_vlu_init,
                                                     res_l2_norm_nm1);
    }))
               .def("get_under_relaxation_factor", (Treal(TclassAitken::*)(Ttime)) &TclassAitken::get_under_relaxation_factor, "")
               .def("get_under_relaxation_factor", (Treal(TclassAitken::*)(Ttime, Titer)) &TclassAitken::get_under_relaxation_factor, "")
               .def("get_residual_L2_Norm", (Treal(TclassAitken::*)(Ttime)) &TclassAitken::get_residual_L2_Norm, "")
               .def("get_residual_L2_Norm", (Treal(TclassAitken::*)(Ttime, Titer)) &TclassAitken::get_residual_L2_Norm, "");
}

void algorithm(py::module &m)
{
#ifdef PYTHON_INT_64
    declare_algorithms<mui::mui_config_1dx>(m);
    declare_algorithms<mui::mui_config_2dx>(m);
    declare_algorithms<mui::mui_config_3dx>(m);
    declare_algorithms<mui::mui_config_1fx>(m);
    declare_algorithms<mui::mui_config_2fx>(m);
    declare_algorithms<mui::mui_config_3fx>(m);
#elif defined PYTHON_INT_32
    declare_algorithms<mui::mui_config_1d>(m);
    declare_algorithms<mui::mui_config_2d>(m);
    declare_algorithms<mui::mui_config_3d>(m);
    declare_algorithms<mui::mui_config_1f>(m);
    declare_algorithms<mui::mui_config_2f>(m);
    declare_algorithms<mui::mui_config_3f>(m);
#else
#error PYTHON_INT_[32|64] not defined.
#endif
}
