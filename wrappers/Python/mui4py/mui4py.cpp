/*
Multiscale Universal Interface Code Coupling Library

Copyright (C) 2017 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis,
                   E. R. Fernandez, S. M. Longshaw, A. Skillen

This software is jointly licensed under the Apache License, Version 2.0
and the GNU General Public License version 3, you may use it according
to either.

** Apache License, version 2.0 **

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

** GNU General Public License, version 3 **

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

** File Details **

Filename: mui4py.cpp
Created: Oct 28, 2018
Author: Eduardo Ramos Fernandez, W. Liu
Description: MUI Python bindings
*/


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../../../mui.h"
#include "../../../uniface.h"
#include <string>
#include <limits>
#include <sstream>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <mpi4py/mpi4py.h>
#include <assert.h> 

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

// Make a FOREACH macro
#define FE_1(WHAT, X) WHAT(X) 
#define FE_2(WHAT, X, ...) WHAT(X)FE_1(WHAT, __VA_ARGS__)
#define FE_3(WHAT, X, ...) WHAT(X)FE_2(WHAT, __VA_ARGS__)
#define FE_4(WHAT, X, ...) WHAT(X)FE_3(WHAT, __VA_ARGS__)
#define FE_5(WHAT, X, ...) WHAT(X)FE_4(WHAT, __VA_ARGS__)
#define FE_6(WHAT, X, ...) WHAT(X)FE_5(WHAT, __VA_ARGS__)

#define GET_MACRO(_1,_2,_3,_4,_5,_6,NAME,...) NAME 
#define FOR_EACH(action,...) \
  GET_MACRO(__VA_ARGS__,FE_6,FE_5,FE_4,FE_3,FE_2,FE_1,)(action,__VA_ARGS__)

using std::string;

#define PUSH_INSTANCE_SINGLE(IO_TYPE)\
    .def("push_" STRINGIFY(IO_TYPE), (void (Tclass::*)(const string&, const mui::point<Treal,Tconfig::D>,\
         const IO_TYPE)) &Tclass::push, "") \
    .def("push_" STRINGIFY(IO_TYPE), (void (Tclass::*)(const string&, const IO_TYPE)) &Tclass::push, "")

#define PUSH_INSTANCE_MANY(IO_TYPE)\
     .def("push_many_" STRINGIFY(IO_TYPE), (void (Tclass::*)(const string& attr, const py::array_t<Treal>& points, const py::array_t<IO_TYPE>& values)) &Tclass::push_many, "") 

//Temporarily disabled the fetch_point function. Bind the variadic template later.
//#define FETCH_POINTS_INSTANCE(IO_TYPE)\
    .def("fetch_points_" STRINGIFY(IO_TYPE), &Tclass::template fetch_points_np<IO_TYPE>) 

//#define DEFINE_MUI_UNIFACE_FETCH_POINTS() \
    FOR_EACH(FETCH_POINTS_INSTANCE, double, float, int64_t, int32_t, string)

#define FETCH_INSTANCE_1ARG(IO_TYPE) \
    .def("fetch_" STRINGIFY(IO_TYPE), (IO_TYPE (Tclass::*)(const string&)) &Tclass::fetch, "")

#define DEFINE_MUI_UNIFACE_FETCH_1ARG() \
    FOR_EACH(FETCH_INSTANCE_1ARG, double, float, int64_t, int32_t, string)

#define DEFINE_MUI_UNIFACE_PUSH() \
    FOR_EACH(PUSH_INSTANCE_SINGLE, double, float, int64_t, int32_t, string)\
    FOR_EACH(PUSH_INSTANCE_MANY, double, float, int64_t, int32_t)

#define FETCH_INSTANCE_SINGLE(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
   .def(replace_str("fetch_" STRINGIFY(IO_TYPE) "_" STRINGIFY(SPATIAL_SAMPLER) "_" STRINGIFY(CHRONO_SAMPLER),"_sampler", "").c_str(),\
        (IO_TYPE (Tclass::*)(const string&, const mui::point<Treal,Tconfig::D>, const Ttime,\
        const mui::SPATIAL_SAMPLER<Tconfig,IO_TYPE,IO_TYPE>&, const mui::CHRONO_SAMPLER<Tconfig>&, bool)) &Tclass::fetch, "") \

#define FETCH_INSTANCE_MANY(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
   .def(replace_str("fetch_many_" STRINGIFY(IO_TYPE) "_" STRINGIFY(SPATIAL_SAMPLER) "_" STRINGIFY(CHRONO_SAMPLER),"_sampler", "").c_str(),\
        (py::array_t<IO_TYPE,py::array::c_style> (Tclass::*) (const string& attr,const py::array_t<Treal,py::array::c_style> points, const Ttime t,\
        const mui::SPATIAL_SAMPLER<Tconfig,IO_TYPE,IO_TYPE>& sampler, const mui::CHRONO_SAMPLER<Tconfig>& t_sampler)) &Tclass::fetch_many, "")

#define FETCH_INSTANCE_SINGLE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
   .def(replace_str("fetch6_" STRINGIFY(IO_TYPE) "_" STRINGIFY(SPATIAL_SAMPLER) "_" STRINGIFY(CHRONO_SAMPLER),"_sampler", "").c_str(),\
        (IO_TYPE (Tclass::*)(const string&, const mui::point<Treal,Tconfig::D>, const Ttime, const Ttime,\
        const mui::SPATIAL_SAMPLER<Tconfig,IO_TYPE,IO_TYPE>&, const mui::CHRONO_SAMPLER<Tconfig>&, bool)) &Tclass::fetch, "")

/* #define FETCH_INSTANCE_MANY6(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
   .def(replace_str("fetch_many6_" STRINGIFY(IO_TYPE) "_" STRINGIFY(SPATIAL_SAMPLER) "_" STRINGIFY(CHRONO_SAMPLER),"_sampler", "").c_str(),\
        (py::array_t<IO_TYPE,py::array::c_style> (Tclass::*) (const string& attr,const py::array_t<Treal,py::array::c_style> points, const Ttime t1, const Ttime t2\
        const mui::SPATIAL_SAMPLER<Tconfig,IO_TYPE,IO_TYPE>& sampler, const mui::CHRONO_SAMPLER<Tconfig>& t_sampler)) &Tclass::fetch_many6, "") */

#define FETCH_INSTANCE(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
    FETCH_INSTANCE_MANY(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
    FETCH_INSTANCE_SINGLE(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE)

#define FETCH_INSTANCE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE) \
    FETCH_INSTANCE_SINGLE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,IO_TYPE)


#define FETCH_SAMPLER_EXACT(CHRONO_SAMPLER) \
       FETCH_INSTANCE(sampler_exact,CHRONO_SAMPLER,double) \
       FETCH_INSTANCE(sampler_exact,CHRONO_SAMPLER,float) \
       FETCH_INSTANCE(sampler_exact,CHRONO_SAMPLER,int64_t) \
       FETCH_INSTANCE(sampler_exact,CHRONO_SAMPLER,int32_t)\
       FETCH_INSTANCE_SINGLE(sampler_exact,CHRONO_SAMPLER,string)
       // FETCH_INSTANCE(sampler_exact,CHRONO_SAMPLER,bool)

#define FETCH_SAMPLER_EXACT6(CHRONO_SAMPLER) \
       FETCH_INSTANCE6(sampler_exact,CHRONO_SAMPLER,double) \
       FETCH_INSTANCE6(sampler_exact,CHRONO_SAMPLER,float) \
       FETCH_INSTANCE6(sampler_exact,CHRONO_SAMPLER,int64_t) \
       FETCH_INSTANCE6(sampler_exact,CHRONO_SAMPLER,int32_t)\
       FETCH_INSTANCE_SINGLE6(sampler_exact,CHRONO_SAMPLER,string)
       // FETCH_INSTANCE6(sampler_exact,CHRONO_SAMPLER,bool)

#define FETCH_SAMPLER_NUMERICAL(SPATIAL_SAMPLER,CHRONO_SAMPLER) \
       FETCH_INSTANCE(SPATIAL_SAMPLER,CHRONO_SAMPLER,double) \
       FETCH_INSTANCE(SPATIAL_SAMPLER,CHRONO_SAMPLER,float) \
       FETCH_INSTANCE(SPATIAL_SAMPLER,CHRONO_SAMPLER,int64_t) \
       FETCH_INSTANCE(SPATIAL_SAMPLER,CHRONO_SAMPLER,int32_t) 

#define FETCH_SAMPLER_NUMERICAL6(SPATIAL_SAMPLER,CHRONO_SAMPLER) \
       FETCH_INSTANCE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,double) \
       FETCH_INSTANCE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,float) \
       FETCH_INSTANCE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,int64_t) \
       FETCH_INSTANCE6(SPATIAL_SAMPLER,CHRONO_SAMPLER,int32_t) 

#ifdef USE_RBF
#define EXPAND_FETCH_NUMERICAL(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_exact,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_moving_average,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_rbf,CHRONO_SAMPLER)
	   
#define EXPAND_FETCH_NUMERICAL6(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_exact,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_moving_average,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_rbf,CHRONO_SAMPLER)

#define EXPAND_FETCH_EXACT(CHRONO_SAMPLER) \
       FETCH_SAMPLER_EXACT(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_moving_average,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_rbf,CHRONO_SAMPLER)

#define EXPAND_FETCH_EXACT6(CHRONO_SAMPLER) \
       FETCH_SAMPLER_EXACT6(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_moving_average,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_rbf,CHRONO_SAMPLER)
#else

#define EXPAND_FETCH_NUMERICAL(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_exact,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_moving_average,CHRONO_SAMPLER)

#define EXPAND_FETCH_NUMERICAL6(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_exact,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_moving_average,CHRONO_SAMPLER)

#define EXPAND_FETCH_EXACT(CHRONO_SAMPLER) \
       FETCH_SAMPLER_EXACT(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL(sampler_moving_average,CHRONO_SAMPLER)

#define EXPAND_FETCH_EXACT6(CHRONO_SAMPLER) \
       FETCH_SAMPLER_EXACT6(CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_gauss,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest2_linear,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_pseudo_nearest_neighbor,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_shepard_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sph_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_sum_quintic,CHRONO_SAMPLER) \
       FETCH_SAMPLER_NUMERICAL6(sampler_moving_average,CHRONO_SAMPLER)

#endif

#define DEFINE_MUI_UNIFACE_FETCH_5ARGS() \
       EXPAND_FETCH_EXACT(chrono_sampler_exact) \
       EXPAND_FETCH_NUMERICAL(chrono_sampler_gauss) \
       EXPAND_FETCH_NUMERICAL(chrono_sampler_mean) \
       EXPAND_FETCH_NUMERICAL(chrono_sampler_sum)

#define DEFINE_MUI_UNIFACE_FETCH_6ARGS() \
       EXPAND_FETCH_EXACT6(chrono_sampler_exact) \
       EXPAND_FETCH_NUMERICAL6(chrono_sampler_gauss) \
       EXPAND_FETCH_NUMERICAL6(chrono_sampler_mean) \
       EXPAND_FETCH_NUMERICAL6(chrono_sampler_sum)


#define DECLARE_MUI_CPP2PY_CLASSES_0ARG(FUNCNAME,CLASSNAME)	\
    DECLARE_MUI_CPP2PY_CLASSES_1ARG(FUNCNAME,CLASSNAME,void)


#ifdef PYTHON_INT_64
#define DECLARE_MUI_CPP2PY_CLASSES_1ARG(FUNCNAME,CLASSNAME,ARG1)	\
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_1dx,ARG1>(m, #CLASSNAME, string("1d_f64_i64"), #ARG1); \
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_2dx,ARG1>(m, #CLASSNAME, string("2d_f64_i64"), #ARG1); \
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_3dx,ARG1>(m, #CLASSNAME, string("3d_f64_i64"), #ARG1); \
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_1fx,ARG1>(m, #CLASSNAME, string("1d_f32_i64"), #ARG1); \
    declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_2fx,ARG1>(m, #CLASSNAME, string("2d_f32_i64"), #ARG1); \
    declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_3fx,ARG1>(m, #CLASSNAME, string("3d_f32_i64"), #ARG1);
#elif defined PYTHON_INT_32
#define DECLARE_MUI_CPP2PY_CLASSES_1ARG(FUNCNAME,CLASSNAME,ARG1)	\
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_1d,ARG1>(m, #CLASSNAME, string("1d_f64_i32"), #ARG1); \
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_2d,ARG1>(m, #CLASSNAME, string("2d_f64_i32"), #ARG1); \
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_3d,ARG1>(m, #CLASSNAME, string("3d_f64_i32"), #ARG1); \
	declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_1f,ARG1>(m, #CLASSNAME, string("1d_f32_i32"), #ARG1); \
    declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_2f,ARG1>(m, #CLASSNAME, string("2d_f32_i32"), #ARG1); \
    declare_##FUNCNAME<mui::CLASSNAME,mui::mui_config_3f,ARG1>(m, #CLASSNAME, string("3d_f32_i32"), #ARG1);
#else
#error PYTHON_INT_[32|64] not defined.
#endif

 #define DECLARE_FUNC_HEADER(NAME) \
    void declare_##NAME(py::module &m, const string &name, const string &typestr, const string &arg1)
 

namespace py = pybind11;

std::string replace_str(std::string data, const std::string &toSearch, const std::string &replaceStr) {
	// Get the first occurrence
	size_t pos = data.find(toSearch);
	// Repeat till end is reached
	while(pos != std::string::npos) {
		data.replace(pos, toSearch.size(), replaceStr);
		pos = data.find(toSearch, pos + toSearch.size());
	}
    return data;
}

const string get_pyclass_name(string name, const string &typestr, const string &arg1 = "") {
   std::size_t pos = name.find("::");
   string sout;
   if (pos != string::npos) {
     sout = name.substr(pos + 2);
     sout[0] = std::toupper(sout[0]);
     sout = "_" + name.substr(0, pos) + "_" + sout + typestr;
   }
   else {
     sout = string(name);
     sout[0] = std::toupper(sout[0]);
     sout = "_" + sout + typestr;
   }
   if (arg1 != "")
       sout = sout + "_" + arg1;
   return sout;
}

// *** GEOMETRY CLASSES *** //

// SHAPE //
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(geometry_shape) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    py::class_<Tclass>(m, pyclass_name.c_str());
    // .def("__add__", (void (Tclass::*)()) &Tclass::forget, "")
    // .def("__add__", [](mui::geometry::any_shape<Tconfig>, mui::geometry::any_shape<Tconfig>));
}

// ANY_SHAPE //
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(geometry_any_shape) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<const mui::geometry::shape<Tconfig>&>());

    // .def("__add__", (void (Tclass::*)()) &Tclass::forget, "")
    // .def("__add__", [](mui::geometry::any_shape<Tconfig>, mui::geometry::any_shape<Tconfig>));
}

//
// OR_SET //
// template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
// DECLARE_FUNC_HEADER(geometry_or_set) {
//     string pyclass_name = get_pyclass_name(name, typestr);
//     using Tclass = TclassTemplate<Tconfig>;
//     py::class_<Tclass, mui::geometry::any_shape<Tconfig>>(m, pyclass_name.c_str())
//     // .def("bbox", 
//     .def(py::init<mui::geometry::any_shape<Tconfig>, mui::geometry::any_shape<Tconfig>&>());
// }
//
// POINT //
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(geometry_point) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    py::class_<Tclass, mui::geometry::shape<Tconfig>>(m, pyclass_name.c_str())
    .def(py::init<const mui::point<typename Tconfig::REAL,Tconfig::D>&>());
}

// SPHERE //
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(geometry_sphere) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    py::class_<Tclass, mui::geometry::shape<Tconfig>>(m, pyclass_name.c_str())
    .def(py::init<const mui::point<typename Tconfig::REAL,Tconfig::D>&, typename Tconfig::REAL>());
}

// BOX //
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(geometry_box) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    using Treal = typename Tconfig::REAL;
    py::class_<Tclass, mui::geometry::shape<Tconfig>>(m, pyclass_name.c_str())
    .def(py::init<const mui::point<Treal,Tconfig::D>&,\
		  const mui::point<Treal,Tconfig::D>&>());
}

// *** UNIFACE CLASS *** //
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(uniface) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    using Treal = typename Tconfig::REAL;
    using Ttime = typename Tconfig::time_type;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def("commit", (int (Tclass::*)(Ttime)) &Tclass::commit, "")
    .def("forecast", (void (Tclass::*)(Ttime)) &Tclass::forecast, "")
    //.def("is_ready", &Tclass::is_ready, "")
    .def("barrier", (void (Tclass::*)(Ttime)) &Tclass::barrier, "")
    .def("barrier", (void (Tclass::*)(Ttime, Ttime)) &Tclass::barrier, "")
    .def("forget", (void (Tclass::*)(Ttime, bool)) &Tclass::forget, "")
    .def("forget", (void (Tclass::*)(Ttime, Ttime, bool)) &Tclass::forget, "")
    .def("set_memory", (void (Tclass::*)(Ttime)) &Tclass::set_memory, "")
    .def("announce_send_span", (void (Tclass::*)(Ttime, Ttime, mui::geometry::any_shape<Tconfig>))\
           &Tclass::announce_send_span,"") 
    .def("announce_recv_span", (void (Tclass::*)(Ttime, Ttime, mui::geometry::any_shape<Tconfig>))\
           &Tclass::announce_recv_span,"")
    DEFINE_MUI_UNIFACE_PUSH()
    DEFINE_MUI_UNIFACE_FETCH_1ARG()
    DEFINE_MUI_UNIFACE_FETCH_5ARGS()
    DEFINE_MUI_UNIFACE_FETCH_6ARGS()
//Temporarily disabled the fetch_point function. Bind the variadic template later.
//    DEFINE_MUI_UNIFACE_FETCH_POINTS()
    .def(py::init<const string &>());
    py::implicitly_convertible<mui::geometry::shape<Tconfig>, mui::geometry::any_shape<Tconfig>>();
}


// [*** CHRONO_SAMPLER CLASSES ***] //


//CHRONO_SAMPLER_EXACT CLASS//
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(chrono_sampler_exact) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    using Ttime = typename Tconfig::time_type;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Ttime>(), py::arg("tol") = Ttime(0.0));
}

//CHRONO_SAMPLER_GAUSS CLASS//
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(chrono_sampler_gauss) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    using Treal = typename Tconfig::REAL;
    using Ttime = typename Tconfig::time_type;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Ttime, Treal>());
}

//CHRONO_SAMPLER_MEAN CLASS//
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(chrono_sampler_mean) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    using Ttime = typename Tconfig::time_type;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Ttime, Ttime>(), py::arg("newleft")=Ttime(0), 
                                   py::arg("newright")=Ttime(0));
}

//CHRONO_SAMPLER_SUM CLASS//
template <template <typename Type> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(chrono_sampler_sum) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Tclass = TclassTemplate<Tconfig>;
    using Ttime = typename Tconfig::time_type;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Ttime, Ttime>(), py::arg("newleft")=Ttime(0), 
                                   py::arg("newright")=Ttime(0));
}

// [*** CHRONO_SAMPLER CLASSES END ***] //



// [*** SPATIAL_SAMPLER CLASSES ***] //

//SPATIAL_SAMPLER_EXACT CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_exact) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    using Treal = typename Tconfig::REAL;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal>(), py::arg("tol") = Treal(std::numeric_limits<Treal>::epsilon()));
}

//SPATIAL_SAMPLER_GAUSS CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_gauss) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal,Treal>());
}

//SPATIAL_SAMPLER_MOVING_AVERAGE CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_moving_average) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<mui::point<Treal,Tconfig::D>>());
}

//SPATIAL_SAMPLER_NEAREST_NEIGHBOR CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_nearest_neighbor) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<>());
}

//SPATIAL_SAMPLER_PSEUDO_NEAREST2_LINEAR CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_pseudo_nearest2_linear) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal>());
}

//SPATIAL_SAMPLER_PSEUDO_NEAREST_NEIGHBOR CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_pseudo_nearest_neighbor) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal>());
}

//SPATIAL_SAMPLER_SHEPARD_QUINTIC CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_shepard_quintic) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal>());
}

//SPATIAL_SAMPLER_SPH_QUINTIC CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_sph_quintic) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal>());
}

//SPATIAL_SAMPLER_SUM_QUINTIC CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_sum_quintic) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal>());
}

#ifdef USE_RBF

//SPATIAL_SAMPLER_RBF CLASS//
template <template <typename Type0, typename Type1, typename Type2> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(sampler_rbf) {
    string pyclass_name = get_pyclass_name(name, typestr, arg1);
    using Treal = typename Tconfig::REAL;
    using Tint = typename Tconfig::INT;
    using Tpoint = typename Tconfig::point_type;
    using Tclass = TclassTemplate<Tconfig,TArg1,TArg1>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<Treal, std::vector<Tpoint> &, int, bool, bool, bool, bool, const std::string&, Treal>());
}

#endif

// [*** SPATIAL_SAMPLER CLASSES END ***] //


// *** POINT CLASS *** //
template <template <typename Type0, uint Type1> class TclassTemplate, typename Tconfig, typename TArg1=void>
DECLARE_FUNC_HEADER(point) {
    string pyclass_name = get_pyclass_name(name, typestr);
    using Treal = typename Tconfig::REAL;
    using Tclass = TclassTemplate<Treal,Tconfig::D>;
    py::class_<Tclass>(m, pyclass_name.c_str())
    .def(py::init<>())
    .def(py::init([](const std::vector<typename Tconfig::REAL>& arr) { 
       return Tclass(arr.data()); 
    }))
    .def("__repr__",
        [](const Tclass &obj) {
            std::stringstream rep_str;
            if (Tconfig::D == 1)
            	rep_str << "( " << obj[0] << " )";
            else if (Tconfig::D == 2)
            	rep_str << "( " << obj[0] << ", " << obj[1] << " )";
            else if (Tconfig::D == 3)
            	rep_str << "( " << obj[0] << ", " << obj[1] << ", " << obj[2] << " )";
            return rep_str.str();
        });
     
}

//TODO: CHange this to general functions
#define DECLARE_MUI_FUNCTION(FUNCNAME,TYPE,TYPE_STR)	\
       m.def("_" STRINGIFY(FUNCNAME) STRINGIFY(TYPE_STR), &create_uniface<mui::mui_config_##TYPE> , "", py::arg("domain"), \
             py::arg("interfaces"), py::arg("world")=NULL);

#ifdef PYTHON_INT_64
#define DECLARE_MUI_CPP2PY_FUNCTIONS(FUNCNAME)	\
    DECLARE_MUI_FUNCTION(FUNCNAME,1dx,1d_f64_i64) \
    DECLARE_MUI_FUNCTION(FUNCNAME,2dx,2d_f64_i64) \
    DECLARE_MUI_FUNCTION(FUNCNAME,3dx,3d_f64_i64) \
    DECLARE_MUI_FUNCTION(FUNCNAME,1fx,1d_f32_i64) \
    DECLARE_MUI_FUNCTION(FUNCNAME,2fx,2d_f32_i64) \
    DECLARE_MUI_FUNCTION(FUNCNAME,3fx,3d_f32_i64)
#elif defined PYTHON_INT_32
#define DECLARE_MUI_CPP2PY_FUNCTIONS(FUNCNAME)	\
    DECLARE_MUI_FUNCTION(FUNCNAME,1d,1d_f64_i32) \
    DECLARE_MUI_FUNCTION(FUNCNAME,2d,2d_f64_i32) \
    DECLARE_MUI_FUNCTION(FUNCNAME,3d,3d_f64_i32) \
    DECLARE_MUI_FUNCTION(FUNCNAME,1f,1d_f32_i32) \
    DECLARE_MUI_FUNCTION(FUNCNAME,2f,2d_f32_i32) \
    DECLARE_MUI_FUNCTION(FUNCNAME,3f,3d_f32_i32)
#else
#error PYTHON_INT_[32|64] not defined.
#endif


template<class CONFIG>
std::vector<std::unique_ptr<mui::uniface<CONFIG>>>
create_uniface(std::string domain, std::vector<std::string> interfaces, py::handle world=NULL) {
    return mui::create_uniface<CONFIG>(domain, interfaces, MPI_COMM_WORLD);
}

py::handle mpi_split_by_app() {
    if (import_mpi4py() < 0) {
        Py_RETURN_NONE;
    }
    return PyMPIComm_New(mui::mpi_split_by_app());
}


string get_mpi_version() {
#ifdef MPI_VERSION_STR
    return MPI_VERSION_STR;
#else
    return "";
#endif
}

string get_compiler_version() {
#ifdef COMPILER_VERSION_STR
    return COMPILER_VERSION_STR;
#else
    return "";
#endif
}

string get_compiler_config() {
#ifdef COMPILER_CONFIG_STR
    return COMPILER_CONFIG_STR;
#else
    return "";
#endif
}

PYBIND11_MODULE(mui4py_mod, m) {
    m.doc() = "MUI bindings for Python."; // optional module docstring

    // Expose numerical limits from C++
    m.attr("numeric_limits_real") = std::numeric_limits<double>::min();
    m.attr("numeric_limits_int") = std::numeric_limits<int>::min();

   // Class bindings
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(geometry_shape,geometry::shape)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(geometry_any_shape,geometry::any_shape)
   // DECLARE_MUI_CPP2PY_CLASSES_0ARG(geometry_or_set,geometry::or_set)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(geometry_point,geometry::point)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(geometry_sphere,geometry::sphere)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(geometry_box,geometry::box)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(uniface,uniface)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(point,point)

   // Spatial samplers
   #define EXPAND_SPATIAL_EXACT_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_exact,sampler_exact,T)
   FOR_EACH(EXPAND_SPATIAL_EXACT_SAMPLER,double,float,int32_t,int64_t,string)
   #define EXPAND_SPATIAL_GAUSS_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_gauss,sampler_gauss,T)
   FOR_EACH(EXPAND_SPATIAL_GAUSS_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_MOVING_AVERAGE_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_moving_average,sampler_moving_average,T)
   FOR_EACH(EXPAND_SPATIAL_MOVING_AVERAGE_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_NEAREST_NEIGHBOR_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_nearest_neighbor,sampler_nearest_neighbor,T)
   FOR_EACH(EXPAND_SPATIAL_NEAREST_NEIGHBOR_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_PSEUDO_NEAREST2_LINEAR_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_pseudo_nearest2_linear,sampler_pseudo_nearest2_linear,T)
   FOR_EACH(EXPAND_SPATIAL_PSEUDO_NEAREST2_LINEAR_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_PSEUDO_NEAREST_NEIGHBOR_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_pseudo_nearest_neighbor,sampler_pseudo_nearest_neighbor,T)
   FOR_EACH(EXPAND_SPATIAL_PSEUDO_NEAREST_NEIGHBOR_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_SHEPARD_QUINTIC_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_shepard_quintic,sampler_shepard_quintic,T)
   FOR_EACH(EXPAND_SPATIAL_SHEPARD_QUINTIC_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_SPH_QUINTIC_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_sph_quintic,sampler_sph_quintic,T)
   FOR_EACH(EXPAND_SPATIAL_SPH_QUINTIC_SAMPLER,double,float,int32_t,int64_t)
   #define EXPAND_SPATIAL_SUM_QUINTIC_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_sum_quintic,sampler_sum_quintic,T)
   FOR_EACH(EXPAND_SPATIAL_SUM_QUINTIC_SAMPLER,double,float,int32_t,int64_t)
#ifdef USE_RBF
    #define EXPAND_SPATIAL_RBF_SAMPLER(T) \
        DECLARE_MUI_CPP2PY_CLASSES_1ARG(sampler_rbf,sampler_rbf,T)
   FOR_EACH(EXPAND_SPATIAL_RBF_SAMPLER,double,float,int32_t,int64_t)
#endif

   // Chrono samplers
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(chrono_sampler_exact,chrono_sampler_exact)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(chrono_sampler_gauss,chrono_sampler_gauss)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(chrono_sampler_mean,chrono_sampler_mean)
   DECLARE_MUI_CPP2PY_CLASSES_0ARG(chrono_sampler_sum,chrono_sampler_sum)

   // Exposed MUI functions
   DECLARE_MUI_CPP2PY_FUNCTIONS(create_uniface)	
   m.def("set_quiet", &mui::set_quiet, "");
   m.def("mpi_split_by_app", &mpi_split_by_app, "");
   m.def("get_mpi_version", &get_mpi_version, "");
   m.def("get_compiler_config", &get_compiler_config, "");
   m.def("get_compiler_version", &get_compiler_version, "");
}
