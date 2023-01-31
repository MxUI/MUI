
#include <mui.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include "config_name.h"

template <typename Tconfig>
void declare_temporal_samplers(py::module &m)
{
    std::string name = "_Temporal_sampler_exact" + config_name<Tconfig>();
    using Ttime = typename Tconfig::time_type;
    py::class_<mui::temporal_sampler_exact<Tconfig>>(m, name.c_str())
        .def(py::init<Ttime>(), py::arg("tol") = Ttime(0.0));

    name = "_Temporal_sampler_gauss" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Ttime = typename Tconfig::time_type;
    py::class_<mui::temporal_sampler_gauss<Tconfig>>(m, name.c_str()).def(py::init<Ttime, Treal>());

    name = "_Temporal_sampler_mean" + config_name<Tconfig>();
    using Ttime = typename Tconfig::time_type;
    py::class_<mui::temporal_sampler_mean<Tconfig>>(m, name.c_str())
        .def(py::init<Ttime, Ttime>(), py::arg("newleft") = Ttime(0),
             py::arg("newright") = Ttime(0));

    name = "_Temporal_sampler_sum" + config_name<Tconfig>();
    using Tclass = mui::temporal_sampler_sum<Tconfig>;
    using Ttime = typename Tconfig::time_type;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<Ttime, Ttime>(), py::arg("newleft") = Ttime(0),
             py::arg("newright") = Ttime(0));
}

void temporal_sampler(py::module &m)
{
#ifdef PYTHON_INT_64
    declare_temporal_samplers<mui::mui_config_1dx>(m);
    declare_temporal_samplers<mui::mui_config_2dx>(m);
    declare_temporal_samplers<mui::mui_config_3dx>(m);
    declare_temporal_samplers<mui::mui_config_1fx>(m);
    declare_temporal_samplers<mui::mui_config_2fx>(m);
    declare_temporal_samplers<mui::mui_config_3fx>(m);
#elif defined PYTHON_INT_32
    declare_temporal_samplers<mui::mui_config_1d>(m);
    declare_temporal_samplers<mui::mui_config_2d>(m);
    declare_temporal_samplers<mui::mui_config_3d>(m);
    declare_temporal_samplers<mui::mui_config_1f>(m);
    declare_temporal_samplers<mui::mui_config_2f>(m);
    declare_temporal_samplers<mui::mui_config_3f>(m);
#else
#error PYTHON_INT_[32|64] not defined.
#endif
}
