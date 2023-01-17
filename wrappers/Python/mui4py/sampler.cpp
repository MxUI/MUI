
#include <mui.h>
#include <pybind11/pybind11.h>
#include <string>
#include "config_name.h"

template <typename Tconfig, typename TArg1>
void declare_sampler_exact(py::module &m)
{
    std::string name = "_sampler_exact_" + config_name<Tconfig>();
    using Tclass = mui::sampler_exact<Tconfig, TArg1, TArg1>;
    using Treal = typename Tconfig::REAL;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<Treal>(),
             py::arg("tol") = Treal(std::numeric_limits<Treal>::epsilon()));
}

// SPATIAL_SAMPLER_GAUSS CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_gauss(py::module &m)
{
    std::string name = "_sampler_gauss_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_gauss<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal, Treal>());
}

// SPATIAL_SAMPLER_MOVING_AVERAGE CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_moving_average(py::module &m)
{
    std::string name = "_sampler_moving_average_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_moving_average<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<mui::point<Treal, Tconfig::D>>());
}

// SPATIAL_SAMPLER_NEAREST_NEIGHBOR CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_nearest_neighbor(py::module &m)
{
    std::string name = "_sampler_nearest_neighbor" + config_name<Tconfig>();
    using Tclass = mui::sampler_nearest_neighbor<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<>());
}

// SPATIAL_sampler_pseudo_n2_linear CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_pseudo_n2_linear(py::module &m)
{
    std::string name = "_sampler_pseudo_n2_linear_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_pseudo_n2_linear<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_PSEUDO_NEAREST_NEIGHBOR CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_pseudo_nearest_neighbor(py::module &m)
{
    std::string name = "_sampler_pseudo_nearest_neighbor" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_pseudo_nearest_neighbor<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_SHEPARD_QUINTIC CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_shepard_quintic(py::module &m)
{
    std::string name = "_sampler_shepard_quintic_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_shepard_quintic<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_SPH_QUINTIC CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_sph_quintic(py::module &m)
{
    std::string name = "_sampler_sph_quintic_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_sph_quintic<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_SUM_QUINTIC CLASS//
template <typename Tconfig, typename TArg1 = void>
void declare_sampler_sum_quintic(py::module &m)
{
    std::string name = "_sampler_sum_quintic_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_sum_quintic<Tconfig, TArg1, TArg1>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

void sampler(py::module &m)
{


    declare_sampler_exact<mui::mui_config_1dx, double>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_1dx, double>(m);
}