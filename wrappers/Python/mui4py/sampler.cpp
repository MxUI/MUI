
#include <mui.h>
#include <pybind11/pybind11.h>
#include <string>
#include "config_name.h"

template <typename Tconfig, typename T>
void declare_sampler_exact_t(py::module &m)
{
    std::string name = "_Sampler_exact" + config_name<Tconfig>() + "_" + type_name<T>();
    using Tclass = mui::sampler_exact<Tconfig, T, T>;
    using Treal = typename Tconfig::REAL;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<Treal>(),
             py::arg("tol") = Treal(std::numeric_limits<Treal>::epsilon()));
}

template <typename Tconfig>
void declare_sampler_exact(py::module &m)
{
    declare_sampler_exact_t<Tconfig, double>(m);
    declare_sampler_exact_t<Tconfig, float>(m);
    declare_sampler_exact_t<Tconfig, std::int32_t>(m);
    declare_sampler_exact_t<Tconfig, std::int64_t>(m);
    declare_sampler_exact_t<Tconfig, std::string>(m);
}

// SPATIAL_SAMPLER_GAUSS CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_gauss_t(py::module &m)
{
    std::string name = "_Sampler_gauss" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_gauss<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal, Treal>());
}

template <typename Tconfig>
void declare_sampler_gauss(py::module &m)
{
    declare_sampler_gauss_t<Tconfig, double>(m);
    declare_sampler_gauss_t<Tconfig, float>(m);
    declare_sampler_gauss_t<Tconfig, std::int32_t>(m);
    declare_sampler_gauss_t<Tconfig, std::int64_t>(m);
    declare_sampler_gauss_t<Tconfig, std::string>(m);
}

// SPATIAL_SAMPLER_MOVING_AVERAGE CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_moving_average_t(py::module &m)
{
    std::string name = "_Sampler_moving_average" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_moving_average<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<mui::point<Treal, Tconfig::D>>());
}

template <typename Tconfig>
void declare_sampler_moving_average(py::module &m)
{
    declare_sampler_moving_average_t<Tconfig, double>(m);
    declare_sampler_moving_average_t<Tconfig, float>(m);
    declare_sampler_moving_average_t<Tconfig, std::int32_t>(m);
    declare_sampler_moving_average_t<Tconfig, std::int64_t>(m);
    declare_sampler_moving_average_t<Tconfig, std::string>(m);
}

// SPATIAL_SAMPLER_NEAREST_NEIGHBOR CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_nearest_neighbor_t(py::module &m)
{
    std::string name = "_Sampler_nearest_neighbor" + config_name<Tconfig>() + "_" + type_name<T>();
    using Tclass = mui::sampler_nearest_neighbor<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<>());
}

template <typename Tconfig>
void declare_sampler_nearest_neighbor(py::module &m)
{
    declare_sampler_nearest_neighbor_t<Tconfig, double>(m);
    declare_sampler_nearest_neighbor_t<Tconfig, float>(m);
    declare_sampler_nearest_neighbor_t<Tconfig, std::int32_t>(m);
    declare_sampler_nearest_neighbor_t<Tconfig, std::int64_t>(m);
    declare_sampler_nearest_neighbor_t<Tconfig, std::string>(m);
}

// SPATIAL_sampler_pseudo_n2_linear CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_pseudo_n2_linear_t(py::module &m)
{
    std::string name = "_Sampler_pseudo_n2_linear" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_pseudo_n2_linear<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_PSEUDO_NEAREST_NEIGHBOR CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_pseudo_nearest_neighbor_t(py::module &m)
{
    std::string name = "_Sampler_pseudo_nearest_neighbor" + config_name<Tconfig>() + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_pseudo_nearest_neighbor<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_SHEPARD_QUINTIC CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_shepard_quintic_t(py::module &m)
{
    std::string name = "_Sampler_shepard_quintic" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_shepard_quintic<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_SPH_QUINTIC CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_sph_quintic(py::module &m)
{
    std::string name = "_Sampler_sph_quintic" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_sph_quintic<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

// SPATIAL_SAMPLER_SUM_QUINTIC CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_sum_quintic_t(py::module &m)
{
    std::string name = "_Sampler_sum_quintic" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_sum_quintic<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

void sampler(py::module &m)
{
#ifdef PYTHON_INT_64

    declare_sampler_exact<mui::mui_config_1dx>(m);
    declare_sampler_exact<mui::mui_config_2dx>(m);
    declare_sampler_exact<mui::mui_config_3dx>(m);
    declare_sampler_exact<mui::mui_config_1fx>(m);
    declare_sampler_exact<mui::mui_config_2fx>(m);
    declare_sampler_exact<mui::mui_config_3fx>(m);

    declare_sampler_gauss<mui::mui_config_1dx>(m);
    declare_sampler_gauss<mui::mui_config_2dx>(m);
    declare_sampler_gauss<mui::mui_config_3dx>(m);
    declare_sampler_gauss<mui::mui_config_1fx>(m);
    declare_sampler_gauss<mui::mui_config_2fx>(m);
    declare_sampler_gauss<mui::mui_config_3fx>(m);

    declare_sampler_moving_average<mui::mui_config_1dx>(m);
    declare_sampler_moving_average<mui::mui_config_2dx>(m);
    declare_sampler_moving_average<mui::mui_config_3dx>(m);
    declare_sampler_moving_average<mui::mui_config_1fx>(m);
    declare_sampler_moving_average<mui::mui_config_2fx>(m);
    declare_sampler_moving_average<mui::mui_config_3fx>(m);

    declare_sampler_nearest_neighbor<mui::mui_config_1dx>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_2dx>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_3dx>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_1fx>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_2fx>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_3fx>(m);
#elif defined PYTHON_INT_32

    declare_sampler_exact<mui::mui_config_1d>(m);
    declare_sampler_gauss<mui::mui_config_1d>(m);
    declare_sampler_moving_average<mui::mui_config_1d>(m);
    declare_sampler_nearest_neighbor<mui::mui_config_1d>(m);

#else
#error PYTHON_INT_[32|64] not defined.
#endif
}