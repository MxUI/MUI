
#include <mui.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include "config_name.h"

template <typename Tconfig>
void declare_point(py::module &m)
{
    std::string name = "_Point" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::point<Treal, Tconfig::D>;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<>())
        .def(py::init([](const std::vector<typename Tconfig::REAL> &arr)
                      { return Tclass(arr.data()); }))
        .def("__repr__", [](const Tclass &obj)
             {
        std::stringstream rep_str;
        if (Tconfig::D == 1)
          rep_str << "( " << obj[0] << " )";
        else if (Tconfig::D == 2)
          rep_str << "( " << obj[0] << ", " << obj[1] << " )";
        else if (Tconfig::D == 3)
          rep_str << "( " << obj[0] << ", " << obj[1] << ", " << obj[2] << " )";
        return rep_str.str(); });
}

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

template <typename Tconfig>
void declare_sampler_pseudo_n2_linear(py::module &m)
{
    declare_sampler_pseudo_n2_linear_t<Tconfig, double>(m);
    declare_sampler_pseudo_n2_linear_t<Tconfig, float>(m);
    declare_sampler_pseudo_n2_linear_t<Tconfig, std::int32_t>(m);
    declare_sampler_pseudo_n2_linear_t<Tconfig, std::int64_t>(m);
}

// SPATIAL_SAMPLER_PSEUDO_NEAREST_NEIGHBOR CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_pseudo_nearest_neighbor_t(py::module &m)
{
    std::string name = "_Sampler_pseudo_nearest_neighbor" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_pseudo_nearest_neighbor<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

template <typename Tconfig>
void declare_sampler_pseudo_nearest_neighbor(py::module &m)
{
    declare_sampler_pseudo_nearest_neighbor_t<Tconfig, double>(m);
    declare_sampler_pseudo_nearest_neighbor_t<Tconfig, float>(m);
    declare_sampler_pseudo_nearest_neighbor_t<Tconfig, std::int32_t>(m);
    declare_sampler_pseudo_nearest_neighbor_t<Tconfig, std::int64_t>(m);
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

template <typename Tconfig>
void declare_sampler_shepard_quintic(py::module &m)
{
    declare_sampler_shepard_quintic_t<Tconfig, double>(m);
    declare_sampler_shepard_quintic_t<Tconfig, float>(m);
    declare_sampler_shepard_quintic_t<Tconfig, std::int32_t>(m);
    declare_sampler_shepard_quintic_t<Tconfig, std::int64_t>(m);
}

// SPATIAL_SAMPLER_SPH_QUINTIC CLASS//
template <typename Tconfig, typename T = void>
void declare_sampler_sph_quintic_t(py::module &m)
{
    std::string name = "_Sampler_sph_quintic" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tclass = mui::sampler_sph_quintic<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str()).def(py::init<Treal>());
}

template <typename Tconfig>
void declare_sampler_sph_quintic(py::module &m)
{
    declare_sampler_sph_quintic_t<Tconfig, double>(m);
    declare_sampler_sph_quintic_t<Tconfig, float>(m);
    declare_sampler_sph_quintic_t<Tconfig, std::int32_t>(m);
    declare_sampler_sph_quintic_t<Tconfig, std::int64_t>(m);
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

template <typename Tconfig>
void declare_sampler_sum_quintic(py::module &m)
{
    declare_sampler_sum_quintic_t<Tconfig, double>(m);
    declare_sampler_sum_quintic_t<Tconfig, float>(m);
    declare_sampler_sum_quintic_t<Tconfig, std::int32_t>(m);
    declare_sampler_sum_quintic_t<Tconfig, std::int64_t>(m);
}

#ifdef USE_RBF
// SPATIAL_SAMPLER_RBF CLASS//
template <typename Tconfig, typename T>
void declare_sampler_rbf_t(py::module &m)
{
    std::string name = "_Sampler_rbf" + config_name<Tconfig>() + "_" + type_name<T>();
    using Treal = typename Tconfig::REAL;
    using Tint = typename Tconfig::INT;
    using Tpoint = typename Tconfig::point_type;
    using Tclass = mui::sampler_rbf<Tconfig, T, T>;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<Treal, std::vector<Tpoint> &, int, bool, bool, bool,
                      bool, const std::string &, Treal, Treal, int, int>());
}

template <typename Tconfig>
void declare_sampler_rbf(py::module &m)
{
    declare_sampler_rbf_t<Tconfig, double>(m);
    declare_sampler_rbf_t<Tconfig, float>(m);
    declare_sampler_rbf_t<Tconfig, std::int32_t>(m);
    declare_sampler_rbf_t<Tconfig, std::int64_t>(m);
}
#endif

template <typename Tconfig>
void declare_samplers(py::module &m)
{
    declare_point<Tconfig>(m);

    declare_sampler_exact<Tconfig>(m);
    declare_sampler_gauss<Tconfig>(m);
    declare_sampler_moving_average<Tconfig>(m);
    declare_sampler_nearest_neighbor<Tconfig>(m);
    declare_sampler_pseudo_n2_linear<Tconfig>(m);
    declare_sampler_pseudo_nearest_neighbor<Tconfig>(m);
    declare_sampler_shepard_quintic<Tconfig>(m);
    declare_sampler_sph_quintic<Tconfig>(m);
    declare_sampler_sum_quintic<Tconfig>(m);
#ifdef USE_RBF
    declare_sampler_rbf<Tconfig>(m);
#endif
}

void sampler(py::module &m)
{
#ifdef PYTHON_INT_64

    declare_samplers<mui::mui_config_1dx>(m);
    declare_samplers<mui::mui_config_2dx>(m);
    declare_samplers<mui::mui_config_3dx>(m);
    declare_samplers<mui::mui_config_1fx>(m);
    declare_samplers<mui::mui_config_2fx>(m);
    declare_samplers<mui::mui_config_3fx>(m);

#elif defined PYTHON_INT_32

    declare_samplers<mui::mui_config_1d>(m);
    declare_samplers<mui::mui_config_2d>(m);
    declare_samplers<mui::mui_config_3d>(m);
    declare_samplers<mui::mui_config_1f>(m);
    declare_samplers<mui::mui_config_2f>(m);
    declare_samplers<mui::mui_config_3f>(m);

#else
#error PYTHON_INT_[32|64] not defined.
#endif
}