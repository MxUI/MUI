
#include <mui.h>
#include <pybind11/chrono.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi4py/mpi4py.h>
#include <string>

#include "config_name.h"

std::string get_mpi_version()
{
#ifdef MPI_VERSION_STR
  return MPI_VERSION_STR;
#else
  return "";
#endif
}

std::string get_compiler_version()
{
#ifdef COMPILER_VERSION_STR
  return COMPILER_VERSION_STR;
#else
  return "";
#endif
}

std::string get_compiler_config()
{
#ifdef COMPILER_CONFIG_STR
  return COMPILER_CONFIG_STR;
#else
  return "";
#endif
}

template <class Tconfig>
void py_create_uniface(py::module &m)
{
  std::string name = "_create_uniface_" + config_name<Tconfig>();
  m.def(name.c_str(), [](std::string domain, std::vector<std::string> interfaces, py::handle world)
        { return mui::create_uniface<Tconfig>(domain, interfaces, MPI_COMM_WORLD); });
}

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

template <typename Tconfig>
void declare_uniface_class(py::module &m)
{
  std::string name = "_uniface_" + config_name<Tconfig>();
  using Tclass = mui::uniface<Tconfig>;
  using Treal = typename Tconfig::REAL;
  using Ttime = typename Tconfig::time_type;
  py::class_<Tclass>(m, name.c_str())
      .def("commit", (int(Tclass::*)(Ttime)) & Tclass::commit, "")
      .def("forecast", (void(Tclass::*)(Ttime)) & Tclass::forecast, "")
      //.def("is_ready", &Tclass::is_ready, "")
      .def("barrier", (void(Tclass::*)(Ttime)) & Tclass::barrier, "")
      .def("barrier", (void(Tclass::*)(Ttime, Ttime)) & Tclass::barrier, "")
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
      // DEFINE_MUI_UNIFACE_PUSH() DEFINE_MUI_UNIFACE_FETCH_1ARG()
      //    DEFINE_MUI_UNIFACE_FETCH_5ARGS() DEFINE_MUI_UNIFACE_FETCH_6ARGS()
      // Temporarily disabled the fetch_point function. Bind the variadic
      // template later.
      //     DEFINE_MUI_UNIFACE_FETCH_POINTS()
      .def(py::init<const std::string &>());
  py::implicitly_convertible<mui::geometry::shape<Tconfig>,
                             mui::geometry::any_shape<Tconfig>>();
}

void chrono_sampler(py::module &m);
void geometry(py::module &m);

PYBIND11_MODULE(mui4py_mod, m)
{
  m.doc() = "MUI bindings for Python."; // optional module docstring

  // Expose numerical limits from C++
  m.attr("numeric_limits_real") = std::numeric_limits<double>::min();
  m.attr("numeric_limits_int") = std::numeric_limits<int>::min();

  geometry(m);

  declare_uniface_class<mui::mui_config_1dx>(m);

  declare_sampler_exact<mui::mui_config_1dx, double>(m);
  declare_sampler_nearest_neighbor<mui::mui_config_1dx, double>(m);

#ifdef PYTHON_INT_64
  py_create_uniface<mui::mui_config_1dx>(m);
  py_create_uniface<mui::mui_config_2dx>(m);
  py_create_uniface<mui::mui_config_3dx>(m);
  py_create_uniface<mui::mui_config_1fx>(m);
  py_create_uniface<mui::mui_config_2fx>(m);
  py_create_uniface<mui::mui_config_3fx>(m);
#elif defined PYTHON_INT_32
  py_create_uniface<mui::mui_config_1d>(m);
  py_create_uniface<mui::mui_config_2d>(m);
  py_create_uniface<mui::mui_config_3d>(m);
  py_create_uniface<mui::mui_config_1f>(m);
  py_create_uniface<mui::mui_config_2f>(m);
  py_create_uniface<mui::mui_config_3f>(m);
#else
#error PYTHON_INT_[32|64] not defined.
#endif

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
