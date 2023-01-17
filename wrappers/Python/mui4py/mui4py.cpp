
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
}

// Declaration of other files
void chrono_sampler(py::module &m);
void sampler(py::module &m);
void geometry(py::module &m);

PYBIND11_MODULE(mui4py_mod, m)
{
  m.doc() = "MUI bindings for Python."; // optional module docstring

  // Expose numerical limits from C++
  m.attr("numeric_limits_real") = std::numeric_limits<double>::min();
  m.attr("numeric_limits_int") = std::numeric_limits<int>::min();

  geometry(m);
  sampler(m);
  chrono_sampler(m);

#ifdef PYTHON_INT_64
  declare_uniface_class<mui::mui_config_1dx>(m);
  declare_uniface_class<mui::mui_config_2dx>(m);
  declare_uniface_class<mui::mui_config_3dx>(m);
  declare_uniface_class<mui::mui_config_1fx>(m);
  declare_uniface_class<mui::mui_config_2fx>(m);
  declare_uniface_class<mui::mui_config_3fx>(m);

  py_create_uniface<mui::mui_config_1dx>(m);
  py_create_uniface<mui::mui_config_2dx>(m);
  py_create_uniface<mui::mui_config_3dx>(m);
  py_create_uniface<mui::mui_config_1fx>(m);
  py_create_uniface<mui::mui_config_2fx>(m);
  py_create_uniface<mui::mui_config_3fx>(m);
#elif defined PYTHON_INT_32
  declare_uniface_class<mui::mui_config_1d>(m);
  declare_uniface_class<mui::mui_config_2d>(m);
  declare_uniface_class<mui::mui_config_3d>(m);
  declare_uniface_class<mui::mui_config_1f>(m);
  declare_uniface_class<mui::mui_config_2f>(m);
  declare_uniface_class<mui::mui_config_3f>(m);

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
