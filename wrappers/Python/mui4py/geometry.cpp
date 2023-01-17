
#include <mui.h>
#include <pybind11/pybind11.h>
#include <string>
#include "config_name.h"

// Geometry

template <typename Tconfig>
void declare_geometry_shape(py::module &m)
{
    std::string name = "_geometry_shape_" + config_name<Tconfig>();
    py::class_<mui::geometry::shape<Tconfig>>(m, name.c_str());
}

template <typename Tconfig>
void declare_geometry_any_shape(py::module &m)
{
    std::string name = "_geometry_any_shape_" + config_name<Tconfig>();
    using Tclass = mui::geometry::any_shape<Tconfig>;
    py::class_<Tclass>(m, name.c_str())
        .def(py::init<const mui::geometry::shape<Tconfig> &>());
}

template <typename Tconfig>
void declare_geometry_point(py::module &m)
{
    std::string name = "_geometry_point_" + config_name<Tconfig>();
    py::class_<mui::geometry::point<Tconfig>, mui::geometry::shape<Tconfig>>(m, name.c_str())
        .def(py::init<const mui::point<typename Tconfig::REAL, Tconfig::D> &>());
}

template <typename Tconfig>
void declare_geometry_sphere(py::module &m)
{
    std::string name = "_geometry_sphere_" + config_name<Tconfig>();
    py::class_<mui::geometry::sphere<Tconfig>, mui::geometry::shape<Tconfig>>(m, name.c_str())
        .def(py::init<const mui::point<typename Tconfig::REAL, Tconfig::D> &,
                      typename Tconfig::REAL>());
}

template <typename Tconfig>
void declare_geometry_box(py::module &m)
{
    std::string name = "_geometry_box_" + config_name<Tconfig>();
    using Treal = typename Tconfig::REAL;
    py::class_<mui::geometry::box<Tconfig>, mui::geometry::shape<Tconfig>>(m, name.c_str())
        .def(py::init<const mui::point<Treal, Tconfig::D> &,
                      const mui::point<Treal, Tconfig::D> &>());
}

void geometry(py::module &m)
{
#ifdef PYTHON_INT_64
    declare_geometry_point<mui::mui_config_1dx>(m);
    declare_geometry_point<mui::mui_config_2dx>(m);
    declare_geometry_point<mui::mui_config_3dx>(m);
    declare_geometry_point<mui::mui_config_1fx>(m);
    declare_geometry_point<mui::mui_config_2fx>(m);
    declare_geometry_point<mui::mui_config_3fx>(m);

    declare_geometry_box<mui::mui_config_1dx>(m);
    declare_geometry_box<mui::mui_config_2dx>(m);
    declare_geometry_box<mui::mui_config_3dx>(m);
    declare_geometry_box<mui::mui_config_1fx>(m);
    declare_geometry_box<mui::mui_config_2fx>(m);
    declare_geometry_box<mui::mui_config_3fx>(m);

    declare_geometry_sphere<mui::mui_config_1dx>(m);
    declare_geometry_sphere<mui::mui_config_2dx>(m);
    declare_geometry_sphere<mui::mui_config_3dx>(m);
    declare_geometry_sphere<mui::mui_config_1fx>(m);
    declare_geometry_sphere<mui::mui_config_2fx>(m);
    declare_geometry_sphere<mui::mui_config_3fx>(m);

    declare_geometry_any_shape<mui::mui_config_1dx>(m);
    declare_geometry_any_shape<mui::mui_config_2dx>(m);
    declare_geometry_any_shape<mui::mui_config_3dx>(m);
    declare_geometry_any_shape<mui::mui_config_1fx>(m);
    declare_geometry_any_shape<mui::mui_config_2fx>(m);
    declare_geometry_any_shape<mui::mui_config_3fx>(m);

    declare_geometry_shape<mui::mui_config_1dx>(m);
    declare_geometry_shape<mui::mui_config_2dx>(m);
    declare_geometry_shape<mui::mui_config_3dx>(m);
    declare_geometry_shape<mui::mui_config_1fx>(m);
    declare_geometry_shape<mui::mui_config_2fx>(m);
    declare_geometry_shape<mui::mui_config_3fx>(m);

#elif defined PYTHON_INT_32

    declare_geometry_point<mui::mui_config_1d>(m);
    declare_geometry_point<mui::mui_config_2d>(m);
    declare_geometry_point<mui::mui_config_3d>(m);
    declare_geometry_point<mui::mui_config_1f>(m);
    declare_geometry_point<mui::mui_config_2f>(m);
    declare_geometry_point<mui::mui_config_3f>(m);

    declare_geometry_box<mui::mui_config_1d>(m);
    declare_geometry_box<mui::mui_config_2d>(m);
    declare_geometry_box<mui::mui_config_3d>(m);
    declare_geometry_box<mui::mui_config_1f>(m);
    declare_geometry_box<mui::mui_config_2f>(m);
    declare_geometry_box<mui::mui_config_3f>(m);

    declare_geometry_sphere<mui::mui_config_1d>(m);
    declare_geometry_sphere<mui::mui_config_2d>(m);
    declare_geometry_sphere<mui::mui_config_3d>(m);
    declare_geometry_sphere<mui::mui_config_1f>(m);
    declare_geometry_sphere<mui::mui_config_2f>(m);
    declare_geometry_sphere<mui::mui_config_3f>(m);

    declare_geometry_any_shape<mui::mui_config_1d>(m);
    declare_geometry_any_shape<mui::mui_config_2d>(m);
    declare_geometry_any_shape<mui::mui_config_3d>(m);
    declare_geometry_any_shape<mui::mui_config_1f>(m);
    declare_geometry_any_shape<mui::mui_config_2f>(m);
    declare_geometry_any_shape<mui::mui_config_3f>(m);

    declare_geometry_shape<mui::mui_config_1d>(m);
    declare_geometry_shape<mui::mui_config_2d>(m);
    declare_geometry_shape<mui::mui_config_3d>(m);
    declare_geometry_shape<mui::mui_config_1f>(m);
    declare_geometry_shape<mui::mui_config_2f>(m);
    declare_geometry_shape<mui::mui_config_3f>(m);

#else
#error PYTHON_INT_[32|64] not defined.
#endif
}