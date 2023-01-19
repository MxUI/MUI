cmake_minimum_required(VERSION 3.18)

project(mui4py)

# Set C++ standard before finding pybind11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Development component which includes both Development.Module and
# Development.Embed is not required for building a Python module.  Correct
# COMPONENT specification Development.Module added only in CMake 3.18 and above.
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

find_package(MPI REQUIRED)

# Check for numpy int size
execute_process(
  COMMAND "${Python3_EXECUTABLE}" "-c"
  "import numpy; print(str(numpy.iinfo(int).dtype)[-2:])"
  OUTPUT_VARIABLE NUMPY_INT_SIZE
  RESULT_VARIABLE NUMPY_COMMAND_RESULT
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "numpy int size = ${NUMPY_INT_SIZE}")

find_package(
  pybind11
  2.7.0
  REQUIRED
  CONFIG
  HINTS
  ${PYBIND11_DIR}
  ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR}
  $ENV{PYBIND11_ROOT}
)

# Create the binding library pybind11 handles its own calls to
# target_link_libraries
pybind11_add_module(mui4py_mod MODULE
  mui4py/cpp/mui4py.cpp
  mui4py/cpp/sampler.cpp
  mui4py/cpp/chrono_sampler.cpp
  mui4py/cpp/geometry.cpp)

if(USE_RBF)
  # Note that -fpermissive is required to compile using the rbf sampler
  # and that it is not supported by clang++
  find_package(Eigen3 REQUIRED NO_MODULE)
  target_include_directories(mui4py_mod PRIVATE ${EIGEN3_INCLUDE_DIR})
  target_compile_definitions(mui4py_mod PRIVATE USE_RBF)
  target_compile_options(mui4py_mod PRIVATE -fpermissive)
endif()

target_include_directories(mui4py_mod PRIVATE ${MPI_INCLUDE_PATH} ../..)
target_compile_definitions(mui4py_mod PRIVATE PYTHON_BINDINGS PYTHON_INT_${NUMPY_INT_SIZE})
target_compile_options(mui4py_mod PRIVATE -Wall)

# In Debug mode override pybind11 symbols visibility Symbols must be visible to
# backtrace_symbols() to produce nice logs
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(mui4py_mod PRIVATE "-fvisibility=default")
endif()

target_link_libraries(mui4py_mod PUBLIC ${MPI_C_LIBRARIES})

# Check for mpi4py
execute_process(
  COMMAND "${Python3_EXECUTABLE}" "-c"
  "import mpi4py; print(mpi4py.get_include())"
  OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
  RESULT_VARIABLE MPI4PY_COMMAND_RESULT
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT MPI4PY_COMMAND_RESULT)
  message(STATUS "Found mpi4py include directory at ${MPI4PY_INCLUDE_DIR}")
  target_include_directories(mui4py_mod PRIVATE ${MPI4PY_INCLUDE_DIR})
else()
  message(FATAL_ERROR "mpi4py could not be found.")
endif()
