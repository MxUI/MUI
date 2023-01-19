cmake_minimum_required(VERSION 3.16)

project(mui4py)

# Set C++ standard before finding pybind11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Development component which includes both Development.Module and
# Development.Embed is not required for building a Python module.  Correct
# COMPONENT specification Development.Module added only in CMake 3.18 and above.
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0")
  find_package(
    Python3
    COMPONENTS Interpreter Development.Module
    REQUIRED
  )
else()
  find_package(
    Python3
    COMPONENTS Interpreter Development
    REQUIRED
  )
endif()

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

get_cmake_property(_variableNames VARIABLES)
list(SORT _variableNames)

foreach(_variableName ${_variableNames})
  message(STATUS "${_variableName}=${${_variableName}}")
endforeach()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I../.. -I${EIGEN3_INCLUDE_DIR} -I${MPI_INCLUDE_PATH} -DPYTHON_BINDINGS -DPYTHON_INT_64")

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
  mui4py/mui4py.cpp
  mui4py/sampler.cpp
  mui4py/chrono_sampler.cpp
  mui4py/geometry.cpp)

# Add strict compiler flags
# include(CheckCXXCompilerFlag)
# check_cxx_compiler_flag("-Wall -Werror -pedantic" HAVE_PEDANTIC)
# if(HAVE_PEDANTIC)
# target_compile_options(cpp PRIVATE -Wall;-Werror;-pedantic)
# endif()

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
