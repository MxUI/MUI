cmake_minimum_required(VERSION 3.18)
project(mui4py)

# Set C++ standard before finding pybind11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force new Python handling policy for CMP0148
cmake_policy(SET CMP0148 NEW)
set(PYBIND11_FINDPYTHON ON)

# Always use the Python interpreter from the active environment
# (Spack, Conda, system, etc.)
find_program(PYTHON_EXECUTABLE python3)
if(NOT PYTHON_EXECUTABLE)
  message(FATAL_ERROR "Could not find 'python3' in PATH. Activate your environment first.")
endif()

# Query Python for its include and site-package paths
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "from sysconfig import get_paths; print(get_paths()['include'])"
  OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Using Python executable: ${PYTHON_EXECUTABLE}")
message(STATUS "Python include dir: ${PYTHON_INCLUDE_DIR}")
message(STATUS "Python site-packages: ${PYTHON_SITE_PACKAGES}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_INCLUDE_DIRECTORIES})
include_directories(${PROJECT_SOURCE_DIR}/../../src)

find_package(MPI REQUIRED)
if(MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
elseif(NOT MPI_FOUND)
    message(SEND_ERROR "MPI not found")
endif(MPI_FOUND)

# Detect Python/Numpy integer size (for PYTHON_INT_32 or PYTHON_INT_64)
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" "-c"
  "import struct; print(struct.calcsize('P') * 8)"
  OUTPUT_VARIABLE PYTHON_BITNESS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Detected Python bitness: ${PYTHON_BITNESS}-bit")
if(PYTHON_BITNESS EQUAL 64)
    set(PYTHON_INT_DEF PYTHON_INT_64)
elseif(PYTHON_BITNESS EQUAL 32)
    set(PYTHON_INT_DEF PYTHON_INT_32)
else()
    message(FATAL_ERROR "Unknown Python bitness: ${PYTHON_BITNESS}")
endif()
message(STATUS "Using Python integer definition: ${PYTHON_INT_DEF}")

# pybind11 detection
find_package(pybind11 2.7.0 REQUIRED CONFIG
  HINTS ${PYTHON_SITE_PACKAGES} ${PYBIND11_DIR} ${PYBIND11_ROOT} $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})

# Create the binding library pybind11 handles its own calls to
# target_link_libraries
pybind11_add_module(mui4py_mod MODULE
  mui4py/cpp/mui4py.cpp
  mui4py/cpp/geometry.cpp
  mui4py/cpp/sampler.cpp
  mui4py/cpp/temporal_sampler.cpp
  mui4py/cpp/algorithm.cpp
  mui4py/cpp/uniface1d.cpp
  mui4py/cpp/uniface2d.cpp
  mui4py/cpp/uniface3d.cpp
  mui4py/cpp/uniface1f.cpp
  mui4py/cpp/uniface2f.cpp
  mui4py/cpp/uniface3f.cpp)

target_include_directories(mui4py_mod PRIVATE ${MPI_INCLUDE_PATH} ../..)
target_compile_definitions(mui4py_mod PRIVATE PYTHON_BINDINGS ${PYTHON_INT_DEF})
target_compile_options(mui4py_mod PRIVATE -Wall)

# In Debug mode override pybind11 symbols visibility Symbols must be visible to
# backtrace_symbols() to produce nice logs
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(mui4py_mod PRIVATE "-fvisibility=default")
endif()

target_link_libraries(mui4py_mod PUBLIC ${MPI_C_LIBRARIES})

# Check for mpi4py
set(ENV{PYTHONPATH} "${PYTHON_SITE_PACKAGES}:$ENV{PYTHONPATH}")
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import mpi4py; print(mpi4py.get_include())"
  OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
  RESULT_VARIABLE MPI4PY_COMMAND_RESULT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT MPI4PY_COMMAND_RESULT EQUAL 0)
  message(FATAL_ERROR "mpi4py could not be found in ${PYTHON_SITE_PACKAGES}. Make sure it is installed.")
endif()
message(STATUS "Found mpi4py include dir: ${MPI4PY_INCLUDE_DIR}")
target_include_directories(mui4py_mod PRIVATE ${MPI4PY_INCLUDE_DIR})
