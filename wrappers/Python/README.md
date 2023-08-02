# Information about the MUI Python wrapper mui4py
### Notes
- Support for `int32`, `int64`, `float`, `double` and `std::string`  is included for  `push()` and `fetch()`   *-like* functions. `push_many()` and `fetch_many()` do not work with `std::string` as they use `numpy` arrays, which  do not play well with strings at the `pybind11` level.

- Only `sampler_exact` and `temporal_sampler_exact` classes work with `std::string` type.
- Interfaces can be configured to be `1d`, `2d` or  `3d` and to use  `float` or `double` for point arithmetics and time stamps.
- A `mui4py-demos` folder has been included in the `mui-demos` repository.

### Extra functionality
- `push_many()` and `fetch_many()` functions. They can push/fetch a list of values at different locations using C looping. They can provide great speedups with respect to python loops. They make use of `numpy.ndarray` .
- `get_mpi_version()`, `get_compiler_version()`and `get_compiler_info()` provide compile-time information about `mui4py`.

# Building

## Installation
- The easiest way to install the Python package is to use:

```
pip3 install .
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In the event that CMake encounter an error related to `mui.h` not being found during compilation, you can resolve this by adding the following flag to the CMAKE_ARGS:

```
-DCMAKE_INCLUDE_DIRECTORIES=/path/to/MUI/src
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In the event that CMake is unable to locate Python3 or is unable to find the correct version of Python3, you can use the following flag to the CMAKE_ARGS to specify the path to the correct Python3 executable:

```
-DPython3_EXECUTABLE=/path/to/python3
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If CMake is unable to locate pybind11 during the compilation process, you can use the following flag to the CMAKE_ARGS to specify the path to the pybind11 directory:

```
-DCMAKE_PREFIX_PATH=/path/to/pybind11
```

- Alternatively `python3 setup.py install` may work.


- The C++ Python bindings can also be built directly with cmake and make, which is useful for testing. An example of built with cmake and make:

```
mkdir build && cd build
cmake ..
make
cd .. && cp build/*.so mui4py/
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note that if built with cmake and make, the path to the MUI Python wrapper should be added to the PYTHONPATH before use:

```
export PYTHONPATH=$PYTHONPATH:/path/to/MUI/wrappers/Python
```

### General

- Bindings tested with:
	- `gcc 7.2.0` and `clang LLVM version 9.1.0`, but it should work with any version of them supporting `c++11`.
	- `Mpich 3.2.1`, `OpenMPI 3.1.0` and `Spectrum MPI 10.02`.

- It is advised to use **conda** or **virtualenv** environments. This will provide isolation and control over the dependencies.

### Dependencies
- mpi4py (3.0.0)
- numpy (>1.13)
- pybind11 (only at compile time)
	> `mpi4py`  package should be compiled with the same compiler and MPI version as `mui4py`. Use `mpi4py.get_config()` and `mui4py.get_compiler_config()`, `mui4py.get_compiler_version()`, `mui4py.get_mpi_version()` to gather compile-time information.

# Contact

Please contact *Eduardo Ramos Fernandez* at eduardo.rf159@gmail.com and *Wendi Liu* at wendi.liu@stfc.ac.uk for issues about overall design or *Chris Richardson (Cambridge University)* for issues relating to packaging.

# ToDo List

- Include boolean OR operation of shapes in the bindings.
- Write exception hierarchy to replace generic exceptions.
- Add 'collide' function to library.
