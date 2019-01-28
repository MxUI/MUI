# Information about the bindings
### Notes
- Support for `int32`, `int64`, `float`, `double` and `std::string`  is included for  `push()` and `fetch()`   *-like* functions. `push_many()` and `fetch_many()` do not work with `std::string` as they use `numpy` arrays, which  do not play well with strings at the `pybind11` level.

- Only `sampler_exact` and `chrono_sampler_exact` classes work with `std::string` type.
- Interfaces can be configured to be `1d`, `2d` or  `3d` and to use  `float` or `double` for point arithmetics and time stamps.
- A `mui4py-demos` folder has been included in the `mui-demos` repository under the `mui4py` branch.
- There has been a few additions/modifications to the previous code base of the library. Some of them were to include useful functionality to the bindings and others to correct *buggy* code (mostly due to inconsistency with typing).

### Extra functionality
- `push_many()` and `fetch_many()` functions. They can push/fetch a list of values at different locations using C looping. They can provide great speedups with respect to python loops. They make use of `numpy.ndarray` .
- `get_mpi_version()`, `get_compiler_version()`and `get_compiler_info()` provide compile-time information about `mui4py`.

### Missing functionality
- Boolean `OR` operations of `shape` classes.
- `or_set` class.
- `collide()`

# Building
### General

- Bindings tested with:
	-  `gcc 7.2.0` and `clang LLVM version 9.1.0`, but it should work with any version of them supporting `c++11`.
	- `Mpich 3.2.1`, `OpenMPI 3.1.0` and `Spectrum MPI 10.02`.

- It is advised to use **conda** or **virtualenv** environments. This will provide isolation and control over the dependencies.

### Dependencies
- mpi4py (3.0.0)
- numpy (>1.13)
- pybind11 (only at compile time)
	> `mpi4py`  package should be compiled with the same compiler and MPI version as `mui4py`. Not keeping this consistency could lead into non-trivial runtime errors. Use `mpi4py.get_config()` and `mui4py.get_compiler_config()`, `mui4py.get_compiler_version()`, `mui4py.get_mpi_version()` to gather compile-time information.

*Note: Indicated versions in brackets are orientative. It might work for earlier versions as well.* 

### Compiling  just the binary module (option 1)

``` bash
make COMPILER=[GCC|CLANG] mui4py_mod
export PYTHONPATH=/pathToMUI/wrappers/Python/mui4py
```

### Create python package (option 2)
``` bash
make COMPILER=[GCC|CLANG] package
make pip-install # To install with pip the created package
```
# Contact

Please contact *Eduardo Ramos Fernandez* at eduardo.rf159@gmail.com if you find any issue or have any enquiry.
