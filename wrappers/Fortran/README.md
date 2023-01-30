# Information about the MUI Fortran wrapper
This creates a Fortran wrapper for the majority of the functionality provided by the C++ MUI library.

# Compilation notes
It relies on C++11 and the f2003 ISO_C_BIND standard and so compliant compilers will be required.

# Wrapper usage
Fortran wrapper usage:
1. Compile application along with 1D/2D/3D module source
2. Link compiled mui_f_wrapper library
3. Enable f2003 standard during compilation

Note: functions that return *arrays* of values perform memory allocation within the C function, therefore it is necessery to
      treat the returned memory locations accordingly within the calling Fortran application. 

# Unit tests
A unit test is provided with the wrapper:
1. unit_test.f90 - Demonstrates creating and using a single instance of a 1D/2D/3D uniface
2. unit_test_multi.f90 - Demonstrates creating and using multiple 1D/2D/3D uniface instances using the create_and_get_uniface_multi_xx_f() helper function 