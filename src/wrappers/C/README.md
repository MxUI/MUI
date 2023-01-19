# Information about the MUI C wrapper
This creates a C wrapper for the majority of the functionality provided by the C++ MUI library.

# Compilation notes
It relies on the C++11 standard and so a compliant compiler will be required.

# Wrapper usage
C wrapper usage:
1. include mui_c_wapper_[1d/2d/3d].h in C source
2. Compile with a C compiler & link against library

# Unit tests
Two unit tests are provided with the wrapper: 
1. unit_test_single.c - Demonstrates creating and using a single instance of a 1D/2D/3D uniface
2. unit_test_multi.c - Demonstrates creating and using multiple 1D/2D/3D uniface instances using the mui::create_uniface<config>() helper function  
