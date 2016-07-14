cmake_minimum_required(VERSION 2.6)

project (MUI CXX)

enable_testing()

# option(BUILD_GTEST ON)
add_subdirectory(tests/googletest/googletest)

find_package (Threads)
add_subdirectory(tests)
