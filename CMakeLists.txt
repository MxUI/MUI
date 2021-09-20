cmake_minimum_required(VERSION 3.0)
project(MUI VERSION 1.1.2 DESCRIPTION "Multiscale Universal Interface" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(USE_RBF "Include RBF support" OFF)
option(C_WRAPPER "Compile C wrapper" OFF)
option(FORTRAN_WRAPPER "Compile Fortran wrapper" OFF) #ToDo: Add Fortran wrapper compilation
option(PYTHON_WRAPPER "Compile Python wrapper" OFF) #ToDo: Add Python wrapper compilation

include(sourcelist.cmake)

find_package(MPI REQUIRED)
if (MPI_FOUND)
	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
	message(SEND_ERROR "MPI not found")
endif (MPI_FOUND)

if(USE_RBF)
	find_package(Eigen3 REQUIRED)
	if( NOT EIGEN3_INCLUDE_DIR )
    		message(FATAL_ERROR "Please set the environment variable EIGEN3_INCLUDE_DIR to the include directory of valid Eigen3 installation")
	endif()
endif(USE_RBF)

add_library(mui INTERFACE)

target_include_directories (mui INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>)

if(USE_RBF)
    target_link_libraries(mui INTERFACE Eigen3::Eigen)
endif(USE_RBF)

target_compile_definitions(mui INTERFACE LIBRARY_HEADER_ONLY)

install(TARGETS mui EXPORT muiTargets INCLUDES DESTINATION include LIBRARY DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include
	FILES_MATCHING PATTERN "*.h"
	PATTERN ".git" EXCLUDE
	PATTERN "build" EXCLUDE)

install(EXPORT muiTargets DESTINATION ${CMAKE_INSTALL_PREFIX}/share/mui/cmake/)
install(FILES mui-config.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/share/mui/cmake)

if(C_WRAPPER)
	add_library(mui_c_wrapper SHARED
		./wrappers/C/mui_c_wrapper_general.cpp
		./wrappers/C/mui_c_wrapper_1d.cpp
		./wrappers/C/mui_c_wrapper_2d.cpp
		./wrappers/C/mui_c_wrapper_3d.cpp)
	set_target_properties(mui_c_wrapper PROPERTIES PUBLIC_HEADER ./wrappers/C/mui_c_wrapper_general.h)
	set_target_properties(mui_c_wrapper PROPERTIES PUBLIC_HEADER ./wrappers/C/mui_c_wrapper_1d.h)
	set_target_properties(mui_c_wrapper PROPERTIES PUBLIC_HEADER ./wrappers/C/mui_c_wrapper_2d.h)
	set_target_properties(mui_c_wrapper PROPERTIES PUBLIC_HEADER ./wrappers/C/mui_c_wrapper_3d.h)
	set_target_properties(mui_c_wrapper PROPERTIES PUBLIC_HEADER ./wrappers/C/config_c_wrapper.h)
	target_include_directories(mui_c_wrapper PRIVATE .)

	if(use_RBF)
    		target_link_libraries(mui_c_wrapper Eigen3::Eigen ${MPI_CXX_LIBRARIES})
	else()
    		target_link_libraries(mui_c_wrapper ${MPI_CXX_LIBRARIES})
	endif(use_RBF)

	install(TARGETS mui_c_wrapper
    		LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/
    		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/wrappers/C/)
endif(C_WRAPPER)

if(FORTRAN_WRAPPER)
endif(FORTRAN_WRAPPER)

if(PYTHON_WRAPPER)
endif(PYTHON_WRAPPER)
