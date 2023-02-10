cmake_minimum_required(VERSION 3.18)
project(MUI VERSION 1.2.4 DESCRIPTION "Multiscale Universal Interface" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Library options
option(USE_RBF "Include RBF support" OFF)
# Wrapper options
option(C_WRAPPER "Compile C wrapper" OFF)
option(FORTRAN_WRAPPER "Compile Fortran wrapper" OFF)
option(PYTHON_WRAPPER "Compile and install Python wrapper using pip" OFF)

include(CMakePackageConfigHelpers)
include(CheckLanguage)

find_package(MPI REQUIRED)
if(MPI_FOUND)
	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
elseif(NOT MPI_FOUND)
	message(SEND_ERROR "MPI not found")
endif(MPI_FOUND)

if(USE_RBF)
	find_package(Eigen3 REQUIRED)
	if( NOT EIGEN3_INCLUDE_DIR )
    		message(FATAL_ERROR "Please set the environment variable EIGEN3_INCLUDE_DIR to a valid Eigen3 installation")
	endif()
endif(USE_RBF)

add_library(MUI INTERFACE)

target_include_directories (MUI INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>)

if(USE_RBF)
	target_link_libraries(MUI INTERFACE Eigen3::Eigen)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_RBF")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_RBF")
endif(USE_RBF)

target_compile_definitions(MUI INTERFACE LIBRARY_HEADER_ONLY)

install(TARGETS MUI EXPORT muiTargets INCLUDES DESTINATION include LIBRARY DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/
	DESTINATION ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/include
	FILES_MATCHING
	PATTERN "*.h"
	PATTERN "*.f90"
	PATTERN ".git" EXCLUDE
	PATTERN "build" EXCLUDE
	PATTERN "wrappers/Python" EXCLUDE
	PATTERN "cmake" EXCLUDE
	PATTERN "doc" EXCLUDE)

install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/lib)

if(C_WRAPPER)
	message("-- MUI C wrapper: Selected")
	
	check_language(C)
    if(CMAKE_C_COMPILER)
	    enable_language(C)
	    
        add_library(MUI_C_wrapper SHARED
		    "./wrappers/C/mui_c_wrapper_general.cpp"
		    "./wrappers/C/mui_c_wrapper_1d.cpp"
		    "./wrappers/C/mui_c_wrapper_2d.cpp"
		    "./wrappers/C/mui_c_wrapper_3d.cpp")
	
        set(MUI_C_wrapper_PUBLIC_HEADERS
            "./wrappers/C/mui_c_wrapper_general.h"
            "./wrappers/C/mui_c_wrapper_1d.h"
            "./wrappers/C/mui_c_wrapper_2d.h"
            "./wrappers/C/mui_c_wrapper_3d.h"
            "./wrappers/C/config_c_wrapper.h")

        set_target_properties(MUI_C_wrapper PROPERTIES PUBLIC_HEADER "${MUI_C_wrapper_PUBLIC_HEADERS}")
	
        target_include_directories(MUI_C_wrapper PRIVATE .)

        add_library(MUI_C_wrapper_static STATIC
            "./wrappers/C/mui_c_wrapper_general.cpp"
            "./wrappers/C/mui_c_wrapper_1d.cpp"
            "./wrappers/C/mui_c_wrapper_2d.cpp"
            "./wrappers/C/mui_c_wrapper_3d.cpp")
 
        set_target_properties(MUI_C_wrapper_static PROPERTIES OUTPUT_NAME MUI_C_wrapper)

        if(USE_RBF)
            target_link_libraries(MUI_C_wrapper Eigen3::Eigen ${MPI_CXX_LIBRARIES})
            target_link_libraries(MUI_C_wrapper_static Eigen3::Eigen ${MPI_CXX_LIBRARIES})
        else()
            target_link_libraries(MUI_C_wrapper ${MPI_CXX_LIBRARIES})
            target_link_libraries(MUI_C_wrapper_static ${MPI_CXX_LIBRARIES})
        endif(USE_RBF)

        install(TARGETS MUI_C_wrapper MUI_C_wrapper_static DESTINATION ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/lib
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/include/wrappers/C)
	else()
        message(STATUS "No C support")
    endif()
endif(C_WRAPPER)

if(FORTRAN_WRAPPER)
	message("-- MUI Fortran wrapper: Selected")
	
	check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
        
        set_source_files_properties(
  		    "./wrappers/Fortran/mui_f_wrapper_general.f90"
		    "./wrappers/Fortran/mui_f_wrapper_1d.f90"
		    "./wrappers/Fortran/mui_f_wrapper_2d.f90"
		    "./wrappers/Fortran/mui_f_wrapper_3d.f90"
  		    PROPERTIES Fortran_PREPROCESS ON)

	    add_library(MUI_Fortran_wrapper SHARED
		    "./wrappers/Fortran/mui_f_wrapper_general.cpp"
		    "./wrappers/Fortran/mui_f_wrapper_1d.cpp"
		    "./wrappers/Fortran/mui_f_wrapper_2d.cpp"
		    "./wrappers/Fortran/mui_f_wrapper_3d.cpp"
		    "./wrappers/Fortran/mui_f_wrapper_general.f90"
            "./wrappers/Fortran/mui_f_wrapper_1d.f90"
            "./wrappers/Fortran/mui_f_wrapper_2d.f90"
            "./wrappers/Fortran/mui_f_wrapper_3d.f90")
            
        set(MUI_Fortran_wrapper_PUBLIC_HEADERS
    	    "./wrappers/Fortran/mui_f_wrapper_general.f90"
		    "./wrappers/Fortran/mui_f_wrapper_1d.f90"
		    "./wrappers/Fortran/mui_f_wrapper_2d.f90"
		    "./wrappers/Fortran/mui_f_wrapper_3d.f90"
		    "./wrappers/Fortran/config_f_wrapper.h")
	
	    set_target_properties(MUI_Fortran_wrapper PROPERTIES PUBLIC_HEADER "${MUI_Fortran_wrapper_PUBLIC_HEADERS}")
	
	    set_target_properties(MUI_Fortran_wrapper PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/lib/mod)
	
	    target_include_directories(MUI_Fortran_wrapper PRIVATE .)

	    add_library(MUI_Fortran_wrapper_static STATIC
    	    "./wrappers/Fortran/mui_f_wrapper_general.cpp"
            "./wrappers/Fortran/mui_f_wrapper_1d.cpp"
            "./wrappers/Fortran/mui_f_wrapper_2d.cpp"
            "./wrappers/Fortran/mui_f_wrapper_3d.cpp")
        
	    set_target_properties(MUI_Fortran_wrapper_static PROPERTIES OUTPUT_NAME MUI_Fortran_wrapper)
	
	    target_include_directories(MUI_Fortran_wrapper_static PRIVATE .)

	    if(USE_RBF)
            target_link_libraries(MUI_Fortran_wrapper Eigen3::Eigen ${MPI_CXX_LIBRARIES})
            target_link_libraries(MUI_Fortran_wrapper_static Eigen3::Eigen ${MPI_CXX_LIBRARIES})
	    else()
            target_link_libraries(MUI_Fortran_wrapper ${MPI_CXX_LIBRARIES})
            target_link_libraries(MUI_Fortran_wrapper_static ${MPI_CXX_LIBRARIES})
        endif(USE_RBF)

        install(TARGETS MUI_Fortran_wrapper MUI_Fortran_wrapper_static DESTINATION ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/lib
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/include/wrappers/Fortran)
    else()
        message(STATUS "No Fortran support")
    endif()
endif(FORTRAN_WRAPPER)

if(PYTHON_WRAPPER)
	message("-- MUI Python wrapper: Selected (Warning - high compilation time and memory use)")
	
	find_package(Python COMPONENTS Interpreter REQUIRED)

	if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})	   
	    set(_pip_args)  
	else()
	    set(_pip_args "--user")	
	endif()
	
	set(_pkg_mui4py "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/Python")

	if(USE_RBF)
	    add_custom_target(mui4py ALL COMMAND bash -c "CMAKE_ARGS=\"-DUSE_RBF=ON -DEigen3_DIR=${Eigen3_DIR}\" ${Python_EXECUTABLE} -m pip install ${_pip_args} ${_pkg_mui4py}/")
	else()
	    add_custom_target(mui4py ALL COMMAND bash -c ${Python_EXECUTABLE} -m pip install ${_pip_args} ${_pkg_mui4py}/)
	endif(USE_RBF)

	set_property(TARGET mui4py APPEND PROPERTY ADDITIONAL_CLEAN_FILES
	             ${_pkg_mui4py}/mui4py.egg-info/
				 ${_pkg_mui4py}/dist/
				 ${_pkg_mui4py}/build/)
endif(PYTHON_WRAPPER)

set(MUI_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/include CACHE STRING "Location of MUI header files")
set(MUI_LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION}/lib CACHE STRING "Location of MUI wrapper libraries")
set(MUI_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/MUI-${PROJECT_VERSION} CACHE STRING "Base location of MUI installation")

configure_package_config_file(cmake/MUIConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MUIConfig.cmake
  INSTALL_DESTINATION ${MUI_INSTALL_DIR}/share/MUI-${PROJECT_VERSION}/cmake
  PATH_VARS MUI_INCLUDE_INSTALL_DIR MUI_LIB_INSTALL_DIR MUI_INSTALL_DIR)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MUIConfig.cmake
        DESTINATION ${MUI_INSTALL_DIR}/share/MUI-${PROJECT_VERSION}/cmake)
