cmake_minimum_required(VERSION 3.9)
project(MUI VERSION 1.0 DESCRIPTION "Multiscale Universal Interface" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(sourcelist.cmake)

find_package(MPI REQUIRED)
if (MPI_FOUND)
	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
	message(SEND_ERROR "MPI not found")
endif (MPI_FOUND)

option(use_RBF OFF)

if(use_RBF)
    find_package(Eigen3 REQUIRED)
endif(use_RBF)

add_library(mui INTERFACE)

target_include_directories (mui INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>
)

if(use_RBF)
    target_link_libraries(mui INTERFACE Eigen3::Eigen)
endif(use_RBF)

target_compile_definitions(mui INTERFACE LIBRARY_HEADER_ONLY)

install( TARGETS mui EXPORT muiTargets INCLUDES DESTINATION include LIBRARY DESTINATION lib )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/ 
	     DESTINATION ${CMAKE_INSTALL_PREFIX}/include 
		 FILES_MATCHING PATTERN "*.h"
	   )

install( EXPORT muiTargets DESTINATION ${CMAKE_INSTALL_PREFIX}/share/mui/cmake/)
install(FILES mui-config.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/share/mui/cmake)


add_library(wrapper_c SHARED ./wrappers/C/mui_3d.cpp)
set_target_properties(wrapper_c PROPERTIES PUBLIC_HEADER ./wrappers/C/mui_3d.h)
target_include_directories(wrapper_c PRIVATE .)
if(use_RBF)
    target_link_libraries(wrapper_c Eigen3::Eigen ${MPI_CXX_LIBRARIES})
else()
    target_link_libraries(wrapper_c ${MPI_CXX_LIBRARIES})
endif(use_RBF)
install(TARGETS wrapper_c
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/wrappers/C/)


add_library(wrapper_f SHARED ./wrappers/Fortran/mui_3df.cpp)
target_include_directories(wrapper_f PRIVATE .)
if(use_RBF)
    target_link_libraries(wrapper_f Eigen3::Eigen)
endif(use_RBF)

install(TARGETS wrapper_f
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
