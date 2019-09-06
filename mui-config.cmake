include(CMakeFindDependencyMacro) 
find_dependency(Eigen3) 
get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include("${SELF_DIR}/muiTargets.cmake")
