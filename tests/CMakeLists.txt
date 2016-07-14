include_directories(${gtest_SOURCE_DIR}/include)
add_executable(run_dummy_tests dummy.cpp)

target_link_libraries(run_dummy_tests ${CMAKE_BINARY_DIR}/tests/googletest/googletest/libgtest_main.a)
target_link_libraries(run_dummy_tests ${CMAKE_BINARY_DIR}/tests/googletest/googletest/libgtest.a)
target_link_libraries(run_dummy_tests ${CMAKE_THREAD_LIBS_INIT})

add_test(dummy_tests run_dummy_tests)
