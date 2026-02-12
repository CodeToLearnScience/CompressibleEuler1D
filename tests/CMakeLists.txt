# Tests CMakeLists.txt

include(GoogleTest)

# Test executable
add_executable(euler1d_tests
    test_eos.cpp
    test_config_parser.cpp
    test_flux.cpp
    test_limiter.cpp
    test_reconstruction.cpp
    test_boundary.cpp
    test_initial_condition.cpp
    test_time_integrator.cpp
    test_solver_integration.cpp
)

target_link_libraries(euler1d_tests
    PRIVATE
        euler1d_lib
        GTest::gtest
        GTest::gtest_main
)

# Discover tests
gtest_discover_tests(euler1d_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
