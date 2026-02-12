cmake_minimum_required(VERSION 3.25)

# Project definition
project(CompressibleEuler1D
    VERSION 1.0.0
    DESCRIPTION "1D Compressible Euler Finite Volume Solver"
    LANGUAGES CXX
)

# -----------------------------------------------------------------------------
# Build Options
# -----------------------------------------------------------------------------

# Precision selection: double (default) or float
set(EULER1D_PRECISION "double" CACHE STRING "Floating point precision (double or float)")
set_property(CACHE EULER1D_PRECISION PROPERTY STRINGS "double" "float")

# Build options
option(EULER1D_BUILD_TESTS "Build unit tests" ON)
option(EULER1D_BUILD_DOCS "Build documentation" OFF)

# Validate precision option
if(NOT EULER1D_PRECISION STREQUAL "double" AND NOT EULER1D_PRECISION STREQUAL "float")
    message(FATAL_ERROR "EULER1D_PRECISION must be 'double' or 'float', got: ${EULER1D_PRECISION}")
endif()

message(STATUS "Euler1D precision: ${EULER1D_PRECISION}")

# -----------------------------------------------------------------------------
# C++ Standard and Compiler Settings
# -----------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# -----------------------------------------------------------------------------
# Compiler Warnings
# -----------------------------------------------------------------------------

# Create interface library for compiler warnings
add_library(euler1d_warnings INTERFACE)
target_compile_options(euler1d_warnings INTERFACE
    $<$<CXX_COMPILER_ID:GNU>:
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wimplicit-fallthrough
    >
    $<$<CXX_COMPILER_ID:Clang>:
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wnull-dereference -Wdouble-promotion -Wformat=2
    >
    $<$<CXX_COMPILER_ID:MSVC>:
        /W4 /permissive-
    >
)

# -----------------------------------------------------------------------------
# Optimization flags
# -----------------------------------------------------------------------------

add_library(euler1d_optimize INTERFACE)
target_compile_options(euler1d_optimize INTERFACE
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Release>>:
        -O3 -march=native -funroll-loops -ffast-math
    >
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CONFIG:Release>>:
        -O3 -march=native -funroll-loops -ffast-math
    >
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:
        /O2 /fp:fast
    >
)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

include(cmake/dependencies.cmake)

# -----------------------------------------------------------------------------
# Main Library
# -----------------------------------------------------------------------------

add_library(euler1d_lib STATIC
    # Config
    src/config/parser.cpp
    # EOS
    src/eos/eos.cpp
    # Mesh
    src/mesh/mesh.cpp
    # Flux
    src/flux/llf.cpp
    src/flux/rusanov.cpp
    src/flux/hll.cpp
    src/flux/hllc.cpp
    # Reconstruction
    src/reconstruction/limiter.cpp
    src/reconstruction/muscl.cpp
    # Boundary
    src/boundary/boundary.cpp
    # Initial conditions
    src/initial/initial_condition.cpp
    # Time integration
    src/time/explicit_euler.cpp
    src/time/ssprk3.cpp
    # Solver
    src/solver/solver.cpp
    src/solver/factory.cpp
    # I/O
    src/io/csv_writer.cpp
    src/io/vtk_writer.cpp
)

target_include_directories(euler1d_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(euler1d_lib
    PUBLIC
        tomlplusplus::tomlplusplus
    PRIVATE
        euler1d_warnings
        euler1d_optimize
)

# Set precision preprocessor definition
target_compile_definitions(euler1d_lib
    PUBLIC
        EULER1D_PRECISION=${EULER1D_PRECISION}
)

# Alias for consistent naming
add_library(euler1d::lib ALIAS euler1d_lib)

# -----------------------------------------------------------------------------
# Main Executable
# -----------------------------------------------------------------------------

add_executable(euler1d apps/main.cpp)

target_link_libraries(euler1d
    PRIVATE
        euler1d_lib
        euler1d_warnings
        euler1d_optimize
)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

if(EULER1D_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------

include(GNUInstallDirs)

install(TARGETS euler1d euler1d_lib
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/euler1d
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

message(STATUS "")
message(STATUS "=== CompressibleEuler1D Configuration ===")
message(STATUS "  Version:     ${PROJECT_VERSION}")
message(STATUS "  Precision:   ${EULER1D_PRECISION}")
message(STATUS "  Build type:  ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Tests:       ${EULER1D_BUILD_TESTS}")
message(STATUS "  Compiler:    ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "==========================================")
message(STATUS "")
