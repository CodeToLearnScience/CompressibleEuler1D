# =============================================================================
# External Dependencies via FetchContent
# =============================================================================

include(FetchContent)

# -----------------------------------------------------------------------------
# toml++ - TOML parser for C++17 and later
# https://github.com/marzer/tomlplusplus
# -----------------------------------------------------------------------------

FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
    GIT_SHALLOW    TRUE
)

# -----------------------------------------------------------------------------
# GoogleTest - Testing framework
# https://github.com/google/googletest
# -----------------------------------------------------------------------------

if(EULER1D_BUILD_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
        GIT_SHALLOW    TRUE
    )
    
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    
    # Disable installation of googletest
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
endif()

# -----------------------------------------------------------------------------
# Make dependencies available
# -----------------------------------------------------------------------------

message(STATUS "Fetching dependencies...")

FetchContent_MakeAvailable(tomlplusplus)

if(EULER1D_BUILD_TESTS)
    FetchContent_MakeAvailable(googletest)
endif()

message(STATUS "Dependencies fetched successfully")
