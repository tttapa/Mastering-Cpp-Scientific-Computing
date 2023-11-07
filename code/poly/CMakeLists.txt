cmake_minimum_required(VERSION 3.21)
project(poly VERSION 1.2.3 LANGUAGES CXX)
# Import CMake testing support
include(CTest)

# Define user options
option(POLY_INIT_NAN
    "Initialize all matrices and vectors to NaN for debugging purposes" Off)

# Include the build script in the src directory
add_subdirectory("src")
# Include the build script in the test (if enabled)
if (BUILD_TESTING)
    add_subdirectory("test")
endif()