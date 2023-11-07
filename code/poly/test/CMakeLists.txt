# Import the Google Test CMake module
include(GoogleTest)
# Find the Google Test library
find_package(GTest REQUIRED)
# Build the test program
add_executable(tests "test-interpolate.cpp")
# The tests depend on the library under test, and on the Google Test framework
target_link_libraries(tests PRIVATE poly GTest::gtest_main)
# Let CMake/CTest discover the available tests
gtest_discover_tests(tests)
