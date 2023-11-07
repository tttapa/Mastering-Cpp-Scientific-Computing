# Find dependencies we rely on
find_package(Eigen3 REQUIRED)
# Define the "poly" library and specify the source files it consists of
add_library(poly
    "include/poly/poly.hpp"
    "include/poly/interpolate.hpp"
    "src/interpolate.cpp"
)
# The library requires the 2020 revision of the C++ standard
target_compile_features(poly PUBLIC cxx_std_20)
# Specify the path of the header files
target_include_directories(poly PUBLIC "include")
# Link the Eigen library into our poly library
target_link_libraries(poly PUBLIC Eigen3::Eigen)
# Add custom compiler options or macros
if (POLY_INIT_NAN)
    target_compile_definitions(poly PUBLIC EIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()
