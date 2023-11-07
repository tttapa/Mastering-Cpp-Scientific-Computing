#include <poly/interpolate.hpp>
#include <poly/poly.hpp>
#include <gtest/gtest.h>

TEST(poly, interpolateCheby) {
    // Generate some arbitrary data vectors
    const poly::index_t N = 5;
    poly::vector_t x      = poly::vector_t::LinSpaced(N, -1.0, 1.0);
    poly::vector_t y(N);
    y << 15.0, 3.5625, 1.0, 0.5625, 3.0;
    // Interpolate the data
    auto p = poly::interpolate_cheby(x, y);
    // Check the coefficients of the interpolant
    poly::ChebyshevPolynomial expected{4.375, -5.0, 4.0, -1.0, 0.625};
    auto error = p.coefficients - expected.coefficients;
    EXPECT_LE(error.norm(), 1e-12);
}