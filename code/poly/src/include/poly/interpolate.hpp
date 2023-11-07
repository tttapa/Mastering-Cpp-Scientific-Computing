#pragma once

#include <poly/poly.hpp>

namespace poly {

/// Given the vectors @f$ x @f$ and @f$ y @f$, interpolate a Chebyshev
/// polynomial @f$ p @f$ such that @f$ p(x_i) = y_i @f$.
/// @param  x   Samples of the independent variable.
/// @param  y   Samples of the dependent variable corresponding to @p x.
/// @pre    @p x and @p y should have the same length.
/// @return The coefficients of polynomial @f$ p @f$ in Chebyshev basis.
ChebyshevPolynomial interpolate_cheby(vector_ref_t x, vector_ref_t y);

} // namespace poly