#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <initializer_list>
#include <utility>

namespace poly {

// Handy type aliases to make the rest of the code more readable.
using matrix_t         = Eigen::MatrixX<double>;
using vector_t         = Eigen::VectorX<double>;
using vector_ref_t     = Eigen::Ref<const vector_t>;
using vector_mut_ref_t = Eigen::Ref<vector_t>;
using index_t          = Eigen::Index;

/// Representation of a polynomial in Chebyshev basis.
struct ChebyshevPolynomial {
    /// Type of the polynomial coefficients.
    using coeff_t = vector_t;

    /// Create an empty polynomial.
    ChebyshevPolynomial() = default;
    /// Create a polynomial with the given coefficients.
    explicit ChebyshevPolynomial(coeff_t init_coeff)
        : coefficients{std::move(init_coeff)} {}
    /// Create a polynomial with all zero coefficients of the given degree.
    explicit ChebyshevPolynomial(index_t degree)
        : coefficients{coeff_t::Zero(degree + 1)} {}
    /// Create a polynomial with the given coefficients.
    ChebyshevPolynomial(std::initializer_list<double> init_coeff)
        : coefficients{init_coeff.size()} {
        std::ranges::copy(init_coeff, std::begin(this->coefficients));
    }
    /// Vector storing the polynomial coefficients in Chebyshev basis.
    coeff_t coefficients;
};

} // namespace poly