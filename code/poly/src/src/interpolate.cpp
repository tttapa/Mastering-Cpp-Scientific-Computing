#include <poly/interpolate.hpp>
#include <Eigen/LU>
#include <cassert>

namespace poly {

namespace {
auto make_chebyshev_basis_matrix(vector_ref_t x, index_t degree) {
    assert(degree >= 0);
    const index_t N = x.size();
    matrix_t V(N, degree + 1);
    V.col(0) = vector_t::Ones(N);
    if (degree >= 1) {
        V.col(1) = x;
        for (index_t i = 0; i < degree - 1; ++i)
            V.col(i + 2) = 2 * V.col(i + 1).cwiseProduct(x) - V.col(i);
    }
    return V;
}
} // namespace

ChebyshevPolynomial interpolate_cheby(vector_ref_t x, vector_ref_t y) {
    assert(x.size() == y.size());
    assert(x.size() > 0);
    // Construct Vandermonde/basis matrix
    auto V = make_chebyshev_basis_matrix(x, x.size() - 1);
    // Scale the system
    const vector_t scaling = V.colwise().norm().cwiseInverse();
    V *= scaling.asDiagonal();
    // Solve the system, undo scaling
    ChebyshevPolynomial solution{V.fullPivLu().solve(y)};
    solution.coefficients.transpose() *= scaling.asDiagonal();
    return solution;
}

} // namespace poly