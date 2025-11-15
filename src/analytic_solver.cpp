#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include "matrix_utils.hpp"
#include "analytic_solver.hpp"

// ----------- analytic solver -----------
std::tuple<arma::vec, arma::mat>
solve_analytic(const int N, double d, double a)
{
    arma::vec eigvals(N);
    arma::mat eigvecs(N, N);    // columns = eigenvectors

    for (int j = 1; j <= N; ++j)
    {
        // j'th eigenvalue
        double arg = j * arma::datum::pi / (N + 1.0);
        eigvals(j - 1) = d + 2.0 * a * std::cos(arg);

        // j'th eigenvector components v_i = sin(i * arg)
        for (int i = 1; i <= N; ++i)
        {
            eigvecs(i - 1, j - 1) = std::sin(i * arg);
        }
        // normalizing j'th eigenvector
        eigvecs.col(j - 1) = arma::normalise(eigvecs.col(j - 1));
    }

    return std::make_tuple(eigvals, eigvecs);
}