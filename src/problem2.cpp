/*
Problem 2
Program that sets up the tridiagonal matrix A for N = 6 and solves
Av = λv using arma::eig_sym. It also checks that the eigenvalues and eigenvectors
from Armadillo agrees with the analytical result for  N = 6 by.
*/

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include "matrix_utils.hpp"
#include "analytic_solver.hpp"


// ----------- numerical solver -----------
std::tuple<arma::vec, arma::mat>
solve_numeric(const arma::mat& A)
{
    arma::vec eigvals;
    arma::mat eigvecs;          // column = eigenvectors
    arma::eig_sym(eigvals, eigvecs, A);

    return std::make_tuple(eigvals, eigvecs);
}


int main()
{
    // parameters
    const int N = 6;                  // (NxN) matrix
    const int n = N + 1;              // discretization steps
    const double h = 1.0 / n;         // step size
    const double a = -1.0 / (h*h);    // sub-, superdiagonal
    const double d = 2.0 / (h*h);     // main diagonal
 
    // create and print A
    arma::mat A = create_symmetric_tridiagonal(N, a, d);
    A.print("A = ");

    // numerical and analytical solutions
    auto [eigvals_num, eigvecs_num] = solve_numeric(A);
    auto [eigvals_ana, eigvecs_ana] = solve_analytic(N, d, a);

    // align signs in analytical and numerical eigenvecs by taking dot product of vectors (columns in matrix)
    for (int j = 1; j < N; ++j) {
            double c = arma::dot(eigvecs_num.col(j), eigvecs_ana.col(j));
            if (c < 0.0) eigvecs_ana.col(j) *= -1.0;   // negative dot product = opposite directions
        }

    // absolute error and max absolute errors
    arma::vec abs_err_vals = arma::abs(eigvals_num - eigvals_ana);
    arma::mat abs_err_vecs = arma::abs(eigvecs_num - eigvecs_ana);
    double max_err_vals = abs_err_vals.max();
    double max_err_vecs = abs_err_vecs.max();

    // printing
    std::cout << std::setprecision(12);
    std::cout << "Numerical eigenvalues:\n"                 << eigvals_num.t();
    std::cout << "Numerical eigenvectors (columns):\n"      << eigvecs_num;
    std::cout << "Analytical eigenvalues:\n"                << eigvals_ana.t();
    std::cout << "Analytical eigenvectors (columns):\n"     << eigvecs_ana;
    std::cout << "Absolute error eigenvalues:\n"            << abs_err_vals.t();
    std::cout << "Absolute error eigenvectors:\n"           << abs_err_vecs;
    std::cout << "Max |λ_num - λ_ana| ="                    << " " << max_err_vals << "\n";
    std::cout << "Max |v_num - v_ana| ="                    << " " << max_err_vecs << "\n";

    return 0;
}