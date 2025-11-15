/*
Program that tests the Jacobi rotation algorithm implemented
in jacobi_solver.cpp with N = 6 (6x6 matrix).
*/

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "matrix_utils.hpp"
#include "jacobi_solver.hpp"
#include "analytic_solver.hpp"

int main()
{
    // define A = tridiag(a,d,a)
    const int N = 6;
    const int n = N + 1;              // discretization steps
    const double h = 1.0 / n;         // step size
    const double a = -1.0 / (h*h);    // sub-, superdiagonal
    const double d = 2.0 / (h*h);     // main diagonal
    
    arma::mat A = create_symmetric_tridiagonal(N, a, d);
    A.print("A = ");

    // parameters for solver
    double eps = 1e-8;
    int maxiter = 10000;
    int iterations;
    bool converged;
    arma::vec eigenvalues(N);
    arma::mat eigenvectors(N, N);

    // perform Jacobi algorithm
    jacobi_eigensolver(A, eps, eigenvalues,eigenvectors, maxiter, iterations, converged);

    // analytical solution
    auto [eigenvalues_ana, eigenvectors_ana] = solve_analytic(N, d, a);

    // align signs in analytical and numerical eigenvecs by taking dot product of vectors (columns in matrix)
    for (int j = 0; j < N; ++j) {
            double c = arma::dot(eigenvectors.col(j), eigenvectors_ana.col(j));
            if (c < 0.0) eigenvectors_ana.col(j) *= -1.0;   // negative dot product = opposite directions
        }

    // absolute error and max absolute errors
    arma::vec abs_err_vals = arma::abs(eigenvalues - eigenvalues_ana);
    arma::mat abs_err_vecs = arma::abs(eigenvectors - eigenvectors_ana);
    double max_err_vals = abs_err_vals.max();
    double max_err_vecs = abs_err_vecs.max();


    // print results
    std::cout << "-----------------------------------------------------------------------------\n";
    std::cout << "Jacobi rotation algorithm results:\n";
    std::cout << "Eigenvalues:\n"                           << eigenvalues.t();
    std::cout << "Eigenvectors (columns:):\n"               << eigenvectors;
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "Analytical results:\n";
    std::cout << "Eigenvalues:\n"                           << eigenvalues_ana.t();
    std::cout << "Eigenvectors (columns):\n"                << eigenvectors_ana;
    std::cout << "Absolute error eigenvalues:\n"            << abs_err_vals.t();
    std::cout << "Absolute error eigenvectors:\n"           << abs_err_vecs;
    std::cout << "Max |λ_jacobi - λ_ana| ="                    << " " << max_err_vals << "\n";
    std::cout << "Max |v_jacobi - v_ana| ="                    << " " << max_err_vecs << "\n";

    return 0;
}