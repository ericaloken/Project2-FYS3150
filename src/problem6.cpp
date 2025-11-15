/*
Problem 6
Solving Av = λv using the Jacobi algorithm implemented in
jacobi_solver.cpp for n = 10 and n = 100 discretization
steps and solving using analytical solver in analytical_solver.cpp.
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
    for (int n : {10, 100}) {

        // define A = tridiag(a,d,a)
        const int N = n - 1;                    // matrix size
        const double h = 1.0 / n;               // step size
        const double a = -1.0 / (h*h);          // sub-, superdiagonal
        const double d = 2.0 / (h*h);           // main diagonal
        
        arma::mat A = create_symmetric_tridiagonal(N, a, d);

        // parameters for solver
        double eps = 1e-8;
        int maxiter = 1000000;
        int iterations;
        bool converged;
        arma::vec eigenvalues(N);
        arma::mat eigenvectors(N, N);
    
        // Jacobi algorithm
        jacobi_eigensolver(A, eps, eigenvalues,eigenvectors, maxiter, iterations, converged);
        
        // analytical solution
        auto [eigenvalues_ana, eigenvectors_ana] = solve_analytic(N, d, a);
    
        // align signs in analytical and numerical eigenvecs by taking dot product of vectors (columns in matrix)
        for (int j = 0; j < N; ++j) {
                double c = arma::dot(eigenvectors.col(j), eigenvectors_ana.col(j));
                if (c < 0.0) eigenvectors_ana.col(j) *= -1.0;   // negative dot product = opposite directions
            }

        // writing first 3 eigenvectors to file with boundary conditions
        std::ofstream ofile("data/problem6_n" + std::to_string(n) + ".txt");
        ofile << std::scientific << std::setprecision(12);
        ofile << "x_i  v1_num  v1_ana  v2_num  v2_ana  v3_num  v3_ana\n" ;

        for (int i = 0; i <= n; ++i) {
            double x = i * h;
            ofile << x;

            // only interested in first 3 columns of matrix
            for (int j = 0; j < 3; ++j) {
                double vi_num;
                double vi_ana;
                if(i == 0 || i == n) {
                    vi_num = 0.0;
                    vi_ana = 0.0;
                } else {
                    vi_num = eigenvectors(i-1, j);
                    vi_ana = eigenvectors_ana(i-1, j);
                }
                ofile << " " << vi_num  << " " << vi_ana;
            }
            ofile << "\n";
        }
        ofile.close();
    }

    return 0;
}