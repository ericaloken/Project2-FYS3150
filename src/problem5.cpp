/*
Problem 5
Program studying how many similarity transformations we need 
to obtain a diagonal matrix. We will run the Jacibi algorithm from
jacobi_solver.cpp with different choices of N and estimate how the 
number of transformations needed scaled with matrix size N.
N and iterations are written to a .txt file "problem5.txt".
*/

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "matrix_utils.hpp"
#include "jacobi_solver.hpp"
#include <fstream>
#include <vector>

int main() 
{
    std::vector<int> Nvals = {10, 20, 40, 60, 80, 100};
    double eps = 1e-8;
    int maxiter = 100000;

    // making file for data
    std::ofstream ofile("data/problem5.txt");
    ofile << std::scientific << std::setprecision(12);
    ofile << "N     Iterations         Converged (1 = yes, 0 = no)\n";

    for (int N : Nvals) {
        const int n = N + 1;              // discretization steps
        const double h = 1.0 / n;
        const double a = -1.0 / (h*h);
        const double d = 2.0 / (h*h);

        arma::mat A = create_symmetric_tridiagonal(N, a, d);
        
        // parameters for solver
        arma::vec eigenvalues(N);
        arma::mat eigenvectors(N, N);
        int iterations;
        bool converged;

        // perform Jacobi algorithm
        jacobi_eigensolver(A, eps, eigenvalues,eigenvectors, maxiter, iterations, converged);

        // writing N and iterations to txt file
        ofile << N << "         " << iterations << "          " << converged << "\n";
    }
    ofile.close();

    return 0;
}