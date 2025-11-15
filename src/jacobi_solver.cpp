/*
Implementation of Jacobi's rotation algorithm to solve eigenvalue problem
Av = λv. A = tridiag(a,d,a), a = -1/h^2, d = 2/h^2
*/

#include <armadillo>
#include <iostream>
#include <iomanip>
#include "matrix_utils.hpp"
#include "jacobi_solver.hpp"
#include <cmath>
#include <cassert>


// -------------- function performing one Jacobi rotation -------------
void jacobi_rotate(arma::mat& A, arma::mat& R, int& k, int& l)
{
    int N = A.n_rows;
    assert(k < l && "We require k < l");
    
    // compute tau
    double tau = (A(l,l) - A(k,k)) / (2*A(k,l));

    // compute tanθ (t), cosθ (c), sinθ (s)
    double t;
    if (tau > 0) {
        t = 1.0 / (tau + std::sqrt(1 + tau*tau));
    } else {
        t = 1.0 / (tau - std::sqrt(1 + tau*tau));
    }

    double c = 1.0 / std::sqrt(1 + t*t);
    double s = c * t;

    // remember old elements in A^(m)
    double a_kk = A(k,k);
    double a_ll = A(l,l);
    double a_kl = A(k,l);

    // update A^(m) -> A^(m+1) starting with k,l elements
    A(k,k) = a_kk * c*c - 2.0 * a_kl * c*s + a_ll * s*s;
    A(l,l) = a_ll * c*c + 2.0 * a_kl * c*s + a_kk * s*s;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    // update other elements
    for (int i = 0; i < N; ++i) {
        if (i != k && i != l) {
            double a_ik = A(i,k);
            double a_il = A(i,l);

            // new elements
            A(i,k) = a_ik * c - a_il * s;
            A(k,i) = A(i,k);
            A(i,l) = a_il * c + a_ik * s;
            A(l,i) = A(i,l);
        }
    }
    // update rotation matrix R^(m) -> R^(m+1)
    for (int i = 0; i < N; ++i) {
        // old elements
        double r_ik = R(i,k);
        double r_il = R(i,l);

        // new elements
        R(i,k) = r_ik * c - r_il * s;
        R(i,l) = r_il * c + r_ik * s;
    }
}


// -------------- function rotating until diagonal within tol -------------
void jacobi_eigensolver(const arma::mat& A_in, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged)
{
    int N = A_in.n_rows;

    // initialize: A^(1) = A_in, R^(1) = I
    arma::mat A = A_in;
    arma::mat R = arma::mat(N, N, arma::fill::eye);

    // indices (k,l) of maximum (abs value) off-diag element
    int k, l;
    double maxval = max_offdiag_symmetric(A, k, l);
    
    // rotate until maxval < eps
    iterations = 0;
    while (maxval > eps && iterations < maxiter) {
        jacobi_rotate(A, R, k, l);
        maxval = max_offdiag_symmetric(A, k, l);
        iterations++;
    }

    // check if solver stopped within tolerance
    converged = (maxval <= eps);
    if (converged){
        std::cout << "Jacobi algorithm converged after " << iterations << " iterations.\n";
    } else {
        std::cout << "Jacobi algorithm did not converge after max iterations.\n";
    }

    // save eigenvalues from diagonal of A in a vector
    eigenvalues = A.diag();

    // save eigenvectors from columns in R in a matrix
    eigenvectors = R;

    // sort eigenvalues and eigenvectors by ascending order
    arma::uvec indices = arma::sort_index(eigenvalues);
    eigenvalues = eigenvalues(indices);
    eigenvectors = eigenvectors.cols(indices);
}