/*
This program has functions that makes matrices and a function that identifies
the largest off-diagonal element. 
*/

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "matrix_utils.hpp"

// --------------- creating tridiag(a, d, a) -----------------
arma::mat create_symmetric_tridiagonal(int N, double a, double d)
{
    assert(N > 1);
    arma::mat A = arma::mat(N, N, arma::fill::zeros);

    // fill main diagonal with d
    A.diag(0).fill(d);

    // fill super- and subdiagonal with a
    A.diag(1).fill(a);
    A.diag(-1).fill(a);
  
  return A;
}


// ------- function fining largest off-diagonal element -----------
double max_offdiag_symmetric(const arma::mat& A, int &k, int &l)
{
    assert(A.is_square());
    int N = A.n_rows;
    assert(N > 1);
    double maxval = 0.0;

    // consistency checks
    assert(A.is_square());     // ensure A is square
    assert(N > 1);             // ensure matrix is bigger than 1x1

    // examine elements in upper triangle (j>i) only (j = column, i = row)
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double elem = std::abs(A(i,j));
            if (elem > maxval) {
                maxval = elem;
                k = i;
                l = j;
            }
        }
    }
    return maxval;
}
