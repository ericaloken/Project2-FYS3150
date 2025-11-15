/*
Problem 3b
Program that tests our function max_offdiag_symmetric on a
symmetric matrix A.
*/

#include <armadillo>
#include <iostream>
#include <iomanip>
#include "matrix_utils.hpp"
#include <cmath>
#include <cassert>

int main()
{
    // building a (4x4) symmetric matrix
    int N = 4;
    arma::mat A(N, N, arma::fill::zeros);
    A.diag().ones();
    A(0,3) = A(3,0) = 0.5;
    A(2,1) = A(1,2) = -0.7;
    A.print("A = ");

    // finding largest off-diagonal element
    int k, l;
    double maxval = max_offdiag_symmetric(A, k, l);

    std::cout << "Largest absolute value of offdiagonal element is " << maxval
              << " with sign A(k,l) = " << A(k,l)
              << " at (" << k << "," << l << ") and (" << l << "," << k << ")\n";
}