#include <armadillo>
#pragma once


// create symmetric tridiagonal (a,d,a) matrix of size (NxN)
arma::mat create_symmetric_tridiagonal(int N, double a, double d);

// find largest off-diaonal element (in abs value)
double max_offdiag_symmetric(const arma::mat& A, int &k, int &l);

