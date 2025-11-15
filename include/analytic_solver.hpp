#include <armadillo>
#include <tuple>
#pragma once

// find eigenvalues and eigenvectors of tridiag matrix analytically
std::tuple<arma::vec, arma::mat>
solve_analytic(const int N, double d, double a);