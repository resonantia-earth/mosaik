#include <Rcpp.h>
using namespace Rcpp;

//' Binary percolation grid (c++)
//'
//' Generate a binary grid where each cell is independently set to 1 with
//' probability \code{p} and 0 otherwise.
//' @param ncol [integer(1)][integer]\cr number of columns.
//' @param nrow [integer(1)][integer]\cr number of rows.
//' @param p [numeric(1)][numeric]\cr probability of a cell being 1.
//' @return An integer vector of length \code{ncol * nrow} with values 0 or 1.
//' @noRd
// [[Rcpp::export]]
IntegerVector percolationCpp(int ncol, int nrow, double p) {

  int n = ncol * nrow;
  IntegerVector out(n);

  for (int i = 0; i < n; i++) {
    out[i] = R::runif(0.0, 1.0) < p ? 1 : 0;
  }

  return out;
}
